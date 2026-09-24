/*!
 * \file CParaviewBinaryFileWriter.cpp
 * \brief Filewriter class for Paraview binary format.
 * \author T. Albring
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/output/filewriter/CParaviewBinaryFileWriter.hpp"
#include "../../../../Common/include/toolboxes/SwapBytes.hpp"
#include <cstdint>
#include <limits>

const string CParaviewBinaryFileWriter::fileExt = ".vtk";

CParaviewBinaryFileWriter::CParaviewBinaryFileWriter(CParallelDataSorter *valDataSorter) :
  CFileWriter(valDataSorter, fileExt){

  /* Check for big endian. We have to swap bytes otherwise.
   * Since size of character is 1 byte when the character pointer
   *  is de-referenced it will contain only first byte of integer. ---*/

  bigEndian = false;
  unsigned int i = 1;
  char *c = (char*)&i;
  bigEndian = *c == 0;
}


CParaviewBinaryFileWriter::~CParaviewBinaryFileWriter()= default;

void CParaviewBinaryFileWriter::WriteData(string val_filename){

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  const vector<string>& fieldNames = dataSorter->GetFieldNames();

  unsigned short iDim = 0, nDim = dataSorter->GetnDim();

  unsigned long iPoint;

  const int MAX_STRING_LENGTH = 255;
  char str_buf[MAX_STRING_LENGTH];

  const int NCOORDS = 3;

  OpenMPIFile(val_filename);

  /*--- The classic (3.0) cell layout stores the cell sizes and node ids of all cells in one Int32 array. When that
   array does not fit in Int32, use the 5.1 layout (VTK >= 9.0) with separate Int64 offsets and connectivity. ---*/

  const unsigned long GlobalCellStorage = dataSorter->GetnConnGlobal() + dataSorter->GetnElemGlobal();
  const bool cellsInt64 = GlobalCellStorage > static_cast<unsigned long>(std::numeric_limits<int32_t>::max());

  string header = string("# vtk DataFile Version ") + (cellsInt64 ? "5.1" : "3.0") + "\n"
                  "vtk output\n"
                  "BINARY\n"
                  "DATASET UNSTRUCTURED_GRID\n";

  WriteMPIString(header, MASTER_NODE);

  /*--- Communicate the number of total points that will be
   written by each rank. After this communication, each proc knows how
   many poinnts will be written before its location in the file and the
   offsets can be correctly set. ---*/

  unsigned long myPoint, GlobalPoint;

  GlobalPoint = dataSorter->GetnPointsGlobal();
  myPoint     = dataSorter->GetnPoints();

  WriteMPIString("POINTS " + std::to_string(GlobalPoint) + " float\n", MASTER_NODE);

  /*--- Load/write the 1D buffer of point coordinates. Note that we
   always have 3 coordinate dimensions, even for 2D problems. ---*/

  vector<float> dataBufferFloat(myPoint*NCOORDS);
  for (iPoint = 0; iPoint < myPoint; iPoint++) {
    for (iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        dataBufferFloat[iPoint*NCOORDS + iDim] = 0.0;
      } else {
        auto val = (float)dataSorter->GetData(iDim, iPoint);
        dataBufferFloat[iPoint*NCOORDS + iDim] = val;
      }
    }
  }

  if (!bigEndian) SwapBytes((char *)dataBufferFloat.data(), sizeof(float), myPoint*NCOORDS);

  /*--- Compute various data sizes --- */

  unsigned long sizeInBytesPerPoint = sizeof(float)*NCOORDS;
  unsigned long sizeInBytesLocal    = sizeInBytesPerPoint*myPoint;
  unsigned long sizeInBytesGlobal   = sizeInBytesPerPoint*GlobalPoint;
  unsigned long offsetInBytes       = sizeInBytesPerPoint*dataSorter->GetnPointCumulative(rank);

  WriteMPIBinaryDataAll(dataBufferFloat.data(), sizeInBytesLocal, sizeInBytesGlobal, offsetInBytes);

  /*--- Compute our local number of elements, the required storage,
   and reduce the total number of elements and storage globally. ---*/

  unsigned long myElem, myElemStorage, GlobalElem, GlobalElemStorage;

  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);

  myElem            = dataSorter->GetnElem();
  myElemStorage     = dataSorter->GetnConn();
  GlobalElem        = dataSorter->GetnElemGlobal();
  GlobalElemStorage = dataSorter->GetnConnGlobal();

  /*--- Loop over the local elements of each type, calling f(type, iElem, nPoints). ---*/

  auto forEachElem = [&](auto f) {
    for (auto type : {LINE, TRIANGLE, QUADRILATERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID}) {
      const auto nPoints = nPointsOfElementType(type);
      for (unsigned long iElem = 0; iElem < dataSorter->GetnElem(type); iElem++) f(type, iElem, nPoints);
    }
  };
  unsigned long iStorage = 0;

  if (!cellsInt64) {

    WriteMPIString("\nCELLS " + std::to_string(GlobalElem) + " " + std::to_string(GlobalCellStorage) + "\n",
                   MASTER_NODE);

    /*--- Load/write the 1D buffer of the number of nodes followed by the node ids of each cell. ---*/

    vector<int32_t> connBuf(myElemStorage + myElem);

    forEachElem([&](GEO_TYPE type, unsigned long iElem, unsigned short nPoints) {
      connBuf[iStorage++] = nPoints;
      for (unsigned short iNode = 0; iNode < nPoints; iNode++)
        connBuf[iStorage++] = static_cast<int32_t>(dataSorter->GetElemConnectivity(type, iElem, iNode) - 1);
    });

    if (!bigEndian) SwapBytes((char *)connBuf.data(), sizeof(int32_t), myElemStorage+myElem);

    sizeInBytesPerPoint = sizeof(int32_t);
    sizeInBytesLocal    = sizeInBytesPerPoint*(myElemStorage + myElem);
    sizeInBytesGlobal   = sizeInBytesPerPoint*GlobalCellStorage;
    offsetInBytes       = sizeInBytesPerPoint*
                          (dataSorter->GetnElemConnCumulative(rank) + dataSorter->GetnElemCumulative(rank));

    WriteMPIBinaryDataAll(connBuf.data(), sizeInBytesLocal, sizeInBytesGlobal, offsetInBytes);

  } else {

    WriteMPIString("\nCELLS " + std::to_string(GlobalElem + 1) + " " + std::to_string(GlobalElemStorage) + "\n",
                   MASTER_NODE);

    /*--- Load the offsets (where each cell ends in the connectivity) and the connectivity. ---*/

    vector<int64_t> offsetBuf(myElem), connBuf(myElemStorage);
    unsigned long iCell = 0;

    forEachElem([&](GEO_TYPE type, unsigned long iElem, unsigned short nPoints) {
      for (unsigned short iNode = 0; iNode < nPoints; iNode++)
        connBuf[iStorage++] = static_cast<int64_t>(dataSorter->GetElemConnectivity(type, iElem, iNode)) - 1;
      offsetBuf[iCell++] = static_cast<int64_t>(iStorage + dataSorter->GetnElemConnCumulative(rank));
    });

    if (!bigEndian) {
      SwapBytes((char *)offsetBuf.data(), sizeof(int64_t), myElem);
      SwapBytes((char *)connBuf.data(), sizeof(int64_t), myElemStorage);
    }

    /*--- The offsets start with a 0, written by the master node. ---*/

    WriteMPIString("OFFSETS vtktypeint64\n", MASTER_NODE);
    const int64_t firstOffset = 0;
    WriteMPIBinaryData(&firstOffset, sizeof(int64_t), MASTER_NODE);
    WriteMPIBinaryDataAll(offsetBuf.data(), sizeof(int64_t)*myElem, sizeof(int64_t)*GlobalElem,
                          sizeof(int64_t)*dataSorter->GetnElemCumulative(rank));

    WriteMPIString("\nCONNECTIVITY vtktypeint64\n", MASTER_NODE);
    WriteMPIBinaryDataAll(connBuf.data(), sizeof(int64_t)*myElemStorage, sizeof(int64_t)*GlobalElemStorage,
                          sizeof(int64_t)*dataSorter->GetnElemConnCumulative(rank));
  }

  WriteMPIString("\nCELL_TYPES " + std::to_string(GlobalElem) + "\n", MASTER_NODE);

  /*--- Load/write the cell type for all elements in the file. ---*/

  vector<int> typeBuf(myElem);
  auto typeIter = typeBuf.begin();

  std::fill(typeIter, typeIter+nParallel_Line, LINE);          typeIter += nParallel_Line;
  std::fill(typeIter, typeIter+nParallel_Tria, TRIANGLE);      typeIter += nParallel_Tria;
  std::fill(typeIter, typeIter+nParallel_Quad, QUADRILATERAL); typeIter += nParallel_Quad;
  std::fill(typeIter, typeIter+nParallel_Tetr, TETRAHEDRON);   typeIter += nParallel_Tetr;
  std::fill(typeIter, typeIter+nParallel_Hexa, HEXAHEDRON);    typeIter += nParallel_Hexa;
  std::fill(typeIter, typeIter+nParallel_Pris, PRISM);         typeIter += nParallel_Pris;
  std::fill(typeIter, typeIter+nParallel_Pyra, PYRAMID);       typeIter += nParallel_Pyra;

  if (!bigEndian) SwapBytes((char *)typeBuf.data(), sizeof(int), myElem);

  /*--- Compute various data sizes --- */

  sizeInBytesPerPoint = sizeof(int);
  sizeInBytesLocal    = sizeInBytesPerPoint*myElem;
  sizeInBytesGlobal   = sizeInBytesPerPoint*GlobalElem;
  offsetInBytes       = sizeInBytesPerPoint*dataSorter->GetnElemCumulative(rank);

  WriteMPIBinaryDataAll(typeBuf.data(), sizeInBytesLocal, sizeInBytesGlobal, offsetInBytes);

  WriteMPIString("\nPOINT_DATA " + std::to_string(GlobalPoint) + "\n", MASTER_NODE);

  /*--- Adjust container start location to avoid point coords. ---*/

  unsigned short varStart = 2;
  if (nDim == 3) varStart++;

  /*--- Loop over all variables that have been registered in the output. ---*/

  unsigned short iField, VarCounter = varStart;
  for (iField = varStart; iField < fieldNames.size(); iField++) {

    string fieldname = fieldNames[iField];
    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'),
                    fieldname.end());

    /*--- Check whether this field is a vector or scalar. ---*/

    bool output_variable = true, isVector = false;
    size_t found = fieldNames[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector        = true;
    }
    found = fieldNames[iField].find("_y");
    if (found!=string::npos) {
      /*--- We have found a vector, so skip the Y component. ---*/
      output_variable = false;
      VarCounter++;
    }
    found = fieldNames[iField].find("_z");
    if (found!=string::npos) {
      /*--- We have found a vector, so skip the Z component. ---*/
      output_variable = false;
      VarCounter++;
    }

    /*--- Write the point data as an <X,Y,Z> vector or a scalar. ---*/

    if (output_variable && isVector) {

      /*--- Adjust the string name to remove the leading "X-" ---*/

      fieldname.erase(fieldname.end()-2,fieldname.end());

      SPRINTF (str_buf, "\nVECTORS %s float\n", fieldname.c_str());
      WriteMPIString(str_buf, MASTER_NODE);

      /*--- Load up the buffer for writing this rank's vector data. ---*/

      float val = 0.0;
      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        for (iDim = 0; iDim < NCOORDS; iDim++) {
          if (nDim == 2 && iDim == 2) {
            dataBufferFloat[iPoint*NCOORDS + iDim] = 0.0;
          } else {
            val = (float)dataSorter->GetData(VarCounter+iDim,iPoint);
            dataBufferFloat[iPoint*NCOORDS + iDim] = val;
          }
        }
      }
      if (!bigEndian)
        SwapBytes((char *)dataBufferFloat.data(), sizeof(float), myPoint*NCOORDS);

      /*--- Compute various data sizes --- */

      sizeInBytesPerPoint = sizeof(float)*NCOORDS;
      sizeInBytesLocal    = sizeInBytesPerPoint*myPoint;
      sizeInBytesGlobal   = sizeInBytesPerPoint*GlobalPoint;
      offsetInBytes       = sizeInBytesPerPoint*dataSorter->GetnPointCumulative(rank);

      WriteMPIBinaryDataAll(dataBufferFloat.data(), sizeInBytesLocal, sizeInBytesGlobal, offsetInBytes);

      VarCounter++;

    } else if (output_variable) {

      SPRINTF (str_buf, "\nSCALARS %s float 1\n", fieldname.c_str());
      WriteMPIString(str_buf, MASTER_NODE);
      WriteMPIString("LOOKUP_TABLE default\n", MASTER_NODE);

      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/

      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        auto val = (float)dataSorter->GetData(VarCounter,iPoint);
        dataBufferFloat[iPoint] = val;
      }

      if (!bigEndian)
        SwapBytes((char *)dataBufferFloat.data(), sizeof(float), myPoint);

      /*--- Compute various data sizes --- */

      sizeInBytesPerPoint = sizeof(float);
      sizeInBytesLocal    = sizeInBytesPerPoint*myPoint;
      sizeInBytesGlobal   = sizeInBytesPerPoint*GlobalPoint;
      offsetInBytes       = sizeInBytesPerPoint*dataSorter->GetnPointCumulative(rank);

      WriteMPIBinaryDataAll(dataBufferFloat.data(), sizeInBytesLocal, sizeInBytesGlobal, offsetInBytes);

      VarCounter++;
    }

  }

  CloseMPIFile();

}
