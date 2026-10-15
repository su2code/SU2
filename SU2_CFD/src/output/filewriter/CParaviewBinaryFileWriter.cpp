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

  WritePoints();

  if (cellsInt64) {
    WriteCellsInt64();
  } else {
    WriteCellsInt32(GlobalCellStorage);
  }

  WriteCellTypes();

  WritePointData();

  CloseMPIFile();
}

void CParaviewBinaryFileWriter::WritePointArray(vector<float>& buffer, unsigned short nComponents) {

  const unsigned long myPoint = dataSorter->GetnPoints();

  if (!bigEndian) SwapBytes((char *)buffer.data(), sizeof(float), myPoint*nComponents);

  const unsigned long sizeInBytesPerPoint = sizeof(float)*nComponents;
  WriteMPIBinaryDataAll(buffer.data(), sizeInBytesPerPoint*myPoint, sizeInBytesPerPoint*dataSorter->GetnPointsGlobal(),
                        sizeInBytesPerPoint*dataSorter->GetnPointCumulative(rank));
}

void CParaviewBinaryFileWriter::WriteVectorArray(unsigned short firstVar) {

  /*--- There are always 3 components, the third one is 0 in 2D. ---*/

  const unsigned short nDim = dataSorter->GetnDim();
  const unsigned long myPoint = dataSorter->GetnPoints();

  vector<float> buffer(myPoint*NCOORDS);
  for (auto iPoint = 0ul; iPoint < myPoint; iPoint++) {
    for (unsigned short iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        buffer[iPoint*NCOORDS + iDim] = 0.0;
      } else {
        buffer[iPoint*NCOORDS + iDim] = (float)dataSorter->GetData(firstVar+iDim, iPoint);
      }
    }
  }
  WritePointArray(buffer, NCOORDS);
}

void CParaviewBinaryFileWriter::WritePoints() {

  WriteMPIString("POINTS " + std::to_string(dataSorter->GetnPointsGlobal()) + " float\n", MASTER_NODE);

  /*--- The coordinates are the first fields of the data sorter. ---*/

  WriteVectorArray(0);
}

void CParaviewBinaryFileWriter::WriteCellsInt32(unsigned long GlobalCellStorage) {

  const unsigned long myElem = dataSorter->GetnElem();
  const unsigned long myElemStorage = dataSorter->GetnConn();

  WriteMPIString("\nCELLS " + std::to_string(dataSorter->GetnElemGlobal()) + " " + std::to_string(GlobalCellStorage) +
                 "\n", MASTER_NODE);

  /*--- Load/write the 1D buffer of the number of nodes followed by the node ids of each cell. ---*/

  vector<int32_t> connBuf(myElemStorage + myElem);
  unsigned long iStorage = 0;

  for (auto type : elemTypes) {
    const auto nPoints = nPointsOfElementType(type);
    for (auto iElem = 0ul; iElem < dataSorter->GetnElem(type); iElem++) {
      connBuf[iStorage++] = nPoints;
      for (unsigned short iNode = 0; iNode < nPoints; iNode++)
        connBuf[iStorage++] = static_cast<int32_t>(dataSorter->GetElemConnectivity(type, iElem, iNode) - 1);
    }
  }

  if (!bigEndian) SwapBytes((char *)connBuf.data(), sizeof(int32_t), myElemStorage+myElem);

  WriteMPIBinaryDataAll(connBuf.data(), sizeof(int32_t)*(myElemStorage + myElem), sizeof(int32_t)*GlobalCellStorage,
                        sizeof(int32_t)*(dataSorter->GetnElemConnCumulative(rank) +
                                         dataSorter->GetnElemCumulative(rank)));
}

void CParaviewBinaryFileWriter::WriteCellsInt64() {

  const unsigned long myElem = dataSorter->GetnElem();
  const unsigned long myElemStorage = dataSorter->GetnConn();
  const unsigned long GlobalElem = dataSorter->GetnElemGlobal();
  const unsigned long GlobalElemStorage = dataSorter->GetnConnGlobal();

  WriteMPIString("\nCELLS " + std::to_string(GlobalElem + 1) + " " + std::to_string(GlobalElemStorage) + "\n",
                 MASTER_NODE);

  /*--- Load the offsets (where each cell ends in the connectivity) and the connectivity. ---*/

  vector<int64_t> offsetBuf(myElem), connBuf(myElemStorage);
  unsigned long iStorage = 0, iCell = 0;

  for (auto type : elemTypes) {
    const auto nPoints = nPointsOfElementType(type);
    for (auto iElem = 0ul; iElem < dataSorter->GetnElem(type); iElem++) {
      for (unsigned short iNode = 0; iNode < nPoints; iNode++)
        connBuf[iStorage++] = static_cast<int64_t>(dataSorter->GetElemConnectivity(type, iElem, iNode)) - 1;
      offsetBuf[iCell++] = static_cast<int64_t>(iStorage + dataSorter->GetnElemConnCumulative(rank));
    }
  }

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

void CParaviewBinaryFileWriter::WriteCellTypes() {

  const unsigned long myElem = dataSorter->GetnElem();
  const unsigned long GlobalElem = dataSorter->GetnElemGlobal();

  WriteMPIString("\nCELL_TYPES " + std::to_string(GlobalElem) + "\n", MASTER_NODE);

  /*--- Load/write the cell type for all elements in the file, in the same order as the cells. ---*/

  vector<int> typeBuf(myElem);
  auto typeIter = typeBuf.begin();
  for (auto type : elemTypes) {
    const auto nElem = dataSorter->GetnElem(type);
    std::fill(typeIter, typeIter+nElem, type);
    typeIter += nElem;
  }

  if (!bigEndian) SwapBytes((char *)typeBuf.data(), sizeof(int), myElem);

  WriteMPIBinaryDataAll(typeBuf.data(), sizeof(int)*myElem, sizeof(int)*GlobalElem,
                        sizeof(int)*dataSorter->GetnElemCumulative(rank));
}

void CParaviewBinaryFileWriter::WritePointData() {

  const vector<string>& fieldNames = dataSorter->GetFieldNames();
  const unsigned long myPoint = dataSorter->GetnPoints();

  WriteMPIString("\nPOINT_DATA " + std::to_string(dataSorter->GetnPointsGlobal()) + "\n", MASTER_NODE);

  /*--- Skip the coordinates, which are the first fields. ---*/

  const unsigned short varStart = (dataSorter->GetnDim() == 3) ? 3 : 2;

  /*--- Loop over all variables that have been registered in the output. A field ending in "_x" starts a
   vector, written with its "_y" (and "_z") components, which are then skipped. ---*/

  unsigned short VarCounter = varStart;
  for (unsigned long iField = varStart; iField < fieldNames.size(); iField++) {

    string fieldname = fieldNames[iField];
    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'), fieldname.end());

    const bool isVector = fieldNames[iField].find("_x") != string::npos;
    const bool isY = fieldNames[iField].find("_y") != string::npos;
    const bool isZ = fieldNames[iField].find("_z") != string::npos;

    if (isY || isZ) {
      VarCounter += isY + isZ;
    } else if (isVector) {

      /*--- Remove the "_x" from the name. ---*/

      fieldname.erase(fieldname.end()-2, fieldname.end());
      WriteMPIString("\nVECTORS " + fieldname + " float\n", MASTER_NODE);
      WriteVectorArray(VarCounter);
      VarCounter++;

    } else {

      WriteMPIString("\nSCALARS " + fieldname + " float 1\n", MASTER_NODE);
      WriteMPIString("LOOKUP_TABLE default\n", MASTER_NODE);

      vector<float> buffer(myPoint);
      for (auto iPoint = 0ul; iPoint < myPoint; iPoint++) buffer[iPoint] = (float)dataSorter->GetData(VarCounter, iPoint);
      WritePointArray(buffer, 1);
      VarCounter++;
    }
  }
}
