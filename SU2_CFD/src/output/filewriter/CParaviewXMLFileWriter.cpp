/*!
 * \file CParaviewXMLFileWriter.cpp
 * \brief Filewriter class for Paraview binary format.
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/output/filewriter/CParaviewXMLFileWriter.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

const string CParaviewXMLFileWriter::fileExt = ".vtu";

CParaviewXMLFileWriter::CParaviewXMLFileWriter(CParallelDataSorter *valDataSorter) :
  CFileWriter(valDataSorter, fileExt){

  /* Check for big endian. We have to swap bytes otherwise.
   * Since size of character is 1 byte when the character pointer
   *  is de-referenced it will contain only first byte of integer. ---*/

  bigEndian = false;
  unsigned int i = 1;
  char *c = (char*)&i;
  bigEndian = *c == 0;

}

CParaviewXMLFileWriter::~CParaviewXMLFileWriter()= default;

void CParaviewXMLFileWriter::WriteData(string val_filename){

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  /*--- We always have 3 coords, independent of the actual value of nDim ---*/

  const int NCOORDS = 3;
  const unsigned short nDim = dataSorter->GetnDim();
  unsigned short iDim = 0;

  /*--- Array containing the field names we want to output ---*/

  const vector<string>& fieldNames = dataSorter->GetFieldNames();

  unsigned long iPoint, iElem;

  char str_buf[255];

  OpenMPIFile(val_filename);

  dataOffset = 0;

  /*--- Communicate the number of total points that will be
   written by each rank. After this communication, each proc knows how
   many poinnts will be written before its location in the file and the
   offsets can be correctly set. ---*/

  unsigned long myPoint, GlobalPoint;

  GlobalPoint = dataSorter->GetnPointsGlobal();
  myPoint     = dataSorter->GetnPoints();

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

  /* Write the ASCII XML header. Note that we use the appended format for the data,
  * which means that all data is appended at the end of the file in one binary blob.
  */

  if (!bigEndian){
    WriteMPIString("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", MASTER_NODE);
  } else {
    WriteMPIString("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"BigEndian\" header_type=\"UInt64\">\n", MASTER_NODE);
  }

  WriteMPIString("<UnstructuredGrid>\n", MASTER_NODE);

  SPRINTF(str_buf, "<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n",
          SU2_TYPE::Int(GlobalPoint), SU2_TYPE::Int(GlobalElem));

  WriteMPIString(std::string(str_buf), MASTER_NODE);
  WriteMPIString("<Points>\n", MASTER_NODE);
  AddDataArray(VTKDatatype::FLOAT32, "", NCOORDS, myPoint*NCOORDS, GlobalPoint*NCOORDS);
  WriteMPIString("</Points>\n", MASTER_NODE);
  WriteMPIString("<Cells>\n", MASTER_NODE);
  AddDataArray(VTKDatatype::INT32, "connectivity", 1, myElemStorage, GlobalElemStorage);
  AddDataArray(VTKDatatype::INT32, "offsets", 1, myElem, GlobalElem);
  AddDataArray(VTKDatatype::UINT8, "types", 1, myElem, GlobalElem);
  WriteMPIString("</Cells>\n", MASTER_NODE);

  WriteMPIString("<PointData>\n", MASTER_NODE);

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

      AddDataArray(VTKDatatype::FLOAT32, fieldname, NCOORDS, myPoint*NCOORDS, GlobalPoint*NCOORDS);

    } else if (output_variable) {

      AddDataArray(VTKDatatype::FLOAT32, fieldname, 1, myPoint, GlobalPoint);

    }

  }
  WriteMPIString("</PointData>\n", MASTER_NODE);
  WriteMPIString("</Piece>\n", MASTER_NODE);
  WriteMPIString("</UnstructuredGrid>\n", MASTER_NODE);

  /*--- Now write all the data we have previously defined into the binary section of the file ---*/

  WriteMPIString("<AppendedData encoding=\"raw\">\n_", MASTER_NODE);

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

  WriteDataArray(dataBufferFloat.data(), VTKDatatype::FLOAT32, NCOORDS*myPoint, GlobalPoint*NCOORDS,
                 dataSorter->GetnPointCumulative(rank)*NCOORDS);

  /*--- Load/write 1D buffers for the connectivity of each element type. ---*/

  vector<int> connBuf(myElemStorage);
  vector<int> offsetBuf(myElem);
  unsigned long iStorage = 0, iElemID = 0;
  unsigned short iNode = 0;

  auto copyToBuffer = [&](GEO_TYPE type, unsigned long nElem, unsigned short nPoints){
    for (iElem = 0; iElem < nElem; iElem++) {
      for (iNode = 0; iNode < nPoints; iNode++){
        connBuf[iStorage+iNode] = int(dataSorter->GetElemConnectivity(type, iElem, iNode)-1);
      }
      iStorage += nPoints;
      offsetBuf[iElemID++] = int(iStorage + dataSorter->GetnElemConnCumulative(rank));
    }
  };

  copyToBuffer(LINE,          nParallel_Line, N_POINTS_LINE);
  copyToBuffer(TRIANGLE,      nParallel_Tria, N_POINTS_TRIANGLE);
  copyToBuffer(QUADRILATERAL, nParallel_Quad, N_POINTS_QUADRILATERAL);
  copyToBuffer(TETRAHEDRON,   nParallel_Tetr, N_POINTS_TETRAHEDRON);
  copyToBuffer(HEXAHEDRON,    nParallel_Hexa, N_POINTS_HEXAHEDRON);
  copyToBuffer(PRISM,         nParallel_Pris, N_POINTS_PRISM);
  copyToBuffer(PYRAMID,       nParallel_Pyra, N_POINTS_PYRAMID);

  WriteDataArray(connBuf.data(), VTKDatatype::INT32, myElemStorage, GlobalElemStorage,
                 dataSorter->GetnElemConnCumulative(rank));
  WriteDataArray(offsetBuf.data(), VTKDatatype::INT32, myElem, GlobalElem, dataSorter->GetnElemCumulative(rank));

  /*--- Load/write the cell type for all elements in the file. ---*/

  vector<uint8_t> typeBuf(myElem);
  auto typeIter = typeBuf.begin();

  std::fill(typeIter, typeIter+nParallel_Line, LINE);          typeIter += nParallel_Line;
  std::fill(typeIter, typeIter+nParallel_Tria, TRIANGLE);      typeIter += nParallel_Tria;
  std::fill(typeIter, typeIter+nParallel_Quad, QUADRILATERAL); typeIter += nParallel_Quad;
  std::fill(typeIter, typeIter+nParallel_Tetr, TETRAHEDRON);   typeIter += nParallel_Tetr;
  std::fill(typeIter, typeIter+nParallel_Hexa, HEXAHEDRON);    typeIter += nParallel_Hexa;
  std::fill(typeIter, typeIter+nParallel_Pris, PRISM);         typeIter += nParallel_Pris;
  std::fill(typeIter, typeIter+nParallel_Pyra, PYRAMID);       typeIter += nParallel_Pyra;

  WriteDataArray(typeBuf.data(), VTKDatatype::UINT8, myElem, GlobalElem, dataSorter->GetnElemCumulative(rank));

  /*--- Loop over all variables that have been registered in the output. ---*/

  VarCounter = varStart;
  for (iField = varStart; iField < fieldNames.size(); iField++) {

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

      WriteDataArray(dataBufferFloat.data(), VTKDatatype::FLOAT32, myPoint*NCOORDS, GlobalPoint*NCOORDS,
                     dataSorter->GetnPointCumulative(rank)*NCOORDS);

      VarCounter++;

    } else if (output_variable) {


      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/

      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        auto val = (float)dataSorter->GetData(VarCounter,iPoint);
        dataBufferFloat[iPoint] = val;
      }

      WriteDataArray(dataBufferFloat.data(), VTKDatatype::FLOAT32, myPoint, GlobalPoint,
                     dataSorter->GetnPointCumulative(rank));

      VarCounter++;
    }

  }

  WriteMPIString("</AppendedData>\n", MASTER_NODE);
  WriteMPIString("</VTKFile>\n", MASTER_NODE);

  CloseMPIFile();

}

void CParaviewXMLFileWriter::WriteDataArray(void* data, VTKDatatype type, unsigned long arraySize,
                                            unsigned long globalSize, unsigned long offset){

  std::string typeStr;
  unsigned long typeSize = 0;

  GetTypeInfo(type, typeStr, typeSize);

  /*--- Compute the size of the data to write in bytes ---*/

  size_t byteSize = arraySize*typeSize;

  /*--- The total data size ---*/
  size_t totalByteSize = globalSize*typeSize;

  /*--- Only the master node writes the total size in bytes as unsigned long in front of the array data ---*/

  if (!WriteMPIBinaryData(&totalByteSize, sizeof(size_t), MASTER_NODE)){
    SU2_MPI::Error("Writing array size failed", CURRENT_FUNCTION);
  }

  /*--- Collectively write all the data ---*/

  if (!WriteMPIBinaryDataAll(data, byteSize, totalByteSize, offset*typeSize)){
    SU2_MPI::Error("Writing data array failed", CURRENT_FUNCTION);
  }
}

void CParaviewXMLFileWriter::AddDataArray(VTKDatatype type, string name,
                                          unsigned short nComponents, unsigned long arraySize, unsigned long globalSize){

  /*--- Add quotation marks around the arguments ---*/

  name = "\"" + name + "\"";

  string nComp ="\"" + PrintingToolbox::to_string(nComponents) + "\"";

  /*--- Full precision ASCII output of offset for 32 bit integer ---*/

  stringstream ss;
  ss.precision(10); ss.setf(std::ios::fixed, std::ios::floatfield);
  ss <<  "\"" << dataOffset <<  "\"";
  string offsetStr = ss.str();

  std::string typeStr;
  unsigned long typeSize = 0;

  GetTypeInfo(type, typeStr, typeSize);

  /*--- Total data size ---*/

  size_t totalByteSize = globalSize*typeSize;

  /*--- Write the ASCII XML header information for this array ---*/

  WriteMPIString(string("<DataArray type=") + typeStr +
                 string(" Name=") + name +
                 string(" NumberOfComponents= ") + nComp +
                 string(" offset=") + offsetStr +
                 string(" format=\"appended\"/>\n"), MASTER_NODE);

  dataOffset += totalByteSize + sizeof(size_t);

}
