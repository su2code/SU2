/*!
 * \file CParaviewBinaryFileWriter.cpp
 * \brief Filewriter class for Paraview binary format.
 * \author T. Albring
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

const string CParaviewBinaryFileWriter::fileExt = ".vtk";

CParaviewBinaryFileWriter::CParaviewBinaryFileWriter(string valFileName, CParallelDataSorter *valDataSorter) :
  CFileWriter(std::move(valFileName), valDataSorter, fileExt){
  
  /* Check for big endian. We have to swap bytes otherwise.
   * Since size of character is 1 byte when the character pointer
   *  is de-referenced it will contain only first byte of integer. ---*/

  bigEndian = false;
  unsigned int i = 1;
  char *c = (char*)&i;
  if (*c) bigEndian = false;
  else bigEndian = true;
}


CParaviewBinaryFileWriter::~CParaviewBinaryFileWriter(){

}

void CParaviewBinaryFileWriter::Write_Data(){

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  const vector<string>& fieldNames = dataSorter->GetFieldNames();

  unsigned short iDim = 0, nDim = dataSorter->GetnDim();

  unsigned long iPoint, iElem;

  const int MAX_STRING_LENGTH = 255;
  char str_buf[MAX_STRING_LENGTH];

  const int NCOORDS = 3;

  OpenMPIFile();
  
  string header = "# vtk DataFile Version 3.0\n"
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


  int *nPoint_Snd = new int[size+1];
  int *nPoint_Cum = new int[size+1];

  nPoint_Snd[0] = 0; nPoint_Cum[0] = 0;
  for (int ii=1; ii < size; ii++) {
    nPoint_Snd[ii] = myPoint; nPoint_Cum[ii] = 0;
  }
  nPoint_Snd[size] = myPoint; nPoint_Cum[size] = 0;

  /*--- Communicate the local counts to all ranks for building offsets. ---*/

  SU2_MPI::Alltoall(&(nPoint_Snd[1]), 1, MPI_INT,
                    &(nPoint_Cum[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Put the counters into cumulative storage format. ---*/

  for (int ii = 0; ii < size; ii++) {
    nPoint_Cum[ii+1] += nPoint_Cum[ii];
  }

  SPRINTF(str_buf, "POINTS %i float\n", SU2_TYPE::Int(GlobalPoint));
  
  WriteMPIString(string(str_buf), MASTER_NODE);

  /*--- Load/write the 1D buffer of point coordinates. Note that we
   always have 3 coordinate dimensions, even for 2D problems. ---*/
  
  vector<float> dataBufferFloat(myPoint*NCOORDS);
  for (iPoint = 0; iPoint < myPoint; iPoint++) {
    for (iDim = 0; iDim < NCOORDS; iDim++) {
      if (nDim == 2 && iDim == 2) {
        dataBufferFloat[iPoint*NCOORDS + iDim] = 0.0;
      } else {
        float val = (float)dataSorter->GetData(iDim, iPoint);
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
  
  unsigned long nTot_Line = dataSorter->GetnElemGlobal(LINE),
                nTot_Tria = dataSorter->GetnElemGlobal(TRIANGLE),
                nTot_Quad = dataSorter->GetnElemGlobal(QUADRILATERAL),
                nTot_Tetr = dataSorter->GetnElemGlobal(TETRAHEDRON),
                nTot_Hexa = dataSorter->GetnElemGlobal(HEXAHEDRON),
                nTot_Pris = dataSorter->GetnElemGlobal(PRISM),
                nTot_Pyra = dataSorter->GetnElemGlobal(PYRAMID);

  myElem        = dataSorter->GetnElem();
  
  myElemStorage = (nParallel_Line*N_POINTS_LINE + 
                   nParallel_Tria*N_POINTS_TRIANGLE + 
                   nParallel_Quad*N_POINTS_QUADRILATERAL +
                   nParallel_Tetr*N_POINTS_TETRAHEDRON +
                   nParallel_Hexa*N_POINTS_HEXAHEDRON +
                   nParallel_Pris*N_POINTS_PRISM +
                   nParallel_Pyra*N_POINTS_PYRAMID);

  GlobalElem        =  dataSorter->GetnElemGlobal();
  GlobalElemStorage = (nTot_Line*N_POINTS_LINE + 
                       nTot_Tria*N_POINTS_TRIANGLE + 
                       nTot_Quad*N_POINTS_QUADRILATERAL +
                       nTot_Tetr*N_POINTS_TETRAHEDRON +
                       nTot_Hexa*N_POINTS_HEXAHEDRON +
                       nTot_Pris*N_POINTS_PRISM +
                       nTot_Pyra*N_POINTS_PYRAMID);

  /*--- Communicate the number of total cells/storage that will be
   written by each rank. After this communication, each proc knows how
   many cells will be written before its location in the file and the
   offsets can be correctly set. ---*/

  int *nElem_Snd = new int[size+1]; int *nElemStorage_Snd = new int[size+1];
  int *nElem_Cum = new int[size+1]; int *nElemStorage_Cum = new int[size+1];

  nElem_Snd[0] = 0; nElemStorage_Snd[0] = 0;
  nElem_Cum[0] = 0; nElemStorage_Cum[0] = 0;
  for (int ii=1; ii < size; ii++) {
    nElem_Snd[ii] = myElem; nElemStorage_Snd[ii] = myElemStorage;
    nElem_Cum[ii] = 0;      nElemStorage_Cum[ii] = 0;
  }
  nElem_Snd[size] = myElem; nElemStorage_Snd[size] = myElemStorage;
  nElem_Cum[size] = 0;      nElemStorage_Cum[size] = 0;

  /*--- Communicate the local counts to all ranks for building offsets. ---*/

  SU2_MPI::Alltoall(&(nElem_Snd[1]), 1, MPI_INT,
                    &(nElem_Cum[1]), 1, MPI_INT, MPI_COMM_WORLD);

  SU2_MPI::Alltoall(&(nElemStorage_Snd[1]), 1, MPI_INT,
                    &(nElemStorage_Cum[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Put the counters into cumulative storage format. ---*/

  for (int ii = 0; ii < size; ii++) {
    nElem_Cum[ii+1]        += nElem_Cum[ii];
    nElemStorage_Cum[ii+1] += nElemStorage_Cum[ii];
  }
  
  SPRINTF(str_buf, "\nCELLS %i %i\n", SU2_TYPE::Int(GlobalElem),
          SU2_TYPE::Int(GlobalElemStorage));
  WriteMPIString(str_buf, MASTER_NODE);
  
  /*--- Load/write 1D buffers for the connectivity of each element type. ---*/

  vector<int> connBuf(myElemStorage);
  unsigned long iStorage = 0, iElemID = 0;
  unsigned short iNode = 0;

  auto copyToBuffer = [&](GEO_TYPE type, unsigned long nElem, unsigned short nPoints){
    for (iElem = 0; iElem < nElem; iElem++) {
      for (iNode = 0; iNode < nPoints; iNode++){
        connBuf[iStorage+iNode] = int(dataSorter->GetElem_Connectivity(type, iElem, iNode)-1);
      }
      iStorage += nPoints;
    }
  };
 
  copyToBuffer(LINE,          nParallel_Line, N_POINTS_LINE);
  copyToBuffer(TRIANGLE,      nParallel_Tria, N_POINTS_TRIANGLE);
  copyToBuffer(QUADRILATERAL, nParallel_Quad, N_POINTS_QUADRILATERAL);
  copyToBuffer(TETRAHEDRON,   nParallel_Tetr, N_POINTS_TETRAHEDRON);
  copyToBuffer(HEXAHEDRON,    nParallel_Hexa, N_POINTS_HEXAHEDRON);
  copyToBuffer(PRISM,         nParallel_Pris, N_POINTS_PRISM);
  copyToBuffer(PYRAMID,       nParallel_Pyra, N_POINTS_PYRAMID);

  if (!bigEndian) SwapBytes((char *)connBuf.data(), sizeof(int), myElemStorage);  
  
  WriteMPIBinaryDataAll(connBuf.data(), 
                        myElemStorage*sizeof(int),
                        GlobalElemStorage*sizeof(int), 
                        nElemStorage_Cum[rank]*sizeof(int));
    
  SPRINTF (str_buf, "\nCELL_TYPES %i\n", SU2_TYPE::Int(GlobalElem));
  WriteMPIString(str_buf, MASTER_NODE);

  /*--- Load/write the cell type for all elements in the file. ---*/

  vector<int> typeBuf(myElem);
  vector<int>::iterator typeIter = typeBuf.begin();

  std::fill(typeIter, typeIter+nParallel_Line, LINE);          typeIter += nParallel_Line;
  std::fill(typeIter, typeIter+nParallel_Tria, TRIANGLE);      typeIter += nParallel_Tria;
  std::fill(typeIter, typeIter+nParallel_Quad, QUADRILATERAL); typeIter += nParallel_Quad;
  std::fill(typeIter, typeIter+nParallel_Tetr, TETRAHEDRON);   typeIter += nParallel_Tetr;
  std::fill(typeIter, typeIter+nParallel_Hexa, HEXAHEDRON);    typeIter += nParallel_Hexa;
  std::fill(typeIter, typeIter+nParallel_Pris, PRISM);         typeIter += nParallel_Pris;
  std::fill(typeIter, typeIter+nParallel_Pyra, PYRAMID);       typeIter += nParallel_Pyra;
  
  if (!bigEndian) SwapBytes((char *)typeBuf.data(), sizeof(int), myElem);

  WriteMPIBinaryDataAll(typeBuf.data(),
                        myElem*sizeof(int),
                        GlobalElem*sizeof(int),
                        nElem_Cum[rank]*sizeof(int));

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

      WriteMPIBinaryDataAll(dataBufferFloat.data(),
                            myPoint*NCOORDS*sizeof(float),
                            GlobalPoint*NCOORDS*sizeof(float), 
                            nPoint_Cum[rank]*NCOORDS*sizeof(float));
      VarCounter++;

    } else if (output_variable) {
      /*--- For now, create a temp 1D buffer to load up the data for writing.
       This will be replaced with a derived data type most likely. ---*/

      for (iPoint = 0; iPoint < myPoint; iPoint++) {
        float val = (float)dataSorter->GetData(VarCounter,iPoint);
        dataBufferFloat[iPoint] = val;
      }
      
      if (!bigEndian)
        SwapBytes((char *)dataBufferFloat.data(), sizeof(float), myPoint);
      
      WriteMPIBinaryDataAll(dataBufferFloat.data(),
                            myPoint*sizeof(float),
                            GlobalPoint*sizeof(float), 
                            nPoint_Cum[rank]*sizeof(float));
      

      VarCounter++;
    }

  }

  CloseMPIFile();

  /*--- Delete the offset counters that we needed for MPI IO. ---*/

  delete [] nElem_Snd;        delete [] nElem_Cum;
  delete [] nElemStorage_Snd; delete [] nElemStorage_Cum;
  delete [] nPoint_Snd;       delete [] nPoint_Cum;

}


/*--- Subroutine to swap bytes, in case we need to convert to
 big endian, which is expected for ParaView binary legacy format. ---*/

void CParaviewBinaryFileWriter::SwapBytes(char *buffer, size_t nBytes, unsigned long nVar) {

  /*--- Store half the number of bytes in kk. ---*/

  const int kk = (int)nBytes/2;

  /*--- Loop over the number of variables in the buffer. ---*/

  for (int j = 0; j < (int)nVar; j++) {

    /*--- Initialize ii and jj, which are used to store the
     indices of the bytes to be swapped. ---*/

    int ii = j*(int)nBytes;
    int jj = ii + (int)nBytes - 1;

    /*--- Swap the bytes. ---*/

    for (int i = 0; i < kk; i++) {
      char tmp   = buffer[jj];
      buffer[jj] = buffer[ii];
      buffer[ii] = tmp;

      ii++;
      jj--;

    }
  }
}
