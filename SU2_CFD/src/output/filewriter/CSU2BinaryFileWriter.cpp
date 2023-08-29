/*!
 * \file CSU2BinaryFileWriter.cpp
 * \brief Filewriter class SU2 native binary format.
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

#include "../../../include/output/filewriter/CSU2BinaryFileWriter.hpp"

const string CSU2BinaryFileWriter::fileExt = ".dat";

CSU2BinaryFileWriter::CSU2BinaryFileWriter(CParallelDataSorter *valDataSorter)  :
  CFileWriter(valDataSorter, fileExt){}


CSU2BinaryFileWriter::~CSU2BinaryFileWriter()= default;

void CSU2BinaryFileWriter::WriteData(string val_filename){

  /*--- Local variables ---*/

  unsigned short iVar;

  const vector<string>& fieldNames = dataSorter->GetFieldNames();
  unsigned short nVar = fieldNames.size();
  unsigned long nParallel_Poin = dataSorter->GetnPoints();
  unsigned long nPoint_Global = dataSorter->GetnPointsGlobal();

  char str_buf[CGNS_STRING_SIZE];

  /*--- Prepare the first ints containing the counts. The first is a
   magic number that we can use to check for binary files (it is the hex
   representation for "SU2"). The second two values are number of variables
   and number of points (DoFs). ---*/

  int var_buf_size = 5;
  int var_buf[5] = {535532, nVar, (int)nPoint_Global, 0, 0};

  /*--- Open the file using MPI I/O ---*/

  OpenMPIFile(val_filename);

  /*--- First, write the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, write the
   variable string names here. Only the master rank writes the header. ---*/

  WriteMPIBinaryData(var_buf, var_buf_size*sizeof(int), MASTER_NODE);

  /*--- Write the variable names to the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    strncpy(str_buf, fieldNames[iVar].c_str(), CGNS_STRING_SIZE);
    WriteMPIBinaryData(str_buf, CGNS_STRING_SIZE*sizeof(char), MASTER_NODE);
  }

  /*--- Compute various data sizes --- */

  unsigned long sizeInBytesPerPoint = sizeof(passivedouble)*nVar;
  unsigned long sizeInBytesLocal    = sizeInBytesPerPoint*nParallel_Poin;
  unsigned long sizeInBytesGlobal   = sizeInBytesPerPoint*nPoint_Global;
  unsigned long offsetInBytes       = sizeInBytesPerPoint*dataSorter->GetnPointCumulative(rank);

  /*--- Collectively write the actual data to file ---*/

  WriteMPIBinaryDataAll(dataSorter->GetData(), sizeInBytesLocal, sizeInBytesGlobal, offsetInBytes);

  /*--- Close the file ---*/

  CloseMPIFile();

}
