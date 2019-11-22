/*!
 * \file CSU2BinaryFileWriter.cpp
 * \brief Filewriter class SU2 native binary format.
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

#include "../../../include/output/filewriter/CSU2BinaryFileWriter.hpp"

const string CSU2BinaryFileWriter::fileExt = ".dat";

CSU2BinaryFileWriter::CSU2BinaryFileWriter(vector<string> fields, unsigned short nDim,
                                           string fileName, CParallelDataSorter *dataSorter)  :
  CFileWriter(std::move(fields), std::move(fileName), dataSorter, fileExt, nDim){}


CSU2BinaryFileWriter::~CSU2BinaryFileWriter(){

}

void CSU2BinaryFileWriter::Write_Data(){

  /*--- Local variables ---*/

  unsigned short iVar;

  unsigned short GlobalField_Counter = fieldnames.size();
  unsigned long nParallel_Poin = dataSorter->GetnPoints();

  ofstream restart_file;
  char str_buf[CGNS_STRING_SIZE], fname[100];

  file_size = 0.0;

  strcpy(fname, fileName.c_str());

  /*--- Prepare the first ints containing the counts. The first is a
   magic number that we can use to check for binary files (it is the hex
   representation for "SU2"). The second two values are number of variables
   and number of points (DoFs). The last two values are for metadata:
   one int for ExtIter and 8 su2doubles. ---*/

  int var_buf_size = 5;
  int var_buf[5] = {535532, GlobalField_Counter, (int)dataSorter->GetnPointsGlobal(), 0, 0};

  /*--- Set a timer for the binary file writing. ---*/

#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif

#ifndef HAVE_MPI

  FILE* fhw;
  fhw = fopen(fname, "wb");

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points. ---*/

  fwrite(var_buf, var_buf_size, sizeof(int), fhw);
  file_size += (su2double)var_buf_size*sizeof(int);

  /*--- Write the variable names to the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. ---*/

  for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
    strncpy(str_buf, fieldnames[iVar].c_str(), CGNS_STRING_SIZE);
    fwrite(str_buf, CGNS_STRING_SIZE, sizeof(char), fhw);
    file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
  }

  /*--- Call to write the entire restart file data in binary in one shot. ---*/

  fwrite(dataSorter->GetData(), nParallel_Poin*GlobalField_Counter, sizeof(passivedouble), fhw);
  file_size += (su2double)nParallel_Poin*GlobalField_Counter*sizeof(passivedouble);

  /*--- Close the file. ---*/

  fclose(fhw);

#else

  /*--- Parallel binary output using MPI I/O. ---*/

  MPI_File fhw;
  SU2_MPI::Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp;
  int ierr;

  /*--- We're writing only su2doubles in the data portion of the file. ---*/

  etype = MPI_DOUBLE;

  /*--- Define a derived datatype for this ranks contiguous chunk of data
   that will be placed in the restart (1D array size = num points * num vars). ---*/

  MPI_Type_contiguous(nParallel_Poin*GlobalField_Counter, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- All ranks open the file using MPI. Here, we try to open the file with
   exclusive so that an error is generated if the file exists. We always want
   to write a fresh restart file, so we delete any existing files and create
   a new one. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                       MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);
  if (ierr != MPI_SUCCESS)  {
    MPI_File_close(&fhw);
    if (rank == 0)
      MPI_File_delete(fname, MPI_INFO_NULL);
    ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                         MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &fhw);
  }

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, write the
   variable string names here. Only the master rank writes the header. ---*/

  if (rank == MASTER_NODE) {
    MPI_File_write(fhw, var_buf, var_buf_size, MPI_INT, MPI_STATUS_IGNORE);
    file_size += (su2double)var_buf_size*sizeof(int);

    /*--- Write the variable names to the file. Note that we are adopting a
     fixed length of 33 for the string length to match with CGNS. This is
     needed for when we read the strings later. ---*/

    for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
      disp = var_buf_size*sizeof(int) + iVar*CGNS_STRING_SIZE*sizeof(char);
      strncpy(str_buf, fieldnames[iVar].c_str(), CGNS_STRING_SIZE);
      MPI_File_write_at(fhw, disp, str_buf, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);
      file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
    }
  }

  /*--- Compute the offset for this rank's linear partition of the data in bytes.
   After the calculations above, we have the partition sizes store in nPoint_Linear
   in cumulative storage format. ---*/

  disp = (var_buf_size*sizeof(int) + GlobalField_Counter*CGNS_STRING_SIZE*sizeof(char) +
          GlobalField_Counter*dataSorter->GetnPointCumulative(rank)*sizeof(passivedouble));

  /*--- Set the view for the MPI file write, i.e., describe the location in
   the file that this rank "sees" for writing its piece of the restart file. ---*/

  MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

  /*--- Collective call for all ranks to write to their view simultaneously. ---*/

  MPI_File_write_all(fhw, dataSorter->GetData(), GlobalField_Counter*nParallel_Poin, MPI_DOUBLE, &status);
  file_size += (su2double)GlobalField_Counter*nParallel_Poin*sizeof(passivedouble);

  /*--- Free the derived datatype. ---*/

  MPI_Type_free(&filetype);

  /*--- Reset the file view before writing the metadata. ---*/

  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL);

  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);

#endif

  /*--- Compute and store the write time. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;

  /*--- Communicate the total file size for the restart ---*/

#ifdef HAVE_MPI
  su2double my_file_size = file_size;
  SU2_MPI::Allreduce(&my_file_size, &file_size, 1,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Compute and store the bandwidth ---*/

  Bandwidth = file_size/(1.0e6)/UsedTime;

}
