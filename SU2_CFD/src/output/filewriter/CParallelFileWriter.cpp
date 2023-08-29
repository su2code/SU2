/*!
 * \file CParallelFileWriter.cpp
 * \brief Filewriter base class.
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

#include <utility>

#include "../../../include/output/filewriter/CFileWriter.hpp"

CFileWriter::CFileWriter(CParallelDataSorter *valDataSorter, string valFileExt):
  fileExt(std::move(valFileExt)),
  dataSorter(valDataSorter){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  fileSize = 0.0;
  bandwidth = 0.0;

}

CFileWriter::CFileWriter(string valFileExt):
  fileExt(std::move(valFileExt)){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  fileSize = 0.0;
  bandwidth = 0.0;

}

CFileWriter::~CFileWriter()= default;

bool CFileWriter::WriteMPIBinaryDataAll(const void *data, unsigned long sizeInBytes,
                                        unsigned long totalSizeInBytes, unsigned long offsetInBytes){

#ifdef HAVE_MPI

  startTime = SU2_MPI::Wtime();

  MPI_Datatype filetype;

  /*--- Prepare to write the actual data ---*/

  MPI_Type_contiguous(int(sizeInBytes), MPI_BYTE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- Set the view for the MPI file write, i.e., describe the
 location in the file that this rank "sees" for writing its
 piece of the file. ---*/

  MPI_File_set_view(fhw, disp + offsetInBytes, MPI_BYTE, filetype,
                    (char*)"native", MPI_INFO_NULL);

  /*--- Collective call for all ranks to write simultaneously. ---*/

  int ierr = MPI_File_write_all(fhw, data, int(sizeInBytes), MPI_BYTE, MPI_STATUS_IGNORE);

  MPI_Type_free(&filetype);

  disp      += totalSizeInBytes;
  fileSize  += sizeInBytes;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (ierr == MPI_SUCCESS);
#else

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten;

  /*--- Write binary data ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);
  fileSize += bytesWritten;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == sizeInBytes);
#endif

}

bool CFileWriter::WriteMPIBinaryData(const void *data, unsigned long sizeInBytes, unsigned short processor){

#ifdef HAVE_MPI

  startTime = SU2_MPI::Wtime();

  int ierr = MPI_SUCCESS;

  /*--- Reset the file view. ---*/

  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                    (char*)"native", MPI_INFO_NULL);

  if (rank == processor)
    ierr = MPI_File_write_at(fhw, disp, data, int(sizeInBytes), MPI_BYTE, MPI_STATUS_IGNORE);

  disp     += sizeInBytes;
  fileSize += sizeInBytes;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (ierr == MPI_SUCCESS);
#else

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten = sizeInBytes;

  /*--- Write the total size in bytes at the beginning of the binary data blob ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == sizeInBytes);

#endif

}

bool CFileWriter::WriteMPIString(const string &str, unsigned short processor){

#ifdef HAVE_MPI

  startTime = SU2_MPI::Wtime();

  int ierr = MPI_SUCCESS;

  /*--- Reset the file view. ---*/

  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE,
                    (char*)"native", MPI_INFO_NULL);

  if (SU2_MPI::GetRank() == processor)
    ierr = MPI_File_write_at(fhw, disp, str.c_str(), str.size(),
                      MPI_CHAR, MPI_STATUS_IGNORE);

  disp += str.size()*sizeof(char);
  fileSize += sizeof(char)*str.size();

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (ierr == MPI_SUCCESS);

#else

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten;
  bytesWritten = fwrite(str.c_str(), sizeof(char), str.size(), fhw);

  fileSize += bytesWritten;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == str.size()*sizeof(char));

#endif

}

bool CFileWriter::OpenMPIFile(string val_filename){

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

#ifdef HAVE_MPI
  int ierr;
  disp     = 0.0;

  /*--- All ranks open the file using MPI. Here, we try to open the file with
   exclusive so that an error is generated if the file exists. We always want
   to write a fresh output file, so we delete any existing files and create
   a new one. ---*/

  ierr = MPI_File_open(SU2_MPI::GetComm(), val_filename.c_str(),
                       MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);
  if (ierr != MPI_SUCCESS)  {
    MPI_File_close(&fhw);
    if (rank == 0)
      MPI_File_delete(val_filename.c_str(), MPI_INFO_NULL);
    ierr = MPI_File_open(SU2_MPI::GetComm(), val_filename.c_str(),
                         MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &fhw);
  }

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open file ") +
                   val_filename, CURRENT_FUNCTION);
  }
#else
  fhw = fopen(val_filename.c_str(), "wb");
  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open file ") +
                   val_filename, CURRENT_FUNCTION);
  }
#endif

  fileSize = 0.0;
  usedTime = 0;

  return true;
}

bool CFileWriter::CloseMPIFile(){

#ifdef HAVE_MPI
  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);
#else
  fclose(fhw);
#endif

  /*--- Communicate the total file size for the restart ---*/

  su2double my_fileSize = fileSize;
  SU2_MPI::Allreduce(&my_fileSize, &fileSize, 1,
                     MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  /*--- Compute and store the bandwidth ---*/

  bandwidth = fileSize/(1.0e6)/usedTime;

  return true;
}

