/*!
 * \file CParallelFileWriter.cpp
 * \brief Filewriter base class.
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

#include <chrono>
#include <numeric>
#include <limits>
#include <utility>
#include <vector>

#include "../../../include/output/filewriter/CFileWriter.hpp"

#ifdef HAVE_MPI
namespace {
void CheckMPIWrite(int error, const MPI_Status& status, unsigned long expected) {
  if (error != MPI_SUCCESS) SU2_MPI::Error("MPI output write failed", CURRENT_FUNCTION);
#if MPI_VERSION >= 3
  MPI_Count written = 0;
  if (MPI_Get_elements_x(&status, MPI_BYTE, &written) != MPI_SUCCESS ||
      written < 0 || static_cast<unsigned long>(written) != expected)
    SU2_MPI::Error("MPI output write failed or was incomplete", CURRENT_FUNCTION);
#else
  // MPI-2 cannot report byte counts above INT_MAX; the final extent is checked at close too.
  int written = 0;
  if (expected <= static_cast<unsigned long>(std::numeric_limits<int>::max()) &&
      (MPI_Get_elements(&status, MPI_BYTE, &written) != MPI_SUCCESS || written < 0 ||
       static_cast<unsigned long>(written) != expected))
    SU2_MPI::Error("MPI output write failed or was incomplete", CURRENT_FUNCTION);
#endif
}
}  // namespace
#endif

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

passivedouble CFileWriter::WallTime() {
  // The no-MPI/no-OpenMP fallback timer measures CPU time, which excludes I/O waits.
  using Clock = std::chrono::steady_clock;
  static const auto origin = Clock::now();
  return std::chrono::duration<passivedouble>(Clock::now() - origin).count();
}

void CFileWriter::UpdateBandwidth(passivedouble elapsed) {
  passivedouble maximum = 0;
  SelectMPIWrapper<passivedouble>::W::Allreduce(&elapsed, &maximum, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
  usedTime = maximum;
  bandwidth = maximum > 0 ? SU2_TYPE::GetValue(fileSize) / (1.0e6 * maximum) : 0.0;
}

bool CFileWriter::WriteMPIBinaryDataAll(const void *data, unsigned long sizeInBytes,
                                        unsigned long totalSizeInBytes, unsigned long offsetInBytes){

#ifdef HAVE_MPI

  startTime = WallTime();

  /*--- Prepare to write the actual data. MPI counts are int, so the local data
   (which can exceed INT_MAX bytes) is described as a number of fixed-size blocks
   plus a remainder, combined in a single datatype that is written with count 1. ---*/

  constexpr unsigned long blockBytes = 1ul << 20;
  CParallelDataSorter::CheckIntCount(sizeInBytes / blockBytes);
  int blockLengths[2] = {static_cast<int>(sizeInBytes / blockBytes), static_cast<int>(sizeInBytes % blockBytes)};
  MPI_Aint displacements[2] = {0, static_cast<MPI_Aint>(sizeInBytes - sizeInBytes % blockBytes)};

  MPI_Datatype blocktype, datatype;
  MPI_Type_contiguous(static_cast<int>(blockBytes), MPI_BYTE, &blocktype);
  MPI_Datatype types[2] = {blocktype, MPI_BYTE};
  MPI_Type_create_struct(2, blockLengths, displacements, types, &datatype);
  MPI_Type_commit(&datatype);

  /*--- Set the view for the MPI file write, i.e., describe the
 location in the file that this rank "sees" for writing its
 piece of the file. The data is contiguous so a byte view is enough. ---*/

  if (MPI_File_set_view(fhw, disp + offsetInBytes, MPI_BYTE, MPI_BYTE,
                        (char*)"native", MPI_INFO_NULL) != MPI_SUCCESS)
    SU2_MPI::Error("Setting MPI output view failed", CURRENT_FUNCTION);

  /*--- Collective call for all ranks to write simultaneously. ---*/

  MPI_Status status{};
  int ierr = MPI_File_write_all(fhw, data, 1, datatype, &status);

  MPI_Type_free(&datatype);
  MPI_Type_free(&blocktype);

  disp      += totalSizeInBytes;
  fileSize  += sizeInBytes;

  stopTime = WallTime();

  usedTime += stopTime - startTime;

  CheckMPIWrite(ierr, status, sizeInBytes);
  return true;
#else

  startTime = WallTime();

  unsigned long bytesWritten;

  /*--- Write binary data ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);
  fileSize += bytesWritten;

  stopTime = WallTime();

  usedTime += stopTime - startTime;

  if (bytesWritten != sizeInBytes) SU2_MPI::Error("Output write failed", CURRENT_FUNCTION);
  return true;
#endif

}

bool CFileWriter::WriteMPIBinaryData(const void *data, unsigned long sizeInBytes, unsigned short processor){
#ifdef HAVE_MPI
  CParallelDataSorter::CheckIntCount(sizeInBytes);
#endif

#ifdef HAVE_MPI

  startTime = WallTime();

  int ierr = MPI_SUCCESS;
  MPI_Status status{};

  /*--- Reset the file view. ---*/

  if (MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL) != MPI_SUCCESS)
    SU2_MPI::Error("Setting MPI output view failed", CURRENT_FUNCTION);

  if (rank == processor)
    ierr = MPI_File_write_at(fhw, disp, data, int(sizeInBytes), MPI_BYTE, &status);

  disp     += sizeInBytes;
  if (rank == processor) fileSize += sizeInBytes;

  stopTime = WallTime();

  usedTime += stopTime - startTime;

  if (rank == processor) CheckMPIWrite(ierr, status, sizeInBytes);
  return true;
#else

  startTime = WallTime();

  unsigned long bytesWritten = sizeInBytes;

  /*--- Write the total size in bytes at the beginning of the binary data blob ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);
  fileSize += bytesWritten;

  stopTime = WallTime();

  usedTime += stopTime - startTime;

  if (bytesWritten != sizeInBytes) SU2_MPI::Error("Output write failed", CURRENT_FUNCTION);
  return true;

#endif

}

bool CFileWriter::WriteMPIString(const string &str, unsigned short processor){
#ifdef HAVE_MPI
  CParallelDataSorter::CheckIntCount(str.size());
#endif

#ifdef HAVE_MPI

  startTime = WallTime();

  int ierr = MPI_SUCCESS;
  MPI_Status status{};

  /*--- Reset the file view. ---*/

  if (MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL) != MPI_SUCCESS)
    SU2_MPI::Error("Setting MPI output view failed", CURRENT_FUNCTION);

  if (SU2_MPI::GetRank() == processor)
    ierr = MPI_File_write_at(fhw, disp, str.c_str(), str.size(),
                      MPI_CHAR, &status);

  disp += str.size()*sizeof(char);
  if (rank == processor) fileSize += sizeof(char)*str.size();

  stopTime = WallTime();

  usedTime += stopTime - startTime;

  if (rank == processor) CheckMPIWrite(ierr, status, str.size());
  return true;

#else

  startTime = WallTime();

  unsigned long bytesWritten;
  bytesWritten = fwrite(str.c_str(), sizeof(char), str.size(), fhw);

  fileSize += bytesWritten;

  stopTime = WallTime();

  usedTime += stopTime - startTime;

  if (bytesWritten != str.size()) SU2_MPI::Error("Output text write failed", CURRENT_FUNCTION);
  return true;

#endif

}

bool CFileWriter::OpenMPIFile(string val_filename, bool append){
  fileStartTime = WallTime();

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

#ifdef HAVE_MPI
  int ierr;
  disp     = 0.0;

  /*--- Continue writing at the end of a file that exists, e.g. to add a zone to a multizone mesh file. ---*/

  if (append) {
    /*--- The file is created if it does not exist yet, for example for the first zone written to a file of its
     own, and the writing starts at its end. ---*/
    ierr = MPI_File_open(SU2_MPI::GetComm(), val_filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL,
                         &fhw);
    if (ierr != MPI_SUCCESS) {
      SU2_MPI::Error(string("Unable to open file ") + val_filename, CURRENT_FUNCTION);
    }
    MPI_Offset fileEnd;
    if (MPI_File_get_size(fhw, &fileEnd) != MPI_SUCCESS)
      SU2_MPI::Error("Reading MPI output size failed", CURRENT_FUNCTION);
    disp = fileEnd;

    fileSize = 0.0;
    usedTime = 0;

    return true;
  }

  /*--- All ranks open the file using MPI. Here, we try to open the file with
   exclusive so that an error is generated if the file exists. We always want
   to write a fresh output file, so we delete any existing files and create
   a new one. ---*/

  ierr = MPI_File_open(SU2_MPI::GetComm(), val_filename.c_str(),
                       MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);
  if (ierr != MPI_SUCCESS)  {
    /*--- The file exists: the failed open gives no valid handle to close, the file is deleted and created again. ---*/
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
  fhw = fopen(val_filename.c_str(), append ? "ab" : "wb");
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

std::pair<unsigned long, unsigned long> CFileWriter::GetRankOffset(unsigned long localCount) const {
  vector<unsigned long> counts(size, localCount);

  SU2_MPI::Allgather(&localCount, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  return std::make_pair(std::accumulate(counts.begin(), counts.begin() + rank, 0ul),
                        std::accumulate(counts.begin(), counts.end(), 0ul));
}

bool CFileWriter::WriteMPIStringAll(const string &str){

  /*--- Each rank writes its own text at the position that follows the text of the ranks before it. ---*/

  const unsigned long sizeInBytes = str.size();
  const auto offsetAndTotal = GetRankOffset(sizeInBytes);

  return WriteMPIBinaryDataAll(str.data(), sizeInBytes, offsetAndTotal.second, offsetAndTotal.first);
}

bool CFileWriter::CloseMPIFile(){

#ifdef HAVE_MPI
  /*--- Check the final extent: some MPI backends report success even after a short write. ---*/
  MPI_Offset expectedEnd = 0, actualEnd = 0;
  MPI_Allreduce(&disp, &expectedEnd, 1, MPI_OFFSET, MPI_MAX, SU2_MPI::GetComm());
  if (MPI_File_get_size(fhw, &actualEnd) != MPI_SUCCESS || actualEnd != expectedEnd)
    SU2_MPI::Error("MPI output file is incomplete", CURRENT_FUNCTION);

  /*--- All ranks close the file after writing. ---*/

  if (MPI_File_close(&fhw) != MPI_SUCCESS)
    SU2_MPI::Error("Closing MPI output failed", CURRENT_FUNCTION);
#else
  if (fclose(fhw) != 0) SU2_MPI::Error("Closing output failed", CURRENT_FUNCTION);
#endif

  const passivedouble elapsed = WallTime() - fileStartTime;
  passivedouble localBytes = SU2_TYPE::GetValue(fileSize), totalBytes = 0;
  SelectMPIWrapper<passivedouble>::W::Allreduce(&localBytes, &totalBytes, 1, MPI_DOUBLE, MPI_SUM,
                                             SU2_MPI::GetComm());
  fileSize = totalBytes;
  UpdateBandwidth(elapsed);

  return true;
}

