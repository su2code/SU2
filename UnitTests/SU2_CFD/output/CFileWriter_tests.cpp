/*!
 * \file CFileWriter_tests.cpp
 * \brief Tests of the size, time and bandwidth accounting of the parallel file writers.
 * \version 8.6.0 "Harrier"
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

#include "catch.hpp"

#include <chrono>
#include <climits>
#include <cstdio>
#include <string>
#include <thread>

#include "../../../SU2_CFD/include/output/filewriter/CFileWriter.hpp"

namespace {

/*--- Gives access to the low level write functions of the base class. ---*/
class CTestFileWriter final : public CFileWriter {
 public:
  CTestFileWriter() : CFileWriter(".bin") {}
  using CFileWriter::CloseMPIFile;
  using CFileWriter::OpenMPIFile;
  using CFileWriter::WriteMPIBinaryData;
  using CFileWriter::WriteMPIBinaryDataAll;
  using CFileWriter::WriteMPIString;
};

}  // namespace

TEST_CASE("Int count limits of the output", "[Output]") {
  /*--- The largest counts that fit, larger ones stop the run. ---*/
  CParallelDataSorter::CheckIntCount(INT_MAX);
  CParallelDataSorter::CheckIntCount(INT_MAX / 8, 8);
  CParallelDataSorter::CheckIntCount(0, 0);
}

TEST_CASE("File writer accounting", "[Output]") {
  const int rank = SU2_MPI::GetRank(), size = SU2_MPI::GetSize();
  const std::string fileName = "file_writer_accounting";

  CTestFileWriter writer;
  writer.OpenMPIFile(fileName);

  /*--- A header written by one rank, then 3 bytes from each rank. ---*/
  const int header = 123;
  writer.WriteMPIString("H", MASTER_NODE);
  writer.WriteMPIBinaryData(&header, sizeof(header), MASTER_NODE);
  const char bytes[] = {'a', 'b', 'c'};
  writer.WriteMPIBinaryDataAll(bytes, 3, 3 * size, 3 * rank);

  /*--- Work between the writes and the close is part of the write time. ---*/
  std::this_thread::sleep_for(std::chrono::milliseconds(80));
  writer.CloseMPIFile();

  const auto fileSize = SU2_TYPE::GetValue(writer.GetFilesize());
  const auto usedTime = SU2_TYPE::GetValue(writer.GetUsedTime());

  /*--- The header is counted once, not once per rank. ---*/
  CHECK(fileSize == 1 + sizeof(header) + 3 * size);
  CHECK(usedTime >= 0.075);
  CHECK(SU2_TYPE::GetValue(writer.GetBandwidth()) == Approx(fileSize / (1.0e6 * usedTime)));

  SU2_MPI::Barrier(SU2_MPI::GetComm());
  if (rank == MASTER_NODE) std::remove((fileName + ".bin").c_str());
}
