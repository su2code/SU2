/*!
 * \file CParaviewPVTUFileWriter_tests.cpp
 * \brief Tests of the field selection of the Paraview writers and of the Paraview parallel (.pvtu) writer.
 * \author A. Rausa
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

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "../../../SU2_CFD/include/output/filewriter/CParaviewFieldLayout.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CParaviewPVTUFileWriter.hpp"

namespace {

/*--- Output data of one hexahedron per rank, field iField of local point iPoint is 1000*rank + 10*iPoint + iField.
 * ---*/
class CHexaDataSorter final : public CParallelDataSorter {
 public:
  explicit CHexaDataSorter(const std::vector<std::string>& names) : CParallelDataSorter(nullptr, names) {
    nDim = 3;
    nPoints = 8;
    nPointsGlobal = nPoints * size;
    linearPartitioner.Initialize(nPointsGlobal, 0);

    dataBuffer = new passivedouble[nPoints * GlobalField_Counter];
    for (auto iPoint = 0ul; iPoint < nPoints; iPoint++)
      for (auto iField = 0u; iField < GlobalField_Counter; iField++)
        dataBuffer[iPoint * GlobalField_Counter + iField] = Value(rank, iPoint, iField);

    Conn_Hexa_Par = new int[nPoints];
    for (auto iPoint = 0ul; iPoint < nPoints; iPoint++) Conn_Hexa_Par[iPoint] = rank * nPoints + iPoint + 1;

    nElem = 1;
    nElemGlobal = size;
    nConn = nPoints;
    nConnGlobal = nPointsGlobal;
    nElemPerType[TypeMap.at(HEXAHEDRON)] = nElem;
    nElemPerTypeGlobal[TypeMap.at(HEXAHEDRON)] = nElemGlobal;
    for (int iRank = 0; iRank <= size; iRank++) {
      nElem_Cum[iRank] = iRank * nElem;
      nElemConn_Cum[iRank] = iRank * nConn;
    }
    connectivitySorted = true;
  }

  static double Value(int iRank, unsigned long iPoint, unsigned short iField) {
    return 1000.0 * iRank + 10.0 * iPoint + iField;
  }
};

/*--- Read the arrays of the appended raw data of a .vtu file, each one is preceded by its size in bytes. ---*/
std::vector<std::vector<char>> ReadAppendedArrays(const std::string& fileName, std::string& header) {
  std::ifstream file(fileName, std::ios::binary);
  const std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

  const std::string start = "<AppendedData encoding=\"raw\">\n_";
  const auto begin = content.find(start);
  REQUIRE(begin != std::string::npos);
  header = content.substr(0, begin);

  const std::string end = "\n</AppendedData>";
  const auto stop = content.rfind(end);
  std::vector<std::vector<char>> arrays;
  for (auto pos = begin + start.size(); pos < stop;) {
    uint64_t nBytes = 0;
    std::memcpy(&nBytes, content.data() + pos, sizeof(nBytes));
    pos += sizeof(nBytes);
    arrays.emplace_back(content.begin() + pos, content.begin() + pos + nBytes);
    pos += nBytes;
  }
  return arrays;
}

template <class T>
std::vector<T> AsValues(const std::vector<char>& bytes) {
  std::vector<T> values(bytes.size() / sizeof(T));
  std::memcpy(values.data(), bytes.data(), bytes.size());
  return values;
}

}  // namespace

TEST_CASE("Paraview field selection", "[Output]") {
  const auto fields3D = GetParaviewFields({"x", "y", "z", "\"Density\"", "Momentum_x", "Momentum_y", "Momentum_z"}, 3);
  REQUIRE(fields3D.size() == 2);
  CHECK(fields3D[0].first == 3);
  CHECK(fields3D[0].components == 1);
  CHECK(fields3D[0].name == "Density");
  CHECK(fields3D[1].first == 4);
  CHECK(fields3D[1].components == 3);
  CHECK(fields3D[1].name == "Momentum");

  const auto fields2D = GetParaviewFields({"x", "y", "Momentum_x", "Momentum_y", "Energy"}, 2);
  REQUIRE(fields2D.size() == 2);
  CHECK(fields2D[0].components == 3);
  CHECK(fields2D[0].name == "Momentum");
  CHECK(fields2D[1].first == 4);
  CHECK(fields2D[1].name == "Energy");

  /*--- Incomplete or unrelated components are scalars. ---*/
  const auto partial = GetParaviewFields({"x", "y", "z", "Momentum_x", "Other_y", "Momentum_z", "Flux_x_extra"}, 3);
  REQUIRE(partial.size() == 4);
  for (const auto& field : partial) CHECK(field.components == 1);
  CHECK(partial[0].name == "Momentum_x");
  CHECK(partial[3].name == "Flux_x_extra");
}

TEST_CASE("Paraview parallel writer", "[Output]") {
  const int rank = SU2_MPI::GetRank();
  const std::string fileName = "pvtu_writer_test";
  const std::vector<std::string> names = {"x", "y", "z", "Density", "Momentum_x", "Momentum_y", "Momentum_z"};

  CHexaDataSorter data(names);
  CParaviewPVTUFileWriter writer(&data, true);
  writer.WriteData(fileName);

  /*--- The index lists the piece of each rank and the fields. ---*/
  std::ifstream indexFile(fileName + ".pvtu");
  const std::string index((std::istreambuf_iterator<char>(indexFile)), std::istreambuf_iterator<char>());
  const auto pieceName = fileName + "_" + std::to_string(rank) + ".vtu";
  CHECK(index.find("<Piece Source=\"" + fileName + "/" + pieceName + "\"/>") != std::string::npos);
  CHECK(index.find("Name=\"Density\" NumberOfComponents=\"1\"") != std::string::npos);
  CHECK(index.find("Name=\"Momentum\" NumberOfComponents=\"3\"") != std::string::npos);

  /*--- The piece has the 8 points and the hexahedron of this rank. ---*/
  std::string header;
  const auto arrays = ReadAppendedArrays(fileName + "/" + pieceName, header);
  CHECK(header.find("NumberOfPoints=\"8\" NumberOfCells=\"1\"") != std::string::npos);
  REQUIRE(arrays.size() == 6);

  const auto coords = AsValues<double>(arrays[0]);
  const auto connectivity = AsValues<int32_t>(arrays[1]);
  const auto offsets = AsValues<int32_t>(arrays[2]);
  const auto types = AsValues<uint8_t>(arrays[3]);
  const auto density = AsValues<double>(arrays[4]);
  const auto momentum = AsValues<double>(arrays[5]);

  REQUIRE(coords.size() == 8 * 3);
  REQUIRE(connectivity.size() == 8);
  REQUIRE(density.size() == 8);
  REQUIRE(momentum.size() == 8 * 3);
  CHECK(offsets == std::vector<int32_t>{8});
  CHECK(types == std::vector<uint8_t>{HEXAHEDRON});

  for (auto iPoint = 0ul; iPoint < 8; iPoint++) {
    CHECK(connectivity[iPoint] == static_cast<int32_t>(iPoint));
    for (auto iDim = 0u; iDim < 3; iDim++) {
      CHECK(coords[3 * iPoint + iDim] == CHexaDataSorter::Value(rank, iPoint, iDim));
      CHECK(momentum[3 * iPoint + iDim] == CHexaDataSorter::Value(rank, iPoint, 4 + iDim));
    }
    CHECK(density[iPoint] == CHexaDataSorter::Value(rank, iPoint, 3));
  }

  SU2_MPI::Barrier(SU2_MPI::GetComm());
  if (rank == MASTER_NODE) {
    std::filesystem::remove_all(fileName);
    std::filesystem::remove(fileName + ".pvtu");
  }
}
