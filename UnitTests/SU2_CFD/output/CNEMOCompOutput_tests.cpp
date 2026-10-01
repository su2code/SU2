/*!
 * \file CNEMOCompOutput_tests.cpp
 * \brief Unit tests for NEMO residual history output.
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

#include "catch.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../SU2_CFD/include/output/CNEMOCompOutput.hpp"
#include "../../../SU2_CFD/include/solvers/CSolver.hpp"
#include "../../../SU2_CFD/include/variables/CVariable.hpp"

namespace {

std::unique_ptr<CConfig> MakeNEMOConfig(bool multizone) {
  std::stringstream options;
  options << "SOLVER= NEMO_EULER\n"
          << "GAS_MODEL= AIR-5\n"
          << "GAS_COMPOSITION= (0.77, 0.23, 0.0, 0.0, 0.0)\n"
          << "FLUID_MODEL= SU2_NONEQ\n"
          << "MATH_PROBLEM= DIRECT\n"
          << "MACH_NUMBER= 5.0\n"
          << "FREESTREAM_PRESSURE= 101325.0\n"
          << "FREESTREAM_TEMPERATURE= 288.15\n"
          << "FREESTREAM_TEMPERATURE_VE= 288.15\n";
  auto config = std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
  config->SetMultizone_Problem(multizone);
  return config;
}

class CTestGeometry final : public CGeometry {
 public:
  explicit CTestGeometry(unsigned short dimension) {
    nDim = dimension;
    nPoint = nPointDomain = 1;
    Global_nPoint = Global_nPointDomain = SU2_MPI::GetSize();
    MGLevel = MESH_0;
  }
};

class CTestSolver final : public CSolver {
 private:
  CVariable variables;

  CVariable* GetBaseClassPointerToNodes() override { return &variables; }

 public:
  CTestSolver(const CConfig* config, unsigned short dimension)
      : variables(0, dimension, config->GetnSpecies() + dimension + 2, config) {
    nDim = dimension;
    nVar = config->GetnSpecies() + dimension + 2;
    nPoint = nPointDomain = 1;

    Residual_RMS.resize(nVar);
    Residual_Max.resize(nVar);
    Residual_BGS.resize(nVar);
    Residual_Max_BGS.resize(nVar);
    Point_Max.resize(nVar);
    Point_Max_BGS.resize(nVar);
    Point_Max_Coord.resize(nVar, nDim);
    Point_Max_Coord_BGS.resize(nVar, nDim);
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      Point_Max[iVar] = Point_Max_BGS[iVar] = rank;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        Point_Max_Coord(iVar, iDim) = -(rank + 1.0) * (iDim + 1.0);
        Point_Max_Coord_BGS(iVar, iDim) = 10.0 * (rank + 1.0) + iDim + 1.0;
      }
    }

    SetBaseClassPointerToNodes();
    SetCFL_Local_Stats(1.0);
    SetResLinSolver(1.0);
  }

  void SeedResidualSums() {
    const double rankAmplitude = rank + 1.0;
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      const auto rmsExponent = static_cast<int>(iVar) + 1;
      const auto maxExponent = static_cast<int>(iVar) + 22;
      const auto bgsExponent = static_cast<int>(iVar) + 12;
      Residual_RMS[iVar] = rankAmplitude * rankAmplitude * std::pow(10.0, -2.0 * rmsExponent);
      Residual_Max[iVar] = rankAmplitude * std::pow(10.0, -maxExponent);
      Residual_BGS[iVar] = rankAmplitude * rankAmplitude * std::pow(10.0, -2.0 * bgsExponent);
      Residual_Max_BGS[iVar] = rankAmplitude * std::pow(10.0, -bgsExponent);
    }

    /* Species 2 is zero everywhere; species 3 has a positive remote maximum in MPI runs. */
    Residual_Max[2] = 0.0;
    if (rank == MASTER_NODE) Residual_Max[3] = 0.0;
  }
};

std::vector<std::string> ExpectedResidualFields(const std::string& prefix, unsigned short dimension,
                                                unsigned short nSpecies) {
  std::vector<std::string> fields;
  for (auto iSpecies = 0u; iSpecies < nSpecies; ++iSpecies) {
    fields.push_back(prefix + "_DENSITY_" + std::to_string(iSpecies));
  }
  fields.push_back(prefix + "_MOMENTUM-X");
  fields.push_back(prefix + "_MOMENTUM-Y");
  if (dimension == 3) fields.push_back(prefix + "_MOMENTUM-Z");
  fields.push_back(prefix + "_ENERGY");
  fields.push_back(prefix + "_ENERGY_VE");
  return fields;
}

std::vector<std::string> ExpectedFieldNames(const std::string& prefix, unsigned short dimension,
                                            unsigned short nSpecies) {
  std::vector<std::string> names;
  for (auto iSpecies = 0u; iSpecies < nSpecies; ++iSpecies) {
    names.push_back(prefix + "[Rho_" + std::to_string(iSpecies) + "]");
  }
  names.push_back(prefix + "[RhoU]");
  names.push_back(prefix + "[RhoV]");
  if (dimension == 3) names.push_back(prefix + "[RhoW]");
  names.push_back(prefix + "[RhoE]");
  names.push_back(prefix + "[RhoEve]");
  return names;
}

void CheckRegisteredGroup(CNEMOCompOutput& output, const std::string& fieldPrefix, const std::string& namePrefix,
                          unsigned short dimension, unsigned short nSpecies) {
  const auto expectedFields = ExpectedResidualFields(fieldPrefix, dimension, nSpecies);
  const auto expectedNames = ExpectedFieldNames(namePrefix, dimension, nSpecies);
  const auto registered = output.GetHistoryGroup(fieldPrefix + "_RES");

  REQUIRE(registered.size() == expectedFields.size());
  const auto& fields = output.GetHistoryFields();
  for (std::size_t i = 0; i < expectedFields.size(); ++i) {
    INFO("history field " << expectedFields[i]);
    REQUIRE(fields.count(expectedFields[i]) == 1);
    CHECK(fields.at(expectedFields[i]).fieldName == expectedNames[i]);
    CHECK(registered[i].fieldName == expectedNames[i]);
  }
}

void CheckLoadedResiduals(unsigned short dimension, bool multizone) {
  auto config = MakeNEMOConfig(multizone);
  REQUIRE(config->GetComm_Level() == COMM_FULL);
  CTestGeometry geometry(dimension);
  CTestSolver flow(config.get(), dimension);
  flow.SeedResidualSums();

  /* Exercise full communication with one local point and distinct residuals on each rank. */
  flow.SetResidual_RMS(&geometry, config.get());
  if (multizone) flow.SetResidual_BGS(&geometry, config.get());

  CNEMOCompOutput output(config.get(), dimension);
  output.SetHistoryOutputFields(config.get());
  const auto nSpecies = config->GetnSpecies();
  const unsigned short nVar = nSpecies + dimension + 2;
  const auto rmsFields = ExpectedResidualFields("RMS", dimension, nSpecies);
  const auto maxFields = ExpectedResidualFields("MAX", dimension, nSpecies);
  const auto bgsFields = ExpectedResidualFields("BGS", dimension, nSpecies);
  REQUIRE(rmsFields.size() == nVar);
  REQUIRE(maxFields.size() == nVar);
  REQUIRE(bgsFields.size() == nVar);

  std::vector<double> initialBGS;
  for (const auto& field : bgsFields) initialBGS.push_back(SU2_TYPE::GetValue(output.GetHistoryFieldValue(field)));

  std::array<CSolver*, MAX_SOLS> solvers{};
  solvers[FLOW_SOL] = &flow;
  output.LoadHistoryData(config.get(), &geometry, solvers.data());

  /* For amplitudes 1,...,N, the mean square is (N+1)(2N+1)/6 and the maximum is N.
   * These expectations do not reuse the production residual getters or indexing helper. */
  const auto nRanks = SU2_MPI::GetSize();
  const double rmsLogShift = 0.5 * std::log10((nRanks + 1.0) * (2.0 * nRanks + 1.0) / 6.0);
  const double maxLogShift = std::log10(nRanks);

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
    const auto rmsValue = SU2_TYPE::GetValue(output.GetHistoryFieldValue(rmsFields[iVar]));
    const auto maxValue = SU2_TYPE::GetValue(output.GetHistoryFieldValue(maxFields[iVar]));
    const auto bgsValue = SU2_TYPE::GetValue(output.GetHistoryFieldValue(bgsFields[iVar]));
    INFO("residual variable " << iVar);
    CHECK(std::isfinite(rmsValue));
    CHECK(rmsValue == Approx(-(static_cast<double>(iVar) + 1.0) + rmsLogShift));
    CHECK(std::isfinite(maxValue));
    const bool zeroMaximum = iVar == 2 || (iVar == 3 && nRanks == 1);
    CHECK(maxValue == Approx(zeroMaximum ? -32.0 : -(static_cast<double>(iVar) + 22.0) + maxLogShift));
    CHECK(std::isfinite(bgsValue));
    if (multizone) {
      CHECK(bgsValue == Approx(-(static_cast<double>(iVar) + 12.0) + rmsLogShift));
      CHECK(std::log10(SU2_TYPE::GetValue(flow.GetRes_Max_BGS(iVar))) ==
            Approx(-(static_cast<double>(iVar) + 12.0) + maxLogShift));
    } else {
      CHECK(bgsValue == initialBGS[iVar]);
    }
  }

  /* In MPI runs, species 3 has a positive maximum on the last rank despite rank 0 having zero. */
  CHECK(flow.GetPoint_Max(3) == static_cast<unsigned long>(nRanks - 1));
  if (multizone) CHECK(flow.GetPoint_Max_BGS(3) == static_cast<unsigned long>(nRanks - 1));
  for (unsigned short iDim = 0; iDim < dimension; ++iDim) {
    CHECK(SU2_TYPE::GetValue(flow.GetPoint_Max_Coord(3)[iDim]) == Approx(-nRanks * (iDim + 1.0)));
    if (multizone) {
      CHECK(SU2_TYPE::GetValue(flow.GetPoint_Max_Coord_BGS(3)[iDim]) == Approx(10.0 * nRanks + iDim + 1.0));
    }
  }
}

}  // namespace

TEST_CASE("NEMO residual history fields use the species-first layout", "[NEMO][Output]") {
  auto config = MakeNEMOConfig(false);

  for (unsigned short dimension : {2, 3}) {
    CAPTURE(dimension);
    CNEMOCompOutput output(config.get(), dimension);
    output.SetHistoryOutputFields(config.get());
    CheckRegisteredGroup(output, "RMS", "rms", dimension, config->GetnSpecies());
    CheckRegisteredGroup(output, "MAX", "max", dimension, config->GetnSpecies());
    CheckRegisteredGroup(output, "BGS", "bgs", dimension, config->GetnSpecies());
    CHECK(output.GetHistoryFields().count("MAX_DENSITY") == 0);
    CHECK(output.GetHistoryFields().count("BGS_DENSITY") == 0);
  }
}

TEST_CASE("NEMO history loads reduced residuals by species-first index", "[NEMO][Output]") {
  for (unsigned short dimension : {2, 3}) {
    CAPTURE(dimension);
    for (bool multizone : {false, true}) {
      CAPTURE(multizone);
      CheckLoadedResiduals(dimension, multizone);
    }
  }
}
