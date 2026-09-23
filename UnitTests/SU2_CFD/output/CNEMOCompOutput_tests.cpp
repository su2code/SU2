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
          << "FREESTREAM_TEMPERATURE_VE= 288.15\n"
          << "COMM_LEVEL= MINIMAL\n";
  auto config = std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
  config->SetMultizone_Problem(multizone);
  return config;
}

class CTestGeometry final : public CGeometry {
 public:
  explicit CTestGeometry(unsigned short dimension) {
    nDim = dimension;
    nPoint = nPointDomain = Global_nPoint = Global_nPointDomain = 1;
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

    SetBaseClassPointerToNodes();
    SetCFL_Local_Stats(1.0);
    SetResLinSolver(1.0);
  }

  void SeedResidualSums() {
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      const auto rmsExponent = static_cast<int>(iVar) + 1;
      const auto maxExponent = static_cast<int>(iVar) + 22;
      const auto bgsExponent = static_cast<int>(iVar) + 12;
      Residual_RMS[iVar] = std::pow(10.0, -2.0 * rmsExponent);
      Residual_Max[iVar] = std::pow(10.0, -maxExponent);
      Residual_BGS[iVar] = std::pow(10.0, -2.0 * bgsExponent);
    }

    /* An exactly-zero species maximum exercises the finite output floor. */
    Residual_Max[2] = 0.0;
  }
};

std::vector<std::string> ExpectedResidualFields(const std::string& prefix, unsigned short dimension) {
  std::vector<std::string> fields;
  for (unsigned short iSpecies = 0; iSpecies < 5; ++iSpecies) {
    fields.push_back(prefix + "_DENSITY_" + std::to_string(iSpecies));
  }
  fields.push_back(prefix + "_MOMENTUM-X");
  fields.push_back(prefix + "_MOMENTUM-Y");
  if (dimension == 3) fields.push_back(prefix + "_MOMENTUM-Z");
  fields.push_back(prefix + "_ENERGY");
  fields.push_back(prefix + "_ENERGY_VE");
  return fields;
}

std::vector<std::string> ExpectedFieldNames(const std::string& prefix, unsigned short dimension) {
  std::vector<std::string> names;
  for (unsigned short iSpecies = 0; iSpecies < 5; ++iSpecies) {
    names.push_back(prefix + "[Rho_" + std::to_string(iSpecies) + "]");
  }
  names.push_back(prefix + "[RhoU]");
  names.push_back(prefix + "[RhoV]");
  if (dimension == 3) names.push_back(prefix + "[RhoW]");
  names.push_back(prefix + "[RhoE]");
  names.push_back(prefix + "[RhoEve]");
  return names;
}

void CheckRegisteredGroup(CNEMOCompOutput& output, const std::string& group, const std::string& prefix,
                          unsigned short dimension) {
  const auto expectedFields = ExpectedResidualFields(group == "MAX_RES" ? "MAX" : "BGS", dimension);
  const auto expectedNames = ExpectedFieldNames(prefix, dimension);
  const auto registered = output.GetHistoryGroup(group);

  REQUIRE(registered.size() == expectedFields.size());
  const auto& fields = output.GetHistoryFields();
  for (std::size_t i = 0; i < expectedFields.size(); ++i) {
    INFO("history field " << expectedFields[i]);
    REQUIRE(fields.count(expectedFields[i]) == 1);
    CHECK(fields.at(expectedFields[i]).fieldName == expectedNames[i]);
    CHECK(registered[i].fieldName == expectedNames[i]);
  }
}

void CheckLoadedResiduals(unsigned short dimension) {
  auto config = MakeNEMOConfig(true);
  CTestGeometry geometry(dimension);
  CTestSolver flow(config.get(), dimension);
  flow.SeedResidualSums();

  /* Exercise the actual reductions used before MAX and multizone BGS output. */
  flow.SetResidual_RMS(&geometry, config.get());
  flow.SetResidual_BGS(&geometry, config.get());

  CNEMOCompOutput output(config.get(), dimension);
  output.SetHistoryOutputFields(config.get());
  std::array<CSolver*, MAX_SOLS> solvers{};
  solvers[FLOW_SOL] = &flow;
  solvers[MESH_SOL] = &flow;
  output.LoadHistoryData(config.get(), &geometry, solvers.data());

  const unsigned short nVar = config->GetnSpecies() + dimension + 2;
  const auto maxFields = ExpectedResidualFields("MAX", dimension);
  const auto bgsFields = ExpectedResidualFields("BGS", dimension);
  REQUIRE(maxFields.size() == nVar);
  REQUIRE(bgsFields.size() == nVar);

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
    const auto maxValue = SU2_TYPE::GetValue(output.GetHistoryFieldValue(maxFields[iVar]));
    const auto bgsValue = SU2_TYPE::GetValue(output.GetHistoryFieldValue(bgsFields[iVar]));
    INFO("residual variable " << iVar);
    CHECK(std::isfinite(maxValue));
    CHECK(maxValue == Approx(iVar == 2 ? -32.0 : -(static_cast<double>(iVar) + 22.0)));
    CHECK(bgsValue == Approx(-(static_cast<double>(iVar) + 12.0)));
  }
}

}  // namespace

TEST_CASE("NEMO MAX and BGS history fields use the species-first layout", "[NEMO][Output]") {
  auto config = MakeNEMOConfig(false);

  SECTION("2D") {
    CNEMOCompOutput output(config.get(), 2);
    output.SetHistoryOutputFields(config.get());
    CheckRegisteredGroup(output, "MAX_RES", "max", 2);
    CheckRegisteredGroup(output, "BGS_RES", "bgs", 2);
    CHECK(output.GetHistoryFields().count("MAX_DENSITY") == 0);
    CHECK(output.GetHistoryFields().count("BGS_DENSITY") == 0);
  }

  SECTION("3D") {
    CNEMOCompOutput output(config.get(), 3);
    output.SetHistoryOutputFields(config.get());
    CheckRegisteredGroup(output, "MAX_RES", "max", 3);
    CheckRegisteredGroup(output, "BGS_RES", "bgs", 3);
    CHECK(output.GetHistoryFields().count("MAX_DENSITY") == 0);
    CHECK(output.GetHistoryFields().count("BGS_DENSITY") == 0);
  }
}

TEST_CASE("NEMO history loads finite MAX and multizone BGS residuals by species-first index", "[NEMO][Output]") {
  SECTION("2D") { CheckLoadedResiduals(2); }
  SECTION("3D") { CheckLoadedResiduals(3); }
}
