/*!
 * \file COutput_convergence_tests.cpp
 * \brief Tests for NEMO convergence field selection and migration.
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
#include <memory>
#include <sstream>
#include <string>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../SU2_CFD/include/output/CFlowCompOutput.hpp"
#include "../../../SU2_CFD/include/output/CMultizoneOutput.hpp"
#include "../../../SU2_CFD/include/output/CNEMOCompOutput.hpp"
#include "../../../SU2_CFD/include/solvers/CSolver.hpp"
#include "../../../SU2_CFD/include/variables/CVariable.hpp"

namespace {

std::unique_ptr<CConfig> MakeConvergenceConfig(const std::string& fields, bool nemo, bool multizone) {
  std::stringstream options;
  options << "SOLVER= " << (nemo ? "NEMO_EULER" : "EULER") << '\n';
  if (nemo) {
    options << "GAS_MODEL= AIR-5\n"
            << "GAS_COMPOSITION= (0.77, 0.23, 0.0, 0.0, 0.0)\n"
            << "FLUID_MODEL= SU2_NONEQ\n"
            << "FREESTREAM_TEMPERATURE_VE= 288.15\n";
  }
  options << "MATH_PROBLEM= DIRECT\n"
          << "MACH_NUMBER= 5.0\n"
          << "FREESTREAM_PRESSURE= 101325.0\n"
          << "FREESTREAM_TEMPERATURE= 288.15\n"
          << "CONV_FIELD= (" << fields << ")\n"
          << "CONV_RESIDUAL_MINVAL= -6\n"
          << "CONV_STARTITER= 1\n"
          << "HISTORY_OUTPUT= ITER\n"
          << "SCREEN_OUTPUT= INNER_ITER\n";
  auto config = std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
  config->SetMultizone_Problem(multizone);
  return config;
}

/* These are already reduced residuals. Reduction and indexing are tested separately. */
class CConvergenceSolver final : public CSolver {
  CVariable variables;
  const unsigned short nSpecies;
  CVariable* GetBaseClassPointerToNodes() override { return &variables; }

 public:
  explicit CConvergenceSolver(const CConfig* config)
      : variables(0, 2, config->GetnSpecies() + 4, config), nSpecies(config->GetnSpecies()) {
    nDim = 2;
    nVar = nSpecies + 4;
    Residual_RMS.resize(nVar);
    Residual_Max.resize(nVar);
    Residual_BGS.resize(nVar);
    SetBaseClassPointerToNodes();
    SetCFL_Local_Stats(1.0);
    SetResLinSolver(1.0);
  }

  void Seed(su2double firstSpecies, su2double lastSpecies, su2double energy) {
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      auto value = su2double(1.0);
      if (iVar == 0) value = firstSpecies;
      if (iVar == nSpecies - 1) value = lastSpecies;
      if (iVar == nSpecies + nDim) value = energy;
      Residual_RMS[iVar] = Residual_Max[iVar] = Residual_BGS[iVar] = value;
    }
  }
};

void CheckNEMOConvergence(const std::string& prefix, bool relative, bool zoned) {
  const bool multizone = zoned || prefix == "BGS";
  auto config = MakeConvergenceConfig("RMS_DENSITY_0", true, multizone);
  const auto fieldPrefix = (relative ? "REL_" : "") + prefix;
  const std::string suffix = zoned ? "[0]" : "";
  const std::array<std::string, 3> fields = {
      fieldPrefix + "_DENSITY_0" + suffix,
      fieldPrefix + "_DENSITY_" + std::to_string(config->GetnSpecies() - 1) + suffix, fieldPrefix + "_ENERGY" + suffix};
  const auto selection = fields[0] + ", " + fields[1] + ", " + fields[2];
  if (!zoned) config = MakeConvergenceConfig(selection, true, multizone);

  CNEMOCompOutput output(config.get(), 2);
  output.PreprocessHistoryOutput(config.get(), false);
  std::array<CConfig*, 1> configs = {config.get()};
  std::array<COutput*, 1> outputs = {&output};
  std::unique_ptr<CConfig> driverConfig;
  std::unique_ptr<CMultizoneOutput> driverOutput;
  COutput* monitor = &output;
  CConfig* monitorConfig = config.get();
  if (zoned) {
    driverConfig = MakeConvergenceConfig(selection, true, true);
    driverOutput = std::make_unique<CMultizoneOutput>(driverConfig.get(), configs.data(), 2);
    driverOutput->PreprocessMultizoneHistoryOutput(outputs.data(), configs.data(), driverConfig.get(), false);
    monitor = driverOutput.get();
    monitorConfig = driverConfig.get();
  }
  const auto retained = monitor->GetResidualConvFields();
  REQUIRE(retained.size() == fields.size());
  for (std::size_t i = 0; i < fields.size(); ++i) CHECK(retained[i].first == fields[i]);

  CConvergenceSolver flow(config.get());
  CGeometry geometry;
  std::array<CSolver*, MAX_SOLS> solvers{};
  solvers[FLOW_SOL] = &flow;
  const auto sample = [&](unsigned long iteration, su2double first, su2double last, su2double energy) {
    flow.Seed(first, last, energy);
    output.SetIteration(0, 0, iteration);
    output.SetHistoryOutput(&geometry, solvers.data(), config.get());
    if (!zoned && !relative) return output.GetConvergence();
    if (zoned) driverOutput->LoadMultizoneHistoryData(outputs.data(), configs.data());
    /* Test the predicate on freshly loaded fields. For unzoned relative residuals,
     * this does not test or change the existing once-per-iteration ordering of
     * convergence monitoring and relative-residual postprocessing. */
    return monitor->ConvergenceMonitoring(monitorConfig, iteration);
  };

  CHECK_FALSE(sample(0, 1.0, 1.0, 1.0));
  CHECK_FALSE(sample(2, 1e-3, 1e-10, 1e-10));
  CHECK_FALSE(sample(3, 1e-10, 1e-3, 1e-10));
  CHECK_FALSE(sample(4, 1e-10, 1e-10, 1e-3));
  CHECK(sample(5, 1e-9, 1e-9, 1e-9));
}

/* Called only by individually launched hidden tests because SU2_MPI::Error aborts. */
void RejectObsoleteNEMOField(const std::string& field, bool zoned, bool mixed = false) {
  const auto selection = field + (mixed ? ", MAX_ENERGY" : "");
  auto config = MakeConvergenceConfig(zoned ? "RMS_DENSITY_0" : selection, true, true);
  CNEMOCompOutput output(config.get(), 2);
  output.PreprocessHistoryOutput(config.get(), false);
  if (zoned) {
    auto driverConfig = MakeConvergenceConfig(selection, true, true);
    std::array<CConfig*, 1> configs = {config.get()};
    std::array<COutput*, 1> outputs = {&output};
    CMultizoneOutput driverOutput(driverConfig.get(), configs.data(), 2);
    driverOutput.PreprocessMultizoneHistoryOutput(outputs.data(), configs.data(), driverConfig.get(), false);
  }
  /* A return is a successful subprocess exit, which the migration runner rejects. */
  SUCCEED("Obsolete convergence field was not rejected");
}

}  // namespace

TEST_CASE("NEMO explicit density fields retain every convergence criterion", "[NEMO][Output][Convergence]") {
  for (const std::string prefix : {"MAX", "BGS"}) {
    for (const bool relative : {false, true}) {
      for (const bool zoned : {false, true}) {
        INFO("prefix=" << prefix << ", relative=" << relative << ", zoned=" << zoned);
        CheckNEMOConvergence(prefix, relative, zoned);
      }
    }
  }
}

TEST_CASE("Compressible density convergence names remain valid", "[Output][Convergence]") {
  const std::array<std::string, 4> fields = {"MAX_DENSITY", "BGS_DENSITY", "REL_MAX_DENSITY", "REL_BGS_DENSITY"};
  const auto selection = fields[0] + ", " + fields[1] + ", " + fields[2] + ", " + fields[3];
  auto config = MakeConvergenceConfig(selection, false, true);
  CFlowCompOutput output(config.get(), 2);
  output.PreprocessHistoryOutput(config.get(), false);
  const auto retained = output.GetResidualConvFields();
  REQUIRE(retained.size() == fields.size());
  for (std::size_t i = 0; i < fields.size(); ++i) CHECK(retained[i].first == fields[i]);

  auto driverConfig =
      MakeConvergenceConfig("MAX_DENSITY[0], BGS_DENSITY[0], REL_MAX_DENSITY[0], REL_BGS_DENSITY[0]", false, true);
  std::array<CConfig*, 1> configs = {config.get()};
  std::array<COutput*, 1> outputs = {&output};
  CMultizoneOutput driverOutput(driverConfig.get(), configs.data(), 2);
  driverOutput.PreprocessMultizoneHistoryOutput(outputs.data(), configs.data(), driverConfig.get(), false);
  const auto zonedFields = driverOutput.GetResidualConvFields();
  REQUIRE(zonedFields.size() == fields.size());
  for (std::size_t i = 0; i < fields.size(); ++i) CHECK(zonedFields[i].first == fields[i] + "[0]");
}

TEST_CASE("NEMO rejects obsolete MAX convergence", "[.NEMOObsoleteMAX]") {
  RejectObsoleteNEMOField("MAX_DENSITY", false);
}
TEST_CASE("NEMO rejects obsolete BGS convergence", "[.NEMOObsoleteBGS]") {
  RejectObsoleteNEMOField("BGS_DENSITY", false);
}
TEST_CASE("NEMO rejects obsolete relative MAX convergence", "[.NEMOObsoleteRelativeMAX]") {
  RejectObsoleteNEMOField("REL_MAX_DENSITY", false);
}
TEST_CASE("NEMO rejects obsolete relative BGS convergence", "[.NEMOObsoleteRelativeBGS]") {
  RejectObsoleteNEMOField("REL_BGS_DENSITY", false);
}
TEST_CASE("NEMO rejects obsolete zoned MAX convergence", "[.NEMOObsoleteZonedMAX]") {
  RejectObsoleteNEMOField("MAX_DENSITY[0]", true);
}
TEST_CASE("NEMO rejects obsolete zoned BGS convergence", "[.NEMOObsoleteZonedBGS]") {
  RejectObsoleteNEMOField("BGS_DENSITY[0]", true);
}
TEST_CASE("NEMO rejects obsolete zoned relative MAX convergence", "[.NEMOObsoleteZonedRelativeMAX]") {
  RejectObsoleteNEMOField("REL_MAX_DENSITY[0]", true);
}
TEST_CASE("NEMO rejects obsolete zoned relative BGS convergence", "[.NEMOObsoleteZonedRelativeBGS]") {
  RejectObsoleteNEMOField("REL_BGS_DENSITY[0]", true);
}
TEST_CASE("NEMO rejects obsolete density in mixed convergence criteria", "[.NEMOObsoleteMixed]") {
  RejectObsoleteNEMOField("MAX_DENSITY", false, true);
}
