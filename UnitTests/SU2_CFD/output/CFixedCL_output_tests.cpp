/*!
 * \file CFixedCL_output_tests.cpp
 * \brief Fixed-CL finite-difference output must not depend on screen frequency.
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */

#include "catch.hpp"

#include <sstream>

#include "../../../Common/include/CConfig.hpp"
#include "../../../SU2_CFD/include/output/CFlowOutput.hpp"

namespace {
class CFixedCLOutputTest : public CFlowOutput {
 public:
  explicit CFixedCLOutputTest(const CConfig* config) : CFlowOutput(config, 2, false) {
    if (config->GetFixed_CL_Mode()) {
      for (const auto* field : {"AOA", "PREV_AOA", "CL_DRIVER_COMMAND"}) {
        AddHistoryOutput(field, field, ScreenOutputFormat::FIXED, "FIXED_CL", "Test field");
      }
    }
  }
  using CFlowOutput::WriteAdditionalFiles;
  using COutput::WriteScreenOutput;
  unsigned long StoredIteration() const { return lastInnerIter; }
  using COutput::SetHistoryOutputValue;
};
}  // namespace

TEST_CASE("Fixed-CL metadata iteration capture is independent of screen frequency", "[output][fixed-cl]") {
  const auto frequency = GENERATE(10, 50, 1000);
  std::stringstream options;
  options << "SOLVER= EULER\n"
          << "MATH_PROBLEM= DIRECT\n"
          << "MACH_NUMBER= 0.8\n"
          << "FIXED_CL_MODE= YES\n"
          << "ITER= 2000\n"
          << "SCREEN_WRT_FREQ_INNER= " << frequency << "\n"
          << "SCREEN_OUTPUT= INNER_ITER\n";
  CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
  CFixedCLOutputTest output(&config);
  output.SetIteration(0, 0, 533);
  output.SetHistoryOutputValue("AOA", 1.141);
  output.SetHistoryOutputValue("PREV_AOA", 1.140);
  output.SetHistoryOutputValue("CL_DRIVER_COMMAND", 0.001);

  // Before the FD phase the ordinary screen frequency still applies.
  CHECK_FALSE(output.WriteScreenOutput(&config));

  config.SetFinite_Difference_Mode(true);
  // The metadata callback captures the iteration even without a screen line.
  output.WriteAdditionalFiles(&config);
  CHECK(output.StoredIteration() == 532);
  CHECK_FALSE(output.WriteScreenOutput(&config));

  output.SetHistoryOutputValue("CL_DRIVER_COMMAND", 0.0);
  // Ordinary intermediate FD iterations are still throttled.
  CHECK_FALSE(output.WriteScreenOutput(&config));

  output.SetIteration(0, 0, 586);
  output.SetHistoryOutputValue("AOA", 1.140);
  // Restoring AoA does not force an off-frequency screen line.
  CHECK_FALSE(output.WriteScreenOutput(&config));

  config.SetFinite_Difference_Mode(false);
  CHECK_FALSE(output.WriteScreenOutput(&config));
  output.SetIteration(0, 0, 1000);
  CHECK(output.WriteScreenOutput(&config));
}

TEST_CASE("Non-fixed-CL output retains its screen cadence", "[output][fixed-cl]") {
  std::stringstream options;
  options << "SOLVER= EULER\n"
          << "FIXED_CL_MODE= NO\n"
          << "ITER= 100\n"
          << "SCREEN_WRT_FREQ_INNER= 10\n";
  CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
  CFixedCLOutputTest output(&config);
  // There are no fixed-CL history fields in this fixture.
  config.SetFinite_Difference_Mode(true);
  output.WriteAdditionalFiles(&config);
  output.SetIteration(0, 0, 3);
  CHECK_FALSE(output.WriteScreenOutput(&config));
  output.SetIteration(0, 0, 10);
  CHECK(output.WriteScreenOutput(&config));
}
