/*!
 * \file CIncIdealGasNASA_tests.cpp
 * \brief Unit tests for the NASA polynomial fluid model class (Single Range).
 * \author Pratyksh Gupta
 * \version 8.4.0 "Harrier"
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
#include <cmath>
#include <sstream>
#include "../../../SU2_CFD/include/fluid/CIncIdealGasNASA.hpp"

// Helper to config an already constructed CConfig object
void ConfigureNASA(CConfig& config) {
  // Nitrogen (N2) Low Range Coefficients (200K - 1000K) in NASA-7 format
  // NASA-9 indexing: a1=0, a2=0 (no inverse T terms), a3-a7 (polynomial), a8-a9 (integration constants)
  su2double nasa_coeffs[9] = {
      0.0,            // a1 (T^-2 term, not used in NASA-7)
      0.0,            // a2 (T^-1 term, not used in NASA-7)
      3.298677,       // a3 (constant term)
      1.4082404e-3,   // a4 (T term)
      -3.963222e-6,   // a5 (T^2 term)
      5.641515e-9,    // a6 (T^3 term)
      -2.444854e-12,  // a7 (T^4 term)
      -1020.8999,     // a8 (H integration constant)
      3.950372        // a9 (S integration constant)
  };

  config.SetCp_NASA_Format(true);
  for (int i = 0; i < 9; ++i) {
    config.SetCp_PolyCoeff(i, nasa_coeffs[i]);
  }

  // Ensure we set the temperature limits that SetCpModel reads
  // Note: These limits are usually set via TEMPERATURE_LIMITS option, here manually.
  config.SetTemperatureLimits(0, 200.0);
  config.SetTemperatureLimits(1, 1000.0);
}

// Helper to populate the stringstream with minimal required options
void SetupConfigStream(std::stringstream& ss) {
  ss << "SOLVER= EULER\n";
  ss << "MATH_PROBLEM= DIRECT\n";
  ss << "KIND_TURB_MODEL= NONE\n";
  ss << "KIND_TRANS_MODEL= NONE\n";
}

TEST_CASE("NASA polynomial Cp calculation", "[CIncIdealGasNASA]") {
  const su2double R_N2 = 296.8;
  const su2double P = 101325.0;
  const su2double T_ref = 298.15;

  std::stringstream ss;
  SetupConfigStream(ss);
  CConfig config(ss, SU2_COMPONENT::SU2_CFD, false);
  ConfigureNASA(config);

  CIncIdealGasNASA<9> nasa_model(R_N2, P, T_ref);
  nasa_model.SetCpModel(&config, T_ref);

  SECTION("Valid temperature range (300K)") {
    nasa_model.SetTDState_T(300.0);
    const su2double cp = nasa_model.GetCp();
    // Cp of N2 at 300K is approx 1040 J/kgK.
    CHECK(cp == Approx(1040.0).epsilon(0.01));
  }
}

TEST_CASE("NASA polynomial enthalpy-temperature inversion", "[CIncIdealGasNASA]") {
  const su2double R_N2 = 296.8;
  const su2double P = 101325.0;
  const su2double T_ref = 298.15;

  std::stringstream ss;
  SetupConfigStream(ss);
  CConfig config(ss, SU2_COMPONENT::SU2_CFD, false);
  ConfigureNASA(config);

  CIncIdealGasNASA<9> nasa_model(R_N2, P, T_ref);
  nasa_model.SetCpModel(&config, T_ref);

  SECTION("Inversion Consistency") {
    su2double T_target = 800.0;
    nasa_model.SetTDState_T(T_target);
    su2double H_target = nasa_model.GetEnthalpy();

    // Reset state to something else
    nasa_model.SetTDState_T(300.0);

    // Invert from H_target
    nasa_model.SetTDState_h(H_target);
    su2double T_recovered = nasa_model.GetTemperature();

    CHECK(T_recovered == Approx(T_target).margin(0.001));
  }
}
