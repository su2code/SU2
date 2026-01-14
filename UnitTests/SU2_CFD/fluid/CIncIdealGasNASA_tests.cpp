/*!
 * \file CIncIdealGasNASA_tests.cpp
 * \brief Unit tests for the NASA polynomial fluid model class.
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
  su2double nasa_low[7] = {3.298677, 1.4082404e-3, -3.963222e-6, 5.641515e-9, -2.444854e-12, -1020.8999, 3.950372};
  su2double nasa_high[7] = {2.92664, 1.4879768e-3, -5.68476e-7, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528};

  for (int i = 0; i < 7; ++i) {
    config.SetNASA_CoeffLowND(i, nasa_low[i]);
    config.SetNASA_CoeffHighND(i, nasa_high[i]);
  }
  config.SetNASA_TempLow(200.0);
  config.SetNASA_TempMid(1000.0);
  config.SetNASA_TempHigh(6000.0);
  // Ensure we set the temperature limits that SetCpModel reads
  config.SetTemperatureLimits(0, 200.0);
  config.SetTemperatureLimits(1, 6000.0);
}

// Helper to populate the stringstream with minimal required options
void SetupConfigStream(std::stringstream& ss) {
  ss << "SOLVER= EULER" << std::endl;
  ss << "MATH_PROBLEM= DIRECT" << std::endl;
  ss << "KIND_TURB_MODEL= NONE" << std::endl;
  ss << "KIND_TRANS_MODEL= NONE" << std::endl;
}

TEST_CASE("NASA polynomial Cp calculation", "[CIncIdealGasNASA]") {
  const su2double R_N2 = 296.8;
  const su2double P = 101325.0;
  const su2double T_ref = 298.15;

  std::stringstream ss;
  SetupConfigStream(ss);
  CConfig config(ss, SU2_COMPONENT::SU2_CFD, false);
  ConfigureNASA(config);

  CIncIdealGasNASA<7> nasa_model(R_N2, P, T_ref);
  nasa_model.SetCpModel(&config, T_ref);

  SECTION("Low temperature range (300K)") {
    nasa_model.SetTDState_T(300.0);
    const su2double cp = nasa_model.GetCp();
    CHECK(cp == Approx(1040.0).epsilon(0.01));
  }

  SECTION("High temperature range (1500K)") {
    nasa_model.SetTDState_T(1500.0);
    const su2double cp = nasa_model.GetCp();
    CHECK(cp == Approx(1240.0).epsilon(0.01));
  }

  SECTION("Temperature range transition") {
    nasa_model.SetTDState_T(999.0);
    const su2double cp_below = nasa_model.GetCp();

    nasa_model.SetTDState_T(1001.0);
    const su2double cp_above = nasa_model.GetCp();

    const su2double discontinuity = std::abs(cp_above - cp_below);
    CHECK(discontinuity < 50.0);
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

  CIncIdealGasNASA<7> nasa_model(R_N2, P, T_ref);
  nasa_model.SetCpModel(&config, T_ref);

  SECTION("Round-trip T->H->T at 500K") {
    const su2double T_target = 500.0;

    nasa_model.SetTDState_T(T_target);
    const su2double H = nasa_model.GetEnthalpy();

    nasa_model.SetTDState_h(H);
    const su2double T_recovered = nasa_model.GetTemperature();

    CHECK(T_recovered == Approx(T_target).epsilon(1e-6));
  }

  SECTION("Multiple temperature points") {
    const su2double test_temps[] = {250.0, 500.0, 1000.0, 2000.0, 4000.0};

    for (const auto T_in : test_temps) {
      nasa_model.SetTDState_T(T_in);
      const su2double H = nasa_model.GetEnthalpy();

      nasa_model.SetTDState_h(H);
      const su2double T_out = nasa_model.GetTemperature();

      CHECK(T_out == Approx(T_in).epsilon(1e-6));
    }
  }
}

TEST_CASE("NASA polynomial boundary conditions", "[CIncIdealGasNASA]") {
  const su2double R_N2 = 296.8;
  const su2double P = 101325.0;
  const su2double T_ref = 298.15;

  std::stringstream ss;
  SetupConfigStream(ss);
  CConfig config(ss, SU2_COMPONENT::SU2_CFD, false);
  ConfigureNASA(config);

  CIncIdealGasNASA<7> nasa_model(R_N2, P, T_ref);
  nasa_model.SetCpModel(&config, T_ref);

  SECTION("Lower temperature bound (200K)") {
    nasa_model.SetTDState_T(200.0);
    const su2double cp = nasa_model.GetCp();

    CHECK(cp > 0.0);
    CHECK(cp < 2000.0);
  }

  SECTION("Upper temperature bound (6000K)") {
    nasa_model.SetTDState_T(6000.0);
    const su2double cp = nasa_model.GetCp();

    CHECK(cp > 0.0);
    CHECK(cp < 3000.0);
  }
}

TEST_CASE("NASA polynomial thermodynamic consistency", "[CIncIdealGasNASA]") {
  const su2double R_N2 = 296.8;
  const su2double P = 101325.0;
  const su2double T_ref = 298.15;

  std::stringstream ss;
  SetupConfigStream(ss);
  CConfig config(ss, SU2_COMPONENT::SU2_CFD, false);
  ConfigureNASA(config);

  CIncIdealGasNASA<7> nasa_model(R_N2, P, T_ref);
  nasa_model.SetCpModel(&config, T_ref);

  SECTION("Density from ideal gas law") {
    const su2double T = 300.0;
    nasa_model.SetTDState_T(T);

    const su2double rho = nasa_model.GetDensity();
    const su2double rho_expected = P / (R_N2 * T);

    CHECK(rho == Approx(rho_expected).epsilon(1e-10));
  }

  SECTION("Cv from Cp (incompressible: gamma = 1)") {
    nasa_model.SetTDState_T(500.0);

    const su2double cp = nasa_model.GetCp();
    const su2double cv = nasa_model.GetCv();

    CHECK(cv == Approx(cp).epsilon(1e-10));
  }
}
