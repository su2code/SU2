/*!
 * \file CFluidCantera_tests.cpp
 * \brief Unit tests for Cantera fluid model.
 * \author C.Morales Ubal
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#if defined(HAVE_CANTERA)
#define USE_CANTERA
#include "../../../SU2_CFD/include/fluid/CFluidCantera.hpp"
#include <cantera/core.h>

using namespace Cantera;
#endif

#ifdef USE_CANTERA
TEST_CASE("Fluid_Cantera", "[Multicomponent_flow]") {
  /*--- Cantera fluid model unit test cases. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig* config = new CConfig("multicomponent_cantera.cfg", val_software, true);
  CFluidCantera* auxFluidModel = nullptr;

  /*--- Create Cantera fluid model. ---*/
  su2double value_pressure_operating = config->GetPressure_Thermodynamic();
  auxFluidModel = new CFluidCantera(value_pressure_operating, config);

  /*--- Get scalar from config file and set temperature. ---*/
  const su2double* scalar = nullptr;
  scalar = config->GetSpecies_Init();
  const su2double Temperature = 300.0;
  /*--- Set state using temperature and scalar. ---*/
  auxFluidModel->SetTDState_T(Temperature, scalar);

  /*--- Check values for density and heat capacity. ---*/

  su2double density = auxFluidModel->GetDensity();
  su2double cp = auxFluidModel->GetCp();
  CHECK(density == Approx(0.924236));
  CHECK(cp == Approx(1277.91));
}
TEST_CASE("Fluid_Cantera_Combustion", "[Reacting_flow]") {
  /*--- Cantera fluid model unit test cases. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig* config = new CConfig("multicomponent_cantera.cfg", val_software, true);
  CFluidCantera* auxFluidModel = nullptr;

  /*--- Create Cantera fluid model. ---*/
  su2double value_pressure_operating = config->GetPressure_Thermodynamic();
  auxFluidModel = new CFluidCantera(value_pressure_operating, config);

  /*--- Get scalar from config file and set temperature. ---*/

  const su2double* scalar = nullptr;
  scalar = config->GetSpecies_Init();
  const su2double Temperature = 1900.0;
  /*--- Set state using temperature and scalar ---*/
  auxFluidModel->SetTDState_T(Temperature, scalar);

  /*--- Compute chemical source terms. ---*/
  auxFluidModel->ComputeChemicalSourceTerm();

  /*--- Check values for source terms. ---*/

  su2double sourceTerm_H2 = auxFluidModel->GetChemicalSourceTerm(0);
  su2double sourceTerm_O2 = auxFluidModel->GetChemicalSourceTerm(1);
  CHECK(sourceTerm_H2 == Approx(-0.13633797171426));
  CHECK(sourceTerm_O2 == Approx(-2.16321066087493));
}
TEST_CASE("Fluid_Cantera_LargeMechanism", "[Multicomponent_flow]") {
  /*--- Three transported species in a mechanism with more species than the transported set. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig config("multicomponent_cantera_gri30.cfg", val_software, true);
  CFluidCantera fluid(config.GetPressure_Thermodynamic(), &config);

  const su2double* scalar = config.GetSpecies_Init();
  fluid.SetTDState_T(300.0, scalar);
  const su2double density = fluid.GetDensity();
  CHECK(density > 1.0);
  CHECK(density < 1.4);
  CHECK(fluid.GetMassDiffusivity(0) > 0.0);
  CHECK(fluid.GetLaminarViscosity() > 0.0);

  /*--- Temperature recovered from the enthalpy of a given state. ---*/

  fluid.SetTDState_T(1500.0, scalar);
  const su2double enthalpy = fluid.GetEnthalpy();
  fluid.SetTDState_h(enthalpy, scalar);
  CHECK(fluid.GetStateFailed() == false);
  CHECK(fluid.GetTemperature() == Approx(1500.0).margin(1e-4));

  /*--- The result does not depend on the temperature guess, which is used only once. ---*/

  fluid.SetTemperatureGuess(1490.0);
  fluid.SetTDState_h(enthalpy, scalar);
  CHECK(fluid.GetTemperature() == Approx(1500.0).margin(1e-4));
  fluid.SetTemperatureGuess(5.0);
  fluid.SetTDState_h(enthalpy, scalar);
  CHECK(fluid.GetStateFailed() == false);
  CHECK(fluid.GetTemperature() == Approx(1500.0).margin(1e-4));

  /*--- An enthalpy outside the range of the thermodynamic data is flagged. ---*/

  fluid.SetTDState_h(1e12, scalar);
  CHECK(fluid.GetStateFailed() == true);
}
TEST_CASE("Fluid_Cantera_SourceJacobian", "[Reacting_flow]") {
  /*--- Diagonal chemical source Jacobian of the one-step methane mechanism. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig config("one_step_ch4.cfg", val_software, true);
  CFluidCantera fluid(config.GetPressure_Thermodynamic(), &config);

  const su2double temperature = 1800.0;
  su2double scalar[4] = {0.2226, 0.05, 0.0446, 0.05};
  fluid.SetTDState_T(temperature, scalar);
  fluid.ComputeChemicalSourceTerm();

  const su2double source_O2 = fluid.GetChemicalSourceTerm(0);
  const su2double source_CH4 = fluid.GetChemicalSourceTerm(2);
  CHECK(source_CH4 < 0.0);

  /*--- Reactants are consumed, so their Jacobian is negative; products have no sink. ---*/

  CHECK(fluid.GetChemicalSourceJacobian(0) < 0.0);
  CHECK(fluid.GetChemicalSourceJacobian(1) == 0.0);
  CHECK(fluid.GetChemicalSourceJacobian(2) < 0.0);
  CHECK(fluid.GetChemicalSourceJacobian(3) == 0.0);

  /*--- The rate is first order in CH4. The finite difference also contains the density change from moving N2
   *    into CH4, which the Jacobian (at fixed density) does not, hence the tolerance. ---*/

  const su2double delta = 1e-6;
  scalar[2] += delta;
  fluid.SetTDState_T(temperature, scalar);
  fluid.ComputeChemicalSourceTerm();
  const su2double slope_CH4 = (fluid.GetChemicalSourceTerm(2) - source_CH4) / delta;
  CHECK(fluid.GetChemicalSourceJacobian(2) == Approx(slope_CH4).epsilon(0.1));
  CHECK(source_O2 < 0.0);

  /*--- Below CANTERA_DC_MIN_TEMP (500 K by default) the chemistry is skipped. ---*/

  scalar[2] -= delta;
  fluid.SetTDState_T(400.0, scalar);
  fluid.ComputeChemicalSourceTerm();
  for (int iVar = 0; iVar < 4; iVar++) {
    CHECK(fluid.GetChemicalSourceTerm(iVar) == 0.0);
    CHECK(fluid.GetChemicalSourceJacobian(iVar) == 0.0);
  }
  CHECK(fluid.GetHeatRelease() == 0.0);

  fluid.SetTDState_T(600.0, scalar);
  fluid.ComputeChemicalSourceTerm();
  CHECK(fluid.GetChemicalSourceTerm(2) != 0.0);
}
#endif
