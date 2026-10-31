/*!
 * \file CFluidCantera_tests.cpp
 * \brief Unit tests for Cantera fluid model.
 * \author C.Morales Ubal
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

#include <cmath>

#ifdef HAVE_CANTERA
#include "../../../SU2_CFD/include/fluid/CFluidCantera.hpp"
#include <cantera/core.h>

using namespace Cantera;

TEST_CASE("Fluid_Cantera", "[Multicomponent_flow]") {
  /*--- Cantera fluid model unit test cases. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig config("src/SU2/UnitTests/SU2_CFD/fluid/multicomponent_cantera.cfg", val_software, true);

  /*--- Create Cantera fluid model. ---*/
  CFluidCantera fluid(config.GetPressure_Thermodynamic(), &config);

  /*--- Set state using the scalars of the config file and the temperature. ---*/
  const su2double* scalar = config.GetSpecies_Init();
  const su2double Temperature = 300.0;
  fluid.SetTDState_T(Temperature, scalar);

  /*--- Check values for density and heat capacity. ---*/

  su2double density = fluid.GetDensity();
  su2double cp = fluid.GetCp();
  CHECK(density == Approx(0.924236));
  CHECK(cp == Approx(1277.91));
}
TEST_CASE("Fluid_Cantera_Combustion", "[Reacting_flow]") {
  /*--- Cantera fluid model unit test cases. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig config("src/SU2/UnitTests/SU2_CFD/fluid/multicomponent_cantera.cfg", val_software, true);

  /*--- Create Cantera fluid model. ---*/
  CFluidCantera fluid(config.GetPressure_Thermodynamic(), &config);

  /*--- Set state using the scalars of the config file and the temperature. ---*/
  const su2double* scalar = config.GetSpecies_Init();
  const su2double Temperature = 1900.0;
  fluid.SetTDState_T(Temperature, scalar);

  /*--- Compute chemical source terms. ---*/
  fluid.ComputeChemicalSourceTerm();

  /*--- Check values for source terms. ---*/

  su2double sourceTerm_H2 = fluid.GetChemicalSourceTerm(0);
  su2double sourceTerm_O2 = fluid.GetChemicalSourceTerm(1);
  CHECK(sourceTerm_H2 == Approx(-0.13633797171426));
  CHECK(sourceTerm_O2 == Approx(-2.16321066087493));
}
TEST_CASE("Fluid_Cantera_LargeMechanism", "[Multicomponent_flow]") {
  /*--- Three transported species in a mechanism with more species than the transported set. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig config("src/SU2/UnitTests/SU2_CFD/fluid/multicomponent_cantera_gri30.cfg", val_software, true);
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
  CHECK(fluid.GetTemperatureIterationFailed() == false);
  CHECK(fluid.GetTemperature() == Approx(1500.0).margin(1e-4));

  /*--- The result does not depend on the temperature guess, which is used only once. ---*/

  fluid.SetTemperatureGuess(1490.0);
  fluid.SetTDState_h(enthalpy, scalar);
  CHECK(fluid.GetTemperature() == Approx(1500.0).margin(1e-4));
  fluid.SetTemperatureGuess(5.0);
  fluid.SetTDState_h(enthalpy, scalar);
  CHECK(fluid.GetTemperatureIterationFailed() == false);
  CHECK(fluid.GetTemperature() == Approx(1500.0).margin(1e-4));

  /*--- An enthalpy outside the range of the thermodynamic data is flagged. ---*/

  fluid.SetTDState_h(1e12, scalar);
  CHECK(fluid.GetTemperatureIterationFailed() == true);
}
TEST_CASE("Fluid_Cantera_SourceJacobian", "[Reacting_flow]") {
  /*--- Diagonal chemical source Jacobian of the one-step methane mechanism. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig config("src/SU2/UnitTests/SU2_CFD/fluid/one_step_ch4.cfg", val_software, true);
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
  CHECK(fluid.GetChemicalSourceJacobian(2) * fluid.GetDensity() == Approx(slope_CH4).epsilon(0.1));
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
TEST_CASE("Fluid_Cantera_Compressible", "[Multicomponent_flow]") {
  /*--- State evaluations of the compressible solver, which set the composition separately. ---*/

  SU2_COMPONENT val_software = SU2_COMPONENT::SU2_CFD;
  CConfig config("src/SU2/UnitTests/SU2_CFD/fluid/multicomponent_cantera.cfg", val_software, true);
  CFluidCantera fluid(config.GetPressure_Thermodynamic(), &config);

  const su2double* scalar = config.GetSpecies_Init();
  const su2double pressure = 2.0e5, temperature = 1500.0;
  fluid.SetComposition(scalar);
  fluid.SetTDState_PT(pressure, temperature);
  const su2double density = fluid.GetDensity();
  const su2double energy = fluid.GetStaticEnergy();
  CHECK(fluid.GetEnthalpy() == Approx(energy + pressure / density));

  /*--- At the operating pressure the state equals the one of the incompressible solver. ---*/

  fluid.SetTDState_T(temperature, scalar);
  const su2double density_operating = fluid.GetDensity();
  fluid.SetTDState_PT(config.GetPressure_Thermodynamic(), temperature);
  CHECK(fluid.GetDensity() == Approx(density_operating));

  /*--- Round trips through the density-energy and pressure-density states. ---*/

  fluid.SetTDState_rhoe(density, energy);
  CHECK(fluid.GetTemperatureIterationFailed() == false);
  CHECK(fluid.GetTemperature() == Approx(temperature).margin(1e-4));
  CHECK(fluid.GetPressure() == Approx(pressure));

  fluid.SetTDState_Prho(pressure, density);
  CHECK(fluid.GetTemperature() == Approx(temperature));
  CHECK(fluid.GetStaticEnergy() == Approx(energy));

  fluid.SetTDState_rhoT(density, temperature);
  CHECK(fluid.GetPressure() == Approx(pressure));

  /*--- Round trips through the entropy-based states of the Riemann boundary condition. ---*/

  const su2double entropy = fluid.GetEntropy();
  const su2double enthalpy = fluid.GetEnthalpy();
  fluid.SetTDState_hs(enthalpy, entropy);
  CHECK(fluid.GetTemperatureIterationFailed() == false);
  CHECK(fluid.GetTemperature() == Approx(temperature).margin(1e-4));
  CHECK(fluid.GetPressure() == Approx(pressure));

  fluid.SetTDState_Ps(pressure, entropy);
  CHECK(fluid.GetTemperature() == Approx(temperature));

  /*--- Derivatives at frozen composition against central differences of the density-energy state. ---*/

  fluid.SetTDState_rhoe(density, energy);
  const su2double dPdrho_e = fluid.GetdPdrho_e();
  const su2double dPde_rho = fluid.GetdPde_rho();
  const su2double dTde_rho = fluid.GetdTde_rho();
  const su2double soundSpeed2 = fluid.GetSoundSpeed2();

  const su2double drho = 1e-4 * density, de = 1e-4 * fabs(energy) + 100.0;
  fluid.SetTDState_rhoe(density + drho, energy);
  const su2double P_rho_plus = fluid.GetPressure();
  fluid.SetTDState_rhoe(density - drho, energy);
  const su2double P_rho_minus = fluid.GetPressure();
  fluid.SetTDState_rhoe(density, energy + de);
  const su2double P_e_plus = fluid.GetPressure(), T_e_plus = fluid.GetTemperature();
  fluid.SetTDState_rhoe(density, energy - de);
  const su2double P_e_minus = fluid.GetPressure(), T_e_minus = fluid.GetTemperature();

  CHECK(dPdrho_e == Approx((P_rho_plus - P_rho_minus) / (2 * drho)).epsilon(1e-5));
  CHECK(dPde_rho == Approx((P_e_plus - P_e_minus) / (2 * de)).epsilon(1e-5));
  CHECK(dTde_rho == Approx((T_e_plus - T_e_minus) / (2 * de)).epsilon(1e-5));
  CHECK(soundSpeed2 == Approx(dPdrho_e + pressure / (density * density) * dPde_rho));

  /*--- The state follows the composition: moving H2 into N2 increases the density. ---*/

  const su2double leaner[2] = {0.5 * scalar[0], scalar[1]};
  fluid.SetComposition(leaner);
  fluid.SetTDState_PT(pressure, temperature);
  CHECK(fluid.GetDensity() > density);
}
#endif
