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
#endif