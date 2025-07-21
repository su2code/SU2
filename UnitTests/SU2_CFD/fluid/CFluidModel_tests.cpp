/*!
 * \file CFluidModel_tests.cpp
 * \brief Unit tests for the fluid model classes.
 * \author E.Bunschoten
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
#include <sstream>
#include "../../../SU2_CFD/include/fluid/CFluidModel.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/fluid/CDataDrivenFluid.hpp"

void FluidModelChecks(CFluidModel* fluid_model, const su2double val_p, const su2double val_T) {
  /*--- Check consistency of reverse look-up ---*/
  {
    fluid_model->SetTDState_PT(val_p, val_T);

    const su2double val_rho_fluidmodel = fluid_model->GetDensity();
    const su2double val_e_fluidmodel = fluid_model->GetStaticEnergy();

    fluid_model->SetTDState_rhoe(val_rho_fluidmodel, val_e_fluidmodel);
    CHECK(Approx(fluid_model->GetPressure()) == val_p);
    CHECK(Approx(fluid_model->GetTemperature()) == val_T);
  }
  /*--- Check internal consistency between primary and derived fluid properties ---*/
  fluid_model->SetTDState_PT(val_p, val_T);
  const su2double val_rho = fluid_model->GetDensity();
  const su2double val_e = fluid_model->GetStaticEnergy();

  const su2double dTdrho_e = fluid_model->GetdTdrho_e();
  const su2double dPdrho_e = fluid_model->GetdPdrho_e();
  const su2double dTde_rho = fluid_model->GetdTde_rho();
  const su2double dPde_rho = fluid_model->GetdPde_rho();

  su2double delta_rho = 1e-2, delta_e = 100.0;
  {
    fluid_model->SetTDState_rhoe(val_rho + delta_rho, val_e);
    const su2double T_plus = fluid_model->GetTemperature();
    const su2double p_plus = fluid_model->GetPressure();

    fluid_model->SetTDState_rhoe(val_rho - delta_rho, val_e);
    const su2double T_minus = fluid_model->GetTemperature();
    const su2double p_minus = fluid_model->GetPressure();
    const su2double dTdrho_e_FD = (T_plus - T_minus) / (2 * delta_rho);
    const su2double dPdrho_e_FD = (p_plus - p_minus) / (2 * delta_rho);

    CHECK(dTdrho_e == Approx(dTdrho_e_FD));
    CHECK(dPdrho_e == Approx(dPdrho_e_FD));
  }
  {
    fluid_model->SetTDState_rhoe(val_rho, val_e + delta_e);
    const su2double T_plus = fluid_model->GetTemperature();
    const su2double p_plus = fluid_model->GetPressure();

    fluid_model->SetTDState_rhoe(val_rho, val_e - delta_e);
    const su2double T_minus = fluid_model->GetTemperature();
    const su2double p_minus = fluid_model->GetPressure();
    const su2double dTde_rho_FD = (T_plus - T_minus) / (2 * delta_e);
    const su2double dPde_rho_FD = (p_plus - p_minus) / (2 * delta_e);

    CHECK(dTde_rho == Approx(dTde_rho_FD));
    CHECK(dPde_rho == Approx(dPde_rho_FD));
  }
}

TEST_CASE("Test case for ideal gas fluid model") {
  CIdealGas* fluid_model = new CIdealGas(1.4, 287.0);

  FluidModelChecks(fluid_model, 101325, 300.0);
  FluidModelChecks(fluid_model, 1e6, 600.0);

  delete fluid_model;
}

TEST_CASE("Test case for data-driven fluid model") {
  std::stringstream config_options;

  config_options << "SOLVER=RANS" << std::endl;
  config_options << "KIND_TURB_MODEL=SA" << std::endl;
  config_options << "SA_OPTIONS= NONE" << std::endl;
  config_options << "REYNOLDS_NUMBER=1e6" << std::endl;
  config_options << "FLUID_MODEL=DATADRIVEN_FLUID" << std::endl;
  config_options << "USE_PINN=YES" << std::endl;
  config_options << "INTERPOLATION_METHOD=MLP" << std::endl;
  config_options << "FILENAMES_INTERPOLATOR=(src/SU2/UnitTests/SU2_CFD/fluid/MLP_PINN.mlp)" << std::endl;
  config_options << "CONV_NUM_METHOD_FLOW=JST" << std::endl;

  /*--- Setup ---*/

  CConfig* config = new CConfig(config_options, SU2_COMPONENT::SU2_CFD, false);

  /*--- Define fluid model ---*/
  CDataDrivenFluid* fluid_model = new CDataDrivenFluid(config, false);

  /*--- Check fluid model consistency for several combinations of pressure-temperature. ---*/
  FluidModelChecks(fluid_model, 1.83e6, 523.0);
  FluidModelChecks(fluid_model, 2e5, 520.0);

  delete config;
  delete fluid_model;
}
