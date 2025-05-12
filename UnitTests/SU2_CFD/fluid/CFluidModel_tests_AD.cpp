/*!
 * \file CFluidModel_tests.cpp
 * \brief AD unit tests for the fluid model classes.
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
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"

void FluidModelADChecks(CFluidModel* fluid_model, su2double val_p, su2double val_T) {
  /*--- Consistency tests for the fluid model secondary variables and AD integration. ---*/

  /*--- Extract fluid density and static energy. ---*/
  fluid_model->SetTDState_PT(val_p, val_T);
  su2double val_rho = fluid_model->GetDensity();
  su2double val_e = fluid_model->GetStaticEnergy();

  /*--- Test 1: check consistency of pressure derivatives. ---*/
  AD::Reset();
  AD::StartRecording();
  AD::RegisterInput(val_rho);
  AD::RegisterInput(val_e);

  fluid_model->SetTDState_rhoe(val_rho, val_e);
  su2double val_p_pred = fluid_model->GetPressure();
  AD::RegisterOutput(val_p_pred);
  AD::StopRecording();

  SU2_TYPE::SetDerivative(val_p_pred, 1.0);
  AD::ComputeAdjoint();

  /*--- Extract pressure derivatives from fluid model and AD. ---*/
  const su2double dpdrho_e_AD = SU2_TYPE::GetDerivative(val_rho);
  const su2double dpde_rho_AD = SU2_TYPE::GetDerivative(val_e);
  const su2double dpdrho_e_pred = fluid_model->GetdPdrho_e();
  const su2double dpde_rho_pred = fluid_model->GetdPde_rho();

  CHECK(SU2_TYPE::GetValue(dpdrho_e_AD) == Approx(SU2_TYPE::GetValue(dpdrho_e_pred)));
  CHECK(SU2_TYPE::GetValue(dpde_rho_AD) == Approx(SU2_TYPE::GetValue(dpde_rho_pred)));

  /*--- Test 2: check consistency of temperature derivatives. ---*/
  AD::Reset();
  AD::StartRecording();
  AD::RegisterInput(val_rho);
  AD::RegisterInput(val_e);

  fluid_model->SetTDState_rhoe(val_rho, val_e);
  su2double val_T_pred = fluid_model->GetTemperature();

  AD::RegisterOutput(val_T_pred);
  AD::StopRecording();

  SU2_TYPE::SetDerivative(val_T_pred, 1.0);
  AD::ComputeAdjoint();

  /*--- Extract temperature derivatives from fluid model and AD. ---*/
  const su2double dTdrho_e_AD = SU2_TYPE::GetDerivative(val_rho);
  const su2double dTde_rho_AD = SU2_TYPE::GetDerivative(val_e);
  const su2double dTdrho_e_pred = fluid_model->GetdTdrho_e();
  const su2double dTde_rho_pred = fluid_model->GetdTde_rho();

  CHECK(SU2_TYPE::GetValue(dTdrho_e_AD) == Approx(SU2_TYPE::GetValue(dTdrho_e_pred)));
  CHECK(SU2_TYPE::GetValue(dTde_rho_AD) == Approx(SU2_TYPE::GetValue(dTde_rho_pred)));

  /*--- Test 3: check consistency of specific heat at constant pressure. ---*/
  const su2double drhode_p = -dpde_rho_AD / dpdrho_e_AD;
  const su2double dTde_p = dTde_rho_AD + dTdrho_e_AD * drhode_p;
  const su2double dhde_rho = 1 + dpde_rho_AD / val_rho;
  const su2double dhdrho_e = -val_p * (1 / pow(val_rho, 2)) + dpdrho_e_AD / val_rho;
  const su2double dhde_p = dhde_rho + drhode_p * dhdrho_e;
  const su2double Cp_AD = dhde_p / dTde_p;
  const su2double Cp_pred = fluid_model->GetCp();

  CHECK(SU2_TYPE::GetValue(Cp_AD) == Approx(SU2_TYPE::GetValue(Cp_pred)));

  /*--- Test 4: check consistency of secondary variables for non-reflecting boundary conditions.---*/
  AD::Reset();
  AD::StartRecording();
  AD::RegisterInput(val_p);
  AD::RegisterInput(val_rho);

  fluid_model->ComputeDerivativeNRBC_Prho(val_p, val_rho);
  su2double val_s_pred = fluid_model->GetEntropy();

  AD::RegisterOutput(val_s_pred);
  AD::StopRecording();

  SU2_TYPE::SetDerivative(val_s_pred, 1.0);
  AD::ComputeAdjoint();
  const su2double dsdp_rho_AD = SU2_TYPE::GetDerivative(val_p);
  const su2double dsdrho_p_AD = SU2_TYPE::GetDerivative(val_rho);
  const su2double dsdp_rho_pred = fluid_model->GetdsdP_rho();
  const su2double dsdrho_p_pred = fluid_model->Getdsdrho_P();

  CHECK(SU2_TYPE::GetValue(dsdp_rho_pred) == Approx(SU2_TYPE::GetValue(dsdp_rho_AD)));
  CHECK(SU2_TYPE::GetValue(dsdrho_p_pred) == Approx(SU2_TYPE::GetValue(dsdrho_p_AD)));
}

TEST_CASE("AD test case for ideal gas fluid model", "[AD tests]") {
  CIdealGas* fluid_model = new CIdealGas(1.4, 287.0, true);
  FluidModelADChecks(fluid_model, 101325, 300.0);

  delete fluid_model;
}

TEST_CASE("AD test case for data-driven fluid model", "[AD tests]") {
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

  FluidModelADChecks(fluid_model, 3e5, 520.0);
  FluidModelADChecks(fluid_model, 1.83e6, 523.0);

  delete config;
  delete fluid_model;
}
