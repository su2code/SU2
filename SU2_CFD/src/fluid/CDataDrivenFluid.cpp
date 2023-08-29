/*!
 * \file CDataDrivenFluid.cpp
 * \brief Source of the data-driven fluid model class
 * \author E.C.Bunschoten M.Mayer A.Capiello
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid/CDataDrivenFluid.hpp"
#if defined(HAVE_MLPCPP)
#include "../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif

CDataDrivenFluid::CDataDrivenFluid(const CConfig* config, bool display) : CFluidModel() {
  rank = SU2_MPI::GetRank();
  Kind_DataDriven_Method = config->GetKind_DataDriven_Method();

  varname_rho = "Density";
  varname_e = "Energy";

  /*--- Set up interpolation algorithm according to data-driven method. Currently only MLP's are supported. ---*/
  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      lookup_mlp = new MLPToolbox::CLookUp_ANN(config->GetNDataDriven_Files(), config->GetDataDriven_FileNames());
      if ((rank == MASTER_NODE) && display) lookup_mlp->DisplayNetworkInfo();
#else
      SU2_MPI::Error("SU2 was not compiled with MLPCpp enabled (-Denable-mlpcpp=true).", CURRENT_FUNCTION);
#endif
      break;
    case ENUM_DATADRIVEN_METHOD::LUT:
      lookup_table = new CLookUpTable(config->GetDataDriven_FileNames()[0], varname_rho, varname_e);
      break;
    default:
      break;
  }

  /*--- Relaxation factor and tolerance for Newton solvers. ---*/
  Newton_Relaxation = config->GetRelaxation_DataDriven();
  Newton_Tolerance = 1e-10;
  MaxIter_Newton = 50;

  /*--- Preprocessing of inputs and outputs for the interpolation method. ---*/
  MapInputs_to_Outputs();

  /*--- Compute approximate ideal gas properties ---*/
  ComputeIdealGasQuantities();
}

CDataDrivenFluid::~CDataDrivenFluid() {
  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      delete iomap_rhoe;
      delete lookup_mlp;
#endif
      break;
    case ENUM_DATADRIVEN_METHOD::LUT:
      delete lookup_table;
      break;
    default:
      break;
  }
}
void CDataDrivenFluid::MapInputs_to_Outputs() {
  /*--- Inputs of the data-driven method are density and internal energy. ---*/
  input_names_rhoe.resize(2);
  idx_rho = 0;
  idx_e = 1;
  input_names_rhoe[idx_rho] = varname_rho;
  input_names_rhoe[idx_e] = varname_e;

  /*--- Required outputs for the interpolation method are entropy and its partial derivatives with respect to energy and
   * density. ---*/
  size_t n_outputs = 6;
  size_t idx_s = 0, idx_dsde_rho = 1, idx_dsdrho_e = 2, idx_d2sde2 = 3, idx_d2sdedrho = 4, idx_d2sdrho2 = 5;

  outputs_rhoe.resize(n_outputs);
  output_names_rhoe.resize(n_outputs);
  output_names_rhoe[idx_s] = "s";
  outputs_rhoe[idx_s] = &Entropy;
  output_names_rhoe[idx_dsde_rho] = "dsde_rho";
  outputs_rhoe[idx_dsde_rho] = &dsde_rho;
  output_names_rhoe[idx_dsdrho_e] = "dsdrho_e";
  outputs_rhoe[idx_dsdrho_e] = &dsdrho_e;
  output_names_rhoe[idx_d2sde2] = "d2sde2";
  outputs_rhoe[idx_d2sde2] = &d2sde2;
  output_names_rhoe[idx_d2sdedrho] = "d2sdedrho";
  outputs_rhoe[idx_d2sdedrho] = &d2sdedrho;
  output_names_rhoe[idx_d2sdrho2] = "d2sdrho2";
  outputs_rhoe[idx_d2sdrho2] = &d2sdrho2;

  /*--- Further preprocessing of input and output variables. ---*/
  if (Kind_DataDriven_Method == ENUM_DATADRIVEN_METHOD::MLP) {
/*--- Map MLP inputs to outputs. ---*/
#ifdef USE_MLPCPP
    iomap_rhoe = new MLPToolbox::CIOMap(input_names_rhoe, output_names_rhoe);
    lookup_mlp->PairVariableswithMLPs(*iomap_rhoe);
    MLP_inputs.resize(2);
#endif
  }
}

void CDataDrivenFluid::SetTDState_rhoe(su2double rho, su2double e) {
  /*--- Compute thermodynamic state based on density and energy. ---*/
  Density = rho;
  StaticEnergy = e;

  /*--- Clip density and energy values to prevent extrapolation. ---*/
  Density = min(rho_max, max(rho_min, Density));
  StaticEnergy = min(e_max, max(e_min, StaticEnergy));

  Evaluate_Dataset(Density, StaticEnergy);

  /*--- Compute speed of sound. ---*/
  auto blue_term = (dsdrho_e * (2 - rho * pow(dsde_rho, -1) * d2sdedrho) + rho * d2sdrho2);
  auto green_term = (-pow(dsde_rho, -1) * d2sde2 * dsdrho_e + d2sdedrho);

  SoundSpeed2 = -rho * pow(dsde_rho, -1) * (blue_term - rho * green_term * (dsdrho_e / dsde_rho));

  /*--- Compute primary flow variables. ---*/
  Temperature = 1.0 / dsde_rho;
  Pressure = -pow(rho, 2) * Temperature * dsdrho_e;
  Density = rho;
  StaticEnergy = e;
  Enthalpy = e + Pressure / rho;

  /*--- Compute secondary flow variables ---*/
  dTde_rho = -pow(dsde_rho, -2) * d2sde2;
  dTdrho_e = -pow(dsde_rho, -2) * d2sdedrho;

  dPde_rho = -pow(rho, 2) * (dTde_rho * dsdrho_e + Temperature * d2sdedrho);
  dPdrho_e = -2 * rho * Temperature * dsdrho_e - pow(rho, 2) * (dTdrho_e * dsdrho_e + Temperature * d2sdrho2);

  /*--- Compute enthalpy and entropy derivatives required for Giles boundary conditions. ---*/
  dhdrho_e = -Pressure * pow(rho, -2) + dPdrho_e / rho;
  dhde_rho = 1 + dPde_rho / rho;

  dhdrho_P = dhdrho_e - dhde_rho * (1 / dPde_rho) * dPdrho_e;
  dhdP_rho = dhde_rho * (1 / dPde_rho);
  dsdrho_P = dsdrho_e - dPdrho_e * (1 / dPde_rho) * dsde_rho;
  dsdP_rho = dsde_rho / dPde_rho;
}

void CDataDrivenFluid::SetTDState_PT(su2double P, su2double T) {

  /*--- Approximate density and static energy with ideal gas law. ---*/
  rho_start = P / (R_idealgas * T);
  e_start = Cv_idealgas * T;
  
  /*--- Run 2D Newton solver for pressure and temperature ---*/
  Run_Newton_Solver(P, T, &Pressure, &Temperature, &dPdrho_e, &dPde_rho, &dTdrho_e, &dTde_rho);
}

void CDataDrivenFluid::SetTDState_Prho(su2double P, su2double rho) {
  /*--- Computing static energy according to pressure and density. ---*/
  SetEnergy_Prho(P, rho);
}

void CDataDrivenFluid::SetEnergy_Prho(su2double P, su2double rho) {
  /*--- Run 1D Newton solver for pressure at constant density. ---*/
  Density = rho;

  /*--- Approximate static energy through ideal gas law. ---*/
  su2double e_idealgas = Cv_idealgas * (P / (R_idealgas * rho));
  StaticEnergy = min(e_max, max(e_idealgas, e_min));

  Run_Newton_Solver(P, &Pressure, &StaticEnergy, &dPde_rho);
}

void CDataDrivenFluid::SetTDState_rhoT(su2double rho, su2double T) {
  /*--- Run 1D Newton solver for temperature at constant density. ---*/
  Density = rho;

  /*--- Approximate static energy through ideal gas law. ---*/
  StaticEnergy = Cv_idealgas * T;

  Run_Newton_Solver(T, &Temperature, &StaticEnergy, &dTde_rho);
}

void CDataDrivenFluid::SetTDState_hs(su2double h, su2double s) {
  /*--- Run 2D Newton solver for enthalpy and entropy. ---*/

  /*--- Approximate density and static energy through ideal gas law under isentropic assumption. ---*/
  su2double T_init = h / Cp_idealgas;
  su2double P_init = P_middle * pow(T_init / T_middle, gamma_idealgas/(gamma_idealgas - 1));

  e_start = h * Cv_idealgas / Cp_idealgas; 
  rho_start = P_init / (R_idealgas * T_init);
  Run_Newton_Solver(h, s, &Enthalpy, &Entropy, &dhdrho_e, &dhde_rho, &dsdrho_e, &dsde_rho);
}

void CDataDrivenFluid::SetTDState_Ps(su2double P, su2double s) {
  /*--- Run 2D Newton solver for pressure and entropy ---*/

  /*--- Approximate initial state through isentropic assumption and ideal gas law. ---*/
  su2double T_init = T_middle * pow(P / P_middle, (gamma_idealgas - 1)/gamma_idealgas);
  e_start = Cv_idealgas * T_init;
  rho_start = P / (R_idealgas * T_init);

  Run_Newton_Solver(P, s, &Pressure, &Entropy, &dPdrho_e, &dPde_rho, &dsdrho_e, &dsde_rho);
}

unsigned long CDataDrivenFluid::Predict_MLP(su2double rho, su2double e) {
  unsigned long exit_code = 0;
/*--- Evaluate MLP collection for the given values for density and energy. ---*/
#ifdef USE_MLPCPP
  MLP_inputs[idx_rho] = rho;
  MLP_inputs[idx_e] = e;
  exit_code = lookup_mlp->PredictANN(iomap_rhoe, MLP_inputs, outputs_rhoe);
#endif
  return exit_code;
}

unsigned long CDataDrivenFluid::Predict_LUT(su2double rho, su2double e) {
  unsigned long exit_code;
  std::vector<std::string> output_names_rhoe_LUT;
  std::vector<su2double*> outputs_LUT;
  output_names_rhoe_LUT.resize(output_names_rhoe.size());
  for (auto iOutput = 0u; iOutput < output_names_rhoe.size(); iOutput++) {
    output_names_rhoe_LUT[iOutput] = output_names_rhoe[iOutput];
  }

  outputs_LUT.resize(outputs_rhoe.size());
  for (auto iOutput = 0u; iOutput < outputs_rhoe.size(); iOutput++) {
    outputs_LUT[iOutput] = outputs_rhoe[iOutput];
  }

  exit_code = lookup_table->LookUp_XY(output_names_rhoe_LUT, outputs_LUT, rho, e);
  return exit_code;
}

void CDataDrivenFluid::Evaluate_Dataset(su2double rho, su2double e) {
  /*--- Evaluate dataset based on regression method. ---*/
  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      outside_dataset = Predict_LUT(rho, e);
      break;
    case ENUM_DATADRIVEN_METHOD::MLP:
      outside_dataset = Predict_MLP(rho, e);
      break;
    default:
      break;
  }
}

void CDataDrivenFluid::Run_Newton_Solver(su2double Y1_target, su2double Y2_target, su2double* Y1, su2double* Y2,
                                         su2double* dY1drho, su2double* dY1de, su2double* dY2drho, su2double* dY2de) {
  /*--- 2D Newton solver, computing the density and internal energy values corresponding to Y1_target and Y2_target.
   * ---*/

  /*--- Setting initial values for density and energy. ---*/
  su2double rho = rho_start, e = e_start;

  bool converged = false;
  unsigned long Iter = 0;

  su2double delta_Y1, delta_Y2, delta_rho, delta_e, determinant;

  /*--- Initiating Newton solver ---*/
  while (!converged && (Iter < MaxIter_Newton)) {
    /*--- Determine thermodynamic state based on current density and energy. ---*/
    SetTDState_rhoe(rho, e);

    /*--- Determine residuals. ---*/
    delta_Y1 = *Y1 - Y1_target;
    delta_Y2 = *Y2 - Y2_target;

    /*--- Continue iterative process if residuals are outside tolerances. ---*/
    if ((abs(delta_Y1 / *Y1) < Newton_Tolerance) && (abs(delta_Y2 / *Y2) < Newton_Tolerance)) {
      converged = true;
    } else {
      /*--- Compute step size for density and energy. ---*/
      determinant = (*dY1drho) * (*dY2de) - (*dY1de) * (*dY2drho);

      delta_rho = (*dY2de * delta_Y1 - *dY1de * delta_Y2) / determinant;
      delta_e = (-*dY2drho * delta_Y1 + *dY1drho * delta_Y2) / determinant;

      /*--- Update density and energy values. ---*/
      rho -= Newton_Relaxation * delta_rho;
      e -= Newton_Relaxation * delta_e;
    }
    Iter++;
  }
  nIter_Newton = Iter;

  /*--- Evaluation of final state. ---*/
  SetTDState_rhoe(rho, e);
}

void CDataDrivenFluid::Run_Newton_Solver(su2double Y_target, su2double* Y, su2double* X, su2double* dYdX) {
  /*--- 1D Newton solver, computing the density or internal energy value corresponding to Y_target. ---*/

  bool converged = false;
  unsigned long Iter = 0;

  su2double delta_Y, delta_X;

  /*--- Initiating Newton solver. ---*/
  while (!converged && (Iter < MaxIter_Newton)) {
    /*--- Determine thermodynamic state based on current density and energy. ---*/
    SetTDState_rhoe(Density, StaticEnergy);

    /*--- Determine residual ---*/
    delta_Y = Y_target - *Y;

    /*--- Continue iterative process if residuals are outside tolerances. ---*/
    if (abs(delta_Y / *Y) < Newton_Tolerance) {
      converged = true;
    } else {
      delta_X = delta_Y / *dYdX;

      /*--- Update energy value ---*/
      *X += Newton_Relaxation * delta_X;
    }
    Iter++;
  }

  /*--- Calculate thermodynamic state based on converged values for density and energy. ---*/
  SetTDState_rhoe(Density, StaticEnergy);

  nIter_Newton = Iter;
}

void CDataDrivenFluid::ComputeIdealGasQuantities() {
  /*--- Compute approximate ideal gas properties from the middle of the reference data set. These properties are used to approximate the initial condition of the Newton solvers using the ideal gas law. ---*/
  su2double rho_average = 1.0, e_average = 1.0;

  /*--- Obtain minimum and maximum density and static energy from data set. ---*/
  switch (Kind_DataDriven_Method)
  {
  case ENUM_DATADRIVEN_METHOD::LUT:
    rho_min = *lookup_table->GetTableLimitsX().first;
    e_min = *lookup_table->GetTableLimitsY().first;
    rho_max = *lookup_table->GetTableLimitsX().second;
    e_max = *lookup_table->GetTableLimitsY().second;
    rho_average = 0.5*(*lookup_table->GetTableLimitsX().first + *lookup_table->GetTableLimitsX().second);
    e_average = 0.5*(*lookup_table->GetTableLimitsY().first + *lookup_table->GetTableLimitsY().second);
    break;
  case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
    rho_min = lookup_mlp->GetInputNorm(iomap_rhoe, idx_rho).first;
    e_min = lookup_mlp->GetInputNorm(iomap_rhoe, idx_e).first;
    rho_max = lookup_mlp->GetInputNorm(iomap_rhoe, idx_rho).second;
    e_max = lookup_mlp->GetInputNorm(iomap_rhoe, idx_e).second;
    rho_average = 0.5*(lookup_mlp->GetInputNorm(iomap_rhoe, idx_rho).first + lookup_mlp->GetInputNorm(iomap_rhoe, idx_rho).second);
    e_average = 0.5*(lookup_mlp->GetInputNorm(iomap_rhoe, idx_e).first + lookup_mlp->GetInputNorm(iomap_rhoe, idx_e).second);
#endif
    break;
  default:
    break;
  }

  /*--- Compute thermodynamic state from middle of data set. ---*/
  SetTDState_rhoe(rho_average, e_average);
  P_middle = Pressure;
  T_middle = Temperature;

  R_idealgas = P_middle / (rho_average * T_middle);
  Cv_idealgas = e_average / T_middle;
  Cp_idealgas = Enthalpy / T_middle;
  gamma_idealgas = (R_idealgas / Cv_idealgas) + 1;
}