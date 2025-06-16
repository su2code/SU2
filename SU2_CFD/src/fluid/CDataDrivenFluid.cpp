/*!
 * \file CDataDrivenFluid.cpp
 * \brief Source of the data-driven fluid model class
 * \author E.C.Bunschoten M.Mayer A.Capiello
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

#include "../../include/fluid/CDataDrivenFluid.hpp"
#if defined(HAVE_MLPCPP)
#include "../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif

CDataDrivenFluid::CDataDrivenFluid(const CConfig* config, bool display) : CFluidModel() {
  rank = SU2_MPI::GetRank();
  DataDrivenFluid_ParsedOptions datadriven_fluid_options = config->GetDataDrivenParsedOptions();

  Kind_DataDriven_Method = datadriven_fluid_options.interp_algorithm_type;

  varname_rho = "Density";
  varname_e = "Energy";

  /*--- Use physics-informed approach ---*/
  use_MLP_derivatives = datadriven_fluid_options.use_PINN;

  /*--- Retrieve initial density and static energy values from config ---*/
  val_custom_init_e = datadriven_fluid_options.rho_init_custom;
  val_custom_init_rho = datadriven_fluid_options.e_init_custom;
  custom_init_e = (val_custom_init_e != -1.0);
  custom_init_rho = (val_custom_init_rho != -1.0);
  

  /*--- Set up interpolation algorithm according to data-driven method. Currently only MLP's are supported. ---*/
  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      lookup_mlp = new MLPToolbox::CLookUp_ANN(datadriven_fluid_options.n_filenames, datadriven_fluid_options.datadriven_filenames);
      if ((rank == MASTER_NODE) && display) lookup_mlp->DisplayNetworkInfo();
#else
      SU2_MPI::Error("SU2 was not compiled with MLPCpp enabled (-Denable-mlpcpp=true).", CURRENT_FUNCTION);
#endif
      break;
    case ENUM_DATADRIVEN_METHOD::LUT:
      if (use_MLP_derivatives && (rank == MASTER_NODE) && display)
        cout << "Physics-informed approach currently only works with MLP-based tabulation." << endl;

      lookup_table = new CLookUpTable(datadriven_fluid_options.datadriven_filenames[0], varname_rho, varname_e);
      break;
    default:
      break;
  }

  /*--- Relaxation factor and tolerance for Newton solvers. ---*/
  Newton_Relaxation = datadriven_fluid_options.Newton_relaxation;
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
  size_t n_outputs, idx_s,idx_dsde_rho = 1, idx_dsdrho_e = 2, idx_d2sde2 = 3, idx_d2sdedrho = 4, idx_d2sdrho2 = 5;
  if (use_MLP_derivatives) {
    n_outputs = 1;
    idx_s = 0;

    outputs_rhoe.resize(n_outputs);
    output_names_rhoe.resize(n_outputs);
    output_names_rhoe[idx_s] = "s";
    outputs_rhoe[idx_s] = &Entropy;

    dsdrhoe.resize(n_outputs);
    d2sdrhoe2.resize(n_outputs);
    dsdrhoe[0].resize(2);
    dsdrhoe[0][idx_rho] = &dsdrho_e;
    dsdrhoe[0][idx_e] = &dsde_rho;

    d2sdrhoe2[0].resize(2);
    d2sdrhoe2[0][idx_rho].resize(2);
    d2sdrhoe2[0][idx_e].resize(2);
    d2sdrhoe2[0][idx_rho][idx_rho] = &d2sdrho2;
    d2sdrhoe2[0][idx_rho][idx_e] = &d2sdedrho;
    d2sdrhoe2[0][idx_e][idx_rho] = &d2sdedrho;
    d2sdrhoe2[0][idx_e][idx_e] = &d2sde2;
  } else {
    n_outputs = 6;
    idx_s = 0;
    idx_dsde_rho = 1, idx_dsdrho_e = 2, idx_d2sde2 = 3, idx_d2sdedrho = 4, idx_d2sdrho2 = 5;

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
  }
  

  /*--- Further preprocessing of input and output variables. ---*/
  if (Kind_DataDriven_Method == ENUM_DATADRIVEN_METHOD::MLP) {
/*--- Map MLP inputs to outputs. ---*/
#ifdef USE_MLPCPP
    iomap_rhoe = new MLPToolbox::CIOMap(input_names_rhoe, output_names_rhoe);
    lookup_mlp->PairVariableswithMLPs(*iomap_rhoe);
    MLP_inputs.resize(2);
#endif
  } else {
    /*--- Retrieve column indices of LUT output variables ---*/
    LUT_idx_s = lookup_table->GetIndexOfVar(output_names_rhoe[idx_s]);
    LUT_idx_dsdrho_e = lookup_table->GetIndexOfVar(output_names_rhoe[idx_dsdrho_e]);
    LUT_idx_dsde_rho = lookup_table->GetIndexOfVar(output_names_rhoe[idx_dsde_rho]);
    LUT_idx_d2sde2 = lookup_table->GetIndexOfVar(output_names_rhoe[idx_d2sde2]);
    LUT_idx_d2sdedrho= lookup_table->GetIndexOfVar(output_names_rhoe[idx_d2sdedrho]);
    LUT_idx_d2sdrho2 = lookup_table->GetIndexOfVar(output_names_rhoe[idx_d2sdrho2]);

    LUT_lookup_indices.push_back(LUT_idx_s);
    LUT_lookup_indices.push_back(LUT_idx_dsde_rho);
    LUT_lookup_indices.push_back(LUT_idx_dsdrho_e);
    LUT_lookup_indices.push_back(LUT_idx_d2sde2);
    LUT_lookup_indices.push_back(LUT_idx_d2sdedrho);
    LUT_lookup_indices.push_back(LUT_idx_d2sdrho2);
  }
}

void CDataDrivenFluid::SetTDState_rhoe(su2double rho, su2double e) {

  Density = max(min(rho, rho_max), rho_min);
  StaticEnergy = max(min(e, e_max), e_min);

  AD::StartPreacc();
  AD::SetPreaccIn(Density);
  AD::SetPreaccIn(StaticEnergy);

  /*--- Compute thermodynamic state based on density and energy. ---*/
  
  Density = max(min(rho, rho_max), rho_min);
  StaticEnergy = max(min(e, e_max), e_min);

  Evaluate_Dataset(Density, StaticEnergy);

  const su2double rho_2 = Density * Density;
  /*--- Compute primary flow variables. ---*/
  Temperature = pow(dsde_rho, -1);
  Pressure = -rho_2 * Temperature * dsdrho_e;
  Enthalpy = StaticEnergy + Pressure / Density;

  /*--- Compute secondary flow variables ---*/
  dTde_rho = -Temperature * Temperature * d2sde2;
  dTdrho_e = -Temperature * Temperature * d2sdedrho;

  /*--- Compute speed of sound. ---*/
  const su2double blue_term = (dsdrho_e * (2 - Density * Temperature * d2sdedrho) + Density * d2sdrho2);
  const su2double green_term = (-Temperature * d2sde2 * dsdrho_e + d2sdedrho);

  SoundSpeed2 = -Density * Temperature * (blue_term - Density * green_term * (dsdrho_e / dsde_rho));

  dPde_rho = -rho_2 * Temperature * (-Temperature * (d2sde2 * dsdrho_e) + d2sdedrho);
  dPdrho_e = - Density * Temperature * (dsdrho_e * (2 - Density * Temperature * d2sdedrho) + Density * d2sdrho2);

  /*--- Compute enthalpy and entropy derivatives required for Giles boundary conditions. ---*/
  dhdrho_e = -Pressure * (1 / rho_2) + dPdrho_e / Density;
  dhde_rho = 1 + dPde_rho / Density;

  /*--- Compute specific heat at constant volume and specific heat at constant pressure. ---*/
  Cv = 1 / dTde_rho;
  dhdrho_P = dhdrho_e - dhde_rho * (1 / dPde_rho) * dPdrho_e;
  dhdP_rho = dhde_rho * (1 / dPde_rho);
  dsdrho_P = dsdrho_e - dPdrho_e * (1 / dPde_rho) * dsde_rho;
  dsdP_rho = dsde_rho / dPde_rho;

  const su2double drhode_p = -dPde_rho/dPdrho_e;
  const su2double dTde_p = dTde_rho + dTdrho_e*drhode_p;
  const su2double dhde_p = dhde_rho + drhode_p*dhdrho_e;
  Cp = dhde_p / dTde_p;

  AD::SetPreaccOut(Temperature);
  AD::SetPreaccOut(SoundSpeed2);
  AD::SetPreaccOut(dPde_rho);
  AD::SetPreaccOut(dPdrho_e);
  AD::SetPreaccOut(dTde_rho);
  AD::SetPreaccOut(dTdrho_e);
  AD::SetPreaccOut(Pressure);
  AD::SetPreaccOut(Entropy);
  AD::SetPreaccOut(Cp);
  AD::SetPreaccOut(Cv);
  AD::EndPreacc();
}

void CDataDrivenFluid::SetTDState_PT(su2double P, su2double T) {

  /*--- Approximate density and static energy with ideal gas law. ---*/
  rho_start = P / (R_idealgas * T);
  e_start = Cv_idealgas * T;
  
  /*--- Run 2D Newton solver for pressure and temperature ---*/
  Run_Newton_Solver(P, T, Pressure, Temperature, dPdrho_e, dPde_rho, dTdrho_e, dTde_rho);
}

void CDataDrivenFluid::SetTDState_Prho(su2double P, su2double rho) {
  /*--- Computing static energy according to pressure and density. ---*/
  SetEnergy_Prho(P, rho);
}

void CDataDrivenFluid::SetEnergy_Prho(su2double P, su2double rho) {
  /*--- Run 1D Newton solver for pressure at constant density. ---*/
  Density = rho;

  /*--- Approximate static energy through ideal gas law. ---*/
  StaticEnergy = Cv_idealgas * (P / (R_idealgas * rho));

  Run_Newton_Solver(P, Pressure, StaticEnergy, dPde_rho);
}

void CDataDrivenFluid::SetTDState_rhoT(su2double rho, su2double T) {
  /*--- Run 1D Newton solver for temperature at constant density. ---*/
  Density = rho;

  /*--- Approximate static energy through ideal gas law. ---*/
  StaticEnergy = Cv_idealgas * T;

  Run_Newton_Solver(T, Temperature, StaticEnergy, dTde_rho);
}

void CDataDrivenFluid::SetTDState_hs(su2double h, su2double s) {
  /*--- Run 2D Newton solver for enthalpy and entropy. ---*/

  e_start = StaticEnergy;
  rho_start = Density;
  Run_Newton_Solver(h, s, Enthalpy, Entropy, dhdrho_e, dhde_rho, dsdrho_e, dsde_rho);
}

void CDataDrivenFluid::SetTDState_Ps(su2double P, su2double s) {
  /*--- Run 2D Newton solver for pressure and entropy ---*/

  e_start = StaticEnergy;
  rho_start = Density;
  Run_Newton_Solver(P, s, Pressure, Entropy, dPdrho_e, dPde_rho, dsdrho_e, dsde_rho);
}


void CDataDrivenFluid::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
  SetTDState_Prho(P, rho);

  dhdrho_P = dhdrho_e - dhde_rho * (1 / dPde_rho) * dPdrho_e;
  dhdP_rho = dhde_rho * (1 / dPde_rho);
  dsdrho_P = dsdrho_e - dPdrho_e * (1 / dPde_rho) * dsde_rho;
  dsdP_rho = dsde_rho / dPde_rho;
}


unsigned long CDataDrivenFluid::Predict_MLP(su2double rho, su2double e) {
  unsigned long exit_code = 0;
/*--- Evaluate MLP collection for the given values for density and energy. ---*/
#ifdef USE_MLPCPP
  MLP_inputs[idx_rho] = rho;
  MLP_inputs[idx_e] = e;
  if (use_MLP_derivatives){
    exit_code = lookup_mlp->PredictANN(iomap_rhoe, MLP_inputs, outputs_rhoe, &dsdrhoe, &d2sdrhoe2);
  } else {
    exit_code = lookup_mlp->PredictANN(iomap_rhoe, MLP_inputs, outputs_rhoe);
  }
#endif
  
  return exit_code;
}

unsigned long CDataDrivenFluid::Predict_LUT(su2double rho, su2double e) {
  bool inside = lookup_table->LookUp_XY(LUT_lookup_indices, outputs_rhoe, rho, e);
  if (inside)
    return 0;
  return 1;
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

void CDataDrivenFluid::Run_Newton_Solver(const su2double Y1_target, const su2double Y2_target, const su2double & Y1, const su2double & Y2, const su2double & dY1drho,
                         const su2double & dY1de, const  su2double & dY2drho, const su2double & dY2de) {
  /*--- 2D Newton solver, computing the density and internal energy values corresponding to Y1_target and Y2_target.
   * ---*/

  AD::StartPreacc();
  AD::SetPreaccIn(Y1_target);
  AD::SetPreaccIn(Y2_target);
  /*--- Setting initial values for density and energy. ---*/
  su2double rho = rho_start, e = e_start;

  bool converged = false;
  unsigned long Iter = 0;

  /*--- Initiating Newton solver ---*/
  while (!converged && (Iter < MaxIter_Newton)) {
    /*--- Determine thermodynamic state based on current density and energy. ---*/
    SetTDState_rhoe(rho, e);

    /*--- Determine residuals. ---*/
    const su2double delta_Y1 = Y1 - Y1_target;
    const su2double delta_Y2 = Y2 - Y2_target;

    /*--- Continue iterative process if residuals are outside tolerances. ---*/
    if ((abs(delta_Y1 / Y1) < Newton_Tolerance) && (abs(delta_Y2 / Y2) < Newton_Tolerance)) {
      converged = true;
    } else {
      /*--- Compute step size for density and energy. ---*/
      const su2double determinant = (dY1drho) * (dY2de) - (dY1de) * (dY2drho);
      const su2double delta_rho = (dY2de * delta_Y1 - dY1de * delta_Y2) / determinant;
      const su2double delta_e = (-dY2drho * delta_Y1 + dY1drho * delta_Y2) / determinant;

      /*--- Update density and energy values. ---*/
      rho -= Newton_Relaxation * delta_rho;
      e -= Newton_Relaxation * delta_e;
    }
    Iter++;
  }
  nIter_Newton = Iter;
  AD::SetPreaccOut(Density);
  AD::SetPreaccOut(StaticEnergy);
  AD::EndPreacc();
  /*--- Evaluation of final state. ---*/
  SetTDState_rhoe(Density, StaticEnergy);
}

void CDataDrivenFluid::Run_Newton_Solver(const su2double Y_target, const su2double & Y, su2double & X, const su2double & dYdX) {
  /*--- 1D Newton solver, computing the density or internal energy value corresponding to Y_target. ---*/

  bool converged = false;
  unsigned long Iter = 0;

  AD::StartPreacc();
  AD::SetPreaccIn(Y_target);
  AD::SetPreaccIn(X);
  /*--- Initiating Newton solver. ---*/
  while (!converged && (Iter < MaxIter_Newton)) {
    /*--- Determine thermodynamic state based on current density and energy. ---*/
    SetTDState_rhoe(Density, StaticEnergy);

    /*--- Determine residual ---*/
    const su2double delta_Y = Y_target - Y;

    /*--- Continue iterative process if residuals are outside tolerances. ---*/
    if (abs(delta_Y / Y) < Newton_Tolerance) {
      converged = true;
    } else {
      const su2double delta_X = delta_Y / dYdX;

      /*--- Update energy value ---*/
      X += Newton_Relaxation * delta_X;
    }
    Iter++;
  }
  AD::SetPreaccOut(Density);
  AD::SetPreaccOut(StaticEnergy);
  AD::EndPreacc();

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
    rho_average = lookup_mlp->GetInputOffset(iomap_rhoe, idx_rho);
    e_average = lookup_mlp->GetInputOffset(iomap_rhoe, idx_e);
#endif
    break;
  default:
    break;
  }
  /*--- Compute thermodynamic state from middle of data set. ---*/
  
  su2double rho_init = custom_init_rho ? val_custom_init_rho : rho_average;
  su2double e_init = custom_init_e ? val_custom_init_e : e_average;
  
  SetTDState_rhoe(rho_init, e_init);
  rho_median = rho_init;
  e_median = e_init;
  P_middle = Pressure;
  T_middle = Temperature;

  R_idealgas = P_middle / (rho_init * T_middle);
  Cv_idealgas = Cv;
  Cp_idealgas = Cp;
  
  gamma_idealgas = (R_idealgas / Cv_idealgas) + 1;
}