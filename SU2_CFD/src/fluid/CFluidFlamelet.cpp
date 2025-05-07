/*!
 * \file CfluidFlamelet.cpp
 * \brief Main subroutines of CFluidFlamelet class
 * \author D. Mayer, T. Economon, N. Beishuizen, E. Bunschoten
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

#include <memory>
#include <string>
#include "../include/fluid/CFluidFlamelet.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"
#if defined(HAVE_MLPCPP)
#include "../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif

CFluidFlamelet::CFluidFlamelet(CConfig* config, su2double value_pressure_operating) : CFluidModel() {
  rank = SU2_MPI::GetRank();
  datadriven_fluid_options = config->GetDataDrivenParsedOptions();
  flamelet_options = config->GetFlameletParsedOptions();

  Kind_DataDriven_Method = datadriven_fluid_options.interp_algorithm_type;

  /* -- number of auxiliary species transport equations, e.g. 1=CO, 2=NOx  --- */
  n_user_scalars = flamelet_options.n_user_scalars;
  n_control_vars = flamelet_options.n_control_vars;
  include_mixture_fraction = (n_control_vars == 3);
  n_scalars = flamelet_options.n_scalars;

  if (rank == MASTER_NODE) {
    cout << "Number of scalars:           " << n_scalars << endl;
    cout << "Number of user scalars:      " << n_user_scalars << endl;
    cout << "Number of control variables: " << n_control_vars << endl;
  }
  scalars_vector.resize(n_scalars);

  table_scalar_names.resize(n_scalars);
  for (auto iCV = 0u; iCV < n_control_vars; iCV++) table_scalar_names[iCV] = flamelet_options.controlling_variable_names[iCV];

  /*--- auxiliary species transport equations---*/
  for (auto i_aux = 0u; i_aux < n_user_scalars; i_aux++) {
    table_scalar_names[n_control_vars + i_aux] = flamelet_options.user_scalar_names[i_aux];
  }

  controlling_variable_names.resize(n_control_vars);
  for (auto iCV = 0u; iCV < n_control_vars; iCV++)
    controlling_variable_names[iCV] =flamelet_options.controlling_variable_names[iCV];

  passive_specie_names.resize(n_user_scalars);
  for (auto i_aux = 0u; i_aux < n_user_scalars; i_aux++) passive_specie_names[i_aux] = flamelet_options.user_scalar_names[i_aux];

  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      if (rank == MASTER_NODE) {
        cout << "*****************************************" << endl;
        cout << "***   initializing the lookup table   ***" << endl;
        cout << "*****************************************" << endl;
      }
      look_up_table = new CLookUpTable(datadriven_fluid_options.datadriven_filenames[0], table_scalar_names[I_PROGVAR],
                                       table_scalar_names[I_ENTH]);
      break;
    default:
      if (rank == MASTER_NODE) {
        cout << "***********************************************" << endl;
        cout << "*** initializing the multi-layer perceptron ***" << endl;
        cout << "***********************************************" << endl;
      }
#ifdef USE_MLPCPP
      lookup_mlp = new MLPToolbox::CLookUp_ANN(datadriven_fluid_options.n_filenames, datadriven_fluid_options.datadriven_filenames);
      if ((rank == MASTER_NODE)) lookup_mlp->DisplayNetworkInfo();
#else
      SU2_MPI::Error("SU2 was not compiled with MLPCpp enabled (-Denable-mlpcpp=true).", CURRENT_FUNCTION);
#endif
      break;
  }

  Pressure = value_pressure_operating;

  PreprocessLookUp(config);

  if (rank == MASTER_NODE) {
    cout << "Preferential diffusion: " << (preferential_diffusion ? "Enabled" : "Disabled") << endl;
  }
}

CFluidFlamelet::~CFluidFlamelet() {
  if (Kind_DataDriven_Method == ENUM_DATADRIVEN_METHOD::LUT)
    delete look_up_table;
#ifdef USE_MLPCPP
  if (Kind_DataDriven_Method == ENUM_DATADRIVEN_METHOD::MLP) {
    delete iomap_TD;
    delete iomap_Sources;
    delete iomap_LookUp;
    delete lookup_mlp;
    if (preferential_diffusion) delete iomap_PD;
    }
#endif
}

void CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars) {
  for (auto iVar = 0u; iVar < n_scalars; iVar++) scalars_vector[iVar] = val_scalars[iVar];

  /*--- Add all quantities and their names to the look up vectors. ---*/
  EvaluateDataSet(scalars_vector, FLAMELET_LOOKUP_OPS::THERMO, val_vars_TD);

  Temperature = val_vars_TD[LOOKUP_TD::TEMPERATURE];
  Cp = val_vars_TD[LOOKUP_TD::HEATCAPACITY];
  Mu = val_vars_TD[LOOKUP_TD::VISCOSITY];
  Kt = val_vars_TD[LOOKUP_TD::CONDUCTIVITY];
  mass_diffusivity = val_vars_TD[LOOKUP_TD::DIFFUSIONCOEFFICIENT];
  switch (density_model) {
    case INC_DENSITYMODEL::FLAMELET:
      Density = val_vars_TD[LOOKUP_TD::MOLARWEIGHT];
      molar_weight = Pressure / (Density * UNIVERSAL_GAS_CONSTANT * Temperature);
      break;
    case INC_DENSITYMODEL::VARIABLE:
      molar_weight = val_vars_TD[LOOKUP_TD::MOLARWEIGHT];
      Density = (molar_weight / 1000) * Pressure / (UNIVERSAL_GAS_CONSTANT * Temperature);
      break;
    default:
      break;
  }
  /*--- Compute Cv from Cp and molar weight of the mixture (ideal gas). ---*/
  Cv = Cp - UNIVERSAL_GAS_CONSTANT / molar_weight;
}

void CFluidFlamelet::PreprocessLookUp(CConfig* config) {
  density_model = config->GetKind_DensityModel();
  /*--- Thermodynamic state variables and names. ---*/
  varnames_TD.resize(LOOKUP_TD::SIZE);
  val_vars_TD.resize(LOOKUP_TD::SIZE);

  /*--- The string in varnames_TD as it appears in the LUT file. ---*/
  varnames_TD[LOOKUP_TD::TEMPERATURE] = "Temperature";
  varnames_TD[LOOKUP_TD::HEATCAPACITY] = "Cp";
  varnames_TD[LOOKUP_TD::VISCOSITY] = "ViscosityDyn";
  varnames_TD[LOOKUP_TD::CONDUCTIVITY] = "Conductivity";
  varnames_TD[LOOKUP_TD::DIFFUSIONCOEFFICIENT] = "DiffusionCoefficient";

  /*--- In case of FLAMELET density model, the density is directly interpolated from the manifold.---*/
  switch (density_model) {
    case INC_DENSITYMODEL::FLAMELET:
      varnames_TD[LOOKUP_TD::MOLARWEIGHT] = "Density";
      break;
    case INC_DENSITYMODEL::VARIABLE:
      varnames_TD[LOOKUP_TD::MOLARWEIGHT] = "MolarWeightMix";
      break;
    default:
      break;
  }
  /*--- Scalar source term variables and names. ---*/
  size_t n_sources = n_control_vars + 2 * n_user_scalars;
  varnames_Sources.resize(n_sources);
  val_vars_Sources.resize(n_sources);
  for (auto iCV = 0u; iCV < n_control_vars; iCV++)
    varnames_Sources[iCV] = flamelet_options.cv_source_names[iCV];
  /*--- No source term for enthalpy ---*/

  /*--- For the auxiliary equations, we use a positive (production) and a negative (consumption) term:
        S_tot = S_PROD + S_CONS * Y ---*/

  for (size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    /*--- Order of the source terms: S_prod_1, S_cons_1, S_prod_2, S_cons_2, ...---*/
    varnames_Sources[n_control_vars + 2 * i_aux] = flamelet_options.user_source_names[2 * i_aux];
    varnames_Sources[n_control_vars + 2 * i_aux + 1] = flamelet_options.user_source_names[2 * i_aux + 1];
  }

  /*--- Passive look-up terms ---*/
  size_t n_lookups = flamelet_options.n_lookups;
  if (n_lookups == 0) {
    varnames_LookUp.resize(1);
    val_vars_LookUp.resize(1);
    varnames_LookUp[0] = "NULL";
  } else {
    varnames_LookUp.resize(n_lookups);
    val_vars_LookUp.resize(n_lookups);
    for (auto iLookup = 0u; iLookup < n_lookups; iLookup++) varnames_LookUp[iLookup] = flamelet_options.lookup_names[iLookup];
  }

  /*--- Preferential diffusion scalars ---*/
  varnames_PD.resize(FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS);
  val_vars_PD.resize(FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS);

  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR] = "Beta_ProgVar";
  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH_THERMAL] = "Beta_Enth_Thermal";
  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH] = "Beta_Enth";
  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC] = "Beta_MixFrac";

  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR] = beta_progvar;
  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH_THERMAL] = beta_enth_thermal;
  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH] = beta_enth;
  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC] = beta_mixfrac;

  preferential_diffusion = flamelet_options.preferential_diffusion;
  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      preferential_diffusion = look_up_table->CheckForVariables(varnames_PD);
      break;
    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      n_betas = 0;
      for (auto iMLP = 0u; iMLP < datadriven_fluid_options.n_filenames; iMLP++) {
        auto outputMap = lookup_mlp->FindVariableIndices(iMLP, varnames_PD, false);
        n_betas += outputMap.size();
      }
      preferential_diffusion = (n_betas == varnames_PD.size());
#endif
      break;
    default:
      break;
  }

  if (!preferential_diffusion && flamelet_options.preferential_diffusion)
    SU2_MPI::Error("Preferential diffusion scalars not included in flamelet manifold.", CURRENT_FUNCTION);

  if (Kind_DataDriven_Method == ENUM_DATADRIVEN_METHOD::MLP) {
#ifdef USE_MLPCPP
    iomap_TD = new MLPToolbox::CIOMap(controlling_variable_names, varnames_TD);
    iomap_Sources = new MLPToolbox::CIOMap(controlling_variable_names, varnames_Sources);
    iomap_LookUp = new MLPToolbox::CIOMap(controlling_variable_names, varnames_LookUp);
    lookup_mlp->PairVariableswithMLPs(*iomap_TD);
    lookup_mlp->PairVariableswithMLPs(*iomap_Sources);
    lookup_mlp->PairVariableswithMLPs(*iomap_LookUp);
    if (preferential_diffusion) {
      iomap_PD = new MLPToolbox::CIOMap(controlling_variable_names, varnames_PD);
      lookup_mlp->PairVariableswithMLPs(*iomap_PD);
    }
#endif
  } else {
    for (auto iVar=0u; iVar < varnames_TD.size(); iVar++) {
      LUT_idx_TD.push_back(look_up_table->GetIndexOfVar(varnames_TD[iVar]));
    }
    for (auto iVar=0u; iVar < varnames_Sources.size(); iVar++) {
      unsigned long LUT_idx;
      if (noSource(varnames_Sources[iVar])) {
        LUT_idx = look_up_table->GetNullIndex();
      } else {
        LUT_idx = look_up_table->GetIndexOfVar(varnames_Sources[iVar]);
      }
      LUT_idx_Sources.push_back(LUT_idx);
    }
    for (auto iVar=0u; iVar < varnames_LookUp.size(); iVar++) {
      unsigned long LUT_idx;
      if (noSource(varnames_LookUp[iVar]))
        LUT_idx = look_up_table->GetNullIndex();
      else 
        LUT_idx = look_up_table->GetIndexOfVar(varnames_LookUp[iVar]);
      LUT_idx_LookUp.push_back(LUT_idx);
    }
    if (preferential_diffusion) {
      for (auto iVar=0u; iVar < varnames_PD.size(); iVar++) {
        LUT_idx_PD.push_back(look_up_table->GetIndexOfVar(varnames_PD[iVar]));
      }
    }
  }
}

unsigned long CFluidFlamelet::EvaluateDataSet(const vector<su2double>& input_scalar, unsigned short lookup_type,
                                              vector<su2double>& output_refs) {
  AD::StartPreacc();
  for (auto iVar = 0u; iVar < input_scalar.size(); iVar++) AD::SetPreaccIn(input_scalar[iVar]);
  
  su2double val_enth = input_scalar[I_ENTH];
  su2double val_prog = input_scalar[I_PROGVAR];
  su2double val_mixfrac = include_mixture_fraction ? input_scalar[I_MIXFRAC] : 0.0;
  vector<su2double> val_vars;
  vector<su2double*> refs_vars;
  vector<unsigned long> LUT_idx;
  switch (lookup_type) {
    case FLAMELET_LOOKUP_OPS::THERMO:
      LUT_idx = LUT_idx_TD;
#ifdef USE_MLPCPP
      iomap_Current = iomap_TD;
#endif
      break;
    case FLAMELET_LOOKUP_OPS::PREFDIF:
      LUT_idx = LUT_idx_PD;
#ifdef USE_MLPCPP
      iomap_Current = iomap_PD;
#endif
      break;
    case FLAMELET_LOOKUP_OPS::SOURCES:
      LUT_idx = LUT_idx_Sources;
#ifdef USE_MLPCPP
      iomap_Current = iomap_Sources;
#endif
      break;
    case FLAMELET_LOOKUP_OPS::LOOKUP:
      LUT_idx = LUT_idx_LookUp;
#ifdef USE_MLPCPP
      iomap_Current = iomap_LookUp;
#endif
      break;
    default:
      break;
  }
  

  /*--- Add all quantities and their names to the look up vectors. ---*/
  bool inside;
  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      if (output_refs.size() != LUT_idx.size())
        SU2_MPI::Error(string("Output vector size incompatible with manifold lookup operation."), CURRENT_FUNCTION);
      if (include_mixture_fraction) {
        inside = look_up_table->LookUp_XYZ(LUT_idx, output_refs, val_prog, val_enth, val_mixfrac);
      } else {
        inside = look_up_table->LookUp_XY(LUT_idx, output_refs, val_prog, val_enth);
      }
      if (inside) extrapolation = 0;
      else extrapolation = 1;
      break;
    case ENUM_DATADRIVEN_METHOD::MLP:
      refs_vars.resize(output_refs.size());
      for (auto iVar = 0u; iVar < output_refs.size(); iVar++) refs_vars[iVar] = &output_refs[iVar];
#ifdef USE_MLPCPP
      extrapolation = lookup_mlp->PredictANN(iomap_Current, input_scalar, refs_vars);
#endif
      break;
    default:
      break;
  }
  for (auto iVar = 0u; iVar < output_refs.size(); iVar++) AD::SetPreaccOut(output_refs[iVar]);
  AD::EndPreacc();
  return extrapolation;
}
