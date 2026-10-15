/*!
 * \file CfluidFlamelet.cpp
 * \brief Main subroutines of CFluidFlamelet class
 * \author D. Mayer, T. Economon, N. Beishuizen, E. Bunschoten
 * \version 8.5.0 "Harrier"
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

#include <memory>
#include <string>
#include "../../include/fluid/CFluidFlamelet.hpp"
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
  for (auto iCV = 0u; iCV < n_control_vars; iCV++)
    table_scalar_names[iCV] = flamelet_options.controlling_variable_names[iCV];

  /*--- auxiliary species transport equations---*/
  for (auto i_aux = 0u; i_aux < n_user_scalars; i_aux++) {
    table_scalar_names[n_control_vars + i_aux] = flamelet_options.user_scalar_names[i_aux];
  }

  controlling_variable_names.resize(n_control_vars);
  for (auto iCV = 0u; iCV < n_control_vars; iCV++)
    controlling_variable_names[iCV] = flamelet_options.controlling_variable_names[iCV];

  passive_specie_names.resize(n_user_scalars);
  for (auto i_aux = 0u; i_aux < n_user_scalars; i_aux++)
    passive_specie_names[i_aux] = flamelet_options.user_scalar_names[i_aux];

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
      lookup_mlp = new MLPToolbox::CLookUp_ANN(datadriven_fluid_options.n_filenames,
                                               datadriven_fluid_options.datadriven_filenames);
      if ((rank == MASTER_NODE)) lookup_mlp->DisplayNetworkInfo();
#else
      SU2_MPI::Error("SU2 was not compiled with MLPCpp enabled (-Denable-mlpcpp=true).", CURRENT_FUNCTION);
#endif
      break;
  }

  Pressure = value_pressure_operating;

  PreprocessLookUp(config);

  if (rank == MASTER_NODE) {
    if (preferential_diffusion) {
      const auto method = flamelet_options.pd_method;
      const char* method_str = (method == FLAMELET_PD_METHOD::BETA_CORRECTION) ? "BETA_CORRECTION"
                             : (method == FLAMELET_PD_METHOD::SOURCE_TERM)     ? "SOURCE_TERM"
                                                                                : "COMBINED";
      cout << "Preferential diffusion: Enabled (" << method_str << ")" << endl;
    } else {
      cout << "Preferential diffusion: Disabled" << endl;
    }
  }
}

CFluidFlamelet::~CFluidFlamelet() {
  if (Kind_DataDriven_Method == ENUM_DATADRIVEN_METHOD::LUT) delete look_up_table;
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

void CFluidFlamelet::SetTDState_h(su2double val_enthalpy, const su2double* val_scalars) {
  /*--- For the fluid flamelet model, the enthalpy (and the temperature) are passive scalars.
  val_scalars contains the enthalpy as the second variable: val_scalars= (Progress_variable, enthalpy,...).
  Consequently, the energy equation is not solved when the fluid flamelet model is used; instead the energy equation for
  enthalpy is solved in the species flamelet solver, and the enthalpy solution is overwritten in the CIncEulerSolver.
  Likewise, The temperature is retrieved from the look up table. Then, the thermodynamics state is fully determined with
  the val_scalars. This is the reason why enthalpy (or temperature) can be passed in either SetTDSTtate_T or
  SetTDState_h without affecting the solution.---*/
  SetTDState_T(val_enthalpy, val_scalars);
}

void CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars) {
  for (auto iVar = 0u; iVar < n_scalars; iVar++) scalars_vector[iVar] = val_scalars[iVar];

  /*--- Add all quantities and their names to the look up vectors. ---*/
  EvaluateDataSet(scalars_vector, FLAMELET_LOOKUP_OPS::THERMO, val_vars_TD);

  Enthalpy = scalars_vector[1];
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
  for (auto iCV = 0u; iCV < n_control_vars; iCV++) varnames_Sources[iCV] = flamelet_options.cv_source_names[iCV];
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
    for (auto iLookup = 0u; iLookup < n_lookups; iLookup++)
      varnames_LookUp[iLookup] = flamelet_options.lookup_names[iLookup];
  }

  /*--- Preferential diffusion scalars: only the terms required by the active PD method are looked up.
   *
   *    BETA_CORRECTION: [Beta_ProgVar, Beta_Enth_Thermal, Beta_Enth, Beta_MixFrac]
   *
   *    SOURCE_TERM (model B2 of Schepers & van Oijen, C&F 280 (2025) 114332), layout:
   *      [0, n_CV)                : S_PV, S_h, (S_Z)     — Eq. (16) closure sources of the
   *                                 NON-major species only, volumetric (rho-weighted) units:
   *                                 [kg/(m^3 s)] for PV and Z, [W/m^3] for h.
   *      [n_CV, n_CV+3)           : Y-H2, Y-H2O, Y-H     — major species mass fractions [-].
   *      [n_CV+3, n_CV+3+3*n_CV)  : D_{cv}_{sp}          — molecular flux coefficients
   *                                 D_{phi_k,i} = c_{phi_k,i} (rho*D_i - lambda/cp) of Eq. (14),
   *                                 cv in {PV, h, Z} (CV order), sp in {H2, H2O, H};
   *                                 units [kg/(m s)] for PV/Z rows, [W/m] for the h row
   *                                 (c_{h,i} = h_i(T) is folded in at generation time).
   *      [n_CV+3+3*n_CV, ...)     : DT_PV, DT_h, (DT_Z)  — combined Soret coefficients
   *                                 DT_{phi_k} = sum_sp c_{phi_k,sp} DT_sp, with DT_sp the
   *                                 Cantera thermal-diffusion coefficient [kg/(m s)] (species
   *                                 flux j_sp = -(rho D_sp grad Y_sp + DT_sp/T grad T)); units
   *                                 [kg/(m s)] for PV/Z, [W/m] for h. The 1/T factor is NOT
   *                                 tabulated; it is applied at runtime with the local CFD T.
   *    The runtime flux applied in CSpeciesFlameletSolver::Viscous_Residual is
   *      J_{phi_k} = - sum_sp D_{cv}_{sp} grad(Y_sp) - (DT_{cv}/T) grad(T),
   *    entering the transport equation as -div(J). ---*/
  const auto pd_method = flamelet_options.pd_method;
  const bool use_beta = (pd_method != FLAMELET_PD_METHOD::SOURCE_TERM);
  const bool use_src  = (pd_method != FLAMELET_PD_METHOD::BETA_CORRECTION);
  const unsigned n_beta_vars = use_beta ? FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS : 0u;
  n_pd_major_species = use_src ? flamelet_options.n_pd_major_species : 0u;
  const unsigned n_majors    = n_pd_major_species;
  const unsigned n_src_vars  = use_src ? (n_control_vars + n_majors + (n_majors + 1) * n_control_vars) : 0u;
  const unsigned n_pd_terms  = n_beta_vars + n_src_vars;
  varnames_PD.resize(n_pd_terms);
  val_vars_PD.resize(n_pd_terms, 0.0);

  if (use_beta) {
    varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR]      = "Beta_ProgVar";
    varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH_THERMAL]  = "Beta_Enth_Thermal";
    varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH]         = "Beta_Enth";
    varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC]      = "Beta_MixFrac";
  }

  if (use_src) {
    const auto* major_names = flamelet_options.pd_major_species_names;
    const auto* cv_names = flamelet_options.controlling_variable_names;

    for (auto iCV = 0u; iCV < n_control_vars; iCV++)
      varnames_PD[n_beta_vars + iCV] = "Res_" + cv_names[iCV];

    unsigned idx = n_beta_vars + n_control_vars;
    for (auto iSp = 0u; iSp < n_majors; iSp++) varnames_PD[idx++] = "Y-" + major_names[iSp];
    for (auto iCV = 0u; iCV < n_control_vars; iCV++)
      for (auto iSp = 0u; iSp < n_majors; iSp++)
        varnames_PD[idx++] = "D_" + cv_names[iCV] + "_" + major_names[iSp];
    for (auto iCV = 0u; iCV < n_control_vars; iCV++) varnames_PD[idx++] = "DT_" + cv_names[iCV];
  }

  preferential_diffusion = flamelet_options.preferential_diffusion;

  if (!preferential_diffusion && flamelet_options.preferential_diffusion)
    SU2_MPI::Error("Preferential diffusion scalars not included in flamelet manifold.", CURRENT_FUNCTION);

  if (Kind_DataDriven_Method == ENUM_DATADRIVEN_METHOD::MLP) {
#ifdef USE_MLPCPP
    iomap_TD = new MLPToolbox::CIOMap(controlling_variable_names, varnames_TD);
    iomap_Sources = new MLPToolbox::CIOMap(controlling_variable_names, varnames_Sources);
    iomap_LookUp = new MLPToolbox::CIOMap(controlling_variable_names, varnames_LookUp);
    lookup_mlp->PairVariableswithMLPs(*iomap_TD);
    lookup_mlp->PairVariableswithMLPs(*iomap_Sources);
    if (n_lookups > 1)
      lookup_mlp->PairVariableswithMLPs(*iomap_LookUp);
    if (preferential_diffusion) {
      iomap_PD = new MLPToolbox::CIOMap(controlling_variable_names, varnames_PD);
      lookup_mlp->PairVariableswithMLPs(*iomap_PD);
    }
#endif
  } else {
    for (auto iVar = 0u; iVar < varnames_TD.size(); iVar++) {
      LUT_idx_TD.push_back(look_up_table->GetIndexOfVar(varnames_TD[iVar]));
    }
    for (auto iVar = 0u; iVar < varnames_Sources.size(); iVar++) {
      unsigned long LUT_idx;
      if (noSource(varnames_Sources[iVar])) {
        LUT_idx = look_up_table->GetNullIndex();
      } else {
        LUT_idx = look_up_table->GetIndexOfVar(varnames_Sources[iVar]);
      }
      LUT_idx_Sources.push_back(LUT_idx);
    }
    for (auto iVar = 0u; iVar < varnames_LookUp.size(); iVar++) {
      unsigned long LUT_idx;
      if (noSource(varnames_LookUp[iVar]))
        LUT_idx = look_up_table->GetNullIndex();
      else
        LUT_idx = look_up_table->GetIndexOfVar(varnames_LookUp[iVar]);
      LUT_idx_LookUp.push_back(LUT_idx);
    }
    if (preferential_diffusion) {
      /*--- Check all preferential diffusion variables in one pass, so that a manifold missing
       several of the required columns (e.g. the Eq. (14) flux coefficients of the SOURCE_TERM
       method) reports the complete list instead of aborting on the first one. ---*/
      std::string missing_vars;
      for (const auto& name : varnames_PD)
        if (!look_up_table->CheckForVariables({name})) missing_vars += "\n  " + name;
      if (!missing_vars.empty())
        SU2_MPI::Error("The manifold is missing the following variables required by the active "
                       "PREFERENTIAL_DIFFUSION_METHOD (see CFluidFlamelet::PreprocessLookUp for "
                       "their definitions and units):" + missing_vars, CURRENT_FUNCTION);

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
  bool inside{true};
  switch (Kind_DataDriven_Method) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      if (output_refs.size() != LUT_idx.size())
        SU2_MPI::Error(string("Output vector size incompatible with manifold lookup operation."), CURRENT_FUNCTION);
      if (include_mixture_fraction) {
        inside = look_up_table->LookUp_XYZ(LUT_idx, output_refs, val_prog, val_enth, val_mixfrac);
      } else {
        inside = look_up_table->LookUp_XY(LUT_idx, output_refs, val_prog, val_enth);
      }
      break;
    case ENUM_DATADRIVEN_METHOD::MLP:
      refs_vars.resize(output_refs.size());
      for (auto iVar = 0u; iVar < output_refs.size(); iVar++) refs_vars[iVar] = &output_refs[iVar];
#ifdef USE_MLPCPP
      inside=lookup_mlp->Predict(*iomap_Current, input_scalar, refs_vars);
#endif
      break;
    default:
      break;
  }
  if (inside) extrapolation = 0;
      else extrapolation = 1;
  for (auto iVar = 0u; iVar < output_refs.size(); iVar++) AD::SetPreaccOut(output_refs[iVar]);
  AD::EndPreacc();
  return extrapolation;
}

void CFluidFlamelet::GetTableCVBounds(su2double& cv1_min, su2double& cv1_max,
                                      su2double& cv2_min, su2double& cv2_max,
                                      su2double& cv3_min, su2double& cv3_max) const {
  if (look_up_table == nullptr) {
    cv1_min = cv2_min = cv3_min = 0.0;
    cv1_max = cv2_max = cv3_max = 1.0;
    return;
  }
  auto lx0 = look_up_table->GetTableLimitsX(0);
  cv1_min = *lx0.first;  cv1_max = *lx0.second;
  auto ly0 = look_up_table->GetTableLimitsY(0);
  cv2_min = *ly0.first;  cv2_max = *ly0.second;
  for (unsigned long l = 1; l < look_up_table->GetNTableLevels(); ++l) {
    auto lx = look_up_table->GetTableLimitsX(l);
    if (*lx.first  < cv1_min) cv1_min = *lx.first;
    if (*lx.second > cv1_max) cv1_max = *lx.second;
    auto ly = look_up_table->GetTableLimitsY(l);
    if (*ly.first  < cv2_min) cv2_min = *ly.first;
    if (*ly.second > cv2_max) cv2_max = *ly.second;
  }
  auto lz = look_up_table->GetTableLimitsZ();
  cv3_min = lz.first;
  cv3_max = lz.second;
}

void CFluidFlamelet::ResetHullMissDistance() {
  if (look_up_table) look_up_table->ResetHullMissDistance();
}

su2double CFluidFlamelet::GetHullMissCV1Dev() const {
  return look_up_table ? look_up_table->GetHullMissCV1Dev() : 0.0;
}

su2double CFluidFlamelet::GetHullMissCV2Dev() const {
  return look_up_table ? look_up_table->GetHullMissCV2Dev() : 0.0;
}

su2double CFluidFlamelet::GetDistanceToNearestZLevel(su2double val_CV3) const {
  return look_up_table ? look_up_table->GetDistanceToNearestZLevel(val_CV3) : 0.0;
}
