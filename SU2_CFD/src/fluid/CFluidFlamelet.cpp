/*!
 * \file CfluidFlamelet.cpp
 * \brief Main subroutines of CFluidFlamelet class
 * \author D. Mayer, T. Economon, N. Beishuizen
 * \version 7.5.1 "Blackbird"
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

#include "../include/fluid/CFluidFlamelet.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"
#if defined(HAVE_MLPCPP)
#include "../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif
#include <chrono>
CFluidFlamelet::CFluidFlamelet(CConfig* config, su2double value_pressure_operating, bool generate_manifold) : CFluidModel() {
  rank = SU2_MPI::GetRank();
  /* -- number of auxiliary species transport equations, e.g. 1=CO, 2=NOx  --- */
  n_user_scalars = config->GetNUserScalars();
  n_control_vars = config->GetNControlVars();
  n_scalars = config->GetNScalars();
  //PreferentialDiffusion = config->GetPreferentialDiffusion();

  if ((rank == MASTER_NODE) && generate_manifold) {
    cout << "Number of scalars:           " << n_scalars << endl;
    cout << "Number of user scalars:      " << n_user_scalars << endl;
    cout << "Number of control variables: " << n_control_vars << endl;
  }

  controlling_variables.resize(n_control_vars);

  n_datadriven_inputs = config->GetNDataDriven_Files();

  include_mixfrac = (n_control_vars > 2);

  controlling_variables[I_PROGVAR] = "ProgressVariable";
  controlling_variables[I_ENTH] = "EnthalpyTot";
  if (include_mixfrac) controlling_variables[I_MIXFRAC] = "MixtureFraction";

  table_scalar_names.resize(n_scalars);
  for (size_t i_CV = 0; i_CV < n_control_vars; i_CV++) table_scalar_names[i_CV] = controlling_variables[i_CV];
  /*--- auxiliary species transport equations---*/
  for (size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    table_scalar_names[n_control_vars + i_aux] = config->GetUserScalarName(i_aux);
  }

  manifold_format = config->GetKind_DataDriven_Method();
  if (generate_manifold) {
    switch (manifold_format) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      if (rank == MASTER_NODE) {
        cout << "*****************************************" << endl;
        cout << "***   initializing the lookup table   ***" << endl;
        cout << "*****************************************" << endl;
      }
      look_up_table = new CLookUpTable(config->GetDataDriven_FileNames()[0], table_scalar_names[I_PROGVAR],
                                       table_scalar_names[I_ENTH]);
      if (look_up_table->GetNDim() != n_control_vars)
        SU2_MPI::Error("Mismatch between table dimension and number of controlling variables.", CURRENT_FUNCTION);
      
      
      break;

    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      if ((rank == MASTER_NODE) && generate_manifold) {
        cout << "***********************************************" << endl;
        cout << "*** initializing the multi-layer perceptron ***" << endl;
        cout << "***********************************************" << endl;
      }
      look_up_ANN = new MLPToolbox::CLookUp_ANN(n_datadriven_inputs, config->GetDataDriven_FileNames());
      if ((rank == MASTER_NODE) && generate_manifold) look_up_ANN->DisplayNetworkInfo();
#else
      SU2_MPI::Error("SU2 was not compiled with MLPCpp enabled (-Denable-mlpcpp=true).", CURRENT_FUNCTION);
#endif
      break;
    default:
      break;
  }

  }

  config->SetLUTScalarNames(table_scalar_names);

  /*--- we currently only need 1 source term from the LUT for the progress variable
        and each auxiliary equations needs 2 source terms ---*/
  n_table_sources = 1 + 2 * n_user_scalars;

  table_source_names.resize(n_table_sources);
  table_sources.resize(n_table_sources);
  table_source_names[I_SRC_TOT_PROGVAR] = "ProdRateTot_PV";
  /*--- No source term for enthalpy ---*/

  /*--- For the auxiliary equations, we use a positive (production) and a negative (consumption) term:
        S_tot = S_PROD + S_CONS * Y ---*/

  for (size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    /*--- Order of the source terms: S_prod_1, S_cons_1, S_prod_2, S_cons_2, ...---*/
    table_source_names[1 + 2 * i_aux] = config->GetUserSourceName(2 * i_aux);
    table_source_names[1 + 2 * i_aux + 1] = config->GetUserSourceName(2 * i_aux + 1);
  }

  config->SetLUTSourceNames(table_source_names);

  n_lookups = config->GetNLookups();
  table_lookup_names.resize(n_lookups);
  for (int i_lookup = 0; i_lookup < n_lookups; ++i_lookup) {
    table_lookup_names[i_lookup] = config->GetLUTLookupName(i_lookup);
  }

  source_scalar.resize(n_scalars);
  lookup_scalar.resize(n_lookups);

  Pressure = value_pressure_operating;
  
  if (generate_manifold){
    PreprocessLookUp();
    config->SetPreferentialDiffusion(PreferentialDiffusion);
    if (rank == MASTER_NODE)
      cout << "Preferential Diffusion: " << (PreferentialDiffusion ? "Enabled" : "Disabled") << endl << endl;
    
  //Measure_Evaluation_Time();
  }
}

CFluidFlamelet::~CFluidFlamelet() {
  switch (manifold_format) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      if (look_up_table != nullptr)
        delete look_up_table;
      break;
    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      if (iomap_TD != nullptr){
        delete iomap_TD;
        if (PreferentialDiffusion) delete iomap_PD;
        delete iomap_Sources;
        delete iomap_LookUp;
        delete look_up_ANN;
      }
#endif
    default:
      break;
  }
}

/*--- do a lookup for the list of variables in table_lookup_names, for visualization purposes ---*/
unsigned long CFluidFlamelet::SetScalarLookups(const su2double* val_scalars) {
  su2double enth = val_scalars[I_ENTH];
  su2double prog = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names[I_ENTH];
  string name_prog = table_scalar_names[I_PROGVAR];

  val_controlling_vars[I_PROGVAR] = prog;
  val_controlling_vars[I_ENTH] = enth;

#ifdef USE_MLPCPP
  iomap_current = iomap_LookUp;
#endif
  /*--- add all quantities and their address to the look up vectors ---*/
  unsigned long exit_code = Evaluate_Dataset(varnames_LookUp, val_vars_LookUp);

  return exit_code;
}

/*--- set the source terms for the transport equations ---*/
unsigned long CFluidFlamelet::SetScalarSources(const su2double* val_scalars) {
  table_sources[0] = 0.0;

  val_controlling_vars[I_PROGVAR] = val_scalars[I_PROGVAR];
  val_controlling_vars[I_ENTH] = val_scalars[I_ENTH];
  if (PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_scalars[I_MIXFRAC];
#ifdef USE_MLPCPP
  iomap_current = iomap_Sources;
#endif
  /*--- add all quantities and their address to the look up vectors ---*/
  unsigned long exit_code = Evaluate_Dataset(varnames_Sources, val_vars_Sources);

  /*--- The source term for progress variable is always positive, we clip from below to makes sure. --- */
  source_scalar[I_PROGVAR] = max(0.0, table_sources[I_SRC_TOT_PROGVAR]-5.0);
  source_scalar[I_ENTH] = 0.0;
  if (PreferentialDiffusion) source_scalar[I_MIXFRAC] = 0.0;

  /*--- Source term for the auxiliary species transport equations ---*/
  for (size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    /*--- The source term for the auxiliary equations consists of a production term and a consumption term:
          S_TOT = S_PROD + S_CONS * Y ---*/
    su2double y_aux = val_scalars[n_control_vars + i_aux];
    su2double source_prod = table_sources[1 + 2 * i_aux];
    su2double source_cons = table_sources[1 + 2 * i_aux + 1];
    source_scalar[n_control_vars + i_aux] = source_prod + source_cons * y_aux;
  }

  return exit_code;
}

void CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars) {
  val_controlling_vars[I_PROGVAR] = val_scalars[I_PROGVAR];
  val_controlling_vars[I_ENTH] = val_scalars[I_ENTH];
  if (PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_scalars[I_MIXFRAC];
#ifdef USE_MLPCPP
  iomap_current = iomap_TD;
#endif
  /*--- add all quantities and their address to the look up vectors ---*/
  Evaluate_Dataset(varnames_TD, val_vars_TD);

  /*--- compute Cv from Cp and molar weight of the mixture (ideal gas) ---*/
  Cv = Cp - UNIVERSAL_GAS_CONSTANT / (molar_weight);

  Density = Pressure * (molar_weight / 1000.0) / (UNIVERSAL_GAS_CONSTANT * Temperature);

  // mass_diffusivity = Kt / (Density * Cp);
}

unsigned long CFluidFlamelet::SetPreferentialDiffusionScalars(su2double* val_scalars) {
  val_controlling_vars[I_PROGVAR] = val_scalars[I_PROGVAR];
  val_controlling_vars[I_ENTH] = val_scalars[I_ENTH];
  if (PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_scalars[I_MIXFRAC];
#ifdef USE_MLPCPP
  iomap_current = iomap_PD;
#endif
  /*--- add all quantities and their address to the look up vectors ---*/
  unsigned long exit_code = Evaluate_Dataset(varnames_PD, val_vars_PD);
  return exit_code;
}

unsigned long CFluidFlamelet::GetEnthFromTemp(su2double& val_enth, const su2double val_prog,
                                              const su2double val_mixfrac, const su2double val_temp,
                                              su2double initial_value) {
  /*--- convergence criterion for temperature in [K], high accuracy needed for restarts. ---*/
  su2double delta_temp_final = 0.001;
  su2double enth_iter = initial_value;
  su2double delta_enth;
  su2double delta_temp_iter = 1e10;

  unsigned long exit_code = 0, counter_limit = 50, counter = 0;

  bool converged = false;
#ifdef USE_MLPCPP
  iomap_current = iomap_TD;
#endif

  val_controlling_vars[I_PROGVAR] = val_prog;
  if (PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_mixfrac;

  while (!converged && (counter++ < counter_limit)) {
    val_controlling_vars[I_ENTH] = enth_iter;
    /*--- look up temperature and heat capacity ---*/
    Evaluate_Dataset(varnames_TD, val_vars_TD);
    /*--- calculate delta_temperature ---*/
    delta_temp_iter = val_temp - Temperature;

    if (abs(delta_temp_iter)/val_temp < 1e-6) {
      converged = true;
    } else {
      /* calculate delta_enthalpy following dh = cp * dT */
      delta_enth = Cp * delta_temp_iter;

      /*--- update enthalpy ---*/
      enth_iter += delta_enth;
    }
  }

  val_enth = enth_iter;

  if (counter >= counter_limit) {
    exit_code = 1;
  }
  return exit_code;
}

su2double CFluidFlamelet::GetBurntProgVar(su2double val_mixfrac) const {
  su2double pv_burnt;
  switch (manifold_format) {
    case ENUM_DATADRIVEN_METHOD::LUT:
    if (PreferentialDiffusion) {
      auto inclusion_levels = look_up_table->FindInclusionLevels(val_mixfrac);
      auto pv_bounds_lower = look_up_table->GetTableLimitsX(inclusion_levels.first);
      auto pv_bounds_upper = look_up_table->GetTableLimitsX(inclusion_levels.second);
      pv_burnt = 0.5*(pv_bounds_lower.second + pv_bounds_upper.second);
    } else {
      auto pv_bounds = look_up_table->GetTableLimitsX();
      pv_burnt = pv_bounds.second;
    }
    break;
    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
    auto pv_bounds = look_up_ANN->GetInputNorm(iomap_TD, I_PROGVAR);
    pv_burnt = pv_bounds.second;
#endif
    break;
  }
  return pv_burnt;
}
void CFluidFlamelet::PreprocessLookUp() {
  /*--- Set lookup names and variables for all relevant lookup processes in the fluid model. ---*/

  val_controlling_vars.resize(n_control_vars);

  /*--- Thermodynamic state variables ---*/
  size_t n_TD = 6;
  varnames_TD.resize(n_TD);
  val_vars_TD.resize(n_TD);

  /*--- The string in varnames_TD as it appears in the LUT file. ---*/
  varnames_TD[0] = "Temperature";
  val_vars_TD[0] = &Temperature;
  varnames_TD[1] = "mean_molar_weight";
  val_vars_TD[1] = &molar_weight;
  varnames_TD[2] = "Cp";
  val_vars_TD[2] = &Cp;
  varnames_TD[3] = "ViscosityDyn";
  val_vars_TD[3] = &Mu;
  varnames_TD[4] = "Conductivity";
  val_vars_TD[4] = &Kt;
  varnames_TD[5] = "DiffusionCoefficient";
  val_vars_TD[5] = &mass_diffusivity;

  /*--- Source term variables ---*/
  varnames_Sources.resize(n_table_sources);
  val_vars_Sources.resize(n_table_sources);

  for (size_t iSource = 0; iSource < n_table_sources; iSource++) {
    varnames_Sources[iSource] = table_source_names[iSource];
    val_vars_Sources[iSource] = &table_sources[iSource];
  }

  /*--- Passive lookups ---*/
  varnames_LookUp.resize(n_lookups);
  val_vars_LookUp.resize(n_lookups);

  for (size_t iLookup = 0; iLookup < n_lookups; iLookup++) {
    varnames_LookUp[iLookup] = table_lookup_names[iLookup];
    val_vars_LookUp[iLookup] = &lookup_scalar[iLookup];
  }

  varnames_CV.resize(n_control_vars);
  val_vars_CV.resize(n_control_vars);
  lookup_CV.resize(n_control_vars);
  for (auto iCV = 0u; iCV < n_control_vars; iCV++) {
    varnames_CV[iCV] = controlling_variables[iCV];
    val_vars_CV[iCV] = &lookup_CV[iCV];
  }

  varnames_PD.resize(FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS);
  val_vars_PD.resize(FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS);

  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR] = "Beta_ProgVar";
  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH_THERMAL] = "Beta_Enth_Thermal";
  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH] = "Beta_Enth";
  varnames_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC] = "Beta_MixFrac";

  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_PROGVAR] = &beta_progvar;
  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH_THERMAL] = &beta_enth_thermal;
  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_ENTH] = &beta_enth;
  val_vars_PD[FLAMELET_PREF_DIFF_SCALARS::I_BETA_MIXFRAC] = &beta_mixfrac;


  size_t n_betas {0};
  PreferentialDiffusion = false;
  switch (manifold_format)
  {
  case ENUM_DATADRIVEN_METHOD::LUT:
    PreferentialDiffusion = look_up_table->CheckForVariables(varnames_PD);
    break;
  case ENUM_DATADRIVEN_METHOD::MLP:
  #ifdef USE_MLPCPP
    for (auto iMLP=0u; iMLP < n_datadriven_inputs; iMLP++) {
      auto outputMap = look_up_ANN->FindVariableIndices(iMLP, varnames_PD, false);
      n_betas += outputMap.size();
    }
    PreferentialDiffusion = (n_betas == varnames_PD.size());
  #endif
    break;
  default:
    break;
  }

  if (PreferentialDiffusion && !include_mixfrac) 
    SU2_MPI::Error("Preferential diffusion can only be used with mixture fraction as a controlling variable.", CURRENT_FUNCTION);


  if (manifold_format == ENUM_DATADRIVEN_METHOD::MLP) {
#ifdef USE_MLPCPP
    iomap_TD = new MLPToolbox::CIOMap(controlling_variables, varnames_TD);
    look_up_ANN->PairVariableswithMLPs(*iomap_TD);
    iomap_Sources = new MLPToolbox::CIOMap(controlling_variables, varnames_Sources);
    look_up_ANN->PairVariableswithMLPs(*iomap_Sources);
    iomap_LookUp = new MLPToolbox::CIOMap(controlling_variables, varnames_LookUp);
    look_up_ANN->PairVariableswithMLPs(*iomap_LookUp);
    if (PreferentialDiffusion){
      iomap_PD = new MLPToolbox::CIOMap(controlling_variables, varnames_PD);
      look_up_ANN->PairVariableswithMLPs(*iomap_PD);
    }
#endif
  }
}

unsigned long CFluidFlamelet::Evaluate_Dataset(vector<string>& varnames, vector<su2double*>& val_vars) {
  unsigned long exit_code = 0;

  vector<string> LUT_varnames;
  vector<su2double*> LUT_val_vars;
  su2matrix<su2double*> gradient_refs;

  su2vector<su2double> CV_LUT;
  CV_LUT.resize(n_control_vars);
  switch (manifold_format) {
    case ENUM_DATADRIVEN_METHOD::LUT:
      LUT_varnames.resize(varnames.size());
      LUT_val_vars.resize(val_vars.size());
      for (auto iVar = 0u; iVar < varnames.size(); iVar++) {
        LUT_varnames[iVar] = varnames[iVar];
        LUT_val_vars[iVar] = val_vars[iVar];
      }
      if (PreferentialDiffusion) {
        exit_code = look_up_table->LookUp_XYZ(LUT_varnames, LUT_val_vars, val_controlling_vars[I_PROGVAR],
                                              val_controlling_vars[I_ENTH], val_controlling_vars[I_MIXFRAC]);
      } else
        exit_code = look_up_table->LookUp_XY(LUT_varnames, LUT_val_vars, val_controlling_vars[I_PROGVAR],
                                             val_controlling_vars[I_ENTH]);

      break;
    case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      exit_code = look_up_ANN->PredictANN(iomap_current, val_controlling_vars, val_vars);
#endif
      break;
    default:
      break;
  }
  return exit_code;
}


void CFluidFlamelet::Measure_Evaluation_Time() {
  // std::string line, word;
  // std::string filename_workers_models = "model_info.txt";
  // std::fstream fileworkers (filename_workers_models, std::ios::in);
  // getline(fileworkers, line);
  // int Worker_TD = stoi(line);
  // getline(fileworkers, line);
  // int Model_TD = stoi(line);
  // getline(fileworkers, line);
  // int Worker_PD = stoi(line);
  // getline(fileworkers, line);
  // int Model_PD = stoi(line);
  // getline(fileworkers, line);
  // int Worker_PPV = stoi(line);
  // getline(fileworkers, line);
  // int Model_PPV = stoi(line);
  // getline(fileworkers, line);
  // int Worker_Sources = stoi(line);
  // getline(fileworkers, line);
  // int Model_Sources = stoi(line);
  //   std::string storage_dir_TD = "/home/ecbunschoten/RACE/fgm_generation_TD/Architecture_Optimization/Architectures/";
  //   std::string storage_dir_PD = "/home/ecbunschoten/RACE/fgm_generation_PD/Architecture_Optimization/Architectures/";
  //   std::string storage_dir_PPV = "/home/ecbunschoten/RACE/fgm_generation/Architecture_Optimization/Architectures/";
  //   std::string storage_dir_Sources = "/home/ecbunschoten/RACE/flameletAI/ArchitectureOptimization/architectures/";

  //   std::string TD_filename = storage_dir_TD + "Worker_"+std::to_string(Worker_TD)+"/Model_"+std::to_string(Model_TD)+"/MLP_TD.mlp";
  //   std::string PD_filename = storage_dir_PD + "Worker_"+std::to_string(Worker_PD)+"/Model_"+std::to_string(Model_PD)+"/MLP_PD.mlp";
  //   std::string PPV_filename = storage_dir_PPV + "Worker_"+std::to_string(Worker_PPV)+"/Model_"+std::to_string(Model_PPV)+"/MLP_PPV.mlp";
  //   std::string Sources_filename = storage_dir_Sources + "Worker_"+std::to_string(Worker_Sources)+"/Model_"+std::to_string(Model_Sources)+"/MLP_Sources.mlp";
  //   size_t n_inputs =4;
  //   std::string input_filenames[] = {TD_filename, PD_filename, PPV_filename, Sources_filename};
    
  // std::vector<std::string> controlling_variables, query_variables;
  //   std::vector<su2double*> output_refs, TD_refs, PD_refs, PPV_refs, Sources_refs;
  //   std::vector<su2double> input_vars, output_vars, output_TD_vars, output_PD_vars, output_PPV_vars, output_Sources_vars;
  // controlling_variables.resize(3);
  //   input_vars.resize(3);
  //   controlling_variables[0] = "ProgressVariable";
  //   controlling_variables[1] = "EnthalpyTot";
  //   controlling_variables[2] = "MixtureFraction";

  //   std::vector<std::string> TD_vars = {"Temperature", "mean_molar_weight", "DiffusionCoefficient", "ViscosityDyn", "Cp", "Conductivity"};
  //   std::vector<std::string> PD_vars = {"Beta_ProgVar", "Beta_Enth_Thermal", "Beta_Enth", "Beta_MixFrac"};
  //   std::vector<std::string> PPV_vars = {"ProdRateTot_PV"};
  //   std::vector<std::string> Sources_vars = {"ProdRateTot_Y-NO", "ProdRateTot_Y-NO2"};

  //   su2double T_interp, M_interp, D_interp, k_interp, mu_interp, cp_interp, beta_pv_interp, beta_h1_interp, beta_h2_interp, beta_z_interp, ppv_interp, pno_interp, pno2_interp;
    
  //   TD_refs.resize(TD_vars.size());
  //   output_TD_vars.resize(TD_vars.size());
  //   TD_refs[0] = &T_interp;
  //   TD_refs[1] = &M_interp;
  //   TD_refs[2] = &D_interp;
  //   TD_refs[3] = &mu_interp;
  //   TD_refs[4] = &cp_interp;
  //   TD_refs[5] = &k_interp;

  //   PD_refs.resize(PD_vars.size());
  //   output_PD_vars.resize(PD_vars.size());
  //   PD_refs[0] = &beta_pv_interp;
  //   PD_refs[1] = &beta_h1_interp;
  //   PD_refs[2] = &beta_h2_interp;
  //   PD_refs[3] = &beta_z_interp;
    
  //   PPV_refs.resize(PPV_vars.size());
  //   output_PPV_vars.resize(PPV_vars.size());
  //   PPV_refs[0] = &ppv_interp;
    
  //   Sources_refs.resize(Sources_vars.size());
  //   output_Sources_vars.resize(Sources_vars.size());
  //   Sources_refs[0] = &pno_interp;
  //   Sources_refs[1] = &pno2_interp;

//if(manifold_format == ENUM_DATADRIVEN_METHOD::MLP){
// #ifdef USE_MLPCPP
//     MLPToolbox::CLookUp_ANN LookUp_MLP = MLPToolbox::CLookUp_ANN(n_inputs, input_filenames);
//     MLPToolbox::CIOMap IOMap_TD = MLPToolbox::CIOMap(controlling_variables, TD_vars);
//     MLPToolbox::CIOMap IOMap_PD = MLPToolbox::CIOMap(controlling_variables, PD_vars);
//     MLPToolbox::CIOMap IOMap_PPV = MLPToolbox::CIOMap(controlling_variables, PPV_vars);
//     MLPToolbox::CIOMap IOMap_Sources = MLPToolbox::CIOMap(controlling_variables, Sources_vars);
//     LookUp_MLP.PairVariableswithMLPs(IOMap_TD);
//     LookUp_MLP.PairVariableswithMLPs(IOMap_PD);
//     LookUp_MLP.PairVariableswithMLPs(IOMap_PPV);
//     LookUp_MLP.PairVariableswithMLPs(IOMap_Sources);
// #endif
//}

    std::string test_filename = "/home/ecbunschoten/RACE/fgm_generation/Flamelet_Generation/FlameletData/FlameletData_phi01-20_totalRR__test.csv";
    std::string line,word;
    std::vector<std::vector<std::string>> content;
    std::vector<std::vector<su2double>> test_data;
	std::vector<su2double> row;
    std::vector<su2double> input_vars;
    input_vars.resize(3);

    std::fstream file (test_filename, std::ios::in);
	if(file.is_open())
	{
        getline(file, line);
		while(getline(file, line))
		{
			row.clear();
 
			std::stringstream str(line);
 
			while(getline(str, word, ','))
				row.push_back(std::stod(word));
			test_data.push_back(row);
		}
	}

    su2double T_max, M_max, D_max, k_max, mu_max, cp_max, beta_pv_max, beta_h1_max, beta_h2_max, beta_z_max, ppv_max, pno_max, pno2_max;
    su2double T_min, M_min, D_min, k_min, mu_min, cp_min, beta_pv_min, beta_h1_min, beta_h2_min, beta_z_min, ppv_min, pno_min, pno2_min;
    size_t i=0;
    T_max = test_data[i][3];
    M_max = test_data[i][4];
    D_max = test_data[i][5];
    k_max = test_data[i][6];
    mu_max = test_data[i][7];
    cp_max =test_data[i][8];
    beta_pv_max = test_data[i][9];
    beta_h1_max = test_data[i][10];
    beta_h2_max = test_data[i][11];
    beta_z_max = test_data[i][12];
    ppv_max = test_data[i][13];
    pno_max = test_data[i][16];
    pno2_max = test_data[i][19];

    T_min = test_data[i][3];
    M_min = test_data[i][4];
    D_min = test_data[i][5];
    k_min = test_data[i][6];
    mu_min = test_data[i][7];
    cp_min =test_data[i][8];
    beta_pv_min = test_data[i][9];
    beta_h1_min = test_data[i][10];
    beta_h2_min = test_data[i][11];
    beta_z_min = test_data[i][12];
    ppv_min = test_data[i][13];
    pno_min = test_data[i][16];
    pno2_min = test_data[i][19];
    
    for(auto i=0u; i<test_data.size(); i++){
        T_max = std::max(T_max, test_data[i][3]);
        T_min = std::min(T_min, test_data[i][ 3]);
        M_max = std::max(M_max, test_data[i][ 4]);
        M_min = std::min(M_min, test_data[i][ 4]);
        D_max = std::max(D_max, test_data[i][ 5]);
        D_min = std::min(D_min, test_data[i][ 5]);
        k_max = std::max(k_max, test_data[i][ 6]);
        k_min = std::min(k_min, test_data[i][ 6]);
        mu_max = std::max(mu_max, test_data[i][ 7]);
        mu_min = std::min(mu_min, test_data[i][ 7]);
        cp_max = std::max(cp_max, test_data[i][ 8]);
        cp_min = std::min(cp_min, test_data[i][ 8]);
        beta_pv_max = std::max(beta_pv_max, test_data[i][ 9]);
        beta_pv_min = std::min(beta_pv_min, test_data[i][ 9]);
        beta_h1_max = std::max(beta_h1_max, test_data[i][ 10]);
        beta_h1_min = std::min(beta_h1_min, test_data[i][ 10]);
        beta_h2_max = std::max(beta_h2_max, test_data[i][ 11]);
        beta_h2_min = std::min(beta_h2_min, test_data[i][ 11]);
        beta_z_max = std::max(beta_z_max, test_data[i][ 12]);
        beta_z_min = std::min(beta_z_min, test_data[i][ 12]);
        ppv_max = std::max(ppv_max, test_data[i][ 13]);
        ppv_min = std::min(ppv_min, test_data[i][ 13]);
        pno_max = std::max(pno_max, test_data[i][ 16]);
        pno_min = std::min(pno_min, test_data[i][ 16]);
        pno2_max = std::max(pno2_max, test_data[i][ 19]);
        pno2_min = std::min(pno2_min, test_data[i][ 19]);
    }
    su2double delta_T = T_max - T_min, delta_M = M_max - M_min, delta_D = D_max - D_min, delta_k = k_max - k_min, delta_mu = mu_max - mu_min, delta_cp = cp_max - cp_min, 
            delta_beta_pv = beta_pv_max - beta_pv_min, delta_beta_h1 = beta_h1_max - beta_h1_min, delta_beta_h2 = beta_h2_max - beta_h2_min, delta_beta_z = beta_z_max - beta_z_min,
            delta_ppv = ppv_max - ppv_min, delta_pno = pno_max - pno_min, delta_pno2 = pno2_max - pno2_min;
    su2double t_average = 0;
    size_t n_eval = 3;
    su2double error_T, error_M, error_D, error_k, error_mu, error_cp, error_beta_pv, error_beta_h1, error_beta_h2, error_beta_z, error_ppv, error_pno, error_pno2;
        
    bool evaluate_accuracy = true;
    unsigned long exit_code;
    ofstream outfile_interp;
    for (auto j=0u; j<n_eval; j++){
        
        size_t i =0;
        su2double T_ref, M_ref, D_ref, k_ref, mu_ref, cp_ref, beta_pv_ref, beta_h1_ref, beta_h2_ref, beta_z_ref, ppv_ref, pno_ref, pno2_ref;
        int nr_of_microseconds = 0;
        if(evaluate_accuracy){
              outfile_interp.open("interpolated_solution.csv");
              outfile_interp << "ProgressVariable,EnthalpyTot,MixtureFraction";
              for(auto iVar=0u; iVar<varnames_TD.size(); iVar++)
                outfile_interp << "," << varnames_TD[iVar];
              for(auto iVar=0u; iVar<varnames_PD.size(); iVar++)
                outfile_interp << "," << varnames_PD[iVar];
              for(auto iVar=0u; iVar<varnames_Sources.size(); iVar++)
                outfile_interp << "," <<varnames_Sources[iVar];
              outfile_interp << endl;
            }
        while (i < test_data.size()) {
            input_vars[0] = test_data[i][0];
            input_vars[1] = test_data[i][1];
            input_vars[2] = test_data[i][2];
            auto start = std::chrono::high_resolution_clock::now();
            
            switch (manifold_format)
            {
            case ENUM_DATADRIVEN_METHOD::LUT:
              look_up_table->LookUp_XYZ(varnames_TD, val_vars_TD, input_vars[0], input_vars[1], input_vars[2]);
              look_up_table->LookUp_XYZ(varnames_PD, val_vars_PD, input_vars[0], input_vars[1], input_vars[2]);
              look_up_table->LookUp_XYZ(varnames_Sources, val_vars_Sources, input_vars[0], input_vars[1], input_vars[2]);
              break;
            case ENUM_DATADRIVEN_METHOD::MLP:
#ifdef USE_MLPCPP
      exit_code = look_up_ANN->PredictANN(iomap_TD, input_vars, val_vars_TD);
      exit_code = look_up_ANN->PredictANN(iomap_PD, input_vars, val_vars_PD);
      exit_code = look_up_ANN->PredictANN(iomap_Sources, input_vars, val_vars_Sources);
#endif
              break;
            default:
              break;
            }
            
            auto end = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            nr_of_microseconds += duration.count();

            if(evaluate_accuracy){
                T_ref = test_data[i][3];
                M_ref = test_data[i][4];
                D_ref = test_data[i][5];
                k_ref = test_data[i][6];
                mu_ref = test_data[i][7];
                cp_ref =test_data[i][8];
                beta_pv_ref = test_data[i][9];
                beta_h1_ref = test_data[i][10];
                beta_h2_ref = test_data[i][11];
                beta_z_ref = test_data[i][12];
                ppv_ref = test_data[i][13];
                pno_ref = test_data[i][16];
                pno2_ref = test_data[i][19];

                // error_T+=(T_interp - T_ref) * (T_interp - T_ref) / (test_data.size()*(delta_T * delta_T));
                // error_mu += (mu_interp - mu_ref) *(mu_interp - mu_ref) / (test_data.size()*delta_mu * delta_mu);
                // error_D += (D_interp - D_ref) * (D_interp - D_ref) / (test_data.size()*delta_D * delta_D);
                // error_k += (k_interp - k_ref) * (k_interp - k_ref) / (test_data.size()*delta_k * delta_k);
                // error_cp += (cp_interp - cp_ref) * (cp_interp - cp_ref) / (test_data.size()*delta_cp * delta_cp);
                // error_M += (M_interp - M_ref) * (M_interp - M_ref) / (test_data.size()*delta_M * delta_M);
                // error_beta_pv += (beta_pv_interp - beta_pv_ref) * (beta_pv_interp - beta_pv_ref) / (test_data.size()*delta_beta_pv * delta_beta_pv);
                // error_beta_h1 += (beta_h1_interp - beta_h1_ref) * (beta_h1_interp - beta_h1_ref) / (test_data.size()*delta_beta_h1 * delta_beta_h1);
                // error_beta_h2 += (beta_h2_interp - beta_h2_ref) * (beta_h2_interp - beta_h2_ref) / (test_data.size()*delta_beta_h2 * delta_beta_h2);
                // error_beta_z += (beta_z_interp - beta_z_ref) * (beta_z_interp - beta_z_ref) / (test_data.size()*delta_beta_z * delta_beta_z);
                // error_ppv += (ppv_interp - ppv_ref) * (ppv_interp - ppv_ref) / (test_data.size()*delta_ppv * delta_ppv);
                // error_pno += (pno_interp - pno_ref) * (pno_interp - pno_ref) / (test_data.size()*delta_pno * delta_pno);
                // error_pno2 += (pno2_interp - pno2_ref) * (pno2_interp - pno2_ref) / (test_data.size()*delta_pno2 * delta_pno2);

                error_T+=(Temperature - T_ref) * (Temperature - T_ref) / (test_data.size()*(delta_T * delta_T));
                error_mu += (Mu- mu_ref) *(Mu - mu_ref) / (test_data.size()*delta_mu * delta_mu);
                error_D += (mass_diffusivity - D_ref) * (mass_diffusivity - D_ref) / (test_data.size()*delta_D * delta_D);
                error_k += (Kt - k_ref) * (Kt - k_ref) / (test_data.size()*delta_k * delta_k);
                error_cp += (Cp - cp_ref) * (Cp - cp_ref) / (test_data.size()*delta_cp * delta_cp);
                error_M += (molar_weight - M_ref) * (molar_weight - M_ref) / (test_data.size()*delta_M * delta_M);
                error_beta_pv += (beta_progvar - beta_pv_ref) * (beta_progvar - beta_pv_ref) / (test_data.size()*delta_beta_pv * delta_beta_pv);
                error_beta_h1 += (beta_enth_thermal - beta_h1_ref) * (beta_enth_thermal - beta_h1_ref) / (test_data.size()*delta_beta_h1 * delta_beta_h1);
                error_beta_h2 += (beta_enth - beta_h2_ref) * (beta_enth - beta_h2_ref) / (test_data.size()*delta_beta_h2 * delta_beta_h2);
                error_beta_z += (beta_mixfrac - beta_z_ref) * (beta_mixfrac - beta_z_ref) / (test_data.size()*delta_beta_z * delta_beta_z);
                error_ppv += (source_scalar[I_PROGVAR] - ppv_ref) * (source_scalar[I_PROGVAR] - ppv_ref) / (test_data.size()*delta_ppv * delta_ppv);
                error_pno += (table_sources[1] - pno_ref) * (table_sources[1] - pno_ref) / (test_data.size()*delta_pno * delta_pno);
                error_pno2 += (table_sources[3] - pno2_ref) * (table_sources[3] - pno2_ref) / (test_data.size()*delta_pno2 * delta_pno2);

                outfile_interp << input_vars[0] << ", " << input_vars[1] << ", " << input_vars[2];
                for(auto iVar=0u; iVar<val_vars_TD.size(); iVar++)
                  outfile_interp << ", " << *val_vars_TD[iVar];
                for (auto iVar=0u; iVar < val_vars_PD.size(); iVar++)
                  outfile_interp << ", " << *val_vars_PD[iVar];
                for (auto iVar=0u; iVar < val_vars_Sources.size(); iVar++)
                  outfile_interp << ", " <<  *val_vars_Sources[iVar];
                outfile_interp << endl;
           }
            i++;
        }
        // if(evaluate_accuracy){
        //     error_T /= test_data.size();
        //     error_mu /= test_data.size();
        //     error_D /= test_data.size();
        //     error_k /= test_data.size();
        //     error_cp /= test_data.size();
        //     error_M /= test_data.size();
        //     error_beta_pv /= test_data.size();
        //     error_beta_h1 /= test_data.size();
        //     error_beta_h2 /= test_data.size();
        //     error_beta_z /= test_data.size();
        //     error_ppv /= test_data.size();
        //     error_pno /= test_data.size();
        //     error_pno2 /= test_data.size();
        // }
        
        evaluate_accuracy = false;
        t_average += nr_of_microseconds;
    }
    outfile_interp.close();
    cout << t_average / n_eval << ", "<< error_T << ", " << error_mu << ", "<< error_D << ", " << error_k << ", "
      << error_cp << ", "<< error_M << ", " << error_beta_pv << ", " << error_beta_h1 << ", "
      << error_beta_h2 << ", "<< error_beta_z << ", " << error_ppv << ", " << error_pno << ", " << error_pno2 << endl;
    ofstream outfile;
    outfile.open("t_errors.csv");
    outfile << t_average / n_eval << ", "<< error_T << ", " << error_mu << ", "<< error_D << ", " << error_k << ", "
      << error_cp << ", "<< error_M << ", " << error_beta_pv << ", " << error_beta_h1 << ", "
      << error_beta_h2 << ", "<< error_beta_z << ", " << error_ppv << ", " << error_pno << ", " << error_pno2 << std::endl;
    outfile.close();
    return;
}