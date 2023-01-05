
/*!
 * \file CfluidFlamelet.cpp
 * \brief Main subroutines of CFluidFlamelet class
 * \author D. Mayer, T. Economon, N. Beishuizen
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../Common/include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"
#include "../../Common/include/toolboxes/multilayer_perceptron/CIOMap.hpp"

CFluidFlamelet::CFluidFlamelet(CConfig* config, su2double value_pressure_operating) : CFluidModel() {
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* -- number of auxiliary species transport equations: 1=CO, 2=NOx --- */
  n_user_scalars = config->GetNUserScalars();
  PreferentialDiffusion = config->GetPreferentialDiffusion();
  if(PreferentialDiffusion){
    n_CV = 3;
  }else
    n_CV = 2;

  n_scalars = n_CV + n_user_scalars;

  if (rank == MASTER_NODE) {
    cout << "n_scalars = " << n_scalars << endl;
    cout << "n_CV = " << n_CV << endl;
    cout << "n_user_scalars = " << n_user_scalars << endl;
  }

  controlling_variables.resize(n_CV);
  config->SetNControllingVars(n_CV);

  controlling_variables[I_PROGVAR] = "ProgressVariable";
  controlling_variables[I_ENTH] = "EnthalpyTot";
  if(PreferentialDiffusion) 
    controlling_variables[I_MIXFRAC] = "MixtureFraction";

  config->SetNScalars(n_scalars);

  table_scalar_names.resize(n_scalars);
  for(size_t i_CV=0; i_CV < n_CV; i_CV++)
    table_scalar_names[i_CV] = controlling_variables[i_CV];
  /*--- auxiliary species transport equations---*/
  for(size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    table_scalar_names[n_CV + i_aux] = config->GetUserScalarName(i_aux);
  }

  manifold_format = config->GetKind_DataDriven_Method();
  switch (manifold_format)
  {
  case ENUM_DATADRIVEN_METHOD::LUT:
    if (rank == MASTER_NODE) {
      cout << "*****************************************" << endl;
      cout << "***   initializing the lookup table   ***" << endl;
      cout << "*****************************************" << endl;
    }
    look_up_table =
      new CLookUpTable(config->GetDataDriven_FileNames()[0], table_scalar_names[I_PROGVAR], table_scalar_names[I_ENTH]);
    break;
  
  case ENUM_DATADRIVEN_METHOD::MLP:
    if (rank == MASTER_NODE) {
      cout << "***********************************************" << endl;
      cout << "*** initializing the multi-layer perceptron ***" << endl;
      cout << "***********************************************" << endl;
    }
    look_up_ANN = new MLPToolbox::CLookUp_ANN(config->GetNDataDriven_Files(), config->GetDataDriven_FileNames());
    break;
  default:
    break;
  }

  config->SetLUTScalarNames(table_scalar_names);

  /*--- we currently only need 1 source term from the LUT for the progress variable
        and each auxiliary equations needs 2 source terms ---*/
  n_table_sources = 1 + 2*n_user_scalars;
  config->SetNLUTSources(n_table_sources);

  table_source_names.resize(n_table_sources);
  table_sources.resize(n_table_sources);
  table_source_names[I_SRC_TOT_PROGVAR] = "ProdRateTot_PV";
  /*--- No source term for enthalpy ---*/

  /*--- For the auxiliary equations, we use a positive (production) and a negative (consumption) term:
        S_tot = S_PROD + S_CONS * Y ---*/ 

  for(size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    /*--- Order of the source terms: S_prod_1, S_cons_1, S_prod_2, S_cons_2, ...---*/
    table_source_names[1 + 2*i_aux]     = config->GetUserSourceName(2*i_aux);
    table_source_names[1 + 2*i_aux + 1] = config->GetUserSourceName(2*i_aux + 1);
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

  PreprocessLookUp();
}

CFluidFlamelet::~CFluidFlamelet() {
  switch (manifold_format)
  {
  case ENUM_DATADRIVEN_METHOD::LUT:
    delete look_up_table;
    break;
  case ENUM_DATADRIVEN_METHOD::MLP:
    delete iomap_TD;
    if(PreferentialDiffusion) delete iomap_PD;
    delete iomap_Sources;
    delete iomap_LookUp;
    delete look_up_ANN;
  default:
    break;
  }
}

/*--- do a lookup for the list of variables in table_lookup_names, for visualization purposes ---*/
unsigned long CFluidFlamelet::SetScalarLookups(su2double* val_scalars) {

  su2double enth = val_scalars[I_ENTH];
  su2double prog = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names[I_ENTH];
  string name_prog = table_scalar_names[I_PROGVAR];

  val_controlling_vars[I_PROGVAR] = prog;
  val_controlling_vars[I_ENTH] = enth;
  /*--- add all quantities and their address to the look up vectors ---*/
  unsigned long exit_code = Evaluate_Dataset(varnames_LookUp, val_vars_LookUp, iomap_LookUp);


  return exit_code;
}

/*--- set the source terms for the transport equations ---*/
unsigned long CFluidFlamelet::SetScalarSources(su2double* val_scalars) {

  table_sources[0] = 0.0;

  val_controlling_vars[I_PROGVAR] = val_scalars[I_PROGVAR];
  val_controlling_vars[I_ENTH] = val_scalars[I_ENTH];
  if(PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_scalars[I_MIXFRAC];
  /*--- add all quantities and their address to the look up vectors ---*/
  unsigned long exit_code = Evaluate_Dataset(varnames_Sources, val_vars_Sources, iomap_Sources);

  /*--- the source term for the progress variable is always positive by construction, but we clip from below  --- */
  source_scalar[I_PROGVAR] = max(EPS, table_sources[I_SRC_TOT_PROGVAR]);
  source_scalar[I_ENTH] = 0.0;
  if(PreferentialDiffusion) source_scalar[I_MIXFRAC] = 0.0;

  /*--- source term for the auxiliary species transport equations---*/
  for(size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    /*--- The source term for the auxiliary equations consists of a production term and a consumption term:
          S_TOT = S_PROD + S_CONS * Y ---*/
    su2double y_aux = val_scalars[n_CV + i_aux];     
    su2double source_prod = table_sources[1 + 2*i_aux];
    su2double source_cons = table_sources[1 + 2*i_aux + 1];
    source_scalar[n_CV + i_aux] = source_prod + source_cons * y_aux;
  }

  return exit_code;
}

void CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars) {

  val_controlling_vars[I_PROGVAR] = val_scalars[I_PROGVAR];
  val_controlling_vars[I_ENTH] = val_scalars[I_ENTH];
  if(PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_scalars[I_MIXFRAC];
  /*--- add all quantities and their address to the look up vectors ---*/
  Evaluate_Dataset(varnames_TD, val_vars_TD, iomap_TD);

  /*--- compute Cv from Cp and molar weight of the mixture (ideal gas) ---*/
  Cv = Cp - UNIVERSAL_GAS_CONSTANT / molar_weight;

  mass_diffusivity = Kt / Cp;
}

unsigned long CFluidFlamelet::SetPreferentialDiffusionScalars(su2double* val_scalars) {
  val_controlling_vars[I_PROGVAR] = val_scalars[I_PROGVAR];
  val_controlling_vars[I_ENTH] = val_scalars[I_ENTH];
  if(PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_scalars[I_MIXFRAC];
  /*--- add all quantities and their address to the look up vectors ---*/
  unsigned long exit_code = Evaluate_Dataset(varnames_PD, val_vars_PD, iomap_PD);
  return exit_code;
}

unsigned long CFluidFlamelet::GetEnthFromTemp(su2double* val_enth, su2double val_prog, su2double val_mixfrac, su2double val_temp, su2double initial_value) {

  string name_prog = table_scalar_names[I_PROGVAR];
  string name_enth = table_scalar_names[I_ENTH];

  su2double   delta_temp_final = 0.01, /* convergence criterion for temperature in [K] */
              relaxation = 0.5,        /* Newton solver relaxation factor. */
              enth_iter = initial_value, /* Initial enthalpy value, default stetting is zero. */
              delta_enth,              /* Enthalpy residual. */
              delta_temp_iter;         /* Temperature residual. */

  unsigned long exit_code = 0,
                counter_limit = 50,
                counter = 0;

  bool converged = false;

  val_controlling_vars[I_PROGVAR] = val_prog;
  if(PreferentialDiffusion) val_controlling_vars[I_MIXFRAC] = val_mixfrac;

  while (!converged && (counter++ < counter_limit)) {
    val_controlling_vars[I_ENTH] = enth_iter;
    /*--- look up temperature and heat capacity ---*/
    Evaluate_Dataset(varnames_TD, val_vars_TD, iomap_TD);
    /*--- calculate delta_temperature ---*/
    delta_temp_iter = val_temp - Temperature;
    if(abs(delta_temp_iter) < delta_temp_final){
      converged = true;
    }else{
      /*--- calculate delta_enthalpy following dh = cp * dT ---*/
      delta_enth = Cp * delta_temp_iter;

      /*--- update enthalpy ---*/
      enth_iter += relaxation * delta_enth;

      counter ++;
    }
  }

  /*--- set enthalpy value ---*/
  *val_enth = enth_iter;

  if (!converged) {
    exit_code = 1;
  }

  return exit_code;
}

void CFluidFlamelet::PreprocessLookUp() {
  /*--- Set lookup names and variables for all relevant lookup processes in the fluid model ---*/

  val_controlling_vars.resize(n_CV);

  /*--- Thermodynamic state variables ---*/
  varnames_TD.resize(5);
  val_vars_TD.resize(5);

  /*--- The string in varnames_TD is the actual string as it appears in the LUT file ---*/
  varnames_TD[0] = "Temperature";
  val_vars_TD[0] = &Temperature;
  varnames_TD[1] = "Density";
  val_vars_TD[1] = &Density;
  varnames_TD[2] = "Cp";
  val_vars_TD[2] = &Cp;
  varnames_TD[3] = "ViscosityDyn";
  val_vars_TD[3] = &Mu;
  varnames_TD[4] = "Conductivity";
  val_vars_TD[4] = &Kt;
  // varnames_TD[5] = "DiffusionCoefficient";
  // val_vars_TD[5] = &mass_diffusivity;
  // varnames_TD[6] = "MolarWeightMix";
  // val_vars_TD[6] = &molar_weight;

  /*--- Source term variables ---*/
  varnames_Sources.resize(n_table_sources);
  val_vars_Sources.resize(n_table_sources);
  
  for(size_t iSource=0; iSource < n_table_sources; iSource++){
    varnames_Sources[iSource] = table_source_names[iSource];
    val_vars_Sources[iSource] = &table_sources[iSource];
  }

  /*--- Passive lookups ---*/
  varnames_LookUp.resize(n_lookups);
  val_vars_LookUp.resize(n_lookups);

  for(size_t iLookup=0; iLookup < n_lookups; iLookup++){
    varnames_LookUp[iLookup] = table_lookup_names[iLookup];
    val_vars_LookUp[iLookup] = &lookup_scalar[iLookup];
  }
  
  if(PreferentialDiffusion){
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
  }
  if(manifold_format == ENUM_DATADRIVEN_METHOD::MLP){
    iomap_TD = new MLPToolbox::CIOMap(look_up_ANN, controlling_variables, varnames_TD);
    iomap_Sources = new MLPToolbox::CIOMap(look_up_ANN, controlling_variables, varnames_Sources);
    iomap_LookUp = new MLPToolbox::CIOMap(look_up_ANN, controlling_variables, varnames_LookUp);
    if(PreferentialDiffusion) iomap_PD = new MLPToolbox::CIOMap(look_up_ANN, controlling_variables, varnames_PD);
  }
}

unsigned long CFluidFlamelet::Evaluate_Dataset(su2vector<string>& varnames, su2vector<su2double*>& val_vars, MLPToolbox::CIOMap* iomap) {
  unsigned long exit_code = 0;

  vector<string> LUT_varnames;
  vector<su2double*> LUT_val_vars;
  switch (manifold_format)
  {
  case ENUM_DATADRIVEN_METHOD::LUT:
    LUT_varnames.resize(varnames.size());
    LUT_val_vars.resize(val_vars.size());
    for(auto iVar=0; iVar<varnames.size(); iVar++){
      LUT_varnames[iVar] = varnames[iVar];
      LUT_val_vars[iVar] = val_vars[iVar];
    }
    if(PreferentialDiffusion){
      exit_code = look_up_table->LookUp_XYZ(LUT_varnames, LUT_val_vars, val_controlling_vars[I_PROGVAR], val_controlling_vars[I_ENTH], val_controlling_vars[I_MIXFRAC]);
    }else
      exit_code = look_up_table->LookUp_XY(LUT_varnames, LUT_val_vars, val_controlling_vars[I_PROGVAR], val_controlling_vars[I_ENTH]);
    break;
  case ENUM_DATADRIVEN_METHOD::MLP:
    exit_code = look_up_ANN->Predict_ANN(iomap,val_controlling_vars, val_vars);
    break;
  default:
    break;
  }
  return exit_code;
}
