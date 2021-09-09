/*!
 * \file CFluidFlamelet.cpp
 * \brief Implementation of flamelet combustion model
 * \author D. Mayer, T. Economon, N. Beishuizen
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid/CFluidFlamelet.hpp"
#include "../../include/numerics/CLookUpTable.hpp"

CFluidFlamelet::CFluidFlamelet(CConfig *config, su2double value_pressure_operating) : CFluidModel() {

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    n_scalars = 4;
    config->SetNScalars(n_scalars);

    table_scalar_names.resize(n_scalars);
    table_scalar_names.at(I_ENTHALPY) = "Enthalpy";
    table_scalar_names.at(I_PROG_VAR) = "ProgVar";
    table_scalar_names.at(I_CO)       = "Y-CO";
    table_scalar_names.at(I_NOX)      = "Y-NOx";

    config->SetScalarNames(table_scalar_names);

    n_table_sources = 3;
    config->SetNTableSources(n_table_sources);

    table_source_names.resize(n_table_sources);
    table_source_names.at(I_SRC_TOT_PROG_VAR) = "ProdRateTot-PV";
    table_source_names.at(I_SRC_TOT_CO)       = "ProdRateTot-CO";
    table_source_names.at(I_SRC_TOT_NOX)      = "ProdRateTot-NOx";


    config->SetTableSourceNames(table_source_names);

    look_up_table = new CLookUpTable(config->GetFileNameLUT(),table_scalar_names.at(I_PROG_VAR),table_scalar_names.at(I_ENTHALPY));

    n_lookups = config->GetNLookups();
    table_lookup_names.resize(n_lookups);
    for (int i_lookup=0; i_lookup < n_lookups; ++i_lookup) {
      table_lookup_names.at(i_lookup) = config->GetLookupName(i_lookup);
    }

    source_scalar.resize(n_scalars);
    lookup_scalar.resize(n_lookups);

    Pressure = value_pressure_operating;
}

CFluidFlamelet::~CFluidFlamelet() {

  if (look_up_table!=NULL) delete   look_up_table;
}

unsigned long CFluidFlamelet::SetScalarLookups(su2double *val_scalars){

  su2double enth   = val_scalars[I_ENTHALPY];
  su2double prog   = val_scalars[I_PROG_VAR];

  string name_enth = table_scalar_names.at(I_ENTHALPY);
  string name_prog = table_scalar_names.at(I_PROG_VAR);

  /* perform table look ups */
  unsigned long exit_code = look_up_table->LookUp_ProgEnth(table_lookup_names, lookup_scalar, prog, enth, name_prog, name_enth);

  return exit_code;
}

unsigned long CFluidFlamelet::SetScalarSources(su2double *val_scalars){

  su2double*         table_sources = new su2double[n_table_sources];
  vector<string>     look_up_tags;
  vector<su2double*> look_up_data;

  su2double enth   = val_scalars[I_ENTHALPY];
  su2double prog   = val_scalars[I_PROG_VAR];

  string name_enth = table_scalar_names.at(I_ENTHALPY);
  string name_prog = table_scalar_names.at(I_PROG_VAR);

  for (int i_source=0; i_source < n_table_sources; ++i_source) {
    look_up_tags.push_back(table_source_names.at(i_source));
    look_up_data.push_back(&table_sources[i_source]);
  }

  /* perform table look ups */
  unsigned long exit_code = look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, prog, enth, name_prog, name_enth);

  source_scalar.at(I_ENTHALPY) = 0;
  source_scalar.at(I_PROG_VAR) = table_sources[I_SRC_TOT_PROG_VAR];
  source_scalar.at(I_CO)       = table_sources[I_SRC_TOT_CO];
  source_scalar.at(I_NOX)      = table_sources[I_SRC_TOT_NOX];


  /*--- we clip at a small positive value --- */
  if (source_scalar.at(I_PROG_VAR)<1.0e-6){
  // cout << "source term < 0!! (c,h)= "<<prog<<", "<<enth << endl;
  /* --- clip negative values of progress variable source term (this should not happen for a good lookup table) ---*/
    source_scalar.at(I_PROG_VAR) = 0.0;
  }

  delete [] table_sources;

  /*--- not used at the moment ---*/ 
  return exit_code;
}

unsigned long CFluidFlamelet::SetTDState_T(su2double val_temperature, su2double *val_scalars){

  su2double val_enth = val_scalars[I_ENTHALPY];
  su2double val_prog = val_scalars[I_PROG_VAR];

  string name_enth = table_scalar_names.at(I_ENTHALPY);
  string name_prog = table_scalar_names.at(I_PROG_VAR);

  vector<string>     look_up_tags;
  vector<su2double*> look_up_data;

  unsigned long exit_code;

  /* add all quantities and their address to the look up vectors */
  // nijso TODO: check if these exist in the lookup table
  look_up_tags.push_back("Temperature");
  look_up_data.push_back(&Temperature);
  look_up_tags.push_back("Density");
  look_up_data.push_back(&Density);
  look_up_tags.push_back("Cp");
  look_up_data.push_back(&Cp);
  look_up_tags.push_back("ViscosityDyn");
  look_up_data.push_back(&Mu);
  look_up_tags.push_back("Conductivity");
  look_up_data.push_back(&Kt);
  look_up_tags.push_back("Diffusivity");
  look_up_data.push_back(&mass_diffusivity);
  look_up_tags.push_back("HeatRelease");
  look_up_data.push_back(&source_energy);

  /* perform table look ups */
  exit_code = look_up_table->LookUp_ProgEnth(look_up_tags,look_up_data, val_prog, val_enth,name_prog,name_enth);

  // nijso: is Cv used somewhere?
  // we could check for the existence of molar_weight_mix in the lookup table, and else we just use gamma
  // default value is 1.4
  Cv = Cp/1.4;
  //Cv = Cp - UNIVERSAL_GAS_CONSTANT / (molar_weight_mix / 1000.);
  return exit_code;
}

unsigned long CFluidFlamelet::GetEnthFromTemp(su2double *val_enth, su2double val_prog, su2double val_temp){

  su2double          delta_temp_final = 0.01 ; /* convergence criterion for temperature in [K] */
  su2double          enth_iter        = 0. ;   /* in CH4/Air flames, 0 is usually a good initial value for the iteration */
  su2double          delta_enth;
  su2double          delta_temp_iter  = 1e10;
  unsigned long      exit_code        = 0;
  vector<string>     look_up_tags;
  vector<su2double*> look_up_data;
  int                counter_limit    = 50 ;
  string name_prog = table_scalar_names.at(I_PROG_VAR);
  string name_enth = table_scalar_names.at(I_ENTHALPY);

  /* set up look up vectors */
  su2double  temp_iter;
  look_up_tags.push_back("Temperature");
  look_up_data.push_back(&temp_iter);

  su2double  cp_iter;
  look_up_tags.push_back("Cp");
  look_up_data.push_back(&cp_iter);

  int counter = 0;
  while ( (abs(delta_temp_iter) > delta_temp_final) && (counter++ < counter_limit) ){

    /* look up temperature and heat capacity */
    look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, val_prog, enth_iter, name_prog, name_enth);

    /* calculate delta_temperature */
    delta_temp_iter = val_temp - temp_iter;

    /* calculate delta_enthalpy following dh = cp * dT */
    delta_enth      = cp_iter * delta_temp_iter;

    /* update enthalpy */
    enth_iter      += delta_enth;
  }

  /* set enthalpy value */
  *val_enth = enth_iter;

  if (counter >= counter_limit) {
    exit_code = 1;

  }

  return exit_code;
}
