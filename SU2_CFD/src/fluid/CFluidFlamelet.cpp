
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

CFluidFlamelet::CFluidFlamelet(CConfig* config, su2double value_pressure_operating) : CFluidModel() {
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE) {
    cout << "n_scalars = " << n_scalars << endl;
  }

  /* -- number of auxiliary species transport equations: 1=CO, 2=NOx --- */

  const unsigned short n_auxiliary = 2;
  n_scalars = 2 + n_auxiliary;
  config->SetNScalars(n_scalars);

  if (rank == MASTER_NODE) {
    cout << "*****************************************" << endl;
    cout << "***   initializing the lookup table   ***" << endl;
    cout << "*****************************************" << endl;
  }

  table_scalar_names.resize(n_scalars);
  table_scalar_names.at(I_ENTH) = "EnthalpyTot";
  table_scalar_names.at(I_PROGVAR) = "ProgressVariable";
  /*--- auxiliary species transport equations---*/
  table_scalar_names.at(I_CO) = "Y-CO";
  table_scalar_names.at(I_NOX) = "Y-NOX";

  config->SetLUTScalarNames(table_scalar_names);

  /*--- we currently only need one source term from the LUT for the progress variable
        and one additional source term for the each of the auxiliary species transport equations ---*/
  n_table_sources = 3;
  config->SetNLUTSources(n_table_sources);

  table_source_names.resize(n_table_sources);
  table_source_names.at(I_SRC_TOT_PROGVAR) = "ProdRateTot-PV";
  /*--- No source term for enthalpy ---*/
  /*--- source terms for auxiliary species transport equations ---*/
  table_source_names.at(I_SRC_TOT_CO) = "ProdRateTot-CO";
  table_source_names.at(I_SRC_TOT_NOX) = "ProdRateTot-NOx";

  config->SetLUTSourceNames(table_source_names);

  look_up_table =
      new CLookUpTable(config->GetFileNameLUT(), table_scalar_names.at(I_PROGVAR), table_scalar_names.at(I_ENTH));

  n_lookups = config->GetNLookups();
  table_lookup_names.resize(n_lookups);
  for (int i_lookup = 0; i_lookup < n_lookups; ++i_lookup) {
    table_lookup_names.at(i_lookup) = config->GetLUTLookupName(i_lookup);
  }

  source_scalar.resize(n_scalars);
  lookup_scalar.resize(n_lookups);

  Pressure = value_pressure_operating;
}

CFluidFlamelet::~CFluidFlamelet() {
  if (look_up_table != NULL) delete look_up_table;
}

/* do a lookup for the list of variables in table_lookup_names, for visualization purposes */
unsigned long CFluidFlamelet::SetScalarLookups(su2double* val_scalars) {
  su2double enth = val_scalars[I_ENTH];
  su2double prog = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names.at(I_ENTH);
  string name_prog = table_scalar_names.at(I_PROGVAR);

  /* perform table look ups */
  unsigned long exit_code =
      look_up_table->LookUp_ProgEnth(table_lookup_names, lookup_scalar, prog, enth, name_prog, name_enth);

  return exit_code;
}

/* set the source terms for the transport equations*/
unsigned long CFluidFlamelet::SetScalarSources(su2double* val_scalars) {
  su2double* table_sources = new su2double[n_table_sources];
  vector<string> look_up_tags;
  vector<su2double*> look_up_data;

  table_sources[0] = 0.0;

  string name_enth = table_scalar_names.at(I_ENTH);
  string name_prog = table_scalar_names.at(I_PROGVAR);

  /*--- value for the progress variable and enthalpy ---*/
  su2double enth = val_scalars[I_ENTH];
  su2double prog = val_scalars[I_PROGVAR];
  /*--- value for the auxiliary species ---*/
  // su2double y_co = val_scalars[I_CO];
  // su2double y_nox = val_scalars[I_NOX];

  for (int i_source = 0; i_source < n_table_sources; ++i_source) {
    look_up_tags.push_back(table_source_names.at(i_source));
    look_up_data.push_back(&table_sources[i_source]);
  }

  /* perform table look ups */
  unsigned long exit_code =
      look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, prog, enth, name_prog, name_enth);

  /*--- source term for the progress variable and enthalpy ---*/
  source_scalar.at(I_PROGVAR) = table_sources[I_SRC_TOT_PROGVAR];
  source_scalar.at(I_ENTH) = 0;
  /*--- source term for the auxiliary species transport equations---*/
  source_scalar.at(I_CO) = table_sources[I_SRC_TOT_CO];
  source_scalar.at(I_NOX) = table_sources[I_SRC_TOT_NOX];

  /*--- the source term for the progress variable is always positive, but we clip it just to be sure --- */
  if (source_scalar.at(I_PROGVAR) < EPS) {
    source_scalar.at(I_PROGVAR) = 0.0;
  }

  delete[] table_sources;

  return exit_code;
}

void CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars) {
  su2double val_enth = val_scalars[I_ENTH];
  su2double val_prog = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names.at(I_ENTH);
  string name_prog = table_scalar_names.at(I_PROGVAR);

  vector<string> look_up_tags;
  vector<su2double*> look_up_data;

  /*--- add all quantities and their address to the look up vectors ---*/
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

  unsigned long exit_code = look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, val_prog, val_enth, name_prog, name_enth);

  /*--- ---*/
  if (exit_code != 0) {
    cout << "lookup was outside of the table" << endl;
  }

  // nijso: is Cv used somewhere?
  // according to cristopher, yes!
  // we could check for the existence of molar_weight_mix in the lookup table, and else we just use gamma
  // default value is 1.4
  Cv = Cp / 1.4;
  // Cv = Cp - UNIVERSAL_GAS_CONSTANT / (molar_weight_mix / 1000.);
}

unsigned long CFluidFlamelet::GetEnthFromTemp(su2double* val_enth, su2double val_prog, su2double val_temp) {
  su2double delta_temp_final = 0.01; /* convergence criterion for temperature in [K] */
  su2double enth_iter = 0.;          /* in CH4/Air flames, 0 is usually a good initial value for the iteration */
  su2double delta_enth;
  su2double delta_temp_iter = 1e10;
  unsigned long exit_code = 0;
  vector<string> look_up_tags;
  vector<su2double*> look_up_data;
  int counter_limit = 50;
  string name_prog = table_scalar_names.at(I_PROGVAR);
  string name_enth = table_scalar_names.at(I_ENTH);

  /* set up look up vectors */
  su2double temp_iter;
  look_up_tags.push_back("Temperature");
  look_up_data.push_back(&temp_iter);

  su2double cp_iter;
  look_up_tags.push_back("Cp");
  look_up_data.push_back(&cp_iter);

  int counter = 0;
  while ((abs(delta_temp_iter) > delta_temp_final) && (counter++ < counter_limit)) {
    /* look up temperature and heat capacity */
    look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, val_prog, enth_iter, name_prog, name_enth);

    /* calculate delta_temperature */
    delta_temp_iter = val_temp - temp_iter;

    /* calculate delta_enthalpy following dh = cp * dT */
    delta_enth = cp_iter * delta_temp_iter;

    /* update enthalpy */
    enth_iter += delta_enth;
  }

  /* set enthalpy value */
  *val_enth = enth_iter;

  if (counter >= counter_limit) {
    exit_code = 1;
  }

  return exit_code;
}
