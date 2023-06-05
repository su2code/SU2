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

CFluidFlamelet::CFluidFlamelet(CConfig* config, su2double value_pressure_operating) : CFluidModel() {
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  /* -- number of auxiliary species transport equations, e.g. 1=CO, 2=NOx  --- */
  n_user_scalars = config->GetNUserScalars();
  n_control_vars = config->GetNControlVars();
  n_scalars = config->GetNScalars();

  if (rank == MASTER_NODE) {
    cout << "Number of scalars:           " << n_scalars << endl;
    cout << "Number of user scalars:      " << n_user_scalars << endl;
    cout << "Number of control variables: " << n_control_vars << endl;
  }

  if (rank == MASTER_NODE) {
    cout << "*****************************************" << endl;
    cout << "***   initializing the lookup table   ***" << endl;
    cout << "*****************************************" << endl;
  }

  table_scalar_names.resize(n_scalars);
  table_scalar_names[I_ENTH] = "EnthalpyTot";
  table_scalar_names[I_PROGVAR] = "ProgressVariable";
  /*--- auxiliary species transport equations---*/
  for (size_t i_aux = 0; i_aux < n_user_scalars; i_aux++) {
    table_scalar_names[n_control_vars + i_aux] = config->GetUserScalarName(i_aux);
  }

  look_up_table = new CLookUpTable(config->GetFileNameLUT(), table_scalar_names[I_PROGVAR], table_scalar_names[I_ENTH]);

  n_lookups = config->GetNLookups();
  table_lookup_names.resize(n_lookups);
  for (int i_lookup = 0; i_lookup < n_lookups; ++i_lookup) {
    table_lookup_names[i_lookup] = config->GetLUTLookupName(i_lookup);
  }

  Pressure = value_pressure_operating;

  /*--- Thermodynamic state variables and names. ---*/
  varnames_TD.resize(LOOKUP_TD::SIZE);
  val_vars_TD.resize(LOOKUP_TD::SIZE);

  /*--- The string in varnames_TD as it appears in the LUT file. ---*/
  varnames_TD[LOOKUP_TD::TEMPERATURE] =  "Temperature";
  varnames_TD[LOOKUP_TD::DENSITY] = "Density";
  varnames_TD[LOOKUP_TD::HEATCAPACITY] = "Cp";
  varnames_TD[LOOKUP_TD::VISCOSITY] = "ViscosityDyn";
  varnames_TD[LOOKUP_TD::CONDUCTIVITY] = "Conductivity";
  varnames_TD[LOOKUP_TD::DIFFUSIONCOEFFICIENT] = "DiffusionCoefficient";
  varnames_TD[LOOKUP_TD::MOLARWEIGHT] = "MolarWeightMix";

}

CFluidFlamelet::~CFluidFlamelet() { delete look_up_table; }

void CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars) {
  su2double val_enth = val_scalars[I_ENTH];
  su2double val_prog = val_scalars[I_PROGVAR];

  /*--- Add all quantities and their names to the look up vectors. ---*/
  look_up_table->LookUp_XY(varnames_TD, val_vars_TD, val_prog, val_enth);

  Temperature = val_vars_TD[LOOKUP_TD::TEMPERATURE];
  Density = val_vars_TD[LOOKUP_TD::DENSITY];
  Cp = val_vars_TD[LOOKUP_TD::HEATCAPACITY];
  Mu = val_vars_TD[LOOKUP_TD::VISCOSITY];
  Kt = val_vars_TD[LOOKUP_TD::CONDUCTIVITY];
  mass_diffusivity = val_vars_TD[LOOKUP_TD::DIFFUSIONCOEFFICIENT];
  molar_weight = val_vars_TD[LOOKUP_TD::MOLARWEIGHT];

  /*--- Compute Cv from Cp and molar weight of the mixture (ideal gas). ---*/
  Cv = Cp - UNIVERSAL_GAS_CONSTANT / molar_weight;
}

/* --- Total enthalpy is the transported variable, but we usually have temperature as a boundary condition,
       so we do a reverse lookup */
unsigned long CFluidFlamelet::GetEnthFromTemp(su2double& val_enth, const su2double val_prog, const su2double val_temp,
                                              const su2double initial_value) {
  /*--- convergence criterion for temperature in [K], high accuracy needed for restarts. ---*/
  su2double delta_temp_final = 0.001;
  su2double enth_iter = initial_value;
  su2double delta_enth;
  su2double delta_temp_iter = 1e10;
  unsigned long exit_code = 0;
  const int counter_limit = 1000;

  int counter = 0;
  while ((abs(delta_temp_iter) > delta_temp_final) && (counter++ < counter_limit)) {
    /*--- Add all quantities and their names to the look up vectors. ---*/
    look_up_table->LookUp_XY(varnames_TD, val_vars_TD, val_prog, enth_iter);
    Temperature = val_vars_TD[LOOKUP_TD::TEMPERATURE];
    Cp = val_vars_TD[LOOKUP_TD::HEATCAPACITY];

    delta_temp_iter = val_temp - Temperature;

    delta_enth = Cp * delta_temp_iter;

    enth_iter += delta_enth;
  }

  val_enth = enth_iter;

  if (counter >= counter_limit) {
    exit_code = 1;
  }

  return exit_code;
}
