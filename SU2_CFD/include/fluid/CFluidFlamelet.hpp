/*!
 * \file CFluidScalar.hpp
 * \brief  Defines the multicomponent incompressible Ideal Gas model for mixtures.
 * \author T. Economon, Mark Heimgartner, Cristopher Morales Ubal
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

#pragma once

#include "../../Common/include/containers/CLookUpTable.hpp"
#include "CFluidModel.hpp"


class CFluidFlamelet final : public CFluidModel {

protected:

  int rank;

  unsigned short n_scalars;
  unsigned short n_lookups;
  unsigned short n_table_sources; 

  vector<string> table_scalar_names;    /*!< \brief vector to store names of scalar variables.   */
  vector<string> table_source_names;    /*!< \brief vector to store names of scalar source variables.   */
  vector<string> table_lookup_names;    /*!< \brief vector to store names of look up variables.   */

  su2double mass_diffusivity;
  //su2double source_energy;
  su2double dDensitydPV;
  su2double dSourcePVdPV;
  su2double dDensitydEnth;

  vector<su2double> source_scalar;
  vector<su2double> lookup_scalar;

  CLookUpTable *look_up_table;

 public:
  CFluidFlamelet(CConfig *config, su2double value_pressure_operating);

  ~CFluidFlamelet();

  void SetTDState_T(su2double val_temperature, const su2double* val_scalars = nullptr) override;

  unsigned long SetScalarSources(su2double *val_scalars);

  unsigned long SetScalarLookups(su2double *val_scalars);

  void SetTDState_prog_enth(su2double val_prog, su2double val_enth);

  unsigned long GetEnthFromTemp(su2double *enthalpy, su2double val_prog, su2double val_temp);

  inline CLookUpTable* GetLookUpTable() {return look_up_table; }

  //inline su2double GetSourceEnergy() { return source_energy; }

  inline su2double GetMassDiffusivity() { return mass_diffusivity; }

  inline su2double GetThermalConductivity() { return Kt; }

  inline su2double GetLaminarViscosity() { return Mu; }

  inline pair<su2double, su2double> GetTableLimitsEnth() { return look_up_table->GetTableLimitsEnth(); }

  inline pair<su2double, su2double> GetTableLimitsProg() { return look_up_table->GetTableLimitsProg(); }

  inline su2double GetdDensitydPV() { return dDensitydPV; }

  inline su2double GetdSourcePVdPV() { return dSourcePVdPV; }

  inline su2double GetdDensitydEnth() { return dDensitydEnth; }

  inline su2double GetScalarSources(int i_scalar) { return source_scalar.at(i_scalar); }

  inline su2double *GetScalarSources(){ return &source_scalar[0]; }

  inline unsigned short GetNScalars() { return n_scalars; }

  inline su2double GetScalarLookups(int i_scalar) { return lookup_scalar.at(i_scalar); }

};
