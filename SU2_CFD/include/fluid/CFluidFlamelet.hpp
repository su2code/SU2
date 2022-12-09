/*!
 * \file CFluidFlamelet.hpp
 * \brief  Defines the flamelet fluid model
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

#pragma once

#include "../../Common/include/containers/CLookUpTable.hpp"
#include "CFluidModel.hpp"

class CFluidFlamelet final : public CFluidModel {
 protected:
  int rank;

  unsigned short n_scalars;
  unsigned short n_lookups;
  unsigned short n_table_sources;
  unsigned short n_reactants; /*!< \brief number of passive reactant species. */
  unsigned short n_CV;        /*!< \brief number of controlling variables. */

  vector<string> table_scalar_names; /*!< \brief vector to store names of scalar variables.   */
  vector<string> table_source_names; /*!< \brief vector to store names of scalar source variables.   */
  vector<string> table_lookup_names; /*!< \brief vector to store names of look up variables.   */

  vector<su2double> table_sources;

  su2double mass_diffusivity; /*!< \brief local mass diffusivity of the mixture */
  su2double molar_weight; /*!< \brief local molar weight of the mixture */

  vector<su2double> source_scalar;
  vector<su2double> lookup_scalar;

  CLookUpTable* look_up_table;

  vector<string> varnames_TD;     // Lookup names for thermodynamic state variables.
  vector<su2double*> val_vars_TD; // References to thermodynamic state variables.

  vector<string> varnames_Sources;      // Lookup names for scalar source terms.
  vector<su2double*> val_vars_Sources;  // References to scalar sources.

  vector<string> varnames_LookUp;     // Lookup names for passive lookup variables.
  vector<su2double*> val_vars_LookUp; // References to lookup variables.

 public:
  CFluidFlamelet(CConfig* config, su2double value_pressure_operating);

  ~CFluidFlamelet();

  /*!
   * \brief Set the thermodynamioc state
   * \param[in] val_temperature - temperature
   * \param[in] val_scalars - pointer to species mass fractions
   */
  void SetTDState_T(su2double val_temperature, const su2double* val_scalars = nullptr) override;

  /*!
   * \brief Set the reaction source terms for the transported species equations
   * \param[in] val_scalars - pointer to species mass fractions
   * \param[out] exit_code = error code
   */
  unsigned long SetScalarSources(su2double* val_scalars);

  /*!
   * \brief Retrieve and set the lookup values for the species
   * \param[in] val_scalars - pointer to species mass fractions
   */
  unsigned long SetScalarLookups(su2double* val_scalars);

  //void SetTDState_prog_enth(su2double val_prog, su2double val_enth);

  /*!
   * \brief Get the total enthalpy from the tabulated temperature and species (inverse lookup)
   * \param[in/out] enthalpy - total enthalpy
   * \param[in] val_prog - progress variable
   * \param[in] val_temp - temperature
   * \param[out] exit_code = error code
   */
  unsigned long GetEnthFromTemp(su2double* enthalpy, su2double val_prog, su2double val_temp, su2double initial_value=0);

  /*!
   * \brief return a pointer to the lookup table
   * \param[out] look_up_table - pointer to lookup table
   */
  inline CLookUpTable* GetLookUpTable() { return look_up_table; }

  // inline su2double GetSourceEnergy() { return source_energy; }

  /*!
   * \brief Get the mass diffusivity of the species
   * \param[in] iVar - index to the species
   * \param[out] mass_diffusivity - value of the mass diffusivity
   */
  inline su2double GetMassDiffusivity(int iVar) final { return mass_diffusivity; }

  /*!
   * \brief Get the thermal conductivity of the species
   * \param[in] iVar - index to the species
   * \param[out] Kt - value of the thermal conductivity
   */
  inline su2double GetThermalConductivity() { return Kt; }

  /*!
   * \brief Get the laminar viscosity of the species
   * \param[in] iVar - index to the species
   * \param[out] Mu - value of the laminar viscosity
   */
  inline su2double GetLaminarViscosity() { return Mu; }

  /*!
   * \brief Get the enthalpy range in the lookup table
   */
  inline pair<su2double, su2double> GetTableLimitsEnth() { return look_up_table->GetTableLimitsEnth(); }

  /*!
   * \brief Get the progress variable range in the lookup table
   */
  inline pair<su2double, su2double> GetTableLimitsProg() { return look_up_table->GetTableLimitsProg(); }

  /*!
   * \brief Get the reaction source term of a species equation
   * \param[in] iVar - index to the species
   */
  inline su2double GetScalarSources(int iVar) { return source_scalar[iVar]; }

  /*!
   * \brief Get the reaction source term of all species equations
   */
  inline su2double* GetScalarSources() { return &source_scalar[0]; }

  /*!
   * \brief Get the value of the looked up variable
   * \param[in] i_scalar - index to the value that we need to retrieve from the lookup table
   */
  inline su2double GetScalarLookups(int i_scalar) { return lookup_scalar[i_scalar]; }

  void PreprocessLookUp();

  inline unsigned short GetNControllingVariables() { return n_CV; }
};
