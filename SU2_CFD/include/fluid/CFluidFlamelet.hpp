/*!
 * \file CFluidFlamelet.hpp
 * \brief  Defines the flamelet fluid model
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

#pragma once

#include "../../Common/include/containers/CLookUpTable.hpp"
#include "CFluidModel.hpp"

class CFluidFlamelet final : public CFluidModel {
 protected:
  int rank;

  unsigned short n_scalars;
  unsigned short n_lookups;
  unsigned short n_table_sources;
  unsigned short n_user_scalars; /*!< \brief number of passive reactant species. */
  unsigned short n_control_vars; /*!< \brief number of controlling variables. */

  vector<string> table_scalar_names; /*!< \brief vector to store names of scalar variables.   */
  vector<string> table_source_names; /*!< \brief vector to store names of scalar source variables.   */
  vector<string> table_lookup_names; /*!< \brief vector to store names of look up variables.   */

  vector<su2double> table_sources;

  su2double mass_diffusivity; /*!< \brief local mass diffusivity of the mixture */
  su2double molar_weight;     /*!< \brief local molar weight of the mixture */

  vector<su2double> source_scalar;

  CLookUpTable* look_up_table;

  vector<string> varnames_TD;     /*!< \brief Lookup names for thermodynamic state variables. */
  vector<su2double> val_vars_TD; /*!< \brief References to thermodynamic state variables. */

  vector<string> varnames_Sources;     /*!< \brief Lookup names for scalar source terms. */
  vector<su2double*> val_vars_Sources; /*!< \brief References to scalar sources. */

 public:
  CFluidFlamelet(CConfig* config, su2double value_pressure_operating);

  ~CFluidFlamelet();

  /*!
   * \brief Set the thermodynamioc state.
   * \param[in] val_temperature - temperature
   * \param[in] val_scalars - pointer to species mass fractions
   */
  void SetTDState_T(su2double val_temperature, const su2double* val_scalars = nullptr) override;

  /*!
   * \brief Set the reaction source terms for the transported species equations.
   * \param[in] val_scalars - pointer to species mass fractions
   * \param[out] exit_code = error code
   */
  unsigned long SetScalarSources(const su2double* val_scalars) override;

  /*!
   * \brief Get the total enthalpy from the tabulated temperature and species (inverse lookup).
   * \param[in/out] enthalpy - total enthalpy
   * \param[in] val_prog - progress variable
   * \param[in] val_temp - temperature
   * \param[in] initial_value - initial value for the iterative lookup procedure.
   * \param[out] exit_code = error code
   */
  unsigned long GetEnthFromTemp(su2double& enthalpy, const su2double val_prog, const su2double val_temp,
                                su2double initial_value = 0);

  /*!
   * \brief Return a pointer to the lookup table
   * \param[out] look_up_table - pointer to lookup table
   */
  inline CLookUpTable* GetLookUpTable() { return look_up_table; }

  /*!
   * \brief Get the mass diffusivity of the species.
   * \param[in] iVar - index to the species
   * \param[out] mass_diffusivity - value of the mass diffusivity
   */
  inline su2double GetMassDiffusivity(int iVar) final { return mass_diffusivity; }

  /*!
   * \brief Get the thermal conductivity of the species.
   * \param[in] iVar - index to the species
   * \param[out] Kt - value of the thermal conductivity
   */
  inline su2double GetThermalConductivity() { return Kt; }

  /*!
   * \brief Get the laminar viscosity of the species.
   * \param[in] iVar - index to the species
   * \param[out] Mu - value of the laminar viscosity
   */
  inline su2double GetLaminarViscosity() { return Mu; }

  /*!
   * \brief Get the reaction source term of all species equations.
   */
  inline const su2double* GetScalarSources() const { return source_scalar.data(); }

  /*!
   * \brief Get the index of a string in a vector.
   */
  inline int GetIdx(vector<string> v,string s) {

    auto it = find(v.begin(), v.end(), s);
  
    if (it != v.end()) 
        return(it - v.begin());
    else 
        return(-1);
  }

  /*!
   * \brief Get the value of the LUT variable using the string name.
   */
  inline su2double GetLUTVal(string varname_TD) {
      return val_vars_TD[GetIdx(varnames_TD,varname_TD)];
    };

 private:

  void PreprocessLookUp();

};
