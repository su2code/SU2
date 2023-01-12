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
#include "../../Common/include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"
#include "../../Common/include/toolboxes/multilayer_perceptron/CIOMap.hpp"
#include "CFluidModel.hpp"

class CFluidFlamelet final : public CFluidModel {
 protected:
  int rank;

  unsigned short n_scalars;
  unsigned short n_lookups;
  unsigned short n_table_sources;
  unsigned short n_user_scalars; /*!< \brief number of passive reactant species. */
  unsigned short n_CV; /*!< \brief number of controlling variables. */

  unsigned short manifold_format = ENUM_DATADRIVEN_METHOD::LUT;

  bool PreferentialDiffusion = false;
  su2vector<string> controlling_variables;
  su2vector<su2double> val_controlling_vars;

  vector<string> table_scalar_names; /*!< \brief vector to store names of scalar variables.   */
  vector<string> table_source_names; /*!< \brief vector to store names of scalar source variables.   */
  vector<string> table_lookup_names; /*!< \brief vector to store names of look up variables.   */

  vector<su2double> table_sources;

  su2double mass_diffusivity; /*!< \brief local mass diffusivity of the mixture */
  su2double molar_weight; /*!< \brief local molar weight of the mixture */
  su2double beta_progvar, 
            beta_enth_thermal, 
            beta_enth, 
            beta_mixfrac;

  vector<su2double> source_scalar;
  vector<su2double> lookup_scalar;

  CLookUpTable* look_up_table;
  MLPToolbox::CLookUp_ANN* look_up_ANN;

  su2vector<string> varnames_TD; /*!< \brief Lookup names for thermodynamic state variables. */
  su2vector<su2double*> val_vars_TD; /*!< \brief References to thermodynamic state variables. */
  su2matrix<su2double> dTD_dCV;
  MLPToolbox::CIOMap * iomap_TD = nullptr;

  su2vector<string> varnames_PD; /*!< \brief Lookup names for preferential diffusion scalars */
  su2vector<su2double*> val_vars_PD; /*!< \brief References to preferential diffusion scalars*/
  su2matrix<su2double> dPD_dCV;
  MLPToolbox::CIOMap * iomap_PD = nullptr;

  su2vector<string> varnames_Sources; /*!< \brief Lookup names for scalar source terms. */
  su2vector<su2double*> val_vars_Sources; /*!< \brief References to scalar sources. */
  su2matrix<su2double> dSources_dCV;
  MLPToolbox::CIOMap * iomap_Sources = nullptr;

  su2vector<string> varnames_LookUp; /*!< \brief Lookup names for passive lookup variables. */
  su2vector<su2double*> val_vars_LookUp; /*!< \brief References to lookup variables. */
  MLPToolbox::CIOMap * iomap_LookUp = nullptr;

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

  /*!
   * \brief Set the preferential diffusion terms for the transported scalar equations.
   * \param[in] val_scalars - pointer to species mass fractions
   */
  unsigned long SetPreferentialDiffusionScalars(su2double* val_scalars);

  /*!
   * \brief Get the total enthalpy from the tabulated temperature and species (inverse lookup)
   * \param[in/out] enthalpy - total enthalpy
   * \param[in] val_prog - progress variable
   * \param[in] val_mixfrac - Mixture fraction
   * \param[in] val_temp - temperature
   * \param[out] exit_code = error code
   */
  unsigned long GetEnthFromTemp(su2double* enthalpy, su2double val_prog, su2double val_mixfrac, su2double val_temp, su2double initial_value=0);

  /*!
   * \brief return a pointer to the lookup table
   * \param[out] look_up_table - pointer to lookup table
   */
  inline CLookUpTable* GetLookUpTable() { return look_up_table; }

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
  inline pair<su2double, su2double> GetTableLimitsEnth() {return look_up_table->GetTableLimitsY(); }

  /*!
   * \brief Get the progress variable range in the lookup table
   */
  inline pair<su2double, su2double> GetTableLimitsProg() { return look_up_table->GetTableLimitsX(); }

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
   * \brief Get the preferential diffusion scalar
   * \param[in] iVar - index to the species
   */
  inline su2double GetPreferentialDiffusionScalar(size_t iVar) const { return *val_vars_PD[iVar]; }

  /*!
   * \brief Get the value of the looked up variable
   * \param[in] i_scalar - index to the value that we need to retrieve from the lookup table
   */
  inline su2double GetScalarLookups(int i_scalar) { return lookup_scalar[i_scalar]; }

  void PreprocessLookUp();

  inline unsigned short GetNControllingVariables() { return n_CV; }

  unsigned long Evaluate_Dataset(su2vector<string>& varnames, su2vector<su2double*>& val_vars, su2matrix<su2double> dOutputs_dInputs, MLPToolbox::CIOMap* iomap=nullptr);
};
