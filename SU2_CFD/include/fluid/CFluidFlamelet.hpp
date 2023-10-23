/*!
 * \file CFluidFlamelet.hpp
 * \brief  Defines the flamelet fluid model
 * \author D. Mayer, T. Economon, N. Beishuizen, E. Bunschoten
 * \version 8.0.0 "Harrier"
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
#if defined(HAVE_MLPCPP)
#define MLP_CUSTOM_TYPE su2double
#include "../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif
#include "CFluidModel.hpp"

class CFluidFlamelet final : public CFluidModel {
 private:
  ENUM_DATADRIVEN_METHOD Kind_DataDriven_Method =
      ENUM_DATADRIVEN_METHOD::LUT; /*!< \brief Interpolation method for data set evaluation. */

  enum LOOKUP_TD { TEMPERATURE, HEATCAPACITY, VISCOSITY, CONDUCTIVITY, DIFFUSIONCOEFFICIENT, MOLARWEIGHT, SIZE };

  int rank;

  bool include_mixture_fraction = false; /*!< \brief include mixture fraction in controlling variables. */

  unsigned short n_scalars, n_lookups, n_user_scalars, /*!< \brief number of passive reactant species. */
      n_control_vars;                                  /*!< \brief number of controlling variables. */

  unsigned long extrapolation;

  INC_DENSITYMODEL density_model;
  vector<string> controlling_variable_names, passive_specie_names;
  vector<string> table_scalar_names; /*!< \brief vector to store names of scalar variables.   */

  su2double mass_diffusivity, /*!< \brief local mass diffusivity of the mixture */
      molar_weight;           /*!< \brief local molar weight of the mixture */

  CLookUpTable* look_up_table;

  /*--- Class variables for the multi-layer perceptron method ---*/
#ifdef USE_MLPCPP
  MLPToolbox::CLookUp_ANN* lookup_mlp; /*!< \brief Multi-layer perceptron collection. */
  MLPToolbox::CIOMap* iomap_TD;        /*!< \brief Input-output map for thermochemical properties. */
  MLPToolbox::CIOMap* iomap_Sources;   /*!< \brief Input-output map for species source terms. */
  MLPToolbox::CIOMap* iomap_LookUp;    /*!< \brief Input-output map for passive look-up terms. */
  MLPToolbox::CIOMap* iomap_Current;
#endif

  vector<su2double> scalars_vector;

  vector<string> varnames_TD, /*!< \brief Lookup names for thermodynamic state variables. */
      varnames_Sources, varnames_LookUp;

  vector<su2double> val_vars_TD, /*!< \brief References to thermodynamic state variables. */
      val_vars_Sources, val_vars_LookUp;

  void PreprocessLookUp(CConfig* config);

 public:
  CFluidFlamelet(CConfig* config, su2double value_pressure_operating);

  ~CFluidFlamelet();

  /*!
   * \brief Set the thermodynamic state.
   * \param[in] val_temperature - temperature
   * \param[in] val_scalars - pointer to species mass fractions
   */
  void SetTDState_T(su2double val_temperature, const su2double* val_scalars = nullptr) override;

  /*!
   * \brief Evaluate the flamelet manifold.
   * \param[in] input_scalar - scalar solution.
   * \param[in] input_varnames - names of variables to evaluate.
   * \param[in] output_refs - output data.
   * \param[out] Extrapolation - scalar solution is within bounds (0) or out of bounds (1).
   */
  inline unsigned long EvaluateDataSet(const vector<su2double>& input_scalar, unsigned short lookup_type,
                                       vector<su2double>& output_refs) override;

  /*!
   * \brief Check for out-of-bounds condition for data set interpolation.
   * \return - within bounds (0) or out of bounds (1).
   */
  unsigned long GetExtrapolation() const override { return extrapolation; }

  /*!
   * \brief Get the mass diffusivity of the species.
   * \param[in] iVar - index to the species
   * \param[out] mass_diffusivity - value of the mass diffusivity
   */
  inline su2double GetMassDiffusivity(int iVar) override { return mass_diffusivity; }

  /*!
   * \brief Get the thermal conductivity of the species.
   * \param[in] iVar - index to the species
   * \param[out] Kt - value of the thermal conductivity
   */
  inline su2double GetThermalConductivity() override { return Kt; }

  /*!
   * \brief Get the laminar viscosity of the species.
   * \param[in] iVar - index to the species
   * \param[out] Mu - value of the laminar viscosity
   */
  inline su2double GetLaminarViscosity() override { return Mu; }
};
