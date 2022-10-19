/*!
 * \file CMLPGas.hpp
 * \brief Defines a template fluid model class using multilayer perceptrons
 *  for theromodynamic state definition
 * \author E.Bunschoten
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

#include "CFluidModel.hpp"
#include "../../../Common/include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"
/*!
 * \class CMLPGas
 * \brief Template class for fluid model definition using multi-layer perceptrons for 
 * fluid dynamic state definition.
 * \author: E.Bunschoten.
 */
class CMLPGas : public CFluidModel {
 protected:
  su2double Gamma{0.0};           /*!< \brief Ratio of Specific Heats. */
  su2double Gamma_Minus_One{0.0}; /*!< \brief Ratio of Specific Heats Minus One. */
  su2double Gas_Constant{0.0};    /*!< \brief Gas Constant. */

  su2double R_u = 8.31451;

  string ann_input_filename;
  MLPToolbox::CLookUp_ANN * lookup_ann;
  vector<string> mlp_input_names, lookup_names;
  vector<su2double*> lookup_data;

  MLPToolbox::CIOMap * iomap_rhoe;
  vector<string> input_names_rhoe, 
                 output_names_rhoe;
  vector<su2double*> outputs_rhoe;
  
  MLPToolbox::CIOMap * iomap_PT;
  vector<string> input_names_PT, 
                 output_names_PT;
  vector<su2double*> outputs_PT;

  MLPToolbox::CIOMap * iomap_Prho;
  vector<string> input_names_Prho, 
                 output_names_Prho;
  vector<su2double*> outputs_Prho;

  void MapInputs_to_Outputs();

 public:
  /*!
   * \brief Constructor of the class.
   */
  CMLPGas(const CConfig* config);

  ~CMLPGas();
  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Pressure  and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  void SetTDState_PT(su2double P, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetTDState_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless Internal Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho(su2double P, su2double rho) override;
};
