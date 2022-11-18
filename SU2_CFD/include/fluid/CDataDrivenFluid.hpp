/*!
 * \file CDataDrivenFluid.hpp
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
 * \class CDataDrivenFluid
 * \brief Template class for fluid model definition using multi-layer perceptrons for 
 * fluid dynamic state definition.
 * \author: E.Bunschoten.
 */
class CDataDrivenFluid : public CFluidModel {
 protected:

  unsigned short Kind_DataDriven_Method = ENUM_DATADRIVEN_METHOD::LUT;
  size_t idx_rho,
         idx_e;

  su2double Newton_Relaxation,
            rho_start,
            e_start;

  string input_filename;

  vector<string> input_names_rhoe,
                 output_names_rhoe;

  vector<su2double*> outputs_rhoe;

  MLPToolbox::CLookUp_ANN * lookup_mlp;
  MLPToolbox::CIOMap * iomap_rhoe;
  vector<su2double> MLP_inputs;

  

  void MapInputs_to_Outputs();

 public:
  /*!
   * \brief Constructor of the class.
   */
  CDataDrivenFluid(const CConfig* config);

  ~CDataDrivenFluid();
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
