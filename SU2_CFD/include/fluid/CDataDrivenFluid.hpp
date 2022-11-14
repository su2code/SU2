/*!
 * \file CDataDrivenFluid.hpp
 * \brief Defines the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 7.3.1 "Blackbird"
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
//#include "../numerics/MultiLayer_Perceptron/LookUp_MLP.hpp"
#include "../../Common/include/toolboxes/multilayer_perceptron/CLookUp_ANN.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"

/*!
 * \class CDataDrivenFluid
 * \brief Child class for defining the ideal gas model.
 * \author: S.Vitale, M.Pini.
 */
class CDataDrivenFluid : public CFluidModel {
 protected:
  size_t idx_rho, 
        idx_e;
  su2double Gamma{0.0};           /*!< \brief Ratio of Specific Heats. */
  su2double Gamma_Minus_One{0.0}; /*!< \brief Ratio of Specific Heats Minus One. */
  su2double Gas_Constant{0.0};    /*!< \brief Gas Constant. */
  bool ComputeEntropy{true},      /*!< \brief Whether or not to compute entropy. */
        failedNewtonSolver{false};
  //LookUp_MLP *ANN;
  MLPToolbox::CLookUp_ANN *ANN;
  MLPToolbox::CIOMap * input_output_map, * input_output_map_lower_rho;

  CLookUpTable *look_up_table;

  vector<string> input_names, input_names_lower_rho;
  vector<su2double> inputs;
  vector<string> output_names;
  
  vector<su2double*> outputs;
  su2double rho_start, e_start, S, ds_de, ds_drho, d2s_dedrho, d2s_de2, d2s_drho2;
  su2double **dOutputs_dInputs;
  
  pair<su2double, su2double> rho_maxmin, e_maxmin;

  unsigned short boundsviolation;
  unsigned long nIter_NewtonSolver{0};
  bool Reevaluate_MLP{true};
  bool LUT = false,
       MLP = false;
  unsigned short Kind_DataDriven_Method = ENUM_DATADRIVEN_METHOD::LUT;
  su2double Newton_Relaxation = 0.05;
  
 public:
  /*!
   * \brief Constructor of the class.
   */
  CDataDrivenFluid(CConfig *config);

  ~CDataDrivenFluid(){
      switch (Kind_DataDriven_Method)
      {
      case ENUM_DATADRIVEN_METHOD::LUT:
            delete look_up_table;
            break;
      case ENUM_DATADRIVEN_METHOD::MLP:
            delete ANN;
            delete input_output_map;
            delete input_output_map_lower_rho;
            break;
      default:
            break;
      }
      for(size_t i=0; i<6; i++){
            delete [] dOutputs_dInputs[i];
      }
      delete [] dOutputs_dInputs;
  };

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

  /*!
   * \brief Set the Dimensionless State using Enthalpy and Entropy
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   *
   */
  void SetTDState_hs(su2double h, su2double s) override;

  /*!
   * \brief Set the Dimensionless State using Density and Temperature
   * \param[in] th1 - first thermodynamic variable (rho).
   * \param[in] th2 - second thermodynamic variable (T).
   *
   */
  void SetTDState_rhoT(su2double rho, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_Ps(su2double P, su2double s) override;


  /*!
   * \brief compute some derivatives of enthalpy and entropy needed for subsonic inflow BC
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   *
   */
  void ComputeDerivativeNRBC_Prho(su2double P, su2double rho) override;

  unsigned short GetClipping() {return boundsviolation;}
  unsigned long GetnIter_NewtonSolver() {return nIter_NewtonSolver;}

  void Predict_MLP(su2double rho, su2double e);
  
  void Predict_LUT(su2double rho, su2double e);

  void Evaluate_Dataset(su2double rho, su2double e);

  void SetDensity(su2double rho) override {rho_start = rho;}
  void SetEnergy(su2double e) override {e_start = e;}
};
