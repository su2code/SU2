/*!
 * \file CDataDrivenFluid.hpp
 * \brief Defines a template fluid model class using multilayer perceptrons
 *  for theromodynamic state definition
 * \author E.C.Bunschoten
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

#include <vector>
#include "../../../Common/include/containers/CLookUpTable.hpp"
#if defined(HAVE_MLPCPP)
#define MLP_CUSTOM_TYPE su2double
#include "../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif
#include "CFluidModel.hpp"

/*!
 * \class CDataDrivenFluid
 * \brief Template class for fluid model definition using multi-layer perceptrons for
 * fluid dynamic state definition.
 * \author: E.C.Bunschoten.
 */
class CDataDrivenFluid final : public CFluidModel {
 protected:
  int rank{MASTER_NODE}; /*!< \brief MPI Rank. */
  ENUM_DATADRIVEN_METHOD Kind_DataDriven_Method =
      ENUM_DATADRIVEN_METHOD::LUT; /*!< \brief Interpolation method for data set evaluation. */

  string varname_rho,  /*!< \brief Controlling variable name for density. */
         varname_e;    /*!< \brief Controlling variable name for static energy. */
         
  size_t idx_rho, /*!< \brief Interpolator index for density input. */
      idx_e;      /*!< \brief Interpolator index for energy input. */

  su2double Newton_Relaxation, /*!< \brief Relaxation factor for Newton solvers. */
      rho_start,               /*!< \brief Starting value for the density in Newton solver processes. */
      e_start,                 /*!< \brief Starting value for the energy in Newton solver processes. */
      Newton_Tolerance,        /*!< \brief Normalized tolerance for Newton solvers. */
      rho_min, rho_max,        /*!< \brief Minimum and maximum density values in data set. */
      e_min, e_max;            /*!< \brief Minimum and maximum energy values in data set. */

  unsigned long MaxIter_Newton; /*!< \brief Maximum number of iterations for Newton solvers. */

  su2double dsde_rho, /*!< \brief Entropy derivative w.r.t. density. */
      dsdrho_e,       /*!< \brief Entropy derivative w.r.t. static energy. */
      d2sde2,         /*!< \brief Entropy second derivative w.r.t. static energy. */
      d2sdedrho,      /*!< \brief Entropy second derivative w.r.t. density and static energy. */
      d2sdrho2;       /*!< \brief Entropy second derivative w.r.t. static density. */

  su2double R_idealgas,     /*!< \brief Approximated ideal gas constant. */
            Cp_idealgas,    /*!< \brief Approximated ideal gas specific heat at constant pressure. */
            gamma_idealgas, /*!< \brief Approximated ideal gas specific heat ratio. */
            Cv_idealgas,    /*!< \brief Approximated ideal gas specific heat at constant volume. */
            P_middle,       /*!< \brief Pressure computed from the centre of the data set. */
            T_middle;       /*!< \brief Temperature computed from the centre of the data set. */

  su2double Enthalpy, /*!< \brief Fluid enthalpy value [J kg^-1] */
      dhdrho_e,       /*!< \brief Enthalpy derivative w.r.t. density. */
      dhde_rho;       /*!< \brief Enthalpy derivative w.r.t. static energy. */

  vector<string> input_names_rhoe, /*!< \brief Data-driven method input variable names of the independent variables
                                      (density, energy). */
      output_names_rhoe; /*!< \brief Output variable names listed in the data-driven method input file name. */

  vector<su2double*> outputs_rhoe; /*!< \brief Pointers to output variables. */

  /*--- Class variables for the multi-layer perceptron method ---*/
#ifdef USE_MLPCPP
  MLPToolbox::CLookUp_ANN* lookup_mlp; /*!< \brief Multi-layer perceptron collection. */
  MLPToolbox::CIOMap* iomap_rhoe;      /*!< \brief Input-output map. */
#endif
  vector<su2double> MLP_inputs; /*!< \brief Inputs for the multi-layer perceptron look-up operation. */

  CLookUpTable* lookup_table; /*!< \brief Look-up table regression object. */

  unsigned long outside_dataset, /*!< \brief Density-energy combination lies outside data set. */
      nIter_Newton;              /*!< \brief Number of Newton solver iterations. */

  /*!
   * \brief Map dataset variables to specific look-up operations.
   */
  void MapInputs_to_Outputs();

  /*!
   * \brief Evaluate dataset through multi-layer perceptron.
   * \param[in] rho - Density value.
   * \param[in] e - Static energy value.
   * \return Query point lies outside MLP normalization range (0 is inside, 1 is outside).
   */
  unsigned long Predict_MLP(su2double rho, su2double e);

  /*!
   * \brief Evaluate dataset through look-up table.
   * \param[in] rho - Density value.
   * \param[in] e - Static energy value.
   * \return Query point lies outside table data range (0 is inside, 1 is outside).
   */
  unsigned long Predict_LUT(su2double rho, su2double e);

  /*!
   * \brief Evaluate the data set.
   * \param[in] rho - Density value.
   * \param[in] e - Static energy value.
   */
  void Evaluate_Dataset(su2double rho, su2double e);

  /*!
   * \brief 2D Newton solver for computing the density and energy corresponding to Y1_target and Y2_target.
   * \param[in] Y1_target - Target value for output quantity 1.
   * \param[in] Y2_target - Target value for output quantity 2.
   * \param[in] Y1 - Pointer to output quantity 1.
   * \param[in] Y2 - Pointer to output quantity 2.
   * \param[in] dY1drho - Pointer to the partial derivative of quantity 1 w.r.t. density at constant energy.
   * \param[in] dY1de - Pointer to the partial derivative of quantity 1 w.r.t. energy at constant density.
   * \param[in] dY2drho - Pointer to the partial derivative of quantity 2 w.r.t. density at constant energy.
   * \param[in] dY2de - Pointer to the partial derivative of quantity 2 w.r.t. energy at constant density.
   */
  void Run_Newton_Solver(su2double Y1_target, su2double Y2_target, su2double* Y1, su2double* Y2, su2double* dY1drho,
                         su2double* dY1de, su2double* dY2drho, su2double* dY2de);

  /*!
   * \brief 1D Newton solver for computing the density or energy corresponding to Y_target.
   * \param[in] Y_target - Target quantity value.
   * \param[in] Y - Pointer to output quantity.
   * \param[in] X - Pointer to controlling variable (density or energy).
   * \param[in] dYdX - Pointer to the partial derivative of target quantity w.r.t. controlling variable.
   */
  void Run_Newton_Solver(su2double Y_target, su2double* Y, su2double* X, su2double* dYdX);

  void ComputeIdealGasQuantities();
 public:
  /*!
   * \brief Constructor of the class.
   */
  CDataDrivenFluid(const CConfig* config, bool display = true);

  ~CDataDrivenFluid();
  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy.
   * \param[in] rho - first thermodynamic variable (density).
   * \param[in] e - second thermodynamic variable (static energy).
   */
  void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Pressure  and Temperature.
   * \param[in] P - first thermodynamic variable (pressure).
   * \param[in] T - second thermodynamic variable (temperature).
   */
  void SetTDState_PT(su2double P, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Density.
   * \param[in] P - first thermodynamic variable (pressure).
   * \param[in] rho - second thermodynamic variable (density).
   */
  void SetTDState_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless Internal Energy using Pressure and Density.
   * \param[in] P - first thermodynamic variable (pressure).
   * \param[in] rho - second thermodynamic variable (density).
   */
  void SetEnergy_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless Internal Energy using Pressure and Density.
   * \param[in] rho - second thermodynamic variable (density).
   */
  void SetTDState_rhoT(su2double rho, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Enthalpy and Entropy.
   * \param[in] h - first thermodynamic variable (h).
   * \param[in] s - second thermodynamic variable (s).
   */
  void SetTDState_hs(su2double h, su2double s) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy.
   * \param[in] P - first thermodynamic variable (P).
   * \param[in] s - second thermodynamic variable (s).
   */
  void SetTDState_Ps(su2double P, su2double s) override;

  /*!
   * \brief Get fluid model extrapolation instance.
   * \return Query point lies outside fluid model data range.
   */
  unsigned long GetExtrapolation() const override { return outside_dataset; }

  /*!
   * \brief Get number of Newton solver iterations.
   * \return Newton solver iteration count at termination.
   */
  unsigned long GetnIter_Newton() override { return nIter_Newton; }
};
