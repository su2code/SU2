/*!
 * \file turb_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../CNumerics.hpp"

/*!
 * \class CAvgGrad_Scalar
 * \brief Template class for computing viscous residual of scalar values
 * \details This class serves as a template for the scalar viscous residual
 *   classes.  The general structure of a viscous residual calculation is the
 *   same for many different  models, which leads to a lot of repeated code.
 *   By using the template design pattern, these sections of repeated code are
 *   moved to a shared base class, and the specifics of each model
 *   are implemented by derived classes.  In order to add a new residual
 *   calculation for a viscous residual, extend this class and implement
 *   the pure virtual functions with model-specific behavior.
 * \ingroup ViscDiscr
 * \author C. Pederson, A. Bueno, and F. Palacios
 */
class CAvgGrad_Scalar : public CNumerics {
protected:
  su2double
  Edge_Vector[MAXNDIM] = {0.0},             /*!< \brief Vector from node i to node j. */
  *Proj_Mean_GradTurbVar_Normal = nullptr,  /*!< \brief Mean_gradTurbVar DOT normal. */
  *Proj_Mean_GradTurbVar_Edge = nullptr,    /*!< \brief Mean_gradTurbVar DOT Edge_Vector. */
  *Proj_Mean_GradTurbVar = nullptr,         /*!< \brief Mean_gradTurbVar DOT normal, corrected if required. */
  dist_ij_2 = 0.0,                          /*!< \brief |Edge_Vector|^2 */
  proj_vector_ij = 0.0,                     /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  *Flux = nullptr,                          /*!< \brief Final result, diffusive flux/residual. */
  **Jacobian_i = nullptr,                   /*!< \brief Flux Jacobian w.r.t. node i. */
  **Jacobian_j = nullptr;                   /*!< \brief Flux Jacobian w.r.t. node j. */

  const bool correct_gradient = false, implicit = false, incompressible = false;

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method, derived classes
   *        should compute the Flux and Jacobians (i/j) inside this method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(const CConfig* config) = 0;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_gradient - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Scalar(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_gradient, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Scalar(void) override;

  /*!
   * \brief Compute the viscous residual using an average of gradients without correction.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGrad_TurbSA final : public CAvgGrad_Scalar {
private:
  const su2double sigma = 2.0/3.0;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_grad, const CConfig* config);
};

/*!
 * \class CAvgGrad_TurbSA_Neg
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author F. Palacios
 */
class CAvgGrad_TurbSA_Neg final : public CAvgGrad_Scalar {
private:
  const su2double sigma = 2.0/3.0;
  const su2double cn1 = 16.0;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                      bool correct_grad, const CConfig* config);
};

/*!
 * \class CAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGrad_TurbSST final : public CAvgGrad_Scalar {
private:
  const su2double
  sigma_k1 = 0.0, /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  sigma_k2 = 0.0,
  sigma_om1 = 0.0,
  sigma_om2 = 0.0;

  su2double F1_i, F1_j; /*!< \brief Menter's first blending function */

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] constants - Constants of the model.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                   const su2double* constants, bool correct_grad, const CConfig* config);

  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) override {
    F1_i = val_F1_i; F1_j = val_F1_j;
  }

};
