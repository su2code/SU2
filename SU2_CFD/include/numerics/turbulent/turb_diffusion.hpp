/*!
 * \file turb_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
 private:

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(su2double *val_residual,
                                  su2double **Jacobian_i,
                                  su2double **Jacobian_j,
                                  CConfig *config) = 0;

 protected:
  bool implicit, incompressible;
  bool correct_gradient;
  unsigned short iVar, iDim;
  su2double **Mean_GradTurbVar;               /*!< \brief Average of gradients at cell face */
  su2double *Edge_Vector,                     /*!< \brief Vector from node i to node j. */
            *Proj_Mean_GradTurbVar_Normal,    /*!< \brief Mean_gradTurbVar DOT normal */
            *Proj_Mean_GradTurbVar_Edge,      /*!< \brief Mean_gradTurbVar DOT Edge_Vector */
            *Proj_Mean_GradTurbVar;           /*!< \brief Mean_gradTurbVar DOT normal, corrected if required*/
  su2double  dist_ij_2,                       /*!< \brief |Edge_Vector|^2 */
             proj_vector_ij;                  /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Scalar(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_gradient, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Scalar(void);

  /*!
   * \brief Compute the viscous residual using an average of gradients without correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                       su2double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGrad_TurbSA : public CAvgGrad_Scalar {
private:

  const su2double sigma;
  su2double nu_i, nu_j, nu_e;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                          su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSA(void);
};

/*!
 * \class CAvgGrad_TurbSA_Neg
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author F. Palacios
 */
class CAvgGrad_TurbSA_Neg : public CAvgGrad_Scalar {
private:

  const su2double sigma;
  const su2double cn1;
  su2double fn, Xi;
  su2double nu_i, nu_j, nu_ij, nu_tilde_ij, nu_e;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                                su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                      bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSA_Neg(void);
};

/*!
 * \class CAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGrad_TurbSST : public CAvgGrad_Scalar {
private:
  su2double sigma_k1, /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  sigma_k2,
  sigma_om1,
  sigma_om2;

  su2double diff_kine,  /*!< \brief Diffusivity for viscous terms of tke eq */
            diff_omega; /*!< \brief Diffusivity for viscous terms of omega eq */

  su2double F1_i, F1_j; /*!< \brief Menter's first blending function */

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                                su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                   const su2double* constants, bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbSST(void);

  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) {
    F1_i = val_F1_i; F1_j = val_F1_j;
  }

};
