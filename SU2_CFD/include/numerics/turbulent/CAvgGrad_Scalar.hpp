/*!
 * \file CAvgGrad_Scalar.hpp
 * \brief Delaration of numerics class CAvgGrad_Scalar, the
 *        implementation is in the CAvgGrad_Scalar.cpp file.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
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
