/*!
 * \file CAvgGrad_AdjFlow.hpp
 * \brief Delaration of numerics class CAvgGrad_AdjFlow, the
 *        implementation is in the CAvgGrad_AdjFlow.cpp file.
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
 * \class CAvgGrad_AdjFlow
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios
 */
class CAvgGrad_AdjFlow : public CNumerics {
private:
  su2double *Velocity_i;  /*!< \brief Auxiliary vector for storing the velocity of point i. */
  su2double *Velocity_j;  /*!< \brief Auxiliary vector for storing the velocity of point j. */
  su2double *Mean_Velocity;
  su2double *Mean_GradPsiE;  /*!< \brief Counter for dimensions of the problem. */
  su2double **Mean_GradPhi;  /*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;    /*!< \brief Vector going from node i to node j. */
  bool implicit;             /*!< \brief Implicit calculus. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjFlow(void);

  /*!
   * \brief Residual computation.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};
