/*!
 * \file CAvgGradCorrected_P1.hpp
 * \brief Declaration and inlines of the class to compute
 *        the viscous residual in the P1 equation using corrected average gradient.
 * \author Ruben Sanchez
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

#include "../CNumericsRadiation.hpp"

class CAvgGradCorrected_P1 : public CNumericsRadiation {
 private:

 protected:

  unsigned short iVar, iDim;
  su2double **Mean_GradP1Var;               /*!< \brief Average of gradients at cell face */
  su2double *Edge_Vector,                   /*!< \brief Vector from node i to node j. */
            *Proj_Mean_GradP1Var_Edge,    /*!< \brief Mean_gradTurbVar DOT normal */
            *Proj_Mean_GradP1Var_Kappa,
            *Proj_Mean_GradP1Var_Corrected;           /*!< \brief Mean_gradTurbVar DOT normal, corrected if required*/
  su2double dist_ij,                        /*!< \brief Length of the edge and face. */
            proj_vector_ij;                 /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  su2double GammaP1;                        /*!< \brief P1 parameter */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_P1(unsigned short val_nDim, unsigned short val_nVar,
                       CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_P1(void);

  /*!
   * \brief Compute the viscous residual of the P1 equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                       su2double **Jacobian_j, CConfig *config);

};
