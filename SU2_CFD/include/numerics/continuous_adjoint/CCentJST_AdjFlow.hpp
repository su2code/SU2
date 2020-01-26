/*!
 * \file CCentJST_AdjFlow.hpp
 * \brief Delaration of numerics class CCentJST_AdjFlow, the
 *        implementation is in the CCentJST_AdjFlow.cpp file.
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
 * \class CCentJST_AdjFlow
 * \brief Class for and adjoint centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CCentJST_AdjFlow : public CNumerics {
private:
  su2double *Diff_Psi, *Diff_Lapl;
  su2double *Velocity_i, *Velocity_j;
  su2double *MeanPhi;
  unsigned short iDim, jDim, iVar, jVar;
  su2double Residual, ProjVelocity_i, ProjVelocity_j, ProjPhi, ProjPhi_Vel, sq_vel, phis1, phis2;
  su2double MeanPsiRho, MeanPsiE, Param_p, Param_Kappa_4, Param_Kappa_2, Local_Lambda_i, Local_Lambda_j, MeanLambda;
  su2double Phi_i, Phi_j, sc4, StretchingFactor, Epsilon_4, Epsilon_2;
  bool implicit, grid_movement;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_AdjFlow(void);

  /*!
   * \brief Compute the adjoint flow residual using a JST method.
   * \param[out] val_resconv_i - Pointer to the convective residual at point i.
   * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
   * \param[out] val_resconv_j - Pointer to the convective residual at point j.
   * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                        su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                        CConfig *config);
};
