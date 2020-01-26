/*!
 * \file CUpwRoe_AdjFlow.hpp
 * \brief Delaration of numerics class CUpwRoe_AdjFlow, the
 *        implementation is in the CUpwRoe_AdjFlow.cpp file.
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
 * \class CUpwRoe_AdjFlow
 * \brief Class for solving an approximate Riemann solver of Roe
 *        for the adjoint flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CUpwRoe_AdjFlow : public CNumerics {
private:
  su2double *Residual_Roe;
  su2double area, Sx, Sy, Sz, rarea, nx, ny, nz, rho_l, u_l, v_l, w_l, h_l, rho_r,
  u_r, v_r, w_r, h_r, psi1, psi2, psi3, psi4, psi5;
  su2double h, u, v, w, c, psi1_l, psi2_l, psi3_l, psi4_l, psi5_l,
  psi1_r, psi2_r, psi3_r, psi4_r, psi5_r, q_l, q_r, Q_l, Q_r, vn,
  rrho_l, weight, rweight1, cc;
  su2double l1psi, l2psi, absQ, absQp, absQm, q2, alpha, beta_u, beta_v, beta_w, Q, l1l2p, l1l2m, eta;
  su2double RoeDensity, RoeSoundSpeed, *RoeVelocity, *Lambda, *Velocity_i, *Velocity_j, **ProjFlux_i, **ProjFlux_j,
  Proj_ModJac_Tensor_ij, **Proj_ModJac_Tensor, Energy_i, Energy_j, **P_Tensor, **invP_Tensor;
  unsigned short iDim, iVar, jVar, kVar;
  bool implicit, grid_movement;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwRoe_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_AdjFlow(void);

  /*!
   * \brief Compute the adjoint Roe's flux between two nodes i and j.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config);
};
