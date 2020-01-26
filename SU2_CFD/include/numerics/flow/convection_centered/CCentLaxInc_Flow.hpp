/*!
 * \file CCentLaxInc_Flow.hpp
 * \brief Delaration of numerics class CCentLaxInc_Flow, the
 *        implementation is in the CCentLaxInc_Flow.cpp file.
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

#include "../../CNumerics.hpp"

/*!
 * \class CCentLaxInc_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme (modified with incompressible preconditioning).
 * \ingroup ConvDiscr
 * \author F. Palacios, T. Economon
 */
class CCentLaxInc_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_V,               /*!< \brief Difference of primitive variables. */
  *Velocity_i, *Velocity_j,        /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity_i,
  ProjVelocity_j,                  /*!< \brief Mean and projected velocities. */
  *ProjFlux,                       /*!< \brief Projected inviscid flux tensor. */
  sq_vel_i, sq_vel_j,              /*!< \brief Modulus of the velocity and the normal vector. */
  Temperature_i, Temperature_j,    /*!< \brief Temperature at node 0 and 1. */
  MeanDensity, MeanPressure,
  MeanBetaInc2, MeanEnthalpy,
  MeanCp, MeanTemperature,         /*!< \brief Mean values of primitive variables. */
  MeandRhodT,                      /*!< \brief Derivative of density w.r.t. temperature (variable density flows). */
  Param_p, Param_Kappa_0,          /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j,
  MeanLambda,                      /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc0,
  StretchingFactor,                /*!< \brief Streching parameters. */
  Epsilon_0;                       /*!< \brief Artificial dissipation values. */
  su2double **Precon;
  bool implicit,                   /*!< \brief Implicit calculation. */
  dynamic_grid,                    /*!< \brief Modification for grid movement. */
  variable_density,                /*!< \brief Variable density incompressible flows. */
  energy;                          /*!< \brief computation with the energy equation. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLaxInc_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLaxInc_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_residual - Pointer to the residual array.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};
