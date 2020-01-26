/*!
 * \file CUpwTurkel_Flow.hpp
 * \brief Delaration of numerics class CUpwTurkel_Flow, the
 *        implementation is in the CUpwTurkel_Flow.cpp file.
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
 * \class CUpwTurkel_Flow
 * \brief Class for solving an approximate Riemann solver of Roe with Turkel Preconditioning for the flow equations.
 * \ingroup ConvDiscr
 * \author A. K. Lonkar
 */
class CUpwTurkel_Flow : public CNumerics {
private:
  bool implicit, dynamic_grid;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **absPeJac, **invRinvPe, **R_Tensor, **Matrix, **Art_Visc;
  su2double sq_vel, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoePressure, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double Beta, Beta_min, Beta_max;
  su2double r_hat, s_hat, t_hat, rhoB2a2, sqr_one_m_Betasqr_Lam1;
  su2double Beta2, one_m_Betasqr, one_p_Betasqr, sqr_two_Beta_c_Area;
  su2double local_Mach;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwTurkel_Flow(void);

  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

};
