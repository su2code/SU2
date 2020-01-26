/*!
 * \file CUpwRoeBase_Flow.hpp
 * \brief Delaration of numerics class CUpwRoeBase_Flow, the
 *        implementation is in the CUpwRoeBase_Flow.cpp file.
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
 * \class CUpwRoeBase_Flow
 * \brief Intermediate base class for Roe schemes on ideal gas.
 * \ingroup ConvDiscr
 * \author A. Bueno, F. Palacios
 */
class CUpwRoeBase_Flow : public CNumerics {
protected:
  bool implicit, dynamic_grid, roe_low_dissipation;
  su2double *Velocity_i, *Velocity_j, *ProjFlux_i, *ProjFlux_j, *Conservatives_i, *Conservatives_j;
  su2double *Diff_U, *Lambda, **P_Tensor, **invP_Tensor;
  su2double *RoeVelocity, RoeDensity, RoeEnthalpy, RoeSoundSpeed, ProjVelocity, RoeSoundSpeed2, kappa;

  /*!
   * \brief Derived classes must specialize this method to add the specifics of the scheme they implement (e.g. low-Mach precond.).
   * \param[out] val_residual - Convective flux.
   * \param[out] val_Jacobian_i - Flux Jacobian wrt node i conservatives (implicit computation).
   * \param[out] val_Jacobian_j - Flux Jacobian wrt node j conservatives (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) = 0;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_low_dissipation - Use a low dissipation formulation.
   */
  CUpwRoeBase_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoeBase_Flow(void);

  /*!
   * \brief Compute the flux from node i to node j, part common to most Roe schemes.
   * \param[out] val_residual - Convective flux.
   * \param[out] val_Jacobian_i - Flux Jacobian wrt node i conservatives (implicit computation).
   * \param[out] val_Jacobian_j - Flux Jacobian wrt node j conservatives (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

};
