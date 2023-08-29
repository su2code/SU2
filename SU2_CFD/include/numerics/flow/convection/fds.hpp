/*!
 * \file fds.hpp
 * \brief Declarations of classes for Flux-Difference-Spliting schemes,
 *        the implementations are in fds.cpp
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CUpwFDSInc_Flow
 * \brief Class for solving a Flux Difference Splitting (FDS) upwind method for the incompressible flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios, T. Economon
 */
class CUpwFDSInc_Flow final : public CNumerics {
private:
  bool implicit,     /*!< \brief Implicit calculation. */
  dynamic_grid,      /*!< \brief Modification for grid movement. */
  variable_density,  /*!< \brief Variable density incompressible flows. */
  energy;            /*!< \brief computation with the energy equation. */
  su2double *Diff_V;
  su2double *Velocity_i, *Velocity_j, *MeanVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **Precon, **invPrecon_A;
  su2double Proj_ModJac_Tensor_ij, Pressure_i,
  Pressure_j, ProjVelocity,
  MeandRhodT, dRhodT_i, dRhodT_j, /*!< \brief Derivative of density w.r.t. temperature (variable density flows). */
  Temperature_i, Temperature_j,   /*!< \brief Temperature at node 0 and 1. */
  MeanDensity, MeanPressure, MeanSoundSpeed, MeanBetaInc2, MeanEnthalpy, MeanCp, MeanTemperature; /*!< \brief Mean values of primitive variables. */
  unsigned short iDim, iVar, jVar, kVar;

  su2double* Flux = nullptr;        /*!< \brief The flux / residual across the edge. */
  su2double** Jacobian_i = nullptr; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j = nullptr; /*!< \brief The Jacobian w.r.t. point j after computation. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwFDSInc_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwFDSInc_Flow(void) override;

  /*!
   * \brief Compute the upwind flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};
