/*!
 * \file pressure_based.hpp
 * \brief Declaration of numerics classes for convective schemes for
 *        the pressure based solver, the implementation is in pressure_based.cpp.
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CCentLinearPB_Flow
 * \brief Class for computing a linear centered scheme.
 * \ingroup ConvDiscr
 * \author T. Aalbers
 */
class CCentLinearPB_Flow : public CNumerics {
private:
  bool implicit, dynamic_grid;
  su2double *Velocity_i, *Velocity_j, *MeanMassFlux;
  su2double MeanDensity;
  unsigned short iDim, iVar, jVar, kVar;
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLinearPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentLinearPB_Flow(void);
  
  /*!
   * \brief Compute the flow residual.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
};

/*!
 * \class CUpwPB_Flow
 * \brief Class for solving an upwind flux for the pressure based incompressible momentum equations.
 * \ingroup ConvDiscr
 * \author A. Koodly
 */
class CUpwPB_Flow : public CNumerics {
private:
  bool implicit, dynamic_grid;
  su2double *Velocity_i, *Velocity_j, *MeanMassFlux, *Velocity_upw;
  su2double Proj_ModJac_Tensor_ij, Pressure_i,
  Pressure_j, MeanDensity, MeanSoundSpeed, MeanPressure, MeanBetaInc2,
  ProjVelocity, FaceVel, Face_Flux;
  unsigned short iDim, iVar, jVar, kVar;
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;
  su2double **Jacobian_upw = nullptr;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwPB_Flow(void);
  
  /*!
   * \brief Compute the upwinded flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
  
  /*!
   * \brief Set the value of face velocity. This is used as a proxy for massflux at a face.
   * \param[in] val_FaceVel.
   */
  inline void SetFaceVel(su2double val_FaceVel) { FaceVel = val_FaceVel; }
};
