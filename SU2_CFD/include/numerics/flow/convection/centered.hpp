/*!
 * \file centered.hpp
 * \brief Declaration of numerics classes for centered schemes,
 *        the implementation is in centered.cpp.
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
 * \class CCentLaxInc_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme (modified with incompressible preconditioning).
 * \ingroup ConvDiscr
 * \author F. Palacios, T. Economon
 */
class CCentLaxInc_Flow final : public CNumerics {
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
  MeandRhodh,                      /*!< \brief Derivative of density w.r.t. enthalpy (variable density flows). */
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

  su2double fix_factor;            /*!< \brief Fix factor for Jacobians. */

  su2double** Jacobian_i = nullptr; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j = nullptr; /*!< \brief The Jacobian w.r.t. point j after computation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLaxInc_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CCentLaxInc_Flow(void) override;

  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CCentJSTInc_Flow
 * \brief Class for centered scheme - modified JST with incompressible preconditioning.
 * \ingroup ConvDiscr
 * \author F. Palacios, T. Economon
 */
class CCentJSTInc_Flow final : public CNumerics {

private:
  unsigned short iDim, iVar, jVar;   /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_V, *Diff_Lapl,     /*!< \brief Diference of primitive variables and undivided laplacians. */
  *Velocity_i, *Velocity_j,          /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity_i,
  ProjVelocity_j,                 /*!< \brief Mean and projected velocities. */
  sq_vel_i, sq_vel_j,             /*!< \brief Modulus of the velocity and the normal vector. */
  Temperature_i, Temperature_j,   /*!< \brief Temperature at node 0 and 1. */
  MeanDensity, MeanPressure,
  MeanBetaInc2, MeanEnthalpy,
  MeanCp, MeanTemperature,        /*!< \brief Mean values of primitive variables. */
  MeandRhodh,                     /*!< \brief Derivative of density w.r.t. enthalpy (variable density flows). */
  Param_p, Param_Kappa_2,
  Param_Kappa_4,                  /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j,
  MeanLambda,                     /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc2, sc4,
  StretchingFactor,               /*!< \brief Streching parameters. */
  *ProjFlux,                      /*!< \brief Projected inviscid flux tensor. */
  Epsilon_2, Epsilon_4;           /*!< \brief Artificial dissipation values. */
  su2double **Precon;
  bool implicit,         /*!< \brief Implicit calculation. */
  dynamic_grid,          /*!< \brief Modification for grid movement. */
  variable_density,      /*!< \brief Variable density incompressible flows. */
  energy;                /*!< \brief computation with the energy equation. */

  su2double fix_factor;  /*!< \brief Fix factor for Jacobians. */

  su2double** Jacobian_i = nullptr; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j = nullptr; /*!< \brief The Jacobian w.r.t. point j after computation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJSTInc_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CCentJSTInc_Flow(void) override;

  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

//TODO: Move this to another file as this is not a centered scheme
/*!
 * \class CUpwPB_Flow
 * \brief Class for solving an upwind flux for the pressure based incompressible momentum equations.
 * \ingroup ConvDiscr
 * \author A. Koodly
 */
class CUpwPB_Flow : public CNumerics {
private:
  bool implicit, dynamic_grid;
  bool gravity;
  su2double Froude, Upw_i, Upw_j;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *MeanVelocity, *Velocity_upw;
  su2double *ProjFlux_i, *ProjFlux_j;
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
   * \brief Compute the Roe's flux between two nodes i and j.
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
  inline void SetFaceVel(su2double val_FaceVel){ FaceVel = val_FaceVel; }
};
