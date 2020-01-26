/*!
 * \file CCentBase_Flow.hpp
 * \brief Delaration of numerics class CCentBase_Flow, the
 *        implementation is in the CCentBase_Flow.cpp file.
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
 * \class CCentBase_Flow
 * \brief Intermediate class to define centered schemes.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CCentBase_Flow : public CNumerics {

protected:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  bool dynamic_grid;               /*!< \brief Consider grid movement. */
  bool implicit;                   /*!< \brief Implicit calculation (compute Jacobians). */
  su2double fix_factor;            /*!< \brief Fix factor for dissipation Jacobians (more diagonal dominance). */

  su2double *Velocity_i, *Velocity_j, *MeanVelocity; /*!< \brief Velocity at nodes i and j and mean. */
  su2double ProjVelocity_i, ProjVelocity_j;          /*!< \brief Velocities in the face normal direction. */
  su2double sq_vel_i,  sq_vel_j;                     /*!< \brief Squared norm of the velocity vectors. */
  su2double Energy_i,  Energy_j,  MeanEnergy;        /*!< \brief Energy at nodes i and j and mean. */
  su2double MeanDensity, MeanPressure, MeanEnthalpy; /*!< \brief Mean density, pressure, and enthalpy. */
  su2double *ProjFlux;                               /*!< \brief Projected inviscid flux. */

  su2double *Diff_U, *Diff_Lapl;                        /*!< \brief Differences of conservatives and undiv. Laplacians. */
  su2double Local_Lambda_i, Local_Lambda_j, MeanLambda; /*!< \brief Local eingenvalues. */
  su2double Param_p, Phi_i, Phi_j, StretchingFactor;    /*!< \brief Streching parameters. */
  su2double cte_0, cte_1;                               /*!< \brief Constants for the scalar dissipation Jacobian. */

  su2double ProjGridVel; /*!< \brief Projected grid velocity. */

  /*!
   * \brief Hook method for derived classes to define preaccumulated variables, optional to implement.
   * \return true if any variable was set as preacc. input, in which case the residual will be output.
   */
  virtual bool SetPreaccInVars(void) {return false;}

  /*!
   * \brief Derived classes must implement this method, called in ComputeResidual after inviscid part.
   * \param[in,out] val_residual - Pointer to the convective flux contribution to the residual.
   * \param[in,out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in,out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   */
  virtual void DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) = 0;

  /*!
   * \brief Add the contribution of a scalar dissipation term to the Jacobians.
   * \param[in,out] val_Jacobian_i - Jacobian of the numerical method at node i.
   * \param[in,out] val_Jacobian_j - Jacobian of the numerical method at node j.
   */
  void ScalarDissipationJacobian(su2double **val_Jacobian_i, su2double **val_Jacobian_j);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentBase_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CCentBase_Flow(void);

  /*!
   * \brief Compute the flow residual using a centered method with artificial dissipation.
   * \param[out] val_residual - Pointer to the convective flux contribution to the residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

};
