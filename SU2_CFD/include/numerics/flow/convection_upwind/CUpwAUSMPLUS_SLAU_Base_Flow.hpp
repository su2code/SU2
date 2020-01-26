/*!
 * \file CUpwAUSMPLUS_SLAU_Base_Flow.hpp
 * \brief Delaration of numerics class CUpwAUSMPLUS_SLAU_Base_Flow, the
 *        implementation is in the CUpwAUSMPLUS_SLAU_Base_Flow.cpp file.
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
 * \class CUpwAUSMPLUS_SLAU_Base_Flow
 * \brief Base class for AUSM+up(2) and SLAU(2) convective schemes.
 * \ingroup ConvDiscr
 * \author Amit Sachdeva
 */
class CUpwAUSMPLUS_SLAU_Base_Flow : public CNumerics {
protected:
  bool implicit;
  bool UseAccurateJacobian;
  bool HasAnalyticalDerivatives;
  su2double FinDiffStep;

  su2double MassFlux, DissFlux, Pressure;
  su2double *Velocity_i, *Velocity_j;
  su2double *psi_i, *psi_j;
  su2double dmdot_dVi[6], dmdot_dVj[6], dpres_dVi[6], dpres_dVj[6];

  /*--- Roe variables (for approximate Jacobian) ---*/
  su2double *Lambda, *Epsilon, *RoeVelocity, **P_Tensor, **invP_Tensor;

  /*!
   * \brief Compute the mass flux and pressure based on Primitives_i/j, derived classes must implement this method.
   * \note See the body of the (empty) default implementation for instructions on how to implement the method.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  virtual void ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure) = 0;

  /*!
   * \brief Compute the flux Jacobians of the Roe scheme to use as an approximation.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   */
  void ApproximateJacobian(su2double **val_Jacobian_i, su2double **val_Jacobian_j);

  /*!
   * \brief Compute the flux Jacobians using a mix of finite differences and manual differentiation.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   */
  void AccurateJacobian(CConfig *config, su2double **val_Jacobian_i, su2double **val_Jacobian_j);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUS_SLAU_Base_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSMPLUS_SLAU_Base_Flow(void);

  /*!
   * \brief Compute the AUSM+ and SLAU family of schemes.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};
