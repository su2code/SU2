/*!
 * \file centered.hpp
 * \brief Delaration of numerics classes for centered schemes,
 *        the implementation is in centered.cpp.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

  su2double Velocity_i[MAXNDIM] = {0.0};             /*!< \brief Velocity at node i. */
  su2double Velocity_j[MAXNDIM] = {0.0};             /*!< \brief Velocity at node j. */
  su2double MeanVelocity[MAXNDIM] = {0.0};           /*!< \brief Mean velocity. */
  su2double ProjVelocity_i, ProjVelocity_j;          /*!< \brief Velocities in the face normal direction. */
  su2double sq_vel_i,  sq_vel_j;                     /*!< \brief Squared norm of the velocity vectors. */
  su2double Energy_i,  Energy_j,  MeanEnergy;        /*!< \brief Energy at nodes i and j and mean. */
  su2double MeanDensity, MeanPressure, MeanEnthalpy; /*!< \brief Mean density, pressure, and enthalpy. */
  su2double *ProjFlux = nullptr;                     /*!< \brief "The" flux. */

  su2double *Diff_U = nullptr, *Diff_Lapl = nullptr;    /*!< \brief Differences of conservatives and undiv. Laplacians. */
  su2double Local_Lambda_i, Local_Lambda_j, MeanLambda; /*!< \brief Local eingenvalues. */
  su2double Param_p, Phi_i, Phi_j, StretchingFactor;    /*!< \brief Streching parameters. */
  su2double cte_0, cte_1;                               /*!< \brief Constants for the scalar dissipation Jacobian. */

  su2double ProjGridVel; /*!< \brief Projected grid velocity. */

  su2double** Jacobian_i = nullptr; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j = nullptr; /*!< \brief The Jacobian w.r.t. point j after computation. */

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
  CCentBase_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CCentBase_Flow(void) override;

  /*!
   * \brief Compute the flow residual using a centered method with artificial dissipation.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;

};

/*!
 * \class CCentLax_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CCentLax_Flow final : public CCentBase_Flow {
private:
  su2double Param_Kappa_0; /*!< \brief Artificial dissipation parameter. */
  su2double sc0;           /*!< \brief Streching parameter. */
  su2double Epsilon_0;     /*!< \brief Artificial dissipation coefficient. */

  /*!
   * \brief Lax-Friedrich first order dissipation term.
   * \param[in,out] val_residual - Pointer to the convective flux contribution to the residual.
   * \param[in,out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in,out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   */
  void DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) override;

  /*!
   * \brief Set input variables for AD preaccumulation.
   * \return true, as we will define inputs.
   */
  bool SetPreaccInVars(void) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

};

/*!
 * \class CCentJST_KE_Flow
 * \brief Class for centered scheme - JST_KE (no 4th dissipation order term).
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CCentJST_KE_Flow final : public CCentBase_Flow {

private:
  su2double Param_Kappa_2; /*!< \brief Artificial dissipation parameter. */
  su2double sc2;           /*!< \brief Streching parameter. */
  su2double Epsilon_2;     /*!< \brief Artificial dissipation coefficient. */

  /*!
   * \brief JST_KE second order dissipation term.
   * \param[in,out] val_residual - Pointer to the convective flux contribution to the residual.
   * \param[in,out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in,out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   */
  void DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) override;

  /*!
   * \brief Set input variables for AD preaccumulation.
   * \return true, as we will define inputs.
   */
  bool SetPreaccInVars(void) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_KE_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

};

/*!
 * \class CCentJST_Flow
 * \brief Class for centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CCentJST_Flow final : public CCentBase_Flow {

private:
  su2double Param_Kappa_2, Param_Kappa_4; /*!< \brief Artificial dissipation parameters. */
  su2double sc2, sc4;                     /*!< \brief Streching parameters. */
  su2double Epsilon_2, Epsilon_4;         /*!< \brief Artificial dissipation coefficients. */

  /*!
   * \brief JST second and forth order dissipation terms.
   * \param[in,out] val_residual - Pointer to the convective flux contribution to the residual.
   * \param[in,out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in,out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   */
  void DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) override;

  /*!
   * \brief Set input variables for AD preaccumulation.
   * \return true, as we will define inputs.
   */
  bool SetPreaccInVars(void) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

};

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
  MeandRhodT,                     /*!< \brief Derivative of density w.r.t. temperature (variable density flows). */
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
