/*!
 * \file roe.hpp
 * \brief Delarations of numerics classes for Roe-type schemes,
 *        implemented in roe.cpp.
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
  virtual ~CUpwRoeBase_Flow(void);

  /*!
   * \brief Compute the flux from node i to node j, part common to most Roe schemes.
   * \param[out] val_residual - Convective flux.
   * \param[out] val_Jacobian_i - Flux Jacobian wrt node i conservatives (implicit computation).
   * \param[out] val_Jacobian_j - Flux Jacobian wrt node j conservatives (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

};

/*!
 * \class CUpwRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author A. Bueno, F. Palacios
 */
class CUpwRoe_Flow : public CUpwRoeBase_Flow {
private:
  /*!
   * \brief Add standard Roe dissipation to the flux.
   * \param[out] val_residual - Convective flux.
   * \param[out] val_Jacobian_i - Flux Jacobian wrt node i conservatives (implicit computation).
   * \param[out] val_Jacobian_j - Flux Jacobian wrt node j conservatives (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_low_dissipation - Use a low dissipation formulation.
   */
  CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation);

};

/*!
 * \class CUpwL2Roe_Flow
 * \brief Class for solving an approximate Riemann solver of L2Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author E. Molina, A. Bueno, F. Palacios
 * \version 7.0.0 "Blackbird"
 */
class CUpwL2Roe_Flow : public CUpwRoeBase_Flow {
private:
  /*!
   * \brief Add L^2 Roe dissipation to the flux (low-Mach scheme).
   * \param[out] val_residual - Convective flux.
   * \param[out] val_Jacobian_i - Flux Jacobian wrt node i conservatives (implicit computation).
   * \param[out] val_Jacobian_j - Flux Jacobian wrt node j conservatives (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwL2Roe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

};

/*!
 * \class CUpwLMRoe_Flow
 * \brief Class for solving an approximate Riemann solver of LMRoe for the flow equations.
 * \ingroup ConvDiscr
 * \author E. Molina, A. Bueno, F. Palacios
 * \version 7.0.0 "Blackbird"
 */
class CUpwLMRoe_Flow : public CUpwRoeBase_Flow {
private:
  /*!
   * \brief Add LMRoe dissipation to the flux (low-Mach scheme).
   * \param[out] val_residual - Convective flux.
   * \param[out] val_Jacobian_i - Flux Jacobian wrt node i conservatives (implicit computation).
   * \param[out] val_Jacobian_j - Flux Jacobian wrt node j conservatives (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwLMRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

};

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

/*!
 * \class CUpwGeneralRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations for a general fluid model.
 * \ingroup ConvDiscr
 * \author S.Vitale, G.Gori, M.Pini
 */
class CUpwGeneralRoe_Flow : public CNumerics {
private:

  bool implicit, dynamic_grid;

  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon, MaxLambda, Delta;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,

  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, RoeSoundSpeed2,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho, kappa;
  unsigned short iDim, iVar, jVar, kVar;


  su2double StaticEnthalpy_i, StaticEnergy_i, StaticEnthalpy_j, StaticEnergy_j, Kappa_i, Kappa_j, Chi_i, Chi_j, Velocity2_i, Velocity2_j;
  su2double RoeKappa, RoeChi;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwGeneralRoe_Flow(void);

  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

  /*!
   * \brief Compute the Average for a general fluid flux between two nodes i and j.
   * Using the approach of Vinokur and Montagne'
   */

  void ComputeRoeAverage();
};
