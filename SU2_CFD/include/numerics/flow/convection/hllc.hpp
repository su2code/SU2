/*!
 * \file hllc.hpp
 * \brief Declaration of HLLC numerics classes, implemented in hllc.cpp.
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
 * \class CUpwHLLC_Flow
 * \brief Class for solving an approximate Riemann HLLC.
 * \ingroup ConvDiscr
 * \author G. Gori, Politecnico di Milano
 * \version 8.0.0 "Harrier"
 */
class CUpwHLLC_Flow final : public CNumerics {
private:
  bool implicit, dynamic_grid;
  unsigned short iDim, jDim, iVar, jVar;

  su2double *IntermediateState;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;

  su2double sq_vel_i, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, ProjVelocity_i;
  su2double sq_vel_j, Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, ProjVelocity_j;

  su2double sq_velRoe, RoeDensity, RoeEnthalpy, RoeSoundSpeed, RoeProjVelocity, ProjInterfaceVel;

  su2double sL, sR, sM, pStar, EStar, rhoSL, rhoSR, Rrho, kappa;

  su2double Omega, RHO, OmegaSM;
  su2double *dSm_dU, *dPI_dU, *drhoStar_dU, *dpStar_dU, *dEStar_dU;

  su2double* Flux;        /*!< \brief The flux accross the face. */
  su2double** Jacobian_i; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j; /*!< \brief The Jacobian w.r.t. point j after computation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwHLLC_Flow(void) override;

  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};

/*!
 * \class CUpwGeneralHLLC_Flow
 * \brief Class for solving an approximate Riemann HLLC.
 * \ingroup ConvDiscr
 * \author G. Gori, Politecnico di Milano
 * \version 8.0.0 "Harrier"
 */
class CUpwGeneralHLLC_Flow final : public CNumerics {
private:
  bool implicit, dynamic_grid;
  unsigned short iDim, jDim, iVar, jVar;

  su2double *IntermediateState;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;

  su2double sq_vel_i, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i, ProjVelocity_i, StaticEnthalpy_i, StaticEnergy_i;
  su2double sq_vel_j, Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, ProjVelocity_j, StaticEnthalpy_j, StaticEnergy_j;

  su2double sq_velRoe, RoeDensity, RoeEnthalpy, RoeSoundSpeed, RoeProjVelocity, ProjInterfaceVel;
  su2double Kappa_i, Kappa_j, Chi_i, Chi_j, RoeKappa, RoeChi, RoeKappaStaticEnthalpy;

  su2double sL, sR, sM, pStar, EStar, rhoSL, rhoSR, Rrho, kappa;

  su2double Omega, RHO, OmegaSM;
  su2double *dSm_dU, *dPI_dU, *drhoStar_dU, *dpStar_dU, *dEStar_dU;

  su2double* Flux;        /*!< \brief The flux accross the face. */
  su2double** Jacobian_i; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j; /*!< \brief The Jacobian w.r.t. point j after computation. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwGeneralHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwGeneralHLLC_Flow(void) override;

  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

  /*!
   * \brief Compute the Average quantities for a general fluid flux between two nodes i and j.
   * Using the approach of Vinokur and Montagne'
   */
  void VinokurMontagne();
};
