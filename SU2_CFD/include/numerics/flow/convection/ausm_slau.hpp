/*!
 * \file ausm_slau.hpp
 * \brief Declaration of numerics classes for the AUSM family of schemes,
 *        including SLAU. The implementation is in ausm.cpp.
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
 * \class CUpwAUSMPLUS_SLAU_Base_Flow
 * \brief Base class for AUSM+up(2) and SLAU(2) convective schemes.
 * \ingroup ConvDiscr
 * \author Amit Sachdeva, P. Gomes
 */
class CUpwAUSMPLUS_SLAU_Base_Flow : public CNumerics {
protected:
  bool implicit;
  bool UseAccurateJacobian;
  bool HasAnalyticalDerivatives;
  su2double FinDiffStep;

  su2double MassFlux, DissFlux, Pressure;
  su2double Velocity_i[MAXNDIM] = {0.0}, Velocity_j[MAXNDIM] = {0.0};
  su2double *psi_i = nullptr, *psi_j = nullptr;
  su2double dmdot_dVi[6], dmdot_dVj[6], dpres_dVi[6], dpres_dVj[6];

  /*--- Roe variables (for approximate Jacobian) ---*/
  su2double *Lambda = nullptr, *Epsilon = nullptr, RoeVelocity[MAXNDIM] = {0.0};
  su2double **P_Tensor = nullptr, **invP_Tensor = nullptr;

  su2double* Flux = nullptr;        /*!< \brief The flux accross the face. */
  su2double** Jacobian_i = nullptr; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j = nullptr; /*!< \brief The Jacobian w.r.t. point j after computation. */

  /*!
   * \brief Compute the mass flux and pressure based on Primitives_i/j, derived classes must implement this method.
   * \note See the body of the (empty) default implementation for instructions on how to implement the method.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  virtual void ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) = 0;

private:
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
  void AccurateJacobian(const CConfig* config, su2double **val_Jacobian_i, su2double **val_Jacobian_j);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUS_SLAU_Base_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSMPLUS_SLAU_Base_Flow(void) override;

  /*!
   * \brief Compute the AUSM+ and SLAU family of schemes.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;

};

/*!
 * \class CUpwAUSMPLUSUP_Flow
 * \brief Class for solving an approximate Riemann AUSM+ -up.
 * \ingroup ConvDiscr
 * \author Amit Sachdeva, P. Gomes
 */
class CUpwAUSMPLUSUP_Flow final : public CUpwAUSMPLUS_SLAU_Base_Flow {
private:
  su2double Kp, Ku, sigma;

  /*!
   * \brief Mass flux and pressure for the AUSM+up scheme.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  void ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) override;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUSUP_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

};

/*!
 * \class CUpwAUSMPLUSUP2_Flow
 * \brief Class for solving an approximate Riemann AUSM+ -up.
 * \ingroup ConvDiscr
 * \author Amit Sachdeva, P. Gomes
 */
class CUpwAUSMPLUSUP2_Flow final : public CUpwAUSMPLUS_SLAU_Base_Flow {
private:
  su2double Kp, sigma;

  /*!
   * \brief Mass flux and pressure for the AUSM+up2 scheme.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  void ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) override;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUSUP2_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

};

/*!
 * \class CUpwSLAU_Flow
 * \brief Class for solving the Low-Dissipation AUSM.
 * \ingroup ConvDiscr
 * \author E. Molina, P. Gomes
 */
class CUpwSLAU_Flow : public CUpwAUSMPLUS_SLAU_Base_Flow {
protected:
  bool slau_low_diss;
  bool slau2;

  /*!
   * \brief Mass flux and pressure for the SLAU and SLAU2 schemes.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  void ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) final;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSLAU_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config, bool val_low_dissipation);

};

/*!
 * \class CUpwSLAU2_Flow
 * \brief Class for solving the Simple Low-Dissipation AUSM 2.
 * \ingroup ConvDiscr
 * \author E. Molina, P. Gomes
 */
class CUpwSLAU2_Flow final : public CUpwSLAU_Flow {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSLAU2_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config, bool val_low_dissipation);

};

/*!
 * \class CUpwAUSM_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CUpwAUSM_Flow final : public CNumerics {
private:
  bool implicit;
  su2double *Diff_U;
  su2double Velocity_i[MAXNDIM], Velocity_j[MAXNDIM], RoeVelocity[MAXNDIM];
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;

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
  CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSM_Flow(void) override;

  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};
