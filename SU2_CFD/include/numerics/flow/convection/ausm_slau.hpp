/*!
 * \file ausm_slau.hpp
 * \brief Declaration of numerics classes for the AUSM family of schemes,
 *        including SLAU. The implementation is in ausm.cpp.
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

/*!
 * \class CUpwAUSMPLUSUP_Flow
 * \brief Class for solving an approximate Riemann AUSM+ -up.
 * \ingroup ConvDiscr
 * \author Amit Sachdeva
 */
class CUpwAUSMPLUSUP_Flow : public CUpwAUSMPLUS_SLAU_Base_Flow {
private:
  su2double Kp, Ku, sigma;

  /*!
   * \brief Mass flux and pressure for the AUSM+up scheme.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  void ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUSUP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

};

/*!
 * \class CUpwAUSMPLUSUP2_Flow
 * \brief Class for solving an approximate Riemann AUSM+ -up.
 * \ingroup ConvDiscr
 * \author Amit Sachdeva
 */
class CUpwAUSMPLUSUP2_Flow : public CUpwAUSMPLUS_SLAU_Base_Flow {
private:
  su2double Kp, sigma;

  /*!
   * \brief Mass flux and pressure for the AUSM+up2 scheme.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  void ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSMPLUSUP2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

};

/*!
 * \class CUpwSLAU_Flow
 * \brief Class for solving the Low-Dissipation AUSM.
 * \ingroup ConvDiscr
 * \author E. Molina
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
  void ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSLAU_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation);

};

/*!
 * \class CUpwSLAU2_Flow
 * \brief Class for solving the Simple Low-Dissipation AUSM 2.
 * \ingroup ConvDiscr
 * \author E. Molina
 */
class CUpwSLAU2_Flow : public CUpwSLAU_Flow {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSLAU2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation);

};

/*!
 * \class CUpwAUSM_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CUpwAUSM_Flow : public CNumerics {
private:
  bool implicit;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *RoeVelocity;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *delta_wave, *delta_vel;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor;
  su2double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
  Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
  ProjVelocity, ProjVelocity_i, ProjVelocity_j;
  unsigned short iDim, iVar, jVar, kVar;
  su2double mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSM_Flow(void);

  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};
