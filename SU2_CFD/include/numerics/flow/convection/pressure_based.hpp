/*!
 * \file pressure_based.hpp
 * \brief Declaration of numerics classes for convective schemes for
 *        the pressure based solver, the implementation is in pressure_based.cpp.
 * \author T. Aalbers
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
 * \class CPBConvection_Base
 * \brief Class for computing a linear centered scheme.
 * \ingroup ConvDiscr
 * \author T. Aalbers
 */
class CPBConvection_Base : public CNumerics {
protected:

  bool implicit, dynamic_grid, energy;

  unsigned short iDim, jDim, iVar, jVar;
  
  su2double *AdvectedVelocity = nullptr;
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;

  su2double MeanDensity;

  /*!
   * \brief Function which defines the velocity that is advected
   */
  void virtual ComputeAdvectedVelocity(void) = 0;

   /*!
   * \brief Function which defines the Jacobian
   */
  void virtual ComputeJacobian(void) = 0;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPBConvection_Base(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CPBConvection_Base(void);
  
  /*!
   * \brief Compute the flow residual.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;
};

/*!
 * \class CPBConvection_Central
 * \brief Class for computing a centered scheme.
 * \ingroup ConvDiscr
 * \author T. Aalbers
 */
class CPBConvection_Central : public CPBConvection_Base {
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPBConvection_Central(unsigned short val_nDim, unsigned short val_nVar, CConfig *config)
    : CPBConvection_Base(val_nDim, val_nVar, config) {}

  /*!
   * \brief Function which defines the velocity that is advected
   */
  void ComputeAdvectedVelocity(void) final;

   /*!
   * \brief Function which defines the Jacobian
   */
  void ComputeJacobian(void) final;
  
};


/*!
 * \class CPBConvection_Upwind
 * \brief Class for computing an upwind scheme.
 * \ingroup ConvDiscr
 * \author T. Aalbers
 */
class CPBConvection_Upwind : public CPBConvection_Base {
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPBConvection_Upwind(unsigned short val_nDim, unsigned short val_nVar, CConfig *config)
    : CPBConvection_Base(val_nDim, val_nVar, config) {}

  /*!
   * \brief Function which defines the velocity that is advected
   */
  void ComputeAdvectedVelocity(void) final;

   /*!
   * \brief Function which defines the Jacobian
   */
  void ComputeJacobian(void) final;
  
};