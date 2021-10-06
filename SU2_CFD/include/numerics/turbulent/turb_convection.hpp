/*!
 * \file turb_convection.hpp
 * \brief Delarations of numerics classes for discretization of
 *        convective fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../scalar/scalar_convection.hpp"

/*!
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaras turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 */
template <class FlowIndices>
class CUpwSca_TurbSA_generic final : public CUpwScalar<FlowIndices> {
private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override;

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

};
using CUpwSca_TurbSA = CUpwSca_TurbSA_generic<CEulerVariable::Indices>;
using CUpwSca_IncTurbSA = CUpwSca_TurbSA_generic<CIncEulerVariable::Indices>;
using CUpwSca_NEMOTurbSA = CUpwSca_TurbSA_generic<CNEMOEulerVariable::Indices>;

/*!
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 */
template <class FlowIndices>
class CUpwSca_TurbSST_generic final : public CUpwScalar<FlowIndices> {
private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override;

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

};
using CUpwSca_TurbSST = CUpwSca_TurbSST_generic<CEulerVariable::Indices>;
using CUpwSca_IncTurbSST = CUpwSca_TurbSST_generic<CIncEulerVariable::Indices>;
using CUpwSca_NEMOTurbSST = CUpwSca_TurbSST_generic<CNEMOEulerVariable::Indices>;