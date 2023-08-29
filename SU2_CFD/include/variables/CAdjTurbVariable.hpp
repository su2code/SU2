/*!
 * \file CAdjTurbVariable.hpp
 * \brief Main class for defining the variables of the adjoint turbulence model.
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

#include "CVariable.hpp"

/*!
 * \class CAdjTurbVariable
 * \brief Main class for defining the variables of the adjoint turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CAdjTurbVariable final : public CVariable {
protected:
  // ToDo: These variables were not being allocated... is this class used?
  MatrixType dmuT_dUTvar;         /*!< \brief Sensitivity of eddy viscosity to mean flow and turbulence vars. */
  CVectorOfMatrix dRTstar_dUTvar;  /*!< \brief Sensitivity of modified turbulence residual (no boundary flux) to mean flow and turbulence vars. */
  CVectorOfMatrix dFT_dUTvar;      /*!< \brief Sensitivity of boundary flux to mean flow and turbulence vars. */
  MatrixType EddyViscSens;        /*!< \brief Eddy Viscosity Sensitivity. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] psinu_inf - Value of the adjoint turbulence variable at the infinity (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjTurbVariable(su2double psinu_inf, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAdjTurbVariable() override = default;

  /*!
   * \brief Set the Eddy Viscosity Sensitivity of the problem.
   * \param[in] val_EddyViscSens - Eddy Viscosity Sensitivity.
   */
  inline void SetEddyViscSens(unsigned long iPoint, const su2double *val_EddyViscSens, unsigned long numTotalVar) override {
    for (unsigned long iVar = 0; iVar < numTotalVar; iVar++)
      EddyViscSens(iPoint,iVar) = val_EddyViscSens[iVar];
  }

  /*!
   * \brief Get the Eddy Viscosity Sensitivity of the problem.
   * \return Pointer to the Eddy Viscosity Sensitivity.
   */
  inline su2double *GetEddyViscSens(unsigned long iPoint) override { return EddyViscSens[iPoint]; }
};
