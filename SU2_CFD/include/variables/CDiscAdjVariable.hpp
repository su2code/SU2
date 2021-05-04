/*!
 * \file CDiscAdjVariable.hpp
 * \brief Main class for defining the variables of the adjoint solver.
 * \author F. Palacios, T. Economon
 * \version 7.1.1 "Blackbird"
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

#include "CVariable.hpp"

/*!
 * \class CDiscAdjVariable
 * \brief Main class for defining the variables of the adjoint solver.
 * \ingroup Discrete_Adjoint
 * \author T. Albring.
 */
class CDiscAdjVariable final : public CVariable {
private:
  MatrixType Sensitivity; /* Vector holding the derivative of target functional with respect to the coordinates at this node*/
  MatrixType Sensitivity_Old; /* Previous time sensitivity holder since inner iterations in FSI problems overwrite sensitivity*/
  MatrixType Solution_Direct;
  MatrixType DualTime_Derivative;
  MatrixType DualTime_Derivative_n;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] sol - Pointer to the adjoint value (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjVariable(const su2double* sol, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjVariable() override = default;

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) override { Sensitivity(iPoint,iDim) = val;}

  /*!
   * \brief Set the previous time sensitivity at the node
   * \param[in] iDim - dimension
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity_Old(unsigned long iPoint, unsigned long iDim, su2double val) override { Sensitivity_Old(iPoint,iDim) = val;}

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity(unsigned long iPoint, unsigned long iDim) const override { return Sensitivity(iPoint,iDim); }

  /*!
   * \brief Get the previous time sensitivity at the node
   * \param[in] iDim - dimension
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity_Old(unsigned long iPoint, unsigned long iDim) const override { return Sensitivity_Old(iPoint,iDim); }

  inline void SetDual_Time_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) override { DualTime_Derivative(iPoint,iVar) = der; }

  inline void SetDual_Time_Derivative_n(unsigned long iPoint, unsigned long iVar, su2double der) override { DualTime_Derivative_n(iPoint,iVar) = der; }

  inline su2double GetDual_Time_Derivative(unsigned long iPoint, unsigned long iVar) const override { return DualTime_Derivative(iPoint,iVar); }

  inline su2double GetDual_Time_Derivative_n(unsigned long iPoint, unsigned long iVar) const override { return DualTime_Derivative_n(iPoint,iVar); }

  inline void SetSolution_Direct(unsigned long iPoint, const su2double *val_solution_direct) override {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution_Direct(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline su2double* GetSolution_Direct(unsigned long iPoint) override { return Solution_Direct[iPoint]; }

};
