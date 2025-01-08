/*!
 * \file CDiscAdjVariable.hpp
 * \brief Main class for defining the variables of the adjoint solver.
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
 * \class CDiscAdjVariable
 * \ingroup DiscAdj
 * \brief Main class for defining the variables of the adjoint solver.
 * \author T. Albring.
 */
class CDiscAdjVariable : public CVariable {
private:
  MatrixType Sensitivity; /*!< \brief Vector holding the derivative of target functional with respect to the coordinates at this node. */
  MatrixType Solution_Direct; /*!< \brief Stores the primal solution of the current timestep in order to be able to reset. */
  MatrixType DualTime_Derivative; /*!< \brief Container holding all/sum-of dual time contributions to the adjoint variable. */
  MatrixType DualTime_Derivative_n; /*!< \brief Container holding dual time contributions to the adjoint variable used in the next timestep. */

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
   * \brief Allocate extra adjoint variables.
   * \param[in] nVarExtra - Number of extra variables.
   */
  void AllocateAdjointSolutionExtra(unsigned long nVarExtra);

  /*!
   * \brief Set the sensitivity at the node
   */
  inline void SetSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) final {
    Sensitivity(iPoint,iDim) = val;
  }

  /*!
   * \brief Get the Sensitivity at the node
   */
  inline su2double GetSensitivity(unsigned long iPoint, unsigned long iDim) const final {
    return Sensitivity(iPoint,iDim);
  }
  inline const MatrixType& GetSensitivity() const final { return Sensitivity; }

  /*!
   * \brief Set/store the dual time contributions to the adjoint variable.
   *        Contains sum of contributions from 2 timesteps for dual time 2nd order.
   */
  inline void SetDual_Time_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) {
    DualTime_Derivative(iPoint,iVar) = der;
  }

  /*!
   * \brief Set/store the dual time contributions to the adjoint variable for upcoming timestep.
   */
  inline void SetDual_Time_Derivative_n(unsigned long iPoint, unsigned long iVar, su2double der) {
    DualTime_Derivative_n(iPoint,iVar) = der;
  }

  /*!
   * \brief Return the dual time contributions to the adjoint variable.
   *        Contains sum of contributions from 2 timesteps for dual time 2nd order.
   */
  inline su2double GetDual_Time_Derivative(unsigned long iPoint, unsigned long iVar) const {
    return DualTime_Derivative(iPoint,iVar);
  }

  /*!
   * \brief Return the dual time contributions to the adjoint variable for upcoming timestep.
   */
  inline su2double GetDual_Time_Derivative_n(unsigned long iPoint, unsigned long iVar) const {
    return DualTime_Derivative_n(iPoint,iVar);
  }

  /*!
   * \brief Set/store the primal solution for all variables of one point.
   */
  inline void SetSolution_Direct(unsigned long iPoint, const su2double *val_solution_direct) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution_Direct(iPoint,iVar) = val_solution_direct[iVar];
  }

  /*!
   * \brief Returns the primal solution for all variables of one point.
   */
  inline su2double* GetSolution_Direct(unsigned long iPoint) final {
    return Solution_Direct[iPoint];
  }

  /*!
   * \brief Set Dual-time derivative contributions to the external.
   */
  void Set_External_To_DualTimeDer() final;

};
