/*!
 * \file CDiscAdjFEAVariable.hpp
 * \brief Main class for defining the variables of the adjoint FEA solver.
 * \author T. Albring, R. Sanchez.
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
 * \class CDiscAdjFEAVariable
 * \brief Main class for defining the variables of the adjoint solver.
 * \ingroup Discrete_Adjoint
 * \author T. Albring, R. Sanchez.
 * \version 7.1.1 "Blackbird"
 */
class CDiscAdjFEAVariable : public CVariable {
protected:
  MatrixType Sensitivity; /* Vector holding the derivative of target functional with respect to the coordinates at this node*/
  MatrixType Sensitivity_Old; /* Previous time sensitivity holder since inner iterations in FSI problems overwrite sensitivity*/
  MatrixType Solution_Direct;

  MatrixType Dynamic_Derivative;
  MatrixType Dynamic_Derivative_n;

  MatrixType Solution_Direct_Vel;
  MatrixType Solution_Direct_Accel;

  /*!
   * \brief Constructor of the class.
   * \param[in] disp - Pointer to the adjoint value (initialization value).
   * \param[in] vel - Pointer to the adjoint value (initialization value).
   * \param[in] accel - Pointer to the adjoint value (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] unsteady - Allocate velocity and acceleration.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFEAVariable(const su2double *disp, const su2double *vel, const su2double *accel,
                      unsigned long npoint, unsigned long ndim, unsigned long nvar, bool unsteady, CConfig *config);

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEAVariable() override = default;

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) final { Sensitivity(iPoint,iDim) = val; }

  /*!
   * \brief Set the previous time sensitivity at the node
   * \param[in] iDim - dimension
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity_Old(unsigned long iPoint, unsigned long iDim, su2double val) final { Sensitivity_Old(iPoint,iDim) = val; }

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity(unsigned long iPoint, unsigned long iDim) const final { return Sensitivity(iPoint,iDim);}

  /*!
   * \brief Get the previous time sensitivity at the node
   * \param[in] iDim - dimension
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity_Old(unsigned long iPoint, unsigned long iDim) const final { return Sensitivity_Old(iPoint,iDim);}

  inline void SetDynamic_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative(iPoint,iVar) = der;
  }

  inline void SetDynamic_Derivative_n(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative_n(iPoint,iVar) = der;
  }

  inline su2double GetDynamic_Derivative(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative(iPoint,iVar);
  }

  inline su2double GetDynamic_Derivative_n(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative_n(iPoint,iVar);
  }

  inline void SetSolution_Direct(unsigned long iPoint, const su2double *val_solution_direct) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Direct(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline void SetSolution_Vel_Direct(unsigned long iPoint, const su2double *val_solution_direct) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Direct_Vel(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline void SetSolution_Accel_Direct(unsigned long iPoint, const su2double *val_solution_direct) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Direct_Accel(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline su2double* GetSolution_Direct(unsigned long iPoint) final { return Solution_Direct[iPoint]; }

  inline su2double* GetSolution_Vel_Direct(unsigned long iPoint) final { return Solution_Direct_Vel[iPoint]; }

  inline su2double* GetSolution_Accel_Direct(unsigned long iPoint) final { return Solution_Direct_Accel[iPoint]; }

};
