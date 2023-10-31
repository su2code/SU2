/*!
 * \file CAdjNSVariable.hpp
 * \brief Main class for defining the variables of the adjoint Navier-Stokes solver.
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

#include "CAdjEulerVariable.hpp"

/*!
 * \class CAdjNSVariable
 * \brief Main class for defining the variables of the adjoint Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios
 */
class CAdjNSVariable final : public CAdjEulerVariable {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] psirho - Value of the adjoint density (initialization value).
   * \param[in] phi - Value of the adjoint velocity (initialization value).
   * \param[in] psie - Value of the adjoint energy (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjNSVariable(su2double psirho, const su2double *phi, su2double psie,
                 unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAdjNSVariable() override = default;

  /*!
   * \brief Set the value of the force projection vector on the old solution vector.
   */
  inline void SetVelSolutionDVector(unsigned long iPoint) override {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Solution(iPoint,iDim+1) = ForceProj_Vector(iPoint,iDim);
  }

};
