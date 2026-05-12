/*!
 * \file CPBIncEulerSolver.hpp
 * \brief Headers of the CPBIncEulerSolver class
 * \author F. Palacios, T. Economon, T. Albring
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

#include "CFVMFlowSolverBase.hpp"
#include "../variables/CPBIncEulerVariable.hpp"

class CPBIncEulerSolver : public CFVMFlowSolverBase<CPBIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE> {
protected:

public:
  /*!
   * \brief Constructor of the class.
   */
  //CPBIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) {};
  
  CPBIncEulerSolver() = delete;

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Grid level.
   * \param[in] navier_stokes - True when the constructor is called by the derived class CPBIncNSSolver.
   */
  CPBIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, const bool navier_stokes = false);

  /*!
   * \brief Destructor of the class.
   */
  ~CPBIncEulerSolver(void) override {}

  /*!
   * \brief Set reference values for pressure, forces, etc.
   */
  void SetReferenceValues(const CConfig& config) final {}

  /*!
   * \brief Print verification error to screen.
   * \param[in] config - Definition of the particular problem.
   */
  void PrintVerificationError(const CConfig* config) const final {}



};
