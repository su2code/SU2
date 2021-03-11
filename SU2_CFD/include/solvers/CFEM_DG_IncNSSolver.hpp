/*!
 * \file CFEM_DG_IncNSSolver.hpp
 * \brief Headers of the CFEM_DG_IncNSSolver class
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFEM_DG_IncEulerSolver.hpp"

/*!
 * \class CFEM_DG_IncNSSolver
 * \brief Main class for defining the incompressible Navier-Stokes Discontinuous Galerkin finite element flow solver.
 * \ingroup Navier_Stokes_Equations
 * \author E. van der Weide, T. Economon, J. Alonso
 * \version 7.1.1 "Blackbird"
 */
class CFEM_DG_IncNSSolver final : public CFEM_DG_IncEulerSolver {

public:

  /*!
   * \brief Constructor of the class, disabled.
   */
  CFEM_DG_IncNSSolver(void) = delete;

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CFEM_DG_IncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEM_DG_IncNSSolver(void) override;

private:

};
