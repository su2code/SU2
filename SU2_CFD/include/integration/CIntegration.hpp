/*!
 * \file CIntegration.hpp
 * \brief Declaration of the main routines to orchestrate space and time integration.
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

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../solvers/CSolver.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/CConfig.hpp"

using namespace std;

/*!
 * \class CIntegration
 * \ingroup Drivers
 * \brief Main class for doing the space integration, time integration, and monitoring
 *        of a system of Partial Differential Equations (PDE).
 * \author F. Palacios
 */
class CIntegration {
protected:
  int rank,      /*!< \brief MPI Rank. */
  size;          /*!< \brief MPI Size. */

  /*!
   * \brief Do the space integration of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Space_Integration(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics,
                         CConfig *config, unsigned short iMesh, unsigned short iRKStep,
                         unsigned short RunTime_EqSystem);

  /*!
   * \brief Do the time integration (explicit or implicit) of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iRKStep, unsigned short RunTime_EqSystem);

public:
  /*!
   * \brief Constructor of the class.
   */
  CIntegration();

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIntegration(void) = default;

  /*!
   * \brief Save the geometry at different time steps.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Mesh solver.
   * \param[in] config - Definition of the particular problem.
   */
  void SetDualTime_Geometry(CGeometry *geometry, CSolver *mesh_solver, const CConfig *config, unsigned short iMesh);

  /*!
   * \brief Save the solution at different time steps, and reset certain fields for the next timestep.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Some solver.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetDualTime_Solver(const CGeometry *geometry, CSolver *solver, const CConfig *config, unsigned short iMesh);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  virtual void MultiGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                   CNumerics ******numerics_container, CConfig **config,
                                   unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) { };

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  virtual void SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                    CNumerics ******numerics_container, CConfig **config,
                                    unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) { };

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  virtual void Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                    CNumerics ******numerics_container, CConfig **config,
                                    unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) { };

};
