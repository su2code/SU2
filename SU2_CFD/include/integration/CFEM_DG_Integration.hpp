/*!
 * \file CFEM_DG_Integration.hpp
 * \brief Declaration of class for integration with the FEM DG solver.
 * \author E. van der Weide, T. Economon
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

#include "CIntegration.hpp"

/*!
 * \class CFEM_DG_Integration
 * \ingroup Drivers
 * \brief Class for integration with the FEM DG solver.
 * \author E. van der Weide, T. Economon
 * \version 8.0.0 "Harrier"
 */
class CFEM_DG_Integration final : public CIntegration {
public:
  /*!
   * \brief Constructor of the class.
   */
  CFEM_DG_Integration();

  /*!
   * \brief Do the numerical integration (implicit) of the turbulence solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                            CNumerics ******numerics_container, CConfig **config,
                            unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) override;

private:
  /*!
   * \brief Perform the spatial integration of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] solver - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iStep - Current step of the Runge-Kutta iteration for the RK schemes
                        and the step in the local time stepping for ADER-DG.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void Space_Integration(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config,
                         unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem);

  /*!
   * \brief Perform the time integration (explicit or implicit) of the numerical system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iStep - Current step of the Runge-Kutta iteration for the RK schemes
                        and the step in the local time stepping for ADER-DG.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                        unsigned short iStep, unsigned short RunTime_EqSystem);
};
