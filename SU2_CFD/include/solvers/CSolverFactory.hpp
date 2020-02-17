/*!
 * \file CSolverFactory.hpp
 * \brief Headers of the CSolverFactory class
 * \author T. Albring
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "CSolver.hpp"

class CSolverFactory {
  
private:
  /*!
   * \brief Create a turbulent solver
   * \param[in] kindTurbModel - Kind of turbulent solver
   * \param[in] solver        - The solver container (used to call preprocessing of the flow solver)
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return              - A pointer to the allocated turbulent solver
   */
  static CSolver* createTurbSolver(ENUM_TURB_MODEL kindTurbModel, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, int adjoint);
  
  /*!
   * \brief Create a heat solver 
   * \param[in] solver        - The solver container 
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return              - A pointer to the allocated heat solver
   */
  static CSolver* createHeatSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint);
  
  /*!
   * \brief Create a mesh solver 
   * \param[in] solver        - The solver container 
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return              - A pointer to the allocated mesh solver
   */
  static CSolver* createMeshSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint);
  
  /*!
   * \brief Create a DG solver 
   * \param[in] kindTurbModel - Kind of DG solver
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \return              - A pointer to the allocated DG solver
   */
  static CSolver* createDGSolver(ENUM_SOLVER kindDGSolver, CGeometry *geometry, CConfig *config, int iMGLevel);
  
  /*!
   * \brief Create a flow solver 
   * \param[in] kindFlowSolver - Kind of flow solver
   * \param[in] solver        - The solver container 
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \return              - A pointer to the allocated flow solver
   */
  static CSolver* createFlowSolver(ENUM_SOLVER kindFlowSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel);
  
  /*!
   * \brief Generic routine to create a solver 
   * \param[in] kindSolver - Kind of solver
   * \param[in] solver        - The solver container 
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \return              - A pointer to the allocated solver
   */
  static CSolver* createSolver(ENUM_SOLVER kindSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel);
  
public:
    
  /*!
   * \brief Create the solver container by allocating the primary solver 
   * and secondary solvers like heat solver, turbulent solver etc
   * \param[in] kindSolver - The kind of primary solver 
   * \param[in] config        - The configuration
   * \param[in] geometry      - The geometry definition
   * \param[in] iMGLevel      - The multigrid level
   * \return                  - Pointer to the allocated solver array
   */
  static CSolver** createSolverContainer(ENUM_SOLVER kindSolver, CConfig *config, CGeometry *geometry, int iMGLevel);

};
