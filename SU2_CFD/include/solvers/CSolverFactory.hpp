/*!
 * \file CSolverFactory.hpp
 * \brief Headers of the CSolverFactory class
 * \author T. Albring
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

#include "../../../Common/include/option_structure.hpp"

/*!
 * \brief Enum of different sub solvers the main solver can use. There is not a 1-to-1 correspondence between the actual classes
 * and the types of sub solvers, as one class can be used for several sub solvers.
 */
enum class SUB_SOLVER_TYPE {
  CONT_ADJ_EULER,          /*!< \brief Continuous Adjoint Euler solver  */
  CONT_ADJ_NAVIER_STOKES,  /*!< \brief Continuous Adjoint Navier Stokes solver  */
  CONT_ADJ_TURB,           /*!< \brief Continuous Adjoint Turbulent solver  */
  BASELINE,                /*!< \brief Baseline solver  */
  TEMPLATE,                /*!< \brief Template solver  */
  BASELINE_FEM,            /*!< \brief Baseline FEM solver */
  DISC_ADJ_FEA,            /*!< \brief Discrete adjoint FEA solver  */
  DISC_ADJ_MESH,           /*!< \brief Discrete adjoint mesh solver */
  DISC_ADJ_FLOW,           /*!< \brief Discrete adjoint flow solver */
  DISC_ADJ_TURB,           /*!< \brief Discrete adjoint turbulence solver */
  DISC_ADJ_SPECIES,        /*!< \brief Discrete adjoint species solver */
  DISC_ADJ_HEAT,           /*!< \brief Discrete adjoint heat solver */
  EULER,                   /*!< \brief Compressible Euler solver */
  NAVIER_STOKES,           /*!< \brief Compressible Navier-Stokes solver */
  NEMO_EULER,              /*!< \brief NEMO Euler solver */
  NEMO_NAVIER_STOKES,      /*!< \brief NEMO Navier-Stokes solver */
  INC_EULER,               /*!< \brief Incompressible Euler solver */
  INC_NAVIER_STOKES,       /*!< \brief Incompressible Navier-stokes solver */
  FEA,                     /*!< \brief Structural Finite-Element solver */
  DG_EULER,                /*!< \brief Higher-order DG Euler solver*/
  DG_NAVIER_STOKES,        /*!< \brief Higher-order DG Navier-Stokes solver*/
  HEAT,                    /*!< \brief Heat solver */
  TRANSITION,              /*!< \brief Transition model solver*/
  TURB_SA,                 /*!< \brief SA turbulence model solver */
  TURB_SST,                /*!< \brief SST turbulence model solver */
  TURB,                    /*!< \brief Turbulence model solver */
  SPECIES,                 /*!< \brief Species model solver */
  MESH,                    /*!< \brief Mesh solver */
  RADIATION,               /*!< \brief Radiation solver */
  DISC_ADJ_RADIATION,      /*!< \brief Discrete adjoint radiation solver */
  NONE
};

enum class INTEGRATION_TYPE{
  MULTIGRID,
  NEWTON,
  SINGLEGRID,
  DEFAULT,
  FEM_DG,
  STRUCTURAL,
  NONE
};

struct SolverMetaData{
  SUB_SOLVER_TYPE  solverType        = SUB_SOLVER_TYPE::NONE;
  INTEGRATION_TYPE integrationType   = INTEGRATION_TYPE::NONE;
};

class CSolver;
class CGeometry;
class CConfig;

class CSolverFactory {

private:

  static std::map<const CSolver*, SolverMetaData> allocatedSolvers;

  /*!
   * \brief Create a turbulent solver
   * \param[in] kindTurbModel - Kind of turbulent solver
   * \param[in] solver        - The solver container (used to call preprocessing of the flow solver)
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return                  - A pointer to the allocated turbulent solver
   */
  static CSolver* CreateTurbSolver(TURB_MODEL kindTurbModel, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, int adjoint);

  /*!
   * \brief Create a transition solver
   * \param[in] kindTransModel - Kind of transition solver
   * \param[in] solver        - The solver container (used to call preprocessing of the flow solver)
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return                  - A pointer to the allocated transition solver
   */
  static CSolver* CreateTransSolver(TURB_TRANS_MODEL kindTransModel , CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, int adjoint);

  /*!
   * \brief Create a species solver
   * \param[in] solver        - The solver container
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return                  - A pointer to the allocated species solver
   */
  static CSolver* CreateSpeciesSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint);

  /*!
   * \brief Create a heat solver
   * \param[in] solver        - The solver container
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return                  - A pointer to the allocated heat solver
   */
  static CSolver* CreateHeatSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint);

  /*!
   * \brief Create a mesh solver
   * \param[in] solver        - The solver container
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \param[in] adjoint       - Boolean indicating whether a primal or adjoint solver should be allocated
   * \return                  - A pointer to the allocated mesh solver
   */
  static CSolver* CreateMeshSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint);

  /*!
   * \brief Create a DG solver
   * \param[in] kindTurbModel - Kind of DG solver
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \return                  - A pointer to the allocated DG solver
   */
  static CSolver* CreateDGSolver(SUB_SOLVER_TYPE kindDGSolver, CGeometry *geometry, CConfig *config, int iMGLevel);

  /*!
   * \brief Create a flow solver
   * \param[in] kindFlowSolver - Kind of flow solver
   * \param[in] solver         - The solver container
   * \param[in] geometry       - The geometry definition
   * \param[in] config         - The configuration
   * \param[in] iMGLevel       - The multigrid level
   * \return                   - A pointer to the allocated flow solver
   */
  static CSolver* CreateFlowSolver(SUB_SOLVER_TYPE kindFlowSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel);

  /*!
   * \brief Generic routine to create a solver
   * \param[in] kindSolver    - Kind of solver
   * \param[in] solver        - The solver container
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \return                  - A pointer to the allocated solver
   */
  static CSolver* CreateSubSolver(SUB_SOLVER_TYPE kindSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel);

public:

  /*!
   * \brief Deleted constructor to avoid creating instances of this class
   */
  CSolverFactory() = delete;

  /*!
   * \brief Create the solver container by allocating the primary solver
   * and secondary solvers like heat solver, turbulent solver etc
   * \param[in] kindSolver    - The kind of primary solver
   * \param[in] config        - The configuration
   * \param[in] geometry      - The geometry definition
   * \param[in] iMGLevel      - The multigrid level
   * \return                  - Pointer to the allocated solver array
   */
  static CSolver** CreateSolverContainer(MAIN_SOLVER kindSolver, CConfig *config, CGeometry *geometry, int iMGLevel);


  /*!
   * \brief Return a sub solver object that contains information about the solver allocated at a specific memory address
   * \param[in] solver - Address of the solver
   * \return sub solver info struct.
   */
  static SolverMetaData GetSolverMeta(const CSolver* solver) { return allocatedSolvers.at(solver); }

  /*!
   * \brief Clear the solver meta data
   */
  static void ClearSolverMeta() { allocatedSolvers.clear(); }

};
