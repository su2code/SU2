/*!
 * \file CMeshSolver.hpp
 * \brief Declaration and inlines of the class to compute the
 *        the discrete adjoint of the linear-elastic mesh solver.
 * \author Ruben Sanchez
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

#include "CSolver.hpp"
#include "../variables/CDiscAdjMeshBoundVariable.hpp"

/*!
 * \class CDiscAdjMeshSolver
 * \ingroup DiscAdj
 * \brief Main class for defining the discrete adjoint solver for mesh deformation problems.
 * \author R. Sanchez
 */
class CDiscAdjMeshSolver final : public CSolver {
private:
  static constexpr size_t MAXNDIM = 3;  /*!< \brief Max number of space dimensions, used in some static arrays. */
  static constexpr size_t MAXNVAR = 3;  /*!< \brief Max number of variables, for static arrays. */

  static constexpr size_t OMP_MAX_SIZE = 1024; /*!< \brief Max chunk size for light point loops. */

  unsigned long omp_chunk_size; /*!< \brief Chunk size used in light point loops. */

  CSolver *direct_solver = nullptr;

  CDiscAdjMeshBoundVariable* nodes = nullptr;   /*!< \brief Variables of the discrete adjoint mesh solver. */

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjMeshSolver() = default;

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
   */
  CDiscAdjMeshSolver(CGeometry *geometry, CConfig *config, CSolver* solver);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjMeshSolver() override;

  /*!
   * \brief Performs the preprocessing of the AD-based mesh adjoint solver.
   *        Registers all necessary variables on the tape. Called while tape is active.
   * \param[in] geometry_container - The geometry container holding all grid levels.
   * \param[in] config_container - The particular config.
   */
  void RegisterSolution(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
   *        after tape has been evaluated.
   * \param[in] geometry - The geometrical definition of the problem.
   * \param[in] config - The particular config.
   */
  void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config, bool CrossTerm) override;

  /*!
   * \brief Extract and set the geometrical sensitivity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] target_solver - The target solver to store the sensitivities.
   */
  void SetSensitivity(CGeometry *geometry, CConfig *config, CSolver* target_solver) override;

  /*!
   * \brief Prepare the solver for a new recording.
   * \param[in] kind_recording - Kind of AD recording.
   */
  void SetRecording(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] reset - If true reset variables to their initial values.
   */
  void RegisterVariables(CGeometry *geometry,
                         CConfig *config,
                         bool reset = false) override;

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry,
                   CSolver ***solver,
                   CConfig *config,
                   int val_iter,
                   bool val_update_geo) override;

};
