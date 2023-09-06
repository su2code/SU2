/*!
 * \file CFEASolver.hpp
 * \brief Base class template for all FEA solvers using the SU2 internal finite elements.
 * \author T. Dick
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

#include <unordered_set>

#include "CSolver.hpp"
#include "../../../Common/include/geometry/elements/CElement.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"

/*!
 * \class CFEASolverBase
 * \ingroup Elasticity_Equations
 * \brief Base class for FEM elasticity solvers.
 */
class CFEASolverBase : public CSolver {
 public:
  enum : size_t {MAXNNODE_2D = 4};
  enum : size_t {MAXNNODE_3D = 8};
  enum : size_t {MAXNDIM = 3};
  enum : size_t {MAXNVAR = 3};
  enum : size_t {OMP_MIN_SIZE = 32};
  enum : size_t {OMP_MAX_SIZE = 512};

 protected:

  CElement*** element_container = nullptr;  /*!< \brief Vector which the define the finite element structure for each problem. */
  unsigned long nElement;                   /*!< \brief Number of elements. */

  vector<unsigned long> ExtraVerticesToEliminate; /*!< Extra vertices for row/column elimination, see CommunicateExtraEliminationVertices. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] mesh_deform_mode - mode of operation for the linear solver.
   */
  CFEASolverBase(LINEAR_SOLVER_MODE mesh_deform_mode = LINEAR_SOLVER_MODE::STANDARD);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] mesh_deform_mode - mode of operation for the linear solver.
   */
  CFEASolverBase(CGeometry *geometry, CConfig *config, LINEAR_SOLVER_MODE mesh_deform_mode = LINEAR_SOLVER_MODE::STANDARD);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEASolverBase();

  /*!
   * \brief Communicate extra vertices for elimination in the linear system.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] myPoints - List of global point indeces to eliminate.
   * \param[in] val_coord - Location (x, y, z) of the max residual point.
   */
  void CommunicateExtraEliminationVertices(const CGeometry* geometry, vector<unsigned long>& myPoints);

};
