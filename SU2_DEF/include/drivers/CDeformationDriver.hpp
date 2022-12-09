/*!
 * \file CDeformationDriver.hpp
 * \brief Headers of the main subroutines for driving the mesh deformation.
 * \author A. Gastaldi, H. Patel
 * \version 7.4.0 "Blackbird"
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

#define ENABLE_MAPS
#include "../../../Common/include/CConfig.hpp"
#undef ENABLE_MAPS

#include "../../../Common/include/drivers/CDriverBase.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../../Common/include/grid_movement/CVolumetricMovement.hpp"
#include "../../../Common/include/parallelization/mpi_structure.hpp"
#include "../../../SU2_CFD/include/numerics/CNumerics.hpp"
#include "../../../SU2_CFD/include/output/COutput.hpp"

class CDeformationDriver : public CDriverBase {
 protected:
  bool haveSurfaceDeformation = false;  // flag used to determine whether surface deformation is available for output

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDeformationDriver(char* confFile, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDeformationDriver(void);

  /*!
   * \brief Preprocess the driver data.
   */
  void Preprocess();

  /*!
   * \brief Launch the driver computation.
   */
  void Run();

  /*!
   * \brief Output the mesh.
   */
  void Output();

  /*!
   * \brief Deallocation routine.
   */
  void Postprocessing();

  /*!
   * \brief Communicate boundary mesh displacements.
   */
  void CommunicateMeshDisplacements(void);

 protected:
  /*!
   * \brief Read in the config and mesh files.
   */
  void Input_Preprocessing();

  /*!
   * \brief Construction of the edge-based data structure.
   */
  void Geometrical_Preprocessing();

  /*!
   * \brief Preprocess the output container.
   */
  void Output_Preprocessing();

  /*!
   * \brief Preprocess the mesh solver container.
   */
  void Solver_Preprocessing();

  /*!
   * \brief Preprocess the numerics container.
   */
  void Numerics_Preprocessing();

  /*!
   * \brief Mesh deformation based on linear elasticity solver (CMeshSolver).
   */
  void Update();

  /*!
   * \brief Mesh deformation based on legacy implementation.
   */
  void Update_Legacy();
};
