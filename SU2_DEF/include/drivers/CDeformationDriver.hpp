 /*!
 * \file CDeformationDriver.hpp
 * \brief Headers of the main subroutines for driving the mesh deformation.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 7.1.1 "Blackbird"
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

#include "../../../Common/include/parallelization/mpi_structure.hpp"

#include "../../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../../Common/include/grid_movement/CVolumetricMovement.hpp"
#include "../../../SU2_CFD/include/output/COutput.hpp"
#include "../../../SU2_CFD/include/numerics/CNumerics.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"

/*!
 * \class CDeformationDriver
 * \brief Class for driving mesh deformation solvers.
 * \author A. Gastaldi, H. Patel
 * \version 7.1.1 "Blackbird"
 */
class CDeformationDriver {
protected:
  char config_file_name[MAX_STRING_SIZE];
  int rank,
      size;
  su2double StartTime,                          /*!< \brief Start point of the timer for performance benchmarking.*/
            StopTime,                           /*!< \brief Stop point of the timer for performance benchmarking.*/
            UsedTimePreproc,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking preprocessing phase.*/
            UsedTimeCompute,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking compute phase.*/
            UsedTime;                           /*!< \brief Elapsed time between Start and Stop point of the timer.*/
  unsigned short iZone, nZone = SINGLE_ZONE;
  CConfig *driver_config;                       /*!< \brief Definition of the driver configuration. */
  CConfig **config_container;                   /*!< \brief Definition of the particular problem. */
  CGeometry **geometry_container;             /*!< \brief Geometrical definition of the problem. */
  CSurfaceMovement **surface_movement;          /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement **grid_movement;         /*!< \brief Volume grid movement classes of the problem. */
  CSolver **solver_container;
  CNumerics ***numerics_container;
  COutput **output_container;                   /*!< \brief Pointer to the COutput class. */

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
   * \brief [Overload] Launch the computation for single-zone problems.
   */
  void Run();

  /*!
   * \brief Output the mesh.
   */
  void Output();

  /*!
   * \brief Deallocation routine
   */
  void Postprocessing();

  /*!
   * \brief Get all the deformable boundary marker tags.
   * \return List of deformable boundary markers tags.
   */
  vector<string> GetAllDeformMeshMarkersTag() const;

  /*!
   * \brief Get all the boundary markers tags with their associated indices.
   * \return List of boundary markers tags with their indices.
   */
  map<string, int> GetAllBoundaryMarkers() const;

  /*!
   * \brief Get all the boundary markers tags with their associated types.
   * \return List of boundary markers tags with their types.
   */
  map<string, string> GetAllBoundaryMarkersType() const;

  /*!
   * \brief Get the number of vertices (halo nodes included) from a specified marker.
   * \param[in] iMarker -  Marker identifier.
   * \return Number of vertices.
   */
  unsigned long GetNumberVertices(unsigned short iMarker) const;

  /*!
   * \brief Get the number of halo vertices from a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \return Number of vertices.
   */
  unsigned long GetNumberHaloVertices(unsigned short iMarker) const;

  /*!
   * \brief Check if a vertex is physical or not (halo node) on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return True if the specified vertex is a halo node.
   */
  bool IsAHaloNode(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the global index of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vertex global index.
   */
  unsigned long GetVertexGlobalIndex(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get undeformed coordinates from the mesh solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x,y,z coordinates of the vertex.
   */
  vector<passivedouble> GetInitialMeshCoord(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the unit normal (vector) at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Unit normal (vector) at the vertex.
   */
  vector<passivedouble> GetVertexNormal(unsigned short iMarker, unsigned long iVertex, bool unitNormal = false) const;

  inline vector<passivedouble> GetVertexUnitNormal(unsigned short iMarker, unsigned long iVertex) const {
    return GetVertexNormal(iMarker, iVertex, true);
  }
  
  /*!
   * \brief Set the mesh displacement for the elasticity mesh solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] DispX - Value of the mesh displacement in the direction X.
   * \param[in] DispY - Value of the mesh displacement in the direction Y.
   * \param[in] DispZ - Value of the mesh displacement in the direction Z.
   */
  void SetMeshDisplacement(unsigned short iMarker, unsigned long iVertex, passivedouble DispX, passivedouble DispY, passivedouble DispZ);

  /*!
   * \brief Communicate the boundary mesh displacements in a python call
   */
  void CommunicateMeshDisplacement(void);

protected:
  /*!
   * \brief Init_Containers
   */
  void SetContainers_Null();

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
