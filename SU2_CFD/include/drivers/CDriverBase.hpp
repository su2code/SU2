/*!
 * \file CDriverBase.hpp
 * \brief Base class for all drivers.
 * \author H. Patel, A. Gastaldi
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/CConfig.hpp"
#include "../numerics/CNumerics.hpp"
#include "../output/COutput.hpp"
#include "../solvers/CSolver.hpp"

/*!
 * \class CDriverBase
 * \ingroup Drivers
 * \brief Base class for all drivers.
 * \author H. Patel, A. Gastaldi
 */

class CDriverBase {
 protected:
  int rank,               /*!< \brief MPI Rank. */
      size;               /*!< \brief MPI Size. */
  char* config_file_name; /*!< \brief Configuration file name of the problem. */

  su2double StartTime, /*!< \brief Start point of the timer for performance benchmarking. */
      StopTime,        /*!< \brief Stop point of the timer for performance benchmarking. */
      UsedTimePreproc, /*!< \brief Elapsed time between Start and Stop point of the timer for tracking preprocessing
                                    phase. */
      UsedTimeCompute, /*!< \brief Elapsed time between Start and Stop point of the timer for tracking compute phase. */
      UsedTime;        /*!< \brief Elapsed time between Start and Stop point of the timer. */

  unsigned long TimeIter;

  unsigned short iMesh,  /*!< \brief Iterator on mesh levels. */
      iZone,             /*!< \brief Iterator on zones. */
      nZone,             /*!< \brief Total number of zones in the problem. */
      nDim,              /*!< \brief Number of dimensions. */
      iInst,             /*!< \brief Iterator on instance levels. */
      *nInst,            /*!< \brief Total number of instances in the problem (per zone). */
      **interface_types; /*!< \brief Type of coupling between the distinct (physical) zones. */

  CConfig* driver_config; /*!< \brief Definition of the driver configuration. */
  COutput* driver_output; /*!< \brief Definition of the driver output. */

  CConfig** config_container;           /*!< \brief Definition of the particular problem. */
  COutput** output_container;           /*!< \brief Pointer to the COutput class. */
  CGeometry**** geometry_container;     /*!< \brief Geometrical definition of the problem. */
  CSolver***** solver_container;        /*!< \brief Container vector with all the solutions. */
  CNumerics****** numerics_container;   /*!< \brief Description of the numerical method (the way in which the equations
                                           are solved). */
  CSurfaceMovement** surface_movement;  /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement*** grid_movement; /*!< \brief Volume grid movement classes of the problem. */
  CFreeFormDefBox*** FFDBox;            /*!< \brief FFD FFDBoxes of the problem. */

  CConfig* main_config;     /*!< \brief Reference to base (i.e. ZONE 0) configuration (used in driver API). */
  CGeometry* main_geometry; /*!< \brief Reference to base (i.e. ZONE, INST, MESH 0) geometry (used in driver API). */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDriverBase(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CDriverBase(void);

  /*!
   * \brief A virtual member.
   */
  virtual void Preprocessing(){};

  /*!
   * \brief A virtual member.
   */
  virtual void Run(){};

  /*!
   * \brief A virtual member.
   */
  virtual void Update(){};

  /*!
   * \brief A virtual member.
   */
  virtual void Update_Legacy(){};

  /*!
   * \brief A virtual member.
   */
  virtual void Output(){};

  /*!
   * \brief A virtual member.
   */
  virtual void Postprocessing(){};

  /*!
   * \brief Get the number of design variables.
   * \return Number of design variables.
   */
  unsigned short GetNumberDesignVariables() const;

  /*!
   * \brief Get the number of FFD boxes.
   * \return Number of FFD boxes.
   */
  unsigned short GetNumberFFDBoxes() const;

  /*!
   * \brief Get the number of markers in the mesh.
   * \return Number of markers.
   */
  unsigned short GetNumberMarkers() const;

  /*!
   * \brief Get all the boundary markers tags with their associated indices.
   * \return List of boundary markers tags with their indices.
   */
  map<string, unsigned short> GetMarkerIndices() const;

  /*!
   * \brief Get all the boundary markers tags with their associated types.
   * \return List of boundary markers tags with their types.
   */
  map<string, string> GetMarkerTypes() const;

  /*!
   * \brief Get all the boundary marker tags.
   * \return List of boundary markers tags.
   */
  vector<string> GetMarkerTags() const;

  /*!
   * \brief Get all the deformable boundary marker tags.
   * \return List of deformable boundary markers tags.
   */
  vector<string> GetDeformableMarkerTags() const;

  /*!
   * \brief Get the number of dimensions of the mesh.
   * \return Number of dimensions.
   */
  unsigned long GetNumberDimensions() const;

  /*!
   * \brief Get the number of elements in the mesh.
   * \return Number of elements.
   */
  unsigned long GetNumberElements() const;

  /*!
   * \brief Get the number of elements in the marker.
   * \param[in] iMarker - Marker index.
   * \return Number of elements.
   */
  unsigned long GetNumberMarkerElements(unsigned short iMarker) const;

  /*!
   * \brief Get the global indices of the mesh elements.
   * \return Global element indices (nElem).
   */
  vector<unsigned long> GetElements() const;

  /*!
   * \brief Get the global index of a mesh element.
   * \param[in] iElem - Mesh element index.
   * \return Global element index.
   */
  unsigned long GetElements(unsigned long iElem) const;

  /*!
   * \brief Get the global indices of the marker elements.
   * \param[in] iMarker - Marker index.
   * \return Global element indices (nElem).
   */
  vector<unsigned long> GetMarkerElements(unsigned short iMarker) const;

  /*!
   * \brief Get the global index of a marker element.
   * \param[in] iMarker - Marker index.
   * \param[in] iElem - Marker element index.
   * \return Global element index.
   */
  unsigned long GetMarkerElements(unsigned short iMarker, unsigned long iElem) const;

  /*!
   * \brief Get the global node indices of the mesh elements.
   * \return Element global node indices (nElem, nNode).
   */
  vector<vector<unsigned long>> GetElementNodes() const;

  /*!
   * \brief Get the global node indices of a mesh element.
   * \param[in] iElem - Mesh element index.
   * \return Element global node indices (nNode).
   */
  vector<unsigned long> GetElementNodes(unsigned long iElem) const;

  /*!
   * \brief Get the global node indices of the marker elements.
   * \param[in] iMarker - Marker index.
   * \return Element global node indices (nElem, nVertex).
   */
  vector<vector<unsigned long>> GetMarkerElementNodes(unsigned short iMarker) const;

  /*!
   * \brief Get the global node indices of a marker element.
   * \param[in] iMarker - Marker index.
   * \param[in] iElem - Marker element index.
   * \return Element global node indices (nVertex).
   */
  vector<unsigned long> GetMarkerElementNodes(unsigned short iMarker, unsigned long iElem) const;

  /*!
   * \brief Get the number of nodes in the mesh.
   * \return Number of nodes.
   */
  unsigned long GetNumberNodes() const;

  /*!
   * \brief Get the number of nodes in the marker.
   * \param[in] iMarker - Marker index.
   * \return Number of nodes.
   */
  unsigned long GetNumberMarkerNodes(unsigned short iMarker) const;

  /*!
   * \brief Get the number of halo nodes in the mesh.
   * \return Number of halo nodes.
   */
  unsigned long GetNumberHaloNodes() const;

  /*!
   * \brief Get the number of halo nodes in the marker.
   * \param[in] iMarker - Marker index.
   * \return Number of halo nodes.
   */
  unsigned long GetNumberMarkerHaloNodes(unsigned short iMarker) const;

  /*!
   * \brief Get the vertices of the marker.
   * \param[in] iMarker - Marker index.
   * \return Marker vertices (nVertex).
   */
  vector<unsigned long> GetMarkerVertices(unsigned short iMarker) const;

  /*!
   * \brief Get the vertex of a marker.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Marker vertex.
   */
  unsigned long GetMarkerVertices(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the global node indices.
   * \return Global node indices (nNode).
   */
  vector<unsigned long> GetNodes() const;

  /*!
   * \brief Get the global node index.
   * \param[in] iPoint - Mesh node index.
   * \return Global node index.
   */
  unsigned long GetNodes(unsigned long iPoint) const;

  /*!
   * \brief Get the global node indices of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \return Global node indices (nVertex).
   */
  vector<unsigned long> GetMarkerNodes(unsigned short iMarker) const;

  /*!
   * \brief Get the global node index of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Global node index.
   */
  unsigned long GetMarkerNodes(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the halo flags of the mesh nodes.
   * \return Node domain flags (nNode).
   */
  vector<bool> GetDomain() const;

  /*!
   * \brief Get the halo flag of a mesh node.
   * \param[in] iPoint - Mesh node index.
   * \return Node domain flag.
   */
  bool GetDomain(unsigned long iPoint) const;

  /*!
   * \brief Get the halo flags of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \return Domain flags (nVertex).
   */
  vector<bool> GetMarkerDomain(unsigned short iMarker) const;

  /*!
   * \brief Get the halo flag of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Domain flag.
   */
  bool GetMarkerDomain(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the initial (un-deformed) coordinates of the mesh nodes.
   * \return Initial node coordinates (nNode, nDim).
   */
  vector<vector<passivedouble>> GetInitialCoordinates() const;

  /*!
   * \brief Get the initial (un-deformed) coordinates of a mesh node.
   * \param[in] iPoint - Mesh node index.
   * \return Initial node coordinates (nDim).
   */
  vector<passivedouble> GetInitialCoordinates(unsigned long iPoint) const;

  /*!
   * \brief Get the initial (un-deformed) coordinates of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \return Initial node coordinates (nVertex, nDim).
   */
  vector<vector<passivedouble>> GetMarkerInitialCoordinates(unsigned short iMarker) const;

  /*!
   * \brief Get the initial (un-deformed) coordinates of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Initial node coordinates (nDim).
   */
  vector<passivedouble> GetMarkerInitialCoordinates(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the coordinates of the mesh nodes.
   * \return Node coordinates (nNode, nDim).
   */
  vector<vector<passivedouble>> GetCoordinates() const;

  /*!
   * \brief Get the coordinates of a mesh node.
   * \param[in] iPoint - Mesh node index.
   * \return Node coordinates (nDim).
   */
  vector<passivedouble> GetCoordinates(unsigned long iPoint) const;

  /*!
   * \brief Get the coordinates of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \return Node coordinates (nVertex, nDim).
   */
  vector<vector<passivedouble>> GetMarkerCoordinates(unsigned short iMarker) const;

  /*!
   * \brief Get the coordinates of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Node coordinates (nDim).
   */
  vector<passivedouble> GetMarkerCoordinates(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Set the coordinates of the mesh nodes.
   * \param[in] values - Node coordinates (nNode, nDim).
   */
  void SetCoordinates(vector<vector<passivedouble>> values);

  /*!
   * \brief Set the coordinates of a mesh node.
   * \param[in] iPoint - Mesh node index.
   * \param[in] values - Node coordinates (nDim).
   */
  void SetCoordinates(unsigned long iPoint, vector<passivedouble> values);

  /*!
   * \brief Set the coordinates of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \param[in] values - Node coordinates (nVertex, nDim).
   */
  void SetMarkerCoordinates(unsigned short iMarker, vector<vector<passivedouble>> values);

  /*!
   * \brief Set the coordinates of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \param[in] values - Node coordinates (nDim).
   */
  void SetMarkerCoordinates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);

  /*!
   * \brief Get the displacements of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \return Node displacements (nVertex, nDim).
   */
  vector<vector<passivedouble>> GetMarkerDisplacements(unsigned short iMarker) const;

  /*!
   * \brief Get the displacements of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Node displacements (nDim).
   */
  vector<passivedouble> GetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Set the displacements of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \param[in] values - Node displacements (nVertex, nDim).
   */
  void SetMarkerDisplacements(unsigned short iMarker, vector<vector<passivedouble>> values);

  /*!
   * \brief Set the displacements of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \param[in] values - Node displacements (nDim).
   */
  void SetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);

  /*!
   * \brief Get the velocities of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \return Node velocities (nVertex, nDim).
   */
  vector<vector<passivedouble>> GetMarkerVelocities(unsigned short iMarker) const;

  /*!
   * \brief Get the velocities of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Node velocities (nDim).
   */
  vector<passivedouble> GetMarkerVelocities(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Set the velocities of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \param[in] values - Node velocities (nVertex, nDim).
   */
  void SetMarkerVelocities(unsigned short iMarker, vector<vector<passivedouble>> values);

  /*!
   * \brief Set the velocities of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \param[in] values - Node velocities (nDim).
   */
  void SetMarkerVelocities(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);

  /*!
   * \brief Get the normal vectors of the marker vertices.
   * \param[in] iMarker - Marker index.
   * \param[in] normalize - If true, the unit (i.e. normalized) normal vector is returned.
   * \return Node normal vectors (nVertex, nDim).
   */
  vector<vector<passivedouble>> GetMarkerVertexNormals(unsigned short iMarker, bool normalize = false) const;

  /*!
   * \brief Get the normal vector of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \param[in] normalize - If true, the unit (i.e. normalized) normal vector is returned.
   * \return Node normal vector (nDim).
   */
  vector<passivedouble> GetMarkerVertexNormals(unsigned short iMarker, unsigned long iVertex,
                                               bool normalize = false) const;

  /*!
   * \brief Communicate the boundary mesh displacements.
   */
  void CommunicateMeshDisplacements(void);

 protected:
  /*!
   * \brief Initialize containers.
   */
  void SetContainers_Null();

  /*!
   * \brief Read in the config and mesh files.
   * \param[in] config - Definition of the particular problem.
   * \param[in] driver_config - Definition of the driver configuration.
   */
  void Input_Preprocessing(CConfig**& config, CConfig*& driver_config);

  /*!
   * \brief Construction of the edge-based data structure and the multi-grid structure.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] dummy - Definition of the dummy driver.
   */
  void Geometrical_Preprocessing(CConfig* config, CGeometry**& geometry, bool dummy);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void Solver_Preprocessing(CConfig* config, CGeometry** geometry, CSolver***& solver);

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   */
  void Numerics_Preprocessing(CConfig* config, CGeometry** geometry, CSolver*** solver, CNumerics****& numerics) const;

  /*!
   * \brief Preprocess the output container.
   * \param[in] config - Definition of the particular problem.
   * \param[in] driver_config - Definition of the driver configuration.
   * \param[in] output_container - Container vector with all the outputs.
   * \param[in] driver_output - Definition of the driver output.
   */
  void Output_Preprocessing(CConfig** config, CConfig* driver_config, COutput**& output_container,
                            COutput*& driver_output);
};
