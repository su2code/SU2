/*!
 * \file CDriverBase.hpp
 * \brief Base class for all drivers.
 * \author H. Patel, A. Gastaldi
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

#include <limits>
#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/containers/CPyWrapperMatrixView.hpp"
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
  unsigned short selected_zone = ZONE_0; /*!< \brief Selected zone for the driver. Defaults to ZONE_0 */
  unsigned short iMesh,  /*!< \brief Iterator on mesh levels. */
      iZone,             /*!< \brief Iterator on zones. */
      nZone,             /*!< \brief Total number of zones in the problem. */
      nDim,              /*!< \brief Number of dimensions. */
      iInst,             /*!< \brief Iterator on instance levels. */
      *nInst,            /*!< \brief Total number of instances in the problem (per zone). */
      **interface_types; /*!< \brief Type of coupling between the distinct (physical) zones. */

  CConfig* driver_config = nullptr; /*!< \brief Definition of the driver configuration. */
  COutput* driver_output = nullptr; /*!< \brief Definition of the driver output. */

  CConfig** config_container;           /*!< \brief Definition of the particular problem. */
  COutput** output_container;           /*!< \brief Pointer to the COutput class. */
  CGeometry**** geometry_container;     /*!< \brief Geometrical definition of the problem. */
  CSolver***** solver_container;        /*!< \brief Container vector with all the solutions. */
  CNumerics****** numerics_container;   /*!< \brief Description of the numerical method (the way in which the equations
                                           are solved). */
  CSurfaceMovement** surface_movement;  /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement*** grid_movement; /*!< \brief Volume grid movement classes of the problem. */

  CConfig* main_config = nullptr;     /*!< \brief Reference to base (i.e. ZONE 0) configuration (used in driver API). */
  CGeometry* main_geometry = nullptr; /*!< \brief Reference to base (i.e. ZONE, INST, MESH 0) geometry (used in driver API). */

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
  virtual void Run(){}

  /*!
   * \brief A virtual member.
   */
  virtual void Finalize(){}

/// \addtogroup PySU2
/// @{

  /*!
   * \brief Get the list of available outputs.
   * \return List of output names.
   */
  inline vector<string> GetOutputNames() const { return output_container[MESH_0]->GetHistoryOutputList(); }

  /*!
   * \brief Get the value of one of the available history outputs.
   * \return Value of the output.
   */
  inline passivedouble GetOutputValue(const std::string& output_name) const {
    return SU2_TYPE::GetValue(output_container[MESH_0]->GetHistoryFieldValue(output_name));
  }

  /*!
   * \brief Get the list of available surface outputs on **both** MARKER_MONITORING and MARKER_ANALYZE.
   * \return List of surface output names.
   */
  inline vector<string> GetMarkerOutputNames() const {
    return output_container[MESH_0]->GetHistoryOutputPerSurfaceList();
  }

  /*!
   * \brief Get the value of one of the available surface outputs at a given MARKER_MONITORING.
   * \return Value of the output.
   */
  inline passivedouble GetMarkerMonitoringOutputValue(const std::string& output_name,
                                                      const std::string& marker_monitoring) const {
    for (auto iMarker = 0u; iMarker < main_config->GetnMarker_Monitoring(); ++iMarker) {
      if (marker_monitoring == main_config->GetMarker_Monitoring_TagBound(iMarker))
        return SU2_TYPE::GetValue(output_container[MESH_0]->GetHistoryFieldValuePerSurface(output_name, iMarker));
    }
    SU2_MPI::Error(marker_monitoring + " is not in MARKER_MONITORING.", CURRENT_FUNCTION);
    return 0;
  }

  /*!
   * \brief Get the value of one of the available surface outputs at a given MARKER_ANALYZE.
   * \return Value of the output.
   */
  inline passivedouble GetMarkerAnalyzeOutputValue(const std::string& output_name,
                                                   const std::string& marker_analyze) const {
    for (auto iMarker = 0u; iMarker < main_config->GetnMarker_Analyze(); ++iMarker) {
      if (marker_analyze == main_config->GetMarker_Analyze_TagBound(iMarker))
        return SU2_TYPE::GetValue(output_container[MESH_0]->GetHistoryFieldValuePerSurface(output_name, iMarker));
    }
    SU2_MPI::Error(marker_analyze + " is not in MARKER_ANALYZE.", CURRENT_FUNCTION);
    return 0;
  }

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
   * \brief Get the global index of a mesh element.
   * \param[in] iElem - Mesh element index.
   * \return Global element index.
   */
  unsigned long GetElementGlobalIndex(unsigned long iElem) const;

  /*!
   * \brief Get the node indices of a mesh element.
   * \param[in] iElem - Mesh element index.
   * \return Element node indices (nNode).
   */
  vector<unsigned long> GetElementNodes(unsigned long iElem) const;

  /*!
   * \brief Get the number of nodes in the mesh (including halos).
   * \return Number of nodes.
   */
  unsigned long GetNumberNodes() const;

  /*!
   * \brief Get the number of halo nodes in the mesh.
   * \return Number of halo nodes.
   */
  unsigned long GetNumberHaloNodes() const;

  /*!
   * \brief Get the global node index.
   * \param[in] iPoint - Mesh node index.
   * \return Global node index.
   */
  unsigned long GetNodeGlobalIndex(unsigned long iPoint) const;

  /*!
   * \brief Get the halo flag of a mesh node.
   * \param[in] iPoint - Mesh node index.
   * \return Node domain flag.
   */
  bool GetNodeDomain(unsigned long iPoint) const;

  /*!
   * \brief Get a read-only view of the initial (undeformed) coordinates of all mesh nodes.
   */
  inline CPyWrapperMatrixView InitialCoordinates() const {
    if (!main_config->GetDeform_Mesh()) {
      SU2_MPI::Error("Initial coordinates are only available with DEFORM_MESH= YES", CURRENT_FUNCTION);
    }
    auto* coords =
        const_cast<su2activematrix*>(solver_container[selected_zone][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord());
    return CPyWrapperMatrixView(*coords, "InitialCoordinates", true);
  }

  /*!
   * \brief Get a read-only view of the initial (undeformed) coordinates of the mesh nodes of a marker.
   */
  inline CPyWrapperMarkerMatrixView MarkerInitialCoordinates(unsigned short iMarker) const {
    if (!main_config->GetDeform_Mesh()) {
      SU2_MPI::Error("Initial coordinates are only available with DEFORM_MESH= YES", CURRENT_FUNCTION);
    }
    if (iMarker >= GetNumberMarkers()) SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);

    auto* coords =
        const_cast<su2activematrix*>(solver_container[selected_zone][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord());
    return CPyWrapperMarkerMatrixView(*coords, main_geometry->vertex[iMarker], main_geometry->GetnVertex(iMarker),
                                      "MarkerInitialCoordinates", true);
  }

  /*!
   * \brief Get a read/write view of the current coordinates of all mesh nodes.
   */
  inline CPyWrapperMatrixView Coordinates() {
    auto& coords = const_cast<su2activematrix&>(main_geometry->nodes->GetCoord());
    return CPyWrapperMatrixView(coords, "Coordinates", false);
  }

  /*!
   * \brief Get a read/write view of the current coordinates of the mesh nodes of a marker.
   */
  inline CPyWrapperMarkerMatrixView MarkerCoordinates(unsigned short iMarker) {
    if (iMarker >= GetNumberMarkers()) SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
    auto& coords = const_cast<su2activematrix&>(main_geometry->nodes->GetCoord());
    return CPyWrapperMarkerMatrixView(coords, main_geometry->vertex[iMarker], main_geometry->GetnVertex(iMarker),
                                      "MarkerCoordinates", false);
  }

  /*!
   * \brief Get the number of markers in the mesh.
   * \return Number of markers.
   */
  unsigned short GetNumberMarkers() const;

  /*!
   * \brief Get all the boundary markers tags with their associated indices.
   * \return Map of boundary markers tags to their indices.
   */
  map<string, unsigned short> GetMarkerIndices() const;

  /*!
   * \brief Get all the boundary markers tags with their associated types.
   * \return Map of boundary markers tags to their types.
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
   * \brief Get all the CHT boundary marker tags.
   * \return List of CHT boundary markers tags.
   */
  vector<string> GetCHTMarkerTags() const;

  /*!
   * \brief Get all the inlet boundary marker tags.
   * \return List of inlet boundary markers tags.
   */
  vector<string> GetInletMarkerTags() const;

  /*!
   * \brief Get the number of elements in the marker.
   * \param[in] iMarker - Marker index.
   * \return Number of elements.
   */
  unsigned long GetNumberMarkerElements(unsigned short iMarker) const;

  /*!
   * \brief Get the global index of a marker element.
   * \param[in] iMarker - Marker index.
   * \param[in] iElem - Marker element index.
   * \return Global element index.
   */
  unsigned long GetMarkerElementGlobalIndex(unsigned short iMarker, unsigned long iElem) const;

  /*!
   * \brief Get the node indices of a marker element.
   * \param[in] iMarker - Marker index.
   * \param[in] iElem - Marker element index.
   * \return Element node indices.
   */
  vector<unsigned long> GetMarkerElementNodes(unsigned short iMarker, unsigned long iElem) const;

  /*!
   * \brief Get the number of nodes in the marker.
   * \param[in] iMarker - Marker index.
   * \return Number of nodes.
   */
  unsigned long GetNumberMarkerNodes(unsigned short iMarker) const;

  /*!
   * \brief Get the node index of a marker.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Marker vertex.
   */
  unsigned long GetMarkerNode(unsigned short iMarker, unsigned long iVertex) const;

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
   * \brief Get the displacements currently imposed of a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Node displacements (nDim).
   */
  inline vector<passivedouble> GetMarkerDisplacement(unsigned short iMarker, unsigned long iVertex) const {
    vector<passivedouble> disp(GetNumberDimensions(), 0.0);
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(MESH_SOL)->GetNodes();

    for (auto iDim = 0u; iDim < GetNumberDimensions(); ++iDim) {
      disp[iDim] = SU2_TYPE::GetValue(nodes->GetBound_Disp(iPoint, iDim));
    }
    return disp;
  }

  /*!
   * \brief Set the mesh displacements of a marker vertex.
   * \note This can be the input of the flow solver in an FSI setting.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \param[in] values - Node displacements (nDim).
   */
  inline void SetMarkerCustomDisplacement(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(MESH_SOL)->GetNodes();

    for (auto iDim = 0u; iDim < GetNumberDimensions(); iDim++) {
      nodes->SetBound_Disp(iPoint, iDim, values[iDim]);
    }
  }

  /*!
   * \brief Get the mesh velocities currently imposed on a marker vertex.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \return Node velocities (nDim).
   */
  inline vector<passivedouble> GetMarkerMeshVelocity(unsigned short iMarker, unsigned long iVertex) const {
    vector<passivedouble> vel(GetNumberDimensions(), 0.0);
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(MESH_SOL)->GetNodes();

    for (auto iDim = 0u; iDim < GetNumberDimensions(); ++iDim) {
      vel[iDim] = SU2_TYPE::GetValue(nodes->GetBound_Vel(iPoint, iDim));
    }
    return vel;
  }

  /*!
   * \brief Set the velocities of a marker vertex.
   * \note This can be the input of the flow solver in an unsteady FSI setting.
   * \param[in] iMarker - Marker index.
   * \param[in] iVertex - Marker vertex index.
   * \param[in] values - Node velocities (nDim).
   */
  inline void SetMarkerCustomMeshVelocity(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values) {
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(MESH_SOL)->GetNodes();

    for (auto iDim = 0u; iDim < GetNumberDimensions(); iDim++) {
      nodes->SetBound_Vel(iPoint, iDim, values[iDim]);
    }
  }

  /*!
   * \brief Communicate the boundary mesh displacements.
   */
  void CommunicateMeshDisplacements(void);

  /*!
   * \brief Get all the active solver names with their associated indices (which can be used to access their data).
   */
  map<string, unsigned short> GetSolverIndices() const;

  /*!
   * \brief Get the structural solver solution variable names with their associated indices.
   * These correspond to the column indices in the matrix returned by e.g. Solution().
   */
  map<string, unsigned short> GetFEASolutionIndices() const;

  /*!
   * \brief Get a read/write view of the current solution on all mesh nodes of a solver.
   */
  inline CPyWrapperMatrixView Solution(unsigned short iSolver) {
    auto* solver = GetSolverAndCheckMarker(iSolver);
    return CPyWrapperMatrixView(solver->GetNodes()->GetSolution(), "Solution of " + solver->GetSolverName(), false);
  }

  /*!
   * \brief Get a read/write view of the current solution on the mesh nodes of a marker.
   */
  inline CPyWrapperMarkerMatrixView MarkerSolution(unsigned short iSolver, unsigned short iMarker) {
    auto* solver = GetSolverAndCheckMarker(iSolver, iMarker);
    return CPyWrapperMarkerMatrixView(
        solver->GetNodes()->GetSolution(), main_geometry->vertex[iMarker], main_geometry->GetnVertex(iMarker),
        "MarkerSolution of " + solver->GetSolverName(), false);
  }

  /*!
   * \brief Get a read/write view of the solution at time N on all mesh nodes of a solver.
   */
  inline CPyWrapperMatrixView SolutionTimeN(unsigned short iSolver) {
    auto* solver = GetSolverAndCheckMarker(iSolver);
    return CPyWrapperMatrixView(
        solver->GetNodes()->GetSolution_time_n(), "SolutionTimeN of " + solver->GetSolverName(), false);
  }

  /*!
   * \brief Get a read/write view of the solution at time N on the mesh nodes of a marker.
   */
  inline CPyWrapperMarkerMatrixView MarkerSolutionTimeN(unsigned short iSolver, unsigned short iMarker) {
    auto* solver = GetSolverAndCheckMarker(iSolver, iMarker);
    return CPyWrapperMarkerMatrixView(
        solver->GetNodes()->GetSolution_time_n(), main_geometry->vertex[iMarker], main_geometry->GetnVertex(iMarker),
        "MarkerSolutionTimeN of " + solver->GetSolverName(), false);
  }

  /*!
   * \brief Get a read/write view of the solution at time N-1 on all mesh nodes of a solver.
   */
  inline CPyWrapperMatrixView SolutionTimeN1(unsigned short iSolver) {
    auto* solver = GetSolverAndCheckMarker(iSolver);
    return CPyWrapperMatrixView(
        solver->GetNodes()->GetSolution_time_n1(), "SolutionTimeN1 of " + solver->GetSolverName(), false);
  }

  /*!
   * \brief Get a read/write view of the solution at time N-1 on the mesh nodes of a marker.
   */
  inline CPyWrapperMarkerMatrixView MarkerSolutionTimeN1(unsigned short iSolver, unsigned short iMarker) {
    auto* solver = GetSolverAndCheckMarker(iSolver, iMarker);
    return CPyWrapperMarkerMatrixView(
        solver->GetNodes()->GetSolution_time_n1(), main_geometry->vertex[iMarker], main_geometry->GetnVertex(iMarker),
        "MarkerSolutionTimeN1 of " + solver->GetSolverName(), false);
  }

  /*!
   * \brief Get the flow solver primitive variable names with their associated indices.
   * These correspond to the column indices in the matrix returned by Primitives.
   */
  map<string, unsigned short> GetPrimitiveIndices() const;

  /*!
   * \brief Get a read/write view of the current primitive variables on all mesh nodes of the flow solver.
   * \warning Primitive variables are only available for flow solvers.
   */
  inline CPyWrapperMatrixView Primitives() {
    auto* solver = GetSolverAndCheckMarker(FLOW_SOL);
    return CPyWrapperMatrixView(const_cast<su2activematrix&>(solver->GetNodes()->GetPrimitive()), "Primitives", false);
  }

  /*!
   * \brief Get a read/write view of the current primitive variables on the mesh nodes of a marker.
   * \warning Primitive variables are only available for flow solvers.
   */
  inline CPyWrapperMarkerMatrixView MarkerPrimitives(unsigned short iMarker) {
    auto* solver = GetSolverAndCheckMarker(FLOW_SOL, iMarker);
    return CPyWrapperMarkerMatrixView(
        const_cast<su2activematrix&>(solver->GetNodes()->GetPrimitive()), main_geometry->vertex[iMarker],
        main_geometry->GetnVertex(iMarker), "MarkerPrimitives", false);
  }

  /*!
   * \brief Get a read-only view of the geometry sensitivity of a discrete adjoint solver.
   */
  inline CPyWrapperMatrixView Sensitivity(unsigned short iSolver) {
    auto* solver = GetSolverAndCheckMarker(iSolver);
    auto& sensitivity = const_cast<su2activematrix&>(solver->GetNodes()->GetSensitivity());
    return CPyWrapperMatrixView(sensitivity, "Sensitivity", true);
  }

  /*!
   * \brief Set the temperature of a vertex on a specified marker (MARKER_PYTHON_CUSTOM).
   * \note This can be the input of a heat or flow solver in a CHT setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] WallTemp - Value of the temperature.
   */
  inline void SetMarkerCustomTemperature(unsigned short iMarker, unsigned long iVertex, passivedouble WallTemp)  {
    main_geometry->SetCustomBoundaryTemperature(iMarker, iVertex, WallTemp);
  }

  /*!
   * \brief Set the wall normal heat flux at a vertex on a specified marker (MARKER_PYTHON_CUSTOM).
   * \note This can be the input of a heat or flow solver in a CHT setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] WallHeatFlux - Value of the normal heat flux.
   */
  inline void SetMarkerCustomNormalHeatFlux(unsigned short iMarker, unsigned long iVertex, passivedouble WallHeatFlux) {
    main_geometry->SetCustomBoundaryHeatFlux(iMarker, iVertex, WallHeatFlux);
  }

  /*!
   * \brief Selects zone to be used for python driver operations.
   * \param[in] iZone - Zone identifier.
   */
  inline void SelectZone(unsigned short iZone) {
    if (iZone >= nZone) SU2_MPI::Error("Zone index out of range", CURRENT_FUNCTION);
    selected_zone = iZone;
    main_geometry = geometry_container[selected_zone][INST_0][MESH_0];
    main_config = config_container[selected_zone];
  }

  /*!
   * \brief Returns the index of the zone selected for python driver operations.
   */
  inline unsigned short SelectedZone() const { return selected_zone; }

  /*!
   * \brief Get the wall normal heat flux at a vertex on a specified marker of the flow or heat solver.
   * \note This can be the output of a heat or flow solver in a CHT setting.
   * \param[in] iSolver - Solver identifier, should be either a flow solver or the heat solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Wall normal component of the heat flux at the vertex.
   */
  inline passivedouble GetMarkerNormalHeatFlux(unsigned short iSolver, unsigned short iMarker, unsigned long iVertex) const {
    if (iSolver != HEAT_SOL && iSolver != FLOW_SOL) {
      SU2_MPI::Error("Normal heat flux is only available for flow or heat solvers.", CURRENT_FUNCTION);
    }
    return SU2_TYPE::GetValue(GetSolverAndCheckMarker(iSolver, iMarker)->GetHeatFlux(iMarker, iVertex));
  }

  /*!
   * \brief Sets the nodal force for the structural solver at a vertex of a marker.
   * \note This can be the input of the FEA solver in an FSI setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] force - Force vector.
   */
  inline void SetMarkerCustomFEALoad(unsigned short iMarker, unsigned long iVertex, std::vector<passivedouble> force) {
    auto* solver = GetSolverAndCheckMarker(FEA_SOL, iMarker);
    std::array<su2double, 3> load{};
    for (auto iDim = 0u; iDim < GetNumberDimensions(); ++iDim) load[iDim] = force[iDim];
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    solver->GetNodes()->Set_FlowTraction(iPoint, load.data());
  }

  /*!
   * \brief Get the fluid force at a vertex of a solid wall marker of the flow solver.
   * \note This can be the output of the flow solver in an FSI setting to then apply it to a structural solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vector of loads.
   */
  inline vector<passivedouble> GetMarkerFlowLoad(unsigned short iMarker, unsigned long iVertex) const {
    vector<passivedouble> FlowLoad(GetNumberDimensions(), 0.0);
    const auto* solver = GetSolverAndCheckMarker(FLOW_SOL, iMarker);

    if (main_config->GetSolid_Wall(iMarker)) {
      for (auto iDim = 0u; iDim < GetNumberDimensions(); ++iDim) {
        FlowLoad[iDim] = SU2_TYPE::GetValue(solver->GetVertexTractions(iMarker, iVertex, iDim));
      }
    }
    return FlowLoad;
  }

  /*!
   * \brief Set the adjoint of the flow tractions of the flow solver.
   * \note This can be the input of the flow solver in an adjoint FSI setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] adjointLoad - Vector of adjoint loads.
   */
  inline void SetMarkerCustomFlowLoadAdjoint(unsigned short iMarker, unsigned long iVertex,
                                             vector<passivedouble> adjointLoad) {
    auto* solver = GetSolverAndCheckMarker(FLOW_SOL, iMarker);
    for (auto iDim = 0u; iDim < GetNumberDimensions(); ++iDim) {
      solver->StoreVertexTractionsAdjoint(iMarker, iVertex, iDim, adjointLoad[iDim]);
    }
  }

  /*!
   * \brief Get the sensitivities of the displacements of the mesh boundary vertices.
   * \note This can be the output of the flow solver in an adjoint FSI setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vector of sensitivities.
   */
  inline vector<passivedouble> GetMarkerDisplacementSensitivity(unsigned short iMarker, unsigned long iVertex) const {
    const auto nDim = GetNumberDimensions();
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(ADJMESH_SOL)->GetNodes();

    vector<passivedouble> sens(nDim, 0.0);
    for (auto iDim = 0u; iDim < nDim; ++iDim) {
      sens[iDim] = SU2_TYPE::GetValue(nodes->GetBoundDisp_Sens(iPoint, iDim));
    }
    return sens;
  }

  /*!
   * \brief Get the sensitivity of the FEA loads of the structural solver (via the adjoint structural solver).
   * \note This can be the output of the FEA solver in an adjoint FSI setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \returns Vector of sensitivities.
   */
  inline vector<passivedouble> GetMarkerFEALoadSensitivity(unsigned short iMarker, unsigned long iVertex) const {
    const auto nDim = GetNumberDimensions();
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(ADJFEA_SOL)->GetNodes();

    vector<passivedouble> sens(nDim, 0.0);
    for (auto iDim = 0u; iDim < nDim; ++iDim) {
      sens[iDim] = SU2_TYPE::GetValue(nodes->GetFlowTractionSensitivity(iPoint, iDim));
    }
    return sens;
  }

  /*!
   * \brief Set the adjoint of the structural displacements.
   * \note This can be the input of the FEA solver in an adjoint FSI setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] adjointDisplacement - Vector of adjoint displacements.
   */
  inline void SetMarkerCustomFEADisplacementAdjoint(unsigned short iMarker, unsigned long iVertex,
                                                    vector<passivedouble> adjointDisplacement) {
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(ADJFEA_SOL)->GetNodes();

    for (auto iDim = 0u; iDim < GetNumberDimensions(); ++iDim) {
      nodes->SetSourceTerm_DispAdjoint(iPoint, iDim, adjointDisplacement[iDim]);
    }
  }

  /*!
   * \brief Set the adjoint of the structural velocities.
   * \note This can be the input of the FEA solver in an unsteady adjoint FSI setting.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] adjointVelocity - Vector of adjoint velocities.
   */
  inline void SetMarkerCustomFEAVelocityAdjoint(unsigned short iMarker, unsigned long iVertex,
                                                vector<passivedouble> adjointVelocity) {
    const auto iPoint = GetMarkerNode(iMarker, iVertex);
    auto* nodes = GetSolverAndCheckMarker(ADJFEA_SOL)->GetNodes();

    for (auto iDim = 0u; iDim < GetNumberDimensions(); ++iDim) {
      nodes->SetSourceTerm_VelAdjoint(iPoint, iDim, adjointVelocity[iDim]);
    }
  }

/// \}

 protected:
  /*!
   * \brief Automates some boilerplate of accessing solution fields for the python wrapper.
   */
  inline CSolver* GetSolverAndCheckMarker(unsigned short iSolver,
                                          unsigned short iMarker = std::numeric_limits<unsigned short>::max()) const {
    if (iMarker < std::numeric_limits<unsigned short>::max() && iMarker > GetNumberMarkers()) {
      SU2_MPI::Error("Marker index exceeds size.", CURRENT_FUNCTION);
    }
    auto* solver = solver_container[selected_zone][INST_0][MESH_0][iSolver];
    if (solver == nullptr) SU2_MPI::Error("The selected solver does not exist.", CURRENT_FUNCTION);
    return solver;
  }

  /*!
   * \brief Initialize containers.
   */
  void InitializeContainers();

  /*!
   * \brief Delete containers.
   */
  void CommonFinalize();

  /*!
   * \brief Read in the config and mesh files.
   * \param[in] config - Definition of the particular problem.
   * \param[in] driver_config - Definition of the driver configuration.
   */
  void InputPreprocessing(CConfig**& config, CConfig*& driver_config);

  /*!
   * \brief Construction of the edge-based data structure and the multi-grid structure.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] dummy - Definition of the dummy driver.
   */
  void InitializeGeometry(CConfig* config, CGeometry**& geometry, bool dummy);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void InitializeSolver(CConfig* config, CGeometry** geometry, CSolver***& solver);

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   */
  void InitializeNumerics(CConfig* config, CGeometry** geometry, CSolver*** solver, CNumerics****& numerics) const;

  /*!
   * \brief Preprocess the output container.
   * \param[in] config - Definition of the particular problem.
   * \param[in] driver_config - Definition of the driver configuration.
   * \param[in] output_container - Container vector with all the outputs.
   * \param[in] driver_output - Definition of the driver output.
   */
  void OutputPreprocessing(CConfig** config, CConfig* driver_config, COutput**& output_container,
                            COutput*& driver_output);
};
