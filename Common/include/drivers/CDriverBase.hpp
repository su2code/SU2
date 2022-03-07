/*!
 * \file CDriverBase.hpp
 * \brief Base class template for all drivers.
 * \author H. Patel, A. Gastaldi
 * \version 7.3.0 "Blackbird"
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

#include "../CConfig.hpp"
#include "../../../SU2_CFD/include/numerics/CNumerics.hpp"
#include "../../../SU2_CFD/include/solvers/CSolver.hpp"
#include "../../../SU2_CFD/include/output/COutput.hpp"

class CDriverBase {
    
protected:
    
    int rank,                                       /*!< \brief MPI Rank. */
    size;                                           /*!< \brief MPI Size. */
    char* config_file_name;                         /*!< \brief Configuration file name of the problem.*/
    
    su2double StartTime,                            /*!< \brief Start point of the timer for performance benchmarking.*/
    StopTime,                                       /*!< \brief Stop point of the timer for performance benchmarking.*/
    UsedTimePreproc,                                /*!< \brief Elapsed time between Start and Stop point of the timer for tracking preprocessing phase.*/
    UsedTimeCompute,                                /*!< \brief Elapsed time between Start and Stop point of the timer for tracking compute phase.*/
    UsedTime;                                       /*!< \brief Elapsed time between Start and Stop point of the timer.*/
    
    unsigned long TimeIter;
    
    unsigned short iMesh,                           /*!< \brief Iterator on mesh levels.*/
    iZone,                                          /*!< \brief Iterator on zones.*/
    nZone,                                          /*!< \brief Total number of zones in the problem. */
    nDim,                                           /*!< \brief Number of dimensions.*/
    iInst,                                          /*!< \brief Iterator on instance levels.*/
    *nInst,                                         /*!< \brief Total number of instances in the problem (per zone). */
    **interface_types;                              /*!< \brief Type of coupling between the distinct (physical) zones.*/
    
    CConfig **config_container;                     /*!< \brief Definition of the particular problem. */
    CConfig *driver_config;                         /*!< \brief Definition of the driver configuration. */
    COutput **output_container;                     /*!< \brief Pointer to the COutput class. */
    COutput *driver_output;                         /*!< \brief Definition of the driver output. */
    CGeometry ****geometry_container;               /*!< \brief Geometrical definition of the problem. */
    CSolver *****solver_container;                  /*!< \brief Container vector with all the solutions. */
    CNumerics ******numerics_container;             /*!< \brief Description of the numerical method (the way in which the equations are solved). */
    CSurfaceMovement **surface_movement;            /*!< \brief Surface movement classes of the problem. */
    CVolumetricMovement ***grid_movement;           /*!< \brief Volume grid movement classes of the problem. */
    CFreeFormDefBox*** FFDBox;                      /*!< \brief FFD FFDBoxes of the problem. */
    
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
    virtual void Preprocessing() {};
    
    /*!
     * \brief A virtual member.
     */
    virtual void Run() {};
    
    /*!
     * \brief A virtual member.
     */
    virtual void Update() {};
    
    /*!
     * \brief A virtual member.
     */
    virtual void Update_Legacy() {};
    
    /*!
     * \brief A virtual member.
     */
    virtual void Output() {};
    
    /*!
     * \brief A virtual member.
     */
    virtual void Postprocessing() {};
    /*!
     * \brief Get the number of markers in the mesh.
     * \return Number of markers.
     */
    unsigned short GetNumberMarkers() const;
    
    /*!
     * \brief Get all the boundary markers tags with their associated indices.
     * \return List of boundary markers tags with their indices.
     */
    map<string, int> GetMarkerIndices() const;
    
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
     * \brief Get the global IDs of the mesh elements.
     * \return Global element IDs.
     */
    vector<unsigned long> GetElementIDs() const;
    
    /*!
     * \brief Get the global ID of a mesh element.
     * \param[in] iElem - Mesh element index.
     * \return Global element ID.
     */
    unsigned long GetElementIDs(unsigned long iElem) const;
    
    /*!
     * \brief Get the global IDs of the marker elements.
     * \param[in] iMarker - Marker index.
     * \return Global element IDs.
     */
    vector<unsigned long> GetMarkerElementIDs(unsigned short iMarker) const;
    
    /*!
     * \brief Get the global IDs of a marker element.
     * \param[in] iMarker - Marker index.
     * \param[in] iBound - Marker element index.
     * \return Global element ID.
     */
    unsigned long GetMarkerElementIDs(unsigned short iMarker, unsigned long iBound) const;
    
    /*!
     * \brief Get the MPI colors for mesh elements.
     * \return Element colors.
     */
    vector<unsigned long> GetElementColors() const;
    
    /*!
     * \brief Get the MPI color for a mesh element.
     * \param[in] iElem - Mesh element index.
     * \return Element color.
     */
    unsigned long GetElementColors(unsigned long iElem) const;
    
    /*!
     * \brief Get the MPI colors for marker elements.
     * \param[in] iMarker - Marker index.
     * \return Element colors.
     */
    vector<unsigned long> GetMarkerElementColors(unsigned short iMarker) const;
     
    /*!
     * \brief Get the MPI color for a marker element.
     * \param[in] iMarker - Marker index.
     * \param[in] iBound - Marker element index.
     * \return Element color.
     */
    unsigned long GetMarkerElementColors(unsigned short iMarker, unsigned long iBound) const;

    /*!
     * \brief Get the table of vertex IDs belonging to the mesh elements.
     * \return Element connectivities (nElem, nNode)
     */
    vector<vector<unsigned long>> GetElementConnectivities() const;
    
    /*!
     * \brief Get the row of vertex IDs belonging to a mesh element.
     * \param[in] iElem - Mesh element index.
     * \return Element connectivity (nNode)
     */
    vector<unsigned long> GetElementConnectivities(unsigned long iElem) const;
    
    /*!
     * \brief Get the table of vertex IDs belonging to the marker elements.
     * \param[in] iMarker - Marker index.
     * \return Element connectivities (nBound, nNode).
     */
    vector<vector<unsigned long>> GetMarkerElementConnectivities(unsigned short iMarker) const;
    
    /*!
     * \brief Get the row of vertex IDs belonging to a marker element.
     * \param[in] iMarker - Marker index.
     * \param[in] iBound - Marker element index.
     * \return Element connectivity (nNode).
     */
    vector<unsigned long> GetMarkerElementConnectivities(unsigned short iMarker, unsigned long iBound) const;
    
    /*!
     * \brief Get the number of vertices in the mesh.
     * \return Number of vertices.
     */
    unsigned long GetNumberVertices() const;
    
    /*!
     * \brief Get the number of vertices in the marker.
     * \param[in] iMarker - Marker index.
     * \return Number of vertices.
     */
    unsigned long GetNumberMarkerVertices(unsigned short iMarker) const;
    
    /*!
     * \brief Get the number of halo vertices in the mesh.
     * \return Number of vertices.
     */
    unsigned long GetNumberHaloVertices() const;
    
    /*!
     * \brief Get the number of halo vertices in the marker.
     * \param[in] iMarker - Marker index.
     * \return Number of vertices.
     */
    unsigned long GetNumberMarkerHaloVertices(unsigned short iMarker) const;
    
    /*!
     * \brief Get the mesh vertex indices of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Mesh vertex indices.
     */
    vector<unsigned long> GetMarkerVertexIndex(unsigned short iMarker) const;
    
    /*!
     * \brief Get the mesh vertex index of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Mesh vertex index.
     */
    unsigned long GetMarkerVertexIndex(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the global IDs of the mesh vertices.
     * \return Global vertex IDs.
     */
    vector<unsigned long> GetVertexIDs() const;

    /*!
     * \brief Get the global ID of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Global vertex ID.
     */
    unsigned long GetVertexIDs(unsigned long iPoint) const;
    
    /*!
     * \brief Get the global IDs of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Global vertex IDs.
     */
    vector<unsigned long> GetMarkerVertexIDs(unsigned short iMarker) const;
    
    /*!
     * \brief Get the global ID of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Global vertex ID.
     */
    unsigned long GetMarkerVertexIDs(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the MPI colors of the mesh vertices.
     * \return Vertex colors.
     */
    vector<unsigned long> GetVertexColors() const;
    
    /*!
     * \brief Get the MPI color of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex color.
     */
    unsigned long GetVertexColors(unsigned long iPoint) const;
    
    /*!
     * \brief Get the MPI colors of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Vertex colors.
     */
    vector<unsigned long> GetMarkerVertexColors(unsigned short iMarker) const;
    
    /*!
     * \brief Get the MPI color of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex color.
     */
    unsigned long GetMarkerVertexColors(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the halo flags of the mesh vertices.
     * \return Vertex domain flags.
     */
    vector<bool> GetDomain() const;
    
    /*!
     * \brief Get the halo flag of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex domain flag.
     */
    bool GetDomain(unsigned long iPoint) const;
    
    /*!
     * \brief Get the halo flags of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Vertex domain flags.
     */
    vector<bool> GetMarkerDomain(unsigned short iMarker) const;
    
    /*!
     * \brief Get the halo flag of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex domain flag.
     */
    bool GetMarkerDomain(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the initial (un-deformed) coordinates of the mesh vertices.
     * \return Initial vertex coordinates (nPoint*nDim).
     */
    vector<passivedouble> GetInitialCoordinates() const;

    /*!
     * \brief Get the initial (un-deformed) coordinates of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Initial vertex coordinates (nDim).
     */
    vector<passivedouble> GetInitialCoordinates(unsigned long iPoint) const;
    
    /*!
     * \brief Get the initial (un-deformed) coordinates of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Initial vertex coordinates (nVertex*nDim).
     */
    vector<passivedouble> GetMarkerInitialCoordinates(unsigned short iMarker) const;
    
    /*!
     * \brief Get the initial (un-deformed) coordinates of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Initial vertex coordinates (nDim).
     */
    vector<passivedouble> GetMarkerInitialCoordinates(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the coordinates of the mesh vertices.
     * \return Vertex coordinates (nPoint*nDim).
     */
    vector<passivedouble> GetCoordinates() const;

    /*!
     * \brief Get the coordinates of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex coordinates (nDim).
     */
    vector<passivedouble> GetCoordinates(unsigned long iPoint) const;
    
    /*!
     * \brief Get the coordinates of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Vertex coordinates (nVertex*nDim).
     */
    vector<passivedouble> GetMarkerCoordinates(unsigned short iMarker) const;
    
    /*!
     * \brief Get the coordinates of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex coordinates (nDim).
     */
    vector<passivedouble> GetMarkerCoordinates(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the coordinates of the mesh vertices.
     * \param[in] values - Vertex coordinates (nPoint*nDim).
     */
    void SetCoordinates(vector<passivedouble> values);
    
    /*!
     * \brief Set the coordinates of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \param[in] values - Vertex coordinates (nDim).
     */
    void SetCoordinates(unsigned long iPoint, vector<passivedouble> values);
    
    /*!
     * \brief Set the coordinates of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \param[in] values - Vertex coordinates (nVertex*nDim).
     */
    void SetMarkerCoordinates(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Set the coordinates of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Vertex coordinates (nDim).
     */
    void SetMarkerCoordinates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);
    
    /*!
     * \brief Get the displacements of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Vertex displacements (nVertex*nDim).
     */
    vector<passivedouble> GetMarkerDisplacements(unsigned short iMarker) const;
    
    /*!
     * \brief Get the displacements of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex displacements (nDim).
     */
    vector<passivedouble> GetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the displacements of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \param[in] values - Vertex displacements (nVertex*nDim).
     */
    void SetMarkerDisplacements(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Set the displacements of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Vertex displacements (nDim).
     */
    void SetMarkerDisplacements(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);
    
    /*!
     * \brief Get the velocities of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \return Vertex velocities (nVertex*nDim).
     */
    vector<passivedouble> GetMarkerVelocities(unsigned short iMarker) const;

    /*!
     * \brief Get the velocities of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex velocities (nVertex*nDim).
     */
    vector<passivedouble> GetMarkerVelocities(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the velocities of the marker vertices.
     * \param[in] iMarker - Marker index.
     * \param[in] values - Vertex velocities (nVertex*nDim).
     */
    void SetMarkerVelocities(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Set the velocities of a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Vertex velocities (nDim).
     */
    void SetMarkerVelocities(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);
    
    /*!
     * \brief Get the normal vectors at the marker vertices.
     * \param[in] iMarker - Marker index.
     * \param[in] normalize - If true, the unit (i.e. normalized) normal vector is returned.
     * \return Normal vector at the vertex (nVertex*nDim).
     */
    vector<passivedouble> GetMarkerVertexNormals(unsigned short iMarker, bool normalize = false) const;
    
    /*!
     * \brief Get the normal vectors at a marker vertex.
     * \param[in] iMarker - Marker index.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] normalize - If true, the unit (i.e. normalized) normal vector is returned.
     * \return Normal vector at the vertex (nDim).
     */
    vector<passivedouble> GetMarkerVertexNormals(unsigned short iMarker, unsigned long iVertex, bool normalize = false) const;
    
    /*!
     * \brief Communicate the boundary mesh displacements in a python call
     */
    void CommunicateMeshDisplacements(void);
    
protected:
    
    /*!
     * \brief Initialize Containers
     */
    void SetContainers_Null();
    
    /*!
     * \brief Read in the config and mesh files.
     */
    void Input_Preprocessing(CConfig **&config, CConfig *&driver_config);
    
    /*!
     * \brief Construction of the edge-based data structure and the multigrid structure.
     */
    void Geometrical_Preprocessing(CConfig *config, CGeometry **&geometry, bool dummy);
    
    /*!
     * \brief Definition and allocation of all solution classes.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void Solver_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***&solver);
    
    /*!
     * \brief Definition and allocation of all solver classes.
     * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
     */
    void Numerics_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CNumerics ****&numerics) const;
    
    /*!
     * \brief Preprocess the output container.
     */
    void Output_Preprocessing(CConfig **config, CConfig *driver_config, COutput **&output_container, COutput *&driver_output);
    
};
