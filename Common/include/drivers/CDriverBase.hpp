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
     * \brief Get all the boundary markers tags with their associated indices.
     * \return List of boundary markers tags with their indices.
     */
    map<string, int> GetBoundaryMarkerIndices() const;
    
    /*!
     * \brief Get all the boundary markers tags with their associated types.
     * \return List of boundary markers tags with their types.
     */
    map<string, string> GetBoundaryMarkerTypes() const;
    
    /*!
     * \brief Get all the deformable boundary marker tags.
     * \return List of deformable boundary markers tags.
     */
    vector<string> GetDeformableMarkerTags() const;
    
    /*!
     * \brief Get the number of mesh dimensions.
     * \return Number of dimensions.
     */
    unsigned long GetNumberDimensions() const;
    
    /*!
     * \brief Get the number of mesh elements.
     * \return Number of elements.
     */
    unsigned long GetNumberElements() const;
    
    /*!
     * \brief Get the number of mesh elements from a specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Number of elements.
     */
    unsigned long GetNumberElementsMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get the number of mesh vertices.
     * \return Number of vertices.
     */
    unsigned long GetNumberVertices() const;
    
    /*!
     * \brief Get the number of mesh vertices from a specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Number of vertices.
     */
    unsigned long GetNumberVerticesMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get the number of halo mesh vertices.
     * \return Number of vertices.
     */
    unsigned long GetNumberHaloVertices() const;
    
    /*!
     * \brief Get the number of halo mesh vertices from a specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Number of vertices.
     */
    unsigned long GetNumberHaloVerticesMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get global IDs of mesh vertices.
     * \return Global vertex IDs.
     */
    vector<unsigned long> GetVertexIDs() const;
    
    /*!
     * \brief Get global IDs of mesh vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Global vertex IDs.
     */
    vector<unsigned long> GetVertexIDsMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get global IDs of mesh elements.
     * \return Global element IDs.
     */
    vector<unsigned long> GetElementIDs() const;
    
    /*!
     * \brief Get global IDs of mesh elements.
     * \param[in] iMarker - Marker identifier.
     * \return Global element IDs.
     */
    vector<unsigned long> GetElementIDsMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get the connected point IDs of mesh elements.
     * \return Element connectivities (nElem, nNode)
     */
    vector<vector<unsigned long>> GetConnectivity() const;
    
    /*!
     * \brief Get the connected point IDs of mesh elements on a specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Element connectivities (nBound, nNode).
     */
    vector<vector<unsigned long>> GetConnectivityMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get halo node stauts of mesh vertices.
     * \return Point domain status.
     */
    vector<bool> GetDomain() const;
    
    /*!
     * \brief Get halo node stauts of mesh marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Point domain status.
     */
    vector<bool> GetDomainMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get the coordinates of the mesh points.
     * \return Point coordinates (nPoint*nDim).
     */
    vector<passivedouble> GetCoordinates() const;
    
    /*!
     * \brief Get the coordinates of the mesh points on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Point coordinates (nVertex*nDim).
     */
    vector<passivedouble> GetCoordinatesMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Set the coordinates of the mesh points.
     * \param[in] values - Point coordinates (nPoint*nDim).
     */
    void SetCoordinates(vector<passivedouble> values);
    
    /*!
     * \brief Set the coordinates of the mesh points on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Point coordinates (nVertex*nDim).
     */
    void SetCoordinatesMarker(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Get the vertex displacements on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex displacements (nVertex*nDim).
     */
    vector<passivedouble> GetDisplacementsMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Set the vertex displacements on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Vertex displacements (nVertex*nDim).
     */
    void SetDisplacementsMarker(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Get the vertex velocities on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex velocities (nVertex*nDim).
     */
    vector<passivedouble> GetVelocitiesMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Set the vertex velocities on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Vertex velocities (nVertex*nDim).
     */
    void SetVelocitiesMarker(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Get undeformed coordinates from mesh solver on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Initial point coordinates (nVertex*nDim).
     */
    vector<passivedouble> GetInitialCoordinatesMarker(unsigned short iMarker) const;
    
    /*!
     * \brief Get the vertex normal vectors on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \param[in] UnitNormal - Boolean to indicate if unit normal vector should be returned.
     * \return Normal vector at the vertex (nVertex*nDim).
     */
    vector<passivedouble> GetVertexNormalsMarker(unsigned short iMarker, bool UnitNormal = false) const;
    
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
