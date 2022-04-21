/*!
 * \file CDiscAdjDeformationDriver.cpp
 * \brief Headers of the main subroutines for driving the projection of sensitivities.
 * \author T. Economon, H. Kline, R. Sanchez
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

#define ENABLE_MAPS
#include "../../../Common/include/CConfig.hpp"
#undef ENABLE_MAPS

#include "../../../Common/include/parallelization/mpi_structure.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/fem/fem_geometry_structure.hpp"
#include "../../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../../Common/include/grid_movement/CVolumetricMovement.hpp"
#include "../../../SU2_CFD/include/output/COutput.hpp"
#include "../../../SU2_CFD/include/output/CBaselineOutput.hpp"
#include "../../../SU2_CFD/include/solvers/CBaselineSolver.hpp"

#include "../../../Common/include/drivers/CDriverBase.hpp"
#include "../../../SU2_CFD/include/solvers/CGradientSmoothingSolver.hpp"
#include "../../../SU2_CFD/include/numerics/CGradSmoothing.hpp"

/*!
 * \class CDiscAdjDeformationDriver
 * \brief Class for driving sensitivity projections.
 * \author A. Gastaldi, H. Patel
 * \version 7.3.0 "Blackbird"
 */
class CDiscAdjDeformationDriver : public CDriverBase {
    
protected:
    su2double** Gradient;
    ofstream Gradient_file; 
    
public:
    /*!
     * \brief Constructor of the class.
     * \param[in] confFile - Configuration file name.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CDiscAdjDeformationDriver(char* confFile, SU2_Comm MPICommunicator);
    
    /*!
     * \brief Destructor of the class.
     */
    ~CDiscAdjDeformationDriver(void);
    
    /*!
     * \brief Pre-process the driver data (includes solution allocation and initialization).
     */
    void Preprocess();
    
    /*!
     * \brief Launch the driver computation.
     */
    void Run();
    
    /*!
     * \brief Deallocation routine
     */
    void Postprocessing();
    
    /*!
     * \brief Get the (total) sensitivity of the objective function w.r.t. the surface coordinates or displacements at the mesh vertices.
     * \return Total sensitivity of the objective function w.r.t. the un-deformed coordinates.
     */
    vector<vector<passivedouble>> GetObjectiveCoordinatesTotalSensitivities() const;
    
    /*!
     * \brief Get the (total) sensitivity of the objective function w.r.t. the surface coordinates or displacements at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Total sensitivity of the objective function w.r.t. the un-deformed coordinates.
     */
    vector<passivedouble> GetObjectiveCoordinatesTotalSensitivities(unsigned long iPoint) const;
    
    /*!
     * \brief Get the (total) sensitivity of the objective function w.r.t. the surface coordinates or displacements at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Total sensitivity of the objective function w.r.t. the surface coordinates.
     */
    vector<vector<passivedouble>> GetMarkerObjectiveCoordinatesTotalSensitivities(unsigned short iMarker) const;
    
    /*!
     * \brief Get the (total) sensitivity of the objective function w.r.t. the surface coordinates or displacements at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Total sensitivity of the objective function w.r.t. the surface coordinates.
     */
    vector<passivedouble> GetMarkerObjectiveCoordinatesTotalSensitivities(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the (total) sensitivity of the objective function w.r.t. a design variable.
     * \return Total sensitivity of the objective function w.r.t. a design variable.
     */
    vector<vector<passivedouble>> GetObjectiveDVsTotalSensitivities() const;
    
    /*!
     * \brief Get the (total) sensitivity of the objective function w.r.t. the design variables.
     * \param[in] iDV - Design Variable index.
     * \return Total sensitivity of the objective function w.r.t. the design variables.
     */
    vector<passivedouble> GetObjectiveDVsTotalSensitivities(unsigned short iDV) const;
    
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
     * \brief Projection of the surface sensitivity using finite differences (FD).
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] surface_movement - Surface movement class of the problem.
     * \param[in] Gradient_file - Output file to store the gradient data.
     */
    
    void SetProjection_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double **Gradient);
    
    /*!
     * \brief Projection of the surface sensitivity using algorithmic differentiation (AD).
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] surface_movement - Surface movement class of the problem.
     * \param[in] Gradient_file - Output file to store the gradient data.
     */
    
    void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double **Gradient);
    
    /*!
     * \brief Prints the gradient information to a file.
     * \param[in] Gradient - The gradient data.
     * \param[in] config - Definition of the particular problem.
     * \param[in] Gradient_file - Output file to store the gradient data.
     */
    
    void OutputGradient(su2double** Gradient, CConfig* config, ofstream& Gradient_file);
    
    /*!
     * \brief Write the sensitivity (including mesh sensitivity) computed with the discrete adjoint method
     *  on the surface and in the volume to a file.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] val_nZone - Number of Zones.
     */
    
    void SetSensitivity_Files(CGeometry ****geometry, CConfig **config, unsigned short val_nZone);
    
    /*!
     * \brief Treatment of derivatives with the Sobolev smoothing solver.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] grid_movement - Volumetric movement class of the problem.
     */
    
    void DerivativeTreatment_MeshSensitivity(CGeometry *geometry, CConfig *config, CVolumetricMovement *grid_movement);
    
    /*!
     * \brief Treatment of derivatives with the Sobolev smoothing solver.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] grid_movement - Volumetric movement class of the problem.
     * \param[in] surface_movement - Surface movement class of the problem.
     * \param[in] Gradient - Output array to store the gradient data.
     */
    
    void DerivativeTreatment_Gradient(CGeometry *geometry, CConfig *config, CVolumetricMovement *grid_movement, CSurfaceMovement *surface_movement, su2double **Gradient);

};
