/*!
 * \file driver_structure.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 7.3.0 "Blackbird"
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

#include "../../../Common/include/parallelization/mpi_structure.hpp"

#include "../integration/CIntegration.hpp"
#include "../solvers/CSolver.hpp"
#include "../interfaces/CInterface.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/drivers/CDriverBase.hpp"

using namespace std;

class COutputLegacy;
class CInterpolator;
class CIteration;
class COutput;

/*!
 * \class CDriver
 * \brief Parent class for driving an iteration of a single or multi-zone problem.
 * \author T. Economon
 */

class CDriver : public CDriverBase {
    
protected:
    
    char runtime_file_name[MAX_STRING_SIZE];
    su2double UsedTimeOutput;                     /*!< \brief Elapsed time between Start and Stop point of the timer for tracking output phase.*/
    
    su2double BandwidthSum = 0.0;                 /*!< \brief Aggregate value of the bandwidth for writing restarts (to be average later).*/
    unsigned long IterCount,                       /*!< \brief Iteration count stored for performance benchmarking.*/
    OutputCount;                                  /*!< \brief Output count stored for performance benchmarking.*/
    unsigned long DOFsPerPoint;                    /*!< \brief Number of unknowns at each vertex, i.e., number of equations solved. */
    su2double Mpoints;                            /*!< \brief Total number of grid points in millions in the calculation (including ghost points).*/
    su2double MpointsDomain;                      /*!< \brief Total number of grid points in millions in the calculation (excluding ghost points).*/
    su2double MDOFs;                              /*!< \brief Total number of DOFs in millions in the calculation (including ghost points).*/
    su2double MDOFsDomain;                        /*!< \brief Total number of DOFs in millions in the calculation (excluding ghost points).*/
    
    ofstream **ConvHist_file;                     /*!< \brief Convergence history file.*/
    ofstream FSIHist_file;                        /*!< \brief FSI convergence history file.*/
    
    bool StopCalc,                                /*!< \brief Stop computation flag.*/
    mixingplane,                                  /*!< \brief mixing-plane simulation flag.*/
    fsi,                                          /*!< \brief FSI simulation flag.*/
    fem_solver;                                   /*!< \brief FEM fluid solver simulation flag. */
    
    CIteration ***iteration_container;            /*!< \brief Container vector with all the iteration methods. */
    CIntegration ****integration_container;       /*!< \brief Container vector with all the integration methods. */
    vector<vector<unique_ptr<CInterpolator> > >
    interpolator_container;                       /*!< \brief Definition of the interpolation method between non-matching discretizations of the interface. */
    CInterface ***interface_container;            /*!< \brief Definition of the interface of information and physics. */
    bool dry_run;                                 /*!< \brief Flag if SU2_CFD was started as dry-run via "SU2_CFD -d <config>.cfg" */
    
public:
    
    /*!
     * \brief Constructor of the class.
     * \param[in] confFile - Configuration file name.
     * \param[in] val_nZone - Total number of zones.
     * \param[in] val_nDim - Number of dimensions.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CDriver(char* confFile,
            unsigned short val_nZone,
            SU2_Comm MPICommunicator, bool dummy_geo);
    
    /*!
     * \brief Destructor of the class.
     */
    virtual ~CDriver(void);
    
    /*!
     * \brief A virtual member.
     */
    virtual void Run() { };
    
protected:
    
    /*!
     * \brief Init_Containers
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
     * \brief Do the geometrical preprocessing for the DG FEM solver.
     */
    void Geometrical_Preprocessing_DGFEM(CConfig *config, CGeometry **&geometry);
    
    /*!
     * \brief Geometrical_Preprocessing_FVM
     */
    void Geometrical_Preprocessing_FVM(CConfig *config, CGeometry **&geometry);
    
    /*!
     * \brief Definition of the physics iteration class or within a single zone.
     * \param[in] iteration_container - Pointer to the iteration container to be instantiated.
     * \param[in] config - Definition of the particular problem.
     * \param[in] iZone - Index of the zone.
     */
    void Iteration_Preprocessing(CConfig *config, CIteration *&iteration) const;
    
    /*!
     * \brief Definition and allocation of all solution classes.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void Solver_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***&solver);
    
    /*!
     * \brief Restart of the solvers from the restart files.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void Solver_Restart(CSolver ***solver, CGeometry **geometry, CConfig *config, bool update_geo);
    
    /*!
     * \brief Definition and allocation of all solution classes.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void Solver_Postprocessing(CSolver ****solver, CGeometry **geometry, CConfig *config, unsigned short val_iInst);
    
    /*!
     * \brief Definition and allocation of all integration classes.
     * \param[in] config - Definition of the particular problem.
     * \param[in] solver - Container vector with all the solutions.
     * \param[out] integration - Container vector with all the integration methods.
     */
    void Integration_Preprocessing(CConfig *config, CSolver **solver, CIntegration **&integration) const;
    
    /*!
     * \brief Definition and allocation of all integration classes.
     * \param[in] integration_container - Container vector with all the integration methods.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void Integration_Postprocessing(CIntegration ***integration, CGeometry **geometry, CConfig *config, unsigned short val_iInst);
    
    /*!
     * \brief Definition and allocation of all interface classes.
     */
    void Interface_Preprocessing(CConfig **config, CSolver *****solver, CGeometry ****geometry,
                                 unsigned short **interface_types, CInterface ***interface,
                                 vector<vector<unique_ptr<CInterpolator> > > &interpolation);
    
    /*!
     * \brief Definition and allocation of all solver classes.
     * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
     */
    void Numerics_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CNumerics ****&numerics) const;
    
    /*!
     * \brief Helper to instantiate turbulence numerics specialized for different flow solvers.
     */
    template <class FlowIndices>
    void InstantiateTurbulentNumerics(unsigned short nVar_Turb, int offset, const CConfig *config,
                                      const CSolver* turb_solver, CNumerics ****&numerics) const;
    
    /*!
     * \brief Helper to instantiate species transport numerics specialized for different flow solvers.
     */
    template <class FlowIndices>
    void InstantiateSpeciesNumerics(unsigned short nVar_Species, int offset, const CConfig *config,
                                    const CSolver* species_solver, CNumerics ****&numerics) const;
    
    /*!
     * \brief Definition and allocation of all solver classes.
     * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void Numerics_Postprocessing(CNumerics *****numerics, CSolver ***solver, CGeometry **geometry, CConfig *config, unsigned short val_iInst);
    
    /*!
     * \brief GridMovement_Preprocessing
     * \param config
     * \param geometry
     * \param solver
     * \param iteration
     * \param grid_movement
     * \param surface_movement
     */
    void DynamicMesh_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CIteration *iteration, CVolumetricMovement *&grid_movement, CSurfaceMovement *&surface_movement) const;
    
    /*!
     * \brief Initialize Python interface functionalities
     */
    void PythonInterface_Preprocessing(CConfig** config, CGeometry**** geometry, CSolver***** solver);
    
    /*!
     * \brief Preprocess the output container.
     */
    void Output_Preprocessing(CConfig **config, CConfig *driver_config, COutput **&output_container, COutput *&driver_output);
    
    /*!
     * \brief Initiate value for static mesh movement such as the gridVel for the ROTATING frame.
     */
    void StaticMesh_Preprocessing(const CConfig *config, CGeometry **geometry);
    
    /*!
     * \brief Initiate value for static mesh movement such as the gridVel for the ROTATING frame.
     */
    void Turbomachinery_Preprocessing(CConfig** config, CGeometry**** geometry, CSolver***** solver,
                                      CInterface*** interface);
    
    /*!
     * \brief A virtual member.
     * \param[in] donorZone - zone in which the displacements will be predicted.
     * \param[in] targetZone - zone which receives the predicted displacements.
     */
    virtual void Predict_Displacements(unsigned short donorZone, unsigned short targetZone) {}
    
    /*!
     * \brief A virtual member.
     * \param[in] donorZone - zone in which the tractions will be predicted.
     * \param[in] targetZone - zone which receives the predicted traction.
     */
    virtual void Predict_Tractions(unsigned short donorZone, unsigned short targetZone) {}
    
    /*!
     * \brief A virtual member.
     * \param[in] donorZone - zone in which the displacements will be transferred.
     * \param[in] targetZone - zone which receives the tractions transferred.
     */
    virtual void Transfer_Displacements(unsigned short donorZone, unsigned short targetZone) {}
    
    /*!
     * \brief A virtual member.
     * \param[in] donorZone - zone from which the tractions will be transferred.
     * \param[in] targetZone - zone which receives the tractions transferred.
     */
    virtual void Transfer_Tractions(unsigned short donorZone, unsigned short targetZone) {}
    
    /*!
     * \brief A virtual member.
     * \param[in] donorZone - origin of the information.
     * \param[in] targetZone - destination of the information.
     * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
     */
    virtual void Relaxation_Displacements(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter) {}
    
    /*!
     * \brief A virtual member.
     * \param[in] donorZone - origin of the information.
     * \param[in] targetZone - destination of the information.
     * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
     */
    virtual void Relaxation_Tractions(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter) {}
    
    /*!
     * \brief A virtual member to run a Block Gauss-Seidel iteration in multizone problems.
     */
    virtual void Run_GaussSeidel(){}
    
    /*!
     * \brief A virtual member to run a Block-Jacobi iteration in multizone problems.
     */
    virtual void Run_Jacobi(){}
    
    /*!
     * \brief A virtual member.
     */
    virtual void Update() {}
    
    /*!
     * \brief Print out the direct residuals.
     * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
     */
    void Print_DirectResidual(RECORDING kind_recording);
    
public:
    
    /*!
     * \brief Launch the computation for all zones and all physics.
     */
    virtual void StartSolver() {}
    
    /*!
     * \brief Deallocation routine
     */
    void Postprocessing();
    
    /*!
     * \brief A virtual member.
     */
    virtual void ResetConvergence();
    
    /*!
     * \brief Perform some pre-processing before an iteration of the physics.
     */
    virtual void Preprocess(unsigned long TimeIter){ }
    
    /*!
     * \brief Monitor the computation.
     */
    virtual bool Monitor(unsigned long TimeIter){ return false; }
    
    /*!
     * \brief Output the solution in solution file.
     */
    virtual void Output(unsigned long TimeIter){ }
    
    /*!
     * \brief Perform a dynamic mesh deformation, including grid velocity computation and update of the multigrid structure.
     */
    virtual void DynamicMeshUpdate(unsigned long TimeIter) { }
    
    /*!
     * \brief Perform a dynamic mesh deformation, including grid velocity computation and update of the multigrid structure.
     */
    virtual void DynamicMeshUpdate(unsigned short val_iZone, unsigned long TimeIter) { }
    
    /*!
     * \brief Perform a mesh deformation as initial condition.
     */
    virtual void SetInitialMesh() { }
    
    /*!
     * \brief Process the boundary conditions and update the multigrid structure.
     */
    void UpdateBoundaryConditions();
    
    /*!
     * \brief Get the total drag force or drag coefficient.
     * \param[in] coefficient - Boolean to indicate if normalized coefficient should be returned.
     * \return Total drag.
     */
    passivedouble GetDrag(bool coefficient = true) const;
    
    /*!
     * \brief Get the total lift force or lift coefficient.
     * \param[in] coefficient - Boolean to indicate if normalized coefficient should be returned.
     * \return Total lift.
     */
    passivedouble GetLift(bool coefficient = true) const;
    
    /*!
     * \brief Get the total roll moment (x-direction) or roll moment coefficient.
     * \param[in] coefficient - Boolean to indicate if normalized coefficient should be returned.
     * \return Total roll moment.
     */
    passivedouble GetRollMoment(bool coefficient = true) const;
    
    /*!
     * \brief Get the total pitch moment (y-direction) or pitch moment coefficient.
     * \param[in] coefficient - Boolean to indicate if normalized coefficient should be returned.
     * \return Total pitch moment.
     */
    passivedouble GetPitchMoment(bool coefficient = true) const;
    
    /*!
     * \brief Get the total yaw moment (z-direction) or yaw moment coefficient.
     * \param[in] coefficient - Boolean to indicate if normalized coefficient should be returned.
     * \return Total yaw moment.
     */
    passivedouble GetYawMoment(bool coefficient = true) const;
    
    /*!
     * \brief Get the objective function.
     * \return Objective function.
     */
    passivedouble GetObjective() const;
    
    /*!
     * \brief Get the number of design variables.
     * \return Number of design variables.
     */
    passivedouble GetNumberDesignVariables() const;
    
    /*!
     * \brief Get the number of FFD boxes.
     * \return Number of FFD boxes.
     */
    passivedouble GetNumberFFDBoxes() const;
    
    /*!
     * \brief Get the number of external iterations.
     * \return Number of external iterations.
     */
    unsigned long GetNumberTimeIterations() const;
    
    /*!
     * \brief Get the current external iteration.
     * \return Current external iteration.
     */
    unsigned long GetCurrentTimeIteration() const;
    
    /*!
     * \brief Get the unsteady time step.
     * \return Unsteady time step.
     */
    passivedouble GetUnsteadyTimeStep() const;
    
    /*!
     * \brief Get the name of the output file for the surface.
     * \return File name for the surface output.
     */
    string GetSurfaceFileName() const;
    
    /*!
     * \brief Get the temperature at the mesh vertices.
     * \return Vertex temperatures (nPoint).
     */
    vector<passivedouble> GetTemperatures() const;
    
    /*!
     * \brief Get the temperature at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex temperature.
     */
    passivedouble GetTemperatures(unsigned long iPoint) const;
    
    /*!
     * \brief Get the temperature at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex temperatures (nVertex).
     */
    vector<passivedouble> GetMarkerTemperatures(unsigned short iMarker) const;
    
    /*!
     * \brief Get the temperature at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex temperature.
     */
    passivedouble GetMarkerTemperatures(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the temperatures at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Vertex temperatures (nVertex).
     */
    void SetMarkerTemperatures(unsigned short iMarker, vector<passivedouble> values);

    /*!
     * \brief Set the temperatures at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Vertex temperature.
     */
    void SetMarkerTemperatures(unsigned short iMarker, unsigned long iVertex, passivedouble value);

    /*!
     * \brief Get the heat fluxes at the mesh vertices.
     * \return Vertex heat fluxes (nPoint).
     */
    vector<vector<passivedouble>> GetHeatFlux() const;
    
    /*!
     * \brief Get the heat fluxes at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex heat flux.
     */
    vector<passivedouble> GetHeatFlux(unsigned long iPoint) const;
    
    /*!
     * \brief Get the heat fluxes at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex heat fluxes (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerHeatFlux(unsigned short iMarker) const;
    
    /*!
     * \brief Get the heat fluxes at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex heat fluxes (nVertex, nDim).
     */
    vector<passivedouble> GetMarkerHeatFlux(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the normal heat fluxes at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex normal heat fluxes (nVertex, nDim).
     */
    vector<passivedouble> GetMarkerNormalHeatFlux(unsigned short iMarker) const;
    
    /*!
     * \brief Get the normal heat fluxes at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex normal heat fluxes (nVertex, nDim).
     */
    passivedouble GetMarkerNormalHeatFlux(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the normal heat fluxes at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Vertex normal heat fluxes (nVertex).
     */
    void SetMarkerNormalHeatFlux(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Set the normal heat flux at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Vertex normal heat flux.
     */
    void SetMarkerNormalHeatFlux(unsigned short iMarker, unsigned long iVertex, passivedouble values);
    
    /*!
     * \brief Get the thermal conductivity at the mesh vertices.
     * \return Vertex thermal conductivities (nPoint).
     */
    vector<passivedouble> GetThermalConductivity() const;

    /*!
     * \brief Get the thermal conductivity at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex thermal conductivity.
     */
    passivedouble GetThermalConductivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the thermal conductivity at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex thermal conductivities (nVertex).
     */
    vector<passivedouble> GetMarkerThermalConductivity(unsigned short iMarker) const;
    
    /*!
     * \brief Get the thermal conductivity at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex thermal conductivity.
     */
    passivedouble GetMarkerThermalConductivity(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the laminar viscosity at the mesh vertices.
     * \return Vertex laminar viscosities (nPoint).
     */
    vector<passivedouble> GetLaminarViscosity() const;
    
    /*!
     * \brief Get the laminar viscosity at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex laminar viscosity.
     */
    passivedouble GetLaminarViscosity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the laminar viscosity at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex laminar viscosities (nVertex).
     */
    vector<passivedouble> GetMarkerLaminarViscosity(unsigned short iMarker) const;
    
    /*!
     * \brief Get the laminar viscosity at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex laminar viscosity.
     */
    passivedouble GetMarkerLaminarViscosity(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the eddy viscosity at the mesh vertices.
     * \return Vertex eddy viscosities (nPoint).
     */
    vector<passivedouble> GetEddyViscosity() const;
    
    /*!
     * \brief Get the eddy viscosity at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Vertex eddy viscosity.
     */
    passivedouble GetEddyViscosity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the eddy viscosity at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Vertex eddy viscosity (nVertex).
     */
    vector<passivedouble> GetMarkerEddyViscosity(unsigned short iMarker) const;
    
    /*!
     * \brief Get eddy viscosity at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Vertex eddy viscosity.
     */
    passivedouble GetMarkerEddyViscosity(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Preprocess the inlets via file input for all solvers.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void Inlet_Preprocessing(CSolver ***solver, CGeometry **geometry, CConfig *config) const;
    
    /*!
     * \brief Get the free-stream Reynolds number.
     * \return Free-stream Reynolds number.
     */
    passivedouble GetReynoldsNumber() const;
    
    /*!
     * \brief Get the free-stream Mach number.
     * \return Free-stream Mach number.
     */
    passivedouble GetMachNumber() const;
    
    /*!
     * \brief Get the free-stream angle of attack (in degrees).
     * \return Free-stream angle of attack.
     */
    passivedouble GetAngleOfAttack() const;
    
    /*!
     * \brief Get the free-stream angle of side-slip (in degrees).
     * \return Free-stream angle of side-slip.
     */
    passivedouble GetAngleOfSideslip() const;
    
    /*!
     * \brief Set the free-stream Reynolds number.
     * \param[in] value - User-defined Reynolds number.
     */
    void SetReynoldsNumber(passivedouble value);
    
    /*!
     * \brief Set the free-stream Mach number.
     * \param[in] value - User-defined Mach number.
     */
    void SetMachNumber(passivedouble value);
    
    /*!
     * \brief Set the free-stream angle of attack (in degrees).
     * \param[in] value - Free-stream angle of attack.
     */
    void SetAngleOfAttack(passivedouble value);
    
    /*!
     * \brief Set the free-stream angle of sideslip (in degrees).
     * \param[in] value - Free-stream angle of sideslip.
     */
    void SetAngleOfSideslip(passivedouble value);
    
    /*!
     * \brief Get the number of conservative state variables.
     * \return Number of conservative state variables.
     */
    unsigned long GetNumberStateVariables() const;
    
    /*!
     * \brief Get the number of primitive state variables.
     * \return Number of primitive state variables.
     */
    unsigned long GetNumberPrimitiveVariables() const;
    
    /*!
     * \brief Get the residuals of the conservative flow variables at the mesh vertices.
     * \return Residuals of the conservative flow variables (nPoint, nVar).
     */
    vector<vector<passivedouble>> GetResiduals() const;
    
    /*!
     * \brief Get the residuals of the conservative flow variables at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Residuals of the conservative flow variables (nVar).
     */
    vector<passivedouble> GetResiduals(unsigned long iPoint) const;
    
    /*!
     * \brief Get the residuals of the conservative flow variables at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Residuals of the conservative flow variables (nVertex, nVar).
     */
    vector<vector<passivedouble>> GetMarkerResiduals(unsigned short iMarker) const;
    
    /*!
     * \brief Get the residuals of the conservative flow variables at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Residuals of the conservative flow variables (nVar).
     */
    vector<passivedouble> GetMarkerResiduals(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the conservative flow states at the mesh vertices.
     * \return Flow states (nPoint, nVar).
     */
    vector<vector<passivedouble>> GetStates() const;
    
    /*!
     * \brief Get the conservative flow states at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Flow states (nVar).
     */
    vector<passivedouble> GetStates(unsigned long iPoint) const;
    
    /*!
     * \brief Get the conservative flow states at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Flow states (nVertex, nVar).
     */
    vector<vector<passivedouble>> GetMarkerStates(unsigned short iMarker) const;
    
    /*!
     * \brief Get the conservative flow states at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Flow states (nVar).
     */
    vector<passivedouble> GetMarkerStates(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the conservative states at the mesh vertices.
     * \param[in] values - Flow states (nPoint, nVar).
     */
    void SetStates(vector<vector<passivedouble>> values);
    
    /*!
     * \brief Set the conservative states at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \param[in] values - Flow states (nVar).
     */
    void SetStates(unsigned long iPoint, vector<passivedouble> values);
    
    /*!
     * \brief Set the conservative states at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Flow states (nVertex, nVar).
     */
    void SetMarkerStates(unsigned short iMarker, vector<vector<passivedouble>> values);
    
    /*!
     * \brief Set the conservative states at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Flow states (nVar).
     */
    void SetMarkerStates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);
    
    /*!
     * \brief Get the primitive state variables (i.e. density, velocities, and pressure) at the mesh vertices.
     * \return Primitive states (nPoint, nPrim).
     */
    vector<vector<passivedouble>> GetPrimitiveStates() const;
    
    /*!
     * \brief Get the primitive state variables (i.e. density, velocities, and pressure) at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Primitive states (nPrim).
     */
    vector<passivedouble> GetPrimitiveStates(unsigned long iPoint) const;
    
    /*!
     * \brief Get the primitive state variables (i.e. density, velocities, and pressure) at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Primitive states (nVertex, nPrim).
     */
    vector<vector<passivedouble>> GetMarkerPrimitiveStates(unsigned short iMarker) const;
    
    /*!
     * \brief Get the primitive state variables (i.e. density, velocities, and pressure) at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Primitive states (nPrim)
     */
    vector<passivedouble> GetMarkerPrimitiveStates(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the local speed of sound at the mesh vertices.
     * \return Speed of sound (nPoint).
     */
    vector<passivedouble> GetSpeedOfSound() const;
    
    /*!
     * \brief Get the local speed of sound at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Speed of sound.
     */
    passivedouble GetSpeedOfSound(unsigned long iPoint) const;
    
    /*!
     * \brief Get the speed of sound on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Speed of sound (nVertex).
     */
    vector<passivedouble> GetMarkerSpeedOfSound(unsigned short iMarker) const;
    
    /*!
     * \brief Get the speed of sound on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Speed of sound.
     */
    passivedouble GetMarkerSpeedOfSound(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the surface forces at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Surface forces (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerForces(unsigned short iMarker) const;
    
    /*!
     * \brief Get the surface forces at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Surface forces (nDim).
     */
    vector<passivedouble> GetMarkerForces(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the adjoint conservative flow states at the mesh vertices.
     * \return Adjoint flow states (nPoint, nVar).
     */
    vector<vector<passivedouble>> GetAdjointStates() const;
    
    /*!
     * \brief Get the adjoint conservative flow states at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Adjoint flow states (nVar).
     */
    vector<passivedouble> GetAdjointStates(unsigned long iPoint) const;
    
    /*!
     * \brief Get the adjoint conservative flow states at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Adjoint flow states (nVertex, nVar).
     */
    vector<vector<passivedouble>> GetMarkerAdjointStates(unsigned short iMarker) const;
    
    /*!
     * \brief Get the adjoint conservative flow states at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Adjoint flow states (nVar).
     */
    vector<passivedouble> GetMarkerAdjointStates(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the adjoint conservative states at the mesh vertices.
     * \param[in] values - Adjoint flow states (nPoint, nVar).
     */
    void SetAdjointStates(vector<vector<passivedouble>> values);
    
    /*!
     * \brief Set the adjoint conservative states at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \param[in] values - Adjoint flow states (nVar).
     */
    void SetAdjointStates(unsigned long iPoint, vector<passivedouble> values);
    
    /*!
     * \brief Set the adjoint conservative states at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Adjoint flow states (nVertex, nVar).
     */
    void SetMarkerAdjointStates(unsigned short iMarker, vector<vector<passivedouble>> values);
    
    /*!
     * \brief Set the adjoint conservative states at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Adjoint flow states (nVar).
     */
    void SetMarkerAdjointStates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);
    
    /*!
     * \brief Get the adjoint flow forces at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Adjoint flow forces (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerAdjointForces(unsigned short iMarker) const;
    
    /*!
     * \brief Get the adjoint flow forces at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Adjoint flow forces (nDim).
     */
    vector<passivedouble> GetMarkerAdjointForces(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the adjoint flow forces at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Adjoint flow forces (nVertex, nDim).
     */
    void SetMarkerAdjointForces(unsigned short iMarker, vector<vector<passivedouble>> values);
    
    /*!
     * \brief Set the adjoint flow forces at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Adjoint flow forces (nDim).
     */
    void SetMarkerAdjointForces(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);
    
    /*!
     * \brief Get the adjoint coordinates of the mesh vertices.
     * \return Adjoint mesh coordinates (nPoint, nDim).
     */
    vector<vector<passivedouble>> GetAdjointCoordinates() const;
    
    /*!
     * \brief Get the adjoint coordinates of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Adjoint mesh coordinates (nDim).
     */
    vector<passivedouble> GetAdjointCoordinates(unsigned long iPoint) const;
    
    /*!
     * \brief Get the adjoint coordinates of the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Adjoint mesh coordinates (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerAdjointCoordinates(unsigned short iMarker) const;
    
    /*!
     * \brief Get the adjoint coordinates of a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Adjoint mesh coordinates (nDim).
     */
    vector<passivedouble> GetMarkerAdjointCoordinates(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Set the adjoint coordinates of the mesh vertices.
     * \param[in] values - Adjoint mesh coordinates (nPoint, nDim).
     */
    void SetAdjointCoordinates(vector<vector<passivedouble>> values);

    /*!
     * \brief Set the adjoint coordinates of a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \param[in] values - Adjoint mesh coordinates (nDim).
     */
    void SetAdjointCoordinates(unsigned long iPoint, vector<passivedouble> values);
    
    /*!
     * \brief Set the adjoint coordinates of the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Adjoint mesh coordinates (nVertex, nDim).
     */
    void SetMarkerAdjointCoordinates(unsigned short iMarker, vector<vector<passivedouble>> values);
    
    /*!
     * \brief Set the adjoint coordinates of a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \param[in] values - Adjoint mesh coordinates (nDim).
     */
    void SetMarkerAdjointCoordinates(unsigned short iMarker, unsigned long iVertex, vector<passivedouble> values);
    
    /*!
     * \brief Get the sensitivity of volume coordinates w.r.t. initial coordinates at the mesh vertices.
     * \return Sensitivity of the deformed mesh coordinates w.r.t. the initial coordinates (nPoint, nDim).
     */
    vector<vector<passivedouble>> GetCoordinatesCoordinatesSensitivity() const;
    
    /*!
     * \brief Get the sensitivity of volume coordinates w.r.t. initial coordinates at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Sensitivity of the deformed mesh coordinates w.r.t. the initial coordinates (nDim).
     */
    vector<passivedouble> GetCoordinatesCoordinatesSensitivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the sensitivity of volume coordinates w.r.t. boundary displacements at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Sensitivity of the deformed mesh coordinates w.r.t. the boundary displacements (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerCoordinatesDisplacementsSensitivity(unsigned short iMarker) const;
    
    /*!
     * \brief Get the sensitivity of volume coordinates w.r.t. boundary displacements at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Sensitivity of the deformed mesh coordinates w.r.t. the boundary displacements (nDim).
     */
    vector<passivedouble> GetMarkerCoordinatesDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the objective function w.r.t. the far-field flow variables (mach, alpha) at the mesh vertices.
     * \return Sensitivity of the objective function w.r.t. the far-field flow variables (nTrim).
     */
    vector<passivedouble> GetObjectiveFarfieldVariablesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow residuals w.r.t. the far-field flow variables (mach, alpha) at a mesh vertex.
     * \return Sensitivity of the flow residuals w.r.t. the far-field flow variables (nTrim).
     */
    vector<passivedouble> GetResidualsFarfieldVariablesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the objective function w.r.t. the conservative flow states at the mesh vertices.
     * \return Sensitivity of the objective function w.r.t. the flow states (nPoint, nVar).
     */
    vector<vector<passivedouble>> GetObjectiveStatesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the objective function w.r.t. the conservative flow states at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Sensitivity of the objective function w.r.t. the flow states (nVar).
     */
    vector<passivedouble> GetObjectiveStatesSensitivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow residuals w.r.t. the conservative flow states at the mesh vertices.
     * \return Sensitivity of the flow residuals w.r.t. the flow states (nPoint, nVar)
     */
    vector<vector<passivedouble>> GetResidualsStatesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow residuals w.r.t. the conservative flow states at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Sensitivity of the flow residuals w.r.t. the flow states (nVar)
     */
    vector<passivedouble> GetResidualsStatesSensitivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow forces w.r.t. the conservative flow states at the mesh vertices.
     * \return Sensitivity of the flow forces w.r.t. the flow states (nPoint, nVar).
     */
    vector<vector<passivedouble>> GetForcesStatesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow forces w.r.t. the conservative flow states at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Sensitivity of the flow forces w.r.t. the flow states (nVar).
     */
    vector<passivedouble> GetForcesStatesSensitivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the objective function w.r.t. the un-deformed coordinates at the mesh vertices.
     * \return Sensitivity of the objective function w.r.t. the flow states (nPoint, nDim).
     */
    vector<vector<passivedouble>> GetObjectiveCoordinatesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the objective function w.r.t. the un-deformed coordinates at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Sensitivity of the objective function w.r.t. the flow states (nDim).
     */
    vector<passivedouble> GetObjectiveCoordinatesSensitivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow residuals w.r.t. the un-deformed mesh coordinates at the mesh vertices.
     * \return Sensitivity of the flow residuals w.r.t. the flow states (nPoint, nDim)
     */
    vector<vector<passivedouble>> GetResidualsCoordinatesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow residuals w.r.t. the un-deformed mesh coordinates at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Sensitivity of the flow residuals w.r.t. the flow states (nDim)
     */
    vector<passivedouble> GetResidualsCoordinatesSensitivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow forces w.r.t. the un-deformed coordinates at the mesh vertices.
     * \return Sensitivity of the flow forces w.r.t. the flow states (nPoint, nDim).
     */
    vector<vector<passivedouble>> GetForcesCoordinatesSensitivity() const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow forces w.r.t. the un-deformed coordinates at a mesh vertex.
     * \param[in] iPoint - Mesh vertex index.
     * \return Sensitivity of the flow forces w.r.t. the flow states (nDim).
     */
    vector<passivedouble> GetForcesCoordinatesSensitivity(unsigned long iPoint) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the objective function w.r.t. the mesh displacements or coordinates at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Sensitivity of the objective function w.r.t. the mesh displacements or coordinates (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerObjectiveDisplacementsSensitivity(unsigned short iMarker) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the objective function w.r.t. the mesh displacements or coordinates at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Sensitivity of the objective function w.r.t. the mesh displacements or coordinates (nDim).
     */
    vector<passivedouble> GetMarkerObjectiveDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow residuals w.r.t. the mesh displacements or coordinates at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Sensitivity of the flow residuals w.r.t. the mesh displacements or coordinates (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerResidualsDisplacementsSensitivity(unsigned short iMarker) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow residuals w.r.t. the mesh displacements or coordinates at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Sensitivity of the flow residuals w.r.t. the mesh displacements or coordinates (nDim).
     */
    vector<passivedouble> GetMarkerResidualsDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow forces w.r.t. the mesh displacements or coordinates at the marker vertices.
     * \param[in] iMarker - Marker identifier.
     * \return Sensitivity of the flow forces w.r.t. the mesh displacements or coordinates (nVertex, nDim).
     */
    vector<vector<passivedouble>> GetMarkerForcesDisplacementsSensitivity(unsigned short iMarker) const;
    
    /*!
     * \brief Get the (partial) sensitivity of the flow forces w.r.t. the mesh displacements or coordinates at a marker vertex.
     * \param[in] iMarker - Marker identifier.
     * \param[in] iVertex - Marker vertex index.
     * \return Sensitivity of the flow forces w.r.t. the mesh displacements or coordinates (nDim).
     */
    vector<passivedouble> GetMarkerForcesDisplacementsSensitivity(unsigned short iMarker, unsigned long iVertex) const;
    
    /*!
     * \brief Get all the flow load boundary marker tags.
     * \return List of flow load boundary markers tags.
     */
    vector<string> GetFluidLoadMarkerTags() const;
    
    /*!
     * \brief Get all the heat transfer boundary markers tags.
     * \return List of heat transfer boundary markers tags.
     */
    vector<string> GetCHTMarkerTags() const;
    
    /*!
     * \brief Get all the inlet boundary marker tags.
     * \return List of inlet boundary markers tags.
     */
    vector<string> GetInletMarkerTags() const;
    
    /*!
     * \brief Get sensitivities of the mesh boundary displacements.
     * \param[in] iMarker - Marker identifier.
     * \return Sensitivity of mesh displacements (nVertex, nDim).
     */
    vector<passivedouble> GetMarkerDisplacementsSensitivity(unsigned short iMarker) const;
    
    /*!
     * \brief Get sensitivities of the flow forces for the structural solver.
     * \param[in] iMarker - Marker identifier.
     * \return Sensitivity of flow forces (nVertex, nDim).
     */
    vector<passivedouble> GetMarkerForcesSensitivity(unsigned short iMarker) const;
    
    /*!
     * \brief Set forces for the structural solver on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - FEA traction components (nVertex, nDim).
     */
    void SetMarkerFEAForces(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Get displacements from the structural solver on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Structural displacements (nVertex, nDim).
     */
    vector<passivedouble> GetMarkerFEADisplacements(unsigned short iMarker) const;
    
    /*!
     * \brief Get velocities from the structural solver on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Structural velocities (nVertex, nDim).
     */
    vector<passivedouble> GetMarkerFEAVelocity(unsigned short iMarker) const;
    
    /*!
     * \brief Get velocities at time n from the structural solver on the specified marker.
     * \param[in] iMarker - Marker identifier.
     * \return Structural velocities at time n (nVertex, nDim).
     */
    vector<passivedouble> GetMarkerCurrentFEAVelocity(unsigned short iMarker) const;
    
    /*!
     * \brief Set the adjoint of the structural displacements (from an outside source).
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Adjoint structural displacements (nVertex, nDim).
     */
    void SetMarkerAdjointDisplacementSourceTerm(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Set the adjoint of the structural velocities (from an outside source).
     * \param[in] iMarker - Marker identifier.
     * \param[in] values - Adjoint structural velocities (nVertex, nDim).
     */
    void SetMarkerAdjointVelocitySourceTerm(unsigned short iMarker, vector<passivedouble> values);
    
    /*!
     * \brief Set the position of the heat source.
     * \param[in] alpha - Angle of rotation respect to Z axis.
     * \param[in] pos_x - Position X.
     * \param[in] pos_y - Position Y.
     * \param[in] pos_z - Position Z.
     */
    void SetHeatSourcePosition(passivedouble alpha, passivedouble pos_x, passivedouble pos_y, passivedouble pos_z);
    
    /*!
     * \brief Set the direction of the inlet.
     * \param[in] iMarker - Marker index.
     * \param[in] alpha - Angle (Zpos).
     */
    void SetInletAngle(unsigned short iMarker, passivedouble alpha);
    
    /*!
     * \brief Sum the number of primal or adjoint variables for all solvers in a given zone.
     * \param[in] iZone - Index of the zone.
     * \param[in] adjoint - True to consider adjoint solvers instead of primal.
     * \return Total number of solution variables.
     */
    unsigned short GetTotalNumberOfVariables(unsigned short iZone, bool adjoint) const {
        unsigned short nVar = 0;
        for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {
            auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
            if (solver && (solver->GetAdjoint() == adjoint)) nVar += solver->GetnVar();
        }
        return nVar;
    }
    
    /*!
     * \brief Set the solution of all solvers (adjoint or primal) in a zone.
     * \param[in] iZone - Index of the zone.
     * \param[in] adjoint - True to consider adjoint solvers instead of primal.
     * \param[in] solution - Solution object with interface (iPoint,iVar).
     * \tparam Old - If true set "old solutions" instead.
     */
    template<class Container, bool Old = false>
    void SetAllSolutions(unsigned short iZone, bool adjoint, const Container& solution) {
        const auto nPoint = geometry_container[iZone][INST_0][MESH_0]->GetnPoint();
        for (auto iSol = 0u, offset = 0u; iSol < MAX_SOLS; ++iSol) {
            auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
            if (!(solver && (solver->GetAdjoint() == adjoint))) continue;
            for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
                for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar)
                    if (!Old) solver->GetNodes()->SetSolution(iPoint, iVar, solution(iPoint,offset+iVar));
                    else solver->GetNodes()->SetSolution_Old(iPoint, iVar, solution(iPoint,offset+iVar));
            offset += solver->GetnVar();
        }
    }
    
    /*!
     * \brief Set the "old solution" of all solvers (adjoint or primal) in a zone.
     */
    template<class Container>
    void SetAllSolutionsOld(unsigned short iZone, bool adjoint, const Container& solution) {
        SetAllSolutions<Container,true>(iZone, adjoint, solution);
    }
    
    /*!
     * \brief Get the solution of all solvers (adjoint or primal) in a zone.
     * \param[in] iZone - Index of the zone.
     * \param[in] adjoint - True to consider adjoint solvers instead of primal.
     * \param[out] solution - Solution object with interface (iPoint,iVar).
     */
    template<class Container>
    void GetAllSolutions(unsigned short iZone, bool adjoint, Container& solution) const {
        const auto nPoint = geometry_container[iZone][INST_0][MESH_0]->GetnPoint();
        for (auto iSol = 0u, offset = 0u; iSol < MAX_SOLS; ++iSol) {
            auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
            if (!(solver && (solver->GetAdjoint() == adjoint))) continue;
            const auto& sol = solver->GetNodes()->GetSolution();
            for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
                for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar)
                    solution(iPoint,offset+iVar) = SU2_TYPE::GetValue(sol(iPoint,iVar));
            offset += solver->GetnVar();
        }
    }
    
};

/*!
 * \class CFluidDriver
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon, G. Gori
 */
class CFluidDriver : public CDriver {
    
protected:
    unsigned long Max_Iter;
    
protected:
    
    /*!
     * \brief Constructor of the class.
     * \param[in] confFile - Configuration file name.
     * \param[in] val_nZone - Total number of zones.
     * \param[in] val_nDim - Number of dimensions.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CFluidDriver(char* confFile,
                 unsigned short val_nZone,
                 SU2_Comm MPICommunicator);
    
public:
    /*!
     * \brief Destructor of the class.
     */
    ~CFluidDriver(void) override;
    
    /*!
     * \brief Launch the computation for all zones and all physics.
     */
    void StartSolver() override;
    
    /*!
     * \brief Run a single iteration of the physics within multiple zones.
     */
    void Run() override;
    
    /*!
     * \brief Update the dual-time solution within multiple zones.
     */
    void Update() override;
    
    /*!
     * \brief Output the solution in solution file.
     */
    void Output(unsigned long InnerIter) override;
    
    /*!
     * \brief Monitor the computation.
     */
    bool Monitor(unsigned long ExtIter) override;
    
    /*!
     * \brief Perform some pre-processing before an iteration of the physics.
     */
    void Preprocess(unsigned long Iter) override;
    
    /*!
     * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure (multiple zone).
     */
    void DynamicMeshUpdate(unsigned long TimeIter) override;
    
    /*!
     * \brief Transfer data among different zones (multiple zone).
     */
    void Transfer_Data(unsigned short donorZone, unsigned short targetZone);
    
};


/*!
 * \class CTurbomachineryDriver
 * \brief Class for driving an iteration for turbomachinery flow analysis.
 * \author S. Vitale
 */
class CTurbomachineryDriver : public CFluidDriver {
private:
    COutputLegacy* output_legacy;
    
public:
    
    /*!
     * \brief Constructor of the class.
     * \param[in] confFile - Configuration file name.
     * \param[in] val_nZone - Total number of zones.
     * \param[in] val_nDim - Number of dimensions.
     * \param[in] val_periodic - Bool for periodic BCs.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CTurbomachineryDriver(char* confFile,
                          unsigned short val_nZone,
                          SU2_Comm MPICommunicator);
    
    /*!
     * \brief Destructor of the class.
     */
    ~CTurbomachineryDriver(void) override;
    
    /*!
     * \brief Run a single iteration of the physics within multiple zones.
     */
    
    void Run() override;
    
    /*!
     * \brief Set Mixing Plane interface within multiple zones.
     */
    void SetMixingPlane(unsigned short iZone);
    
    /*!
     * \brief Set Mixing Plane interface within multiple zones.
     */
    void SetTurboPerformance(unsigned short targetZone);
    
    /*!
     * \brief Monitor the computation.
     */
    bool Monitor(unsigned long TimeIter) override;
    
    
    
};

/*!
 * \class CHBDriver
 * \brief Class for driving an iteration of Harmonic Balance (HB) method problem using multiple time zones.
 * \author T. Economon
 */
class CHBDriver : public CFluidDriver {
    
private:
    COutputLegacy* output_legacy;
    unsigned short nInstHB;
    su2double **D; /*!< \brief Harmonic Balance operator. */
    
public:
    
    /*!
     * \brief Constructor of the class.
     * \param[in] confFile - Configuration file name.
     * \param[in] val_nZone - Total number of zones.
     * \param[in] val_nDim - Number of dimensions.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CHBDriver(char* confFile,
              unsigned short val_nZone,
              SU2_Comm MPICommunicator);
    
    /*!
     * \brief Destructor of the class.
     */
    ~CHBDriver(void) override;
    
    /*!
     * \brief Run a single iteration of a Harmonic Balance problem.
     */
    void Run() override;
    
    /*!
     * \brief Computation and storage of the Harmonic Balance method source terms.
     * \author T. Economon, K. Naik
     * \param[in] iZone - Current zone number.
     */
    void SetHarmonicBalance(unsigned short iZone);
    
    /*!
     * \brief Precondition Harmonic Balance source term for stability
     * \author J. Howison
     */
    void StabilizeHarmonicBalance();
    
    /*!
     * \brief Computation of the Harmonic Balance operator matrix for harmonic balance.
     * \author A. Rubino, S. Nimmagadda
     */
    void ComputeHB_Operator();
    
    /*!
     * \brief Update the solution for the Harmonic Balance.
     */
    void Update() override;
    
    /*!
     * \brief Reset the convergence flag (set to false) of the solver for the Harmonic Balance.
     */
    void ResetConvergence() override;
};
