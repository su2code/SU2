/*!
 * \file driver_structure.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/mpi_structure.hpp"
#include "../iteration_structure.hpp"
#include "../solver_structure.hpp"
#include "../integration_structure.hpp"

#include "../numerics_structure.hpp"
/*--- Transfer includes ---*/
#include "../interfaces/CInterface.hpp"
#include "../interfaces/cfd/CConservativeVarsInterface.hpp"
#include "../interfaces/cfd/CMixingPlaneInterface.hpp"
#include "../interfaces/cfd/CSlidingInterface.hpp"
#include "../interfaces/cht/CConjugateHeatInterface.hpp"
#include "../interfaces/fsi/CDisplacementsInterface.hpp"
#include "../interfaces/fsi/CFlowTractionInterface.hpp"
#include "../interfaces/fsi/CDiscAdjFlowTractionInterface.hpp"
#include "../interfaces/fsi/CDisplacementsInterfaceLegacy.hpp"
#include "../interfaces/fsi/CDiscAdjDisplacementsInterfaceLegacy.hpp"
#include "../solvers/CDiscAdjMeshSolver.hpp"
#include "../solvers/CMeshSolver.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/grid_movement_structure.hpp"
#include "../../../Common/include/config_structure.hpp"
#include "../../../Common/include/interpolation_structure.hpp"

#include "../output/COutputLegacy.hpp"

#include "../output/COutput.hpp"
#include "../output/CMultizoneOutput.hpp"
#include "../output/CElasticityOutput.hpp"
#include "../output/CAdjElasticityOutput.hpp"
#include "../output/CFlowCompOutput.hpp"
#include "../output/CAdjFlowOutput.hpp"
#include "../output/CFlowCompFEMOutput.hpp"
#include "../output/CFlowIncOutput.hpp"
#include "../output/CAdjFlowIncOutput.hpp"
#include "../output/CHeatOutput.hpp"
#include "../output/CAdjHeatOutput.hpp"

using namespace std;

/*!
 * \class CDriver
 * \brief Parent class for driving an iteration of a single or multi-zone problem.
 * \author T. Economon
 */
class CDriver {
protected:
  int rank,   /*!< \brief MPI Rank. */
  size;         /*!< \brief MPI Size. */
  char* config_file_name;                       /*!< \brief Configuration file name of the problem.*/
  char runtime_file_name[MAX_STRING_SIZE];
  su2double StartTime,                          /*!< \brief Start point of the timer for performance benchmarking.*/
            StopTime,                           /*!< \brief Stop point of the timer for performance benchmarking.*/
            UsedTimePreproc,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking preprocessing phase.*/
            UsedTimeCompute,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking compute phase.*/
            UsedTimeOutput,                     /*!< \brief Elapsed time between Start and Stop point of the timer for tracking output phase.*/
            UsedTime;                           /*!< \brief Elapsed time between Start and Stop point of the timer.*/
  su2double BandwidthSum;                       /*!< \brief Aggregate value of the bandwidth for writing restarts (to be average later).*/
  unsigned long IterCount,                      /*!< \brief Iteration count stored for performance benchmarking.*/
  OutputCount;                                  /*!< \brief Output count stored for performance benchmarking.*/
  unsigned long DOFsPerPoint;                   /*!< \brief Number of unknowns at each vertex, i.e., number of equations solved. */
  su2double Mpoints;                            /*!< \brief Total number of grid points in millions in the calculation (including ghost points).*/
  su2double MpointsDomain;                      /*!< \brief Total number of grid points in millions in the calculation (excluding ghost points).*/
  su2double MDOFs;                              /*!< \brief Total number of DOFs in millions in the calculation (including ghost points).*/
  su2double MDOFsDomain;                        /*!< \brief Total number of DOFs in millions in the calculation (excluding ghost points).*/
  unsigned long TimeIter;                       /*!< \brief External iteration.*/
  ofstream **ConvHist_file;                     /*!< \brief Convergence history file.*/
  ofstream FSIHist_file;                        /*!< \brief FSI convergence history file.*/
  unsigned short iMesh,                         /*!< \brief Iterator on mesh levels.*/
                iZone,                          /*!< \brief Iterator on zones.*/
                nZone,                          /*!< \brief Total number of zones in the problem. */
                nDim,                           /*!< \brief Number of dimensions.*/
                iInst,                          /*!< \brief Iterator on instance levels.*/
                *nInst,                         /*!< \brief Total number of instances in the problem (per zone). */
                **interface_types;              /*!< \brief Type of coupling between the distinct (physical) zones.*/
  bool StopCalc,                                /*!< \brief Stop computation flag.*/
       mixingplane,                             /*!< \brief mixing-plane simulation flag.*/
       fsi,                                     /*!< \brief FSI simulation flag.*/
       fem_solver;                              /*!< \brief FEM fluid solver simulation flag. */
  CIteration ***iteration_container;            /*!< \brief Container vector with all the iteration methods. */
  COutput **output_container;                   /*!< \brief Pointer to the COutput class. */
  CIntegration ****integration_container;       /*!< \brief Container vector with all the integration methods. */
  CGeometry ****geometry_container;             /*!< \brief Geometrical definition of the problem. */
  CSolver *****solver_container;                /*!< \brief Container vector with all the solutions. */
  CNumerics ******numerics_container;           /*!< \brief Description of the numerical method (the way in which the equations are solved). */
  CConfig **config_container;                   /*!< \brief Definition of the particular problem. */
  CConfig *driver_config;                       /*!< \brief Definition of the driver configuration. */
  COutput *driver_output;                       /*!< \brief Definition of the driver output. */
  CSurfaceMovement **surface_movement;          /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement ***grid_movement;         /*!< \brief Volume grid movement classes of the problem. */
  CFreeFormDefBox*** FFDBox;                    /*!< \brief FFD FFDBoxes of the problem. */
  CInterpolator ***interpolator_container;      /*!< \brief Definition of the interpolation method between non-matching discretizations of the interface. */
  CInterface ***interface_container;            /*!< \brief Definition of the interface of information and physics. */
  su2double PyWrapVarCoord[3],                  /*!< \brief This is used to store the VarCoord of each vertex. */
            PyWrapNodalForce[3],                /*!< \brief This is used to store the force at each vertex. */
            PyWrapNodalForceDensity[3],         /*!< \brief This is used to store the force density at each vertex. */
            PyWrapNodalHeatFlux[3];             /*!< \brief This is used to store the heat flux at each vertex. */
  bool dummy_geometry;

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
  void Iteration_Preprocessing(CConfig *config, CIteration *&iteration);

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
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Preprocessing(CConfig *config, CIntegration **&integration);

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
                               unsigned short **interface_types, CInterface ***&interface,
                               CInterpolator ***&interpolation);

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CNumerics ****&numerics);

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
  void DynamicMesh_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CIteration *iteration, CVolumetricMovement *&grid_movement, CSurfaceMovement *&surface_movement);

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
  void StaticMesh_Preprocessing(CConfig *config, CGeometry **geometry, CSurfaceMovement *surface_movement);

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
   * \brief Perform a static mesh deformation, without considering grid velocity.
   */
  virtual void StaticMeshUpdate() { }

  /*!
   * \brief Perform a mesh deformation as initial condition.
   */
  virtual void SetInitialMesh() { }

  /*!
   * \brief Process the boundary conditions and update the multigrid structure.
   */
  virtual void BoundaryConditionsUpdate() { }

  /*!
   * \brief Get the total drag.
   * \return Total drag.
   */
  passivedouble Get_Drag();

  /*!
   * \brief Get the total lift.
   * \return Total lift.
   */
  passivedouble Get_Lift();

  /*!
   * \brief Get the total x moment.
   * \return Total x moment.
   */
  passivedouble Get_Mx();

  /*!
   * \brief Get the total y moment.
   * \return Total y moment.
   */
  passivedouble Get_My();

  /*!
   * \brief Get the total z moment.
   * \return Total z moment.
   */
  passivedouble Get_Mz();

  /*!
   * \brief Get the total drag coefficient.
   * \return Total drag coefficient.
   */
  passivedouble Get_DragCoeff();

  /*!
   * \brief Get the total lift coefficient.
   * \return Total lift coefficient.
   */
  passivedouble Get_LiftCoeff();

  /*!
   * \brief Get the moving marker identifier.
   * \return Moving marker identifier.
   */
  unsigned short GetMovingMarker();

  /*!
   * \brief Get the number of vertices (halo nodes included) from a specified marker.
   * \param[in] iMarker -  Marker identifier.
   * \return Number of vertices.
   */
  unsigned long GetNumberVertices(unsigned short iMarker);

  /*!
   * \brief Get the number of halo vertices from a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \return Number of vertices.
   */
  unsigned long GetNumberHaloVertices(unsigned short iMarker);

  /*!
   * \brief Check if a vertex is physical or not (halo node) on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return True if the specified vertex is a halo node.
   */
  bool IsAHaloNode(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the number of external iterations.
   * \return Number of external iterations.
   */
  unsigned long GetnTimeIter();

  /*!
   * \brief Get the current external iteration.
   * \return Current external iteration.
   */
  unsigned long GetTime_Iter();

  /*!
   * \brief Get the unsteady time step.
   * \return Unsteady time step.
   */
  passivedouble GetUnsteady_TimeStep();

  /*!
   * \brief Get the global index of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vertex global index.
   */
  unsigned long GetVertexGlobalIndex(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the x coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x coordinate of the vertex.
   */
  passivedouble GetVertexCoordX(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the y coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y coordinate of the vertex.
   */
  passivedouble GetVertexCoordY(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the z coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z coordinate of the vertex.
   */
  passivedouble GetVertexCoordZ(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Compute the total force (pressure and shear stress) at a vertex on a specified marker (3 components).
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return True if the vertex is a halo node (non physical force).
   */
  bool ComputeVertexForces(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the x component of the force at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x component of the force at the vertex.
   */
  passivedouble GetVertexForceX(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the y component of the force at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y component of the force at the vertex.
   */
  passivedouble GetVertexForceY(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the z component of the force at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z component of the force at the vertex.
   */
  passivedouble GetVertexForceZ(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the x component of the force density at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x component of the force density at the vertex.
   */
  passivedouble GetVertexForceDensityX(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the y component of the force density at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y component of the force density at the vertex.
   */
  passivedouble GetVertexForceDensityY(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the z component of the force density at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z component of the force density at the vertex.
   */
  passivedouble GetVertexForceDensityZ(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Set the x coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] newPosX - New x coordinate of the vertex.
   */
  void SetVertexCoordX(unsigned short iMarker, unsigned long iVertex, passivedouble newPosX);

  /*!
   * \brief Set the y coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] newPosY - New y coordinate of the vertex.
   */
  void SetVertexCoordY(unsigned short iMarker, unsigned long iVertex, passivedouble newPosY);

  /*!
   * \brief Set the z coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] newPosZ - New z coordinate of the vertex.
   */
  void SetVertexCoordZ(unsigned short iMarker, unsigned long iVertex, passivedouble newPosZ);

  /*!
   * \brief Set the VarCoord of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Norm of the VarCoord.
   */
  passivedouble SetVertexVarCoord(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the temperature at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Temperature of the vertex.
   */
  passivedouble GetVertexTemperature(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Set the temperature of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_WallTemp - Value of the temperature.
   */
  void SetVertexTemperature(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallTemp);

  /*!
   * \brief Compute the heat flux at a vertex on a specified marker (3 components).
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return True if the vertex is a halo node.
   */
  bool ComputeVertexHeatFluxes(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the x component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x component of the heat flux at the vertex.
   */
  passivedouble GetVertexHeatFluxX(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the y component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y component of the heat flux at the vertex.
   */
  passivedouble GetVertexHeatFluxY(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the z component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z component of the heat flux at the vertex.
   */
  passivedouble GetVertexHeatFluxZ(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the wall normal component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Wall normal component of the heat flux at the vertex.
   */
  passivedouble GetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Set the wall normal component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_WallHeatFlux - Value of the normal heat flux.
   */
  void SetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallHeatFlux);

  /*!
   * \brief Get the thermal conductivity at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Thermal conductivity at the vertex.
   */
  passivedouble GetThermalConductivity(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Preprocess the inlets via file input for all solvers.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Inlet_Preprocessing(CSolver ***solver, CGeometry **geometry,
                                    CConfig *config);

  /*!
   * \brief Get the unit normal (vector) at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Unit normal (vector) at the vertex.
   */
  vector<passivedouble> GetVertexUnitNormal(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get all the boundary markers tags.
   * \return List of boundary markers tags.
   */
  vector<string> GetAllBoundaryMarkersTag();

  /*!
   * \brief Get all the moving boundary markers tags.
   * \return List of moving boundary markers tags.
   */
  vector<string> GetAllMovingMarkersTag();

  /*!
   * \brief Get all the deformable boundary marker tags.
   * \return List of deformable boundary markers tags.
   */
  vector<string> GetAllDeformMeshMarkersTag();

  /*!
   * \brief Get all the fluid load boundary marker tags.
   * \return List of fluid load boundary markers tags.
   */
  vector<string> GetAllFluidLoadMarkersTag();

  /*!
   * \brief Get all the heat transfer boundary markers tags.
   * \return List of heat transfer boundary markers tags.
   */
  vector<string> GetAllCHTMarkersTag();

  /*!
   * \brief Get all the (subsonic) inlet boundary markers tags.
   * \return List of inlet boundary markers tags.
   */
  vector<string> GetAllInletMarkersTag();

  /*!
   * \brief Get all the boundary markers tags with their associated indices.
   * \return List of boundary markers tags with their indices.
   */
  map<string, int> GetAllBoundaryMarkers();

  /*!
   * \brief Get all the boundary markers tags with their associated types.
   * \return List of boundary markers tags with their types.
   */
  map<string, string> GetAllBoundaryMarkersType();

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

  /*!
   * \brief Return the sensitivities of the mesh boundary vertices.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vector of sensitivities.
   */
  vector<passivedouble> GetMeshDisp_Sensitivity(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Set the load in X direction for the structural solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] LoadX - Value of the load in the direction X.
   * \param[in] LoadX - Value of the load in the direction Y.
   * \param[in] LoadX - Value of the load in the direction Z.
   */
  void SetFEA_Loads(unsigned short iMarker, unsigned long iVertex, passivedouble LoadX,
                    passivedouble LoadY, passivedouble LoadZ);

  /*!
   * \brief Return the displacements from the FEA solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vector of displacements.
   */
  vector<passivedouble> GetFEA_Displacements(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Return the velocities from the FEA Solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vector of velocities.
   */
  vector<passivedouble> GetFEA_Velocity(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Return the velocities from the FEA Solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vector of velocities at time n.
   */
  vector<passivedouble> GetFEA_Velocity_n(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the sensitivity of the flow loads for the structural solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] LoadX - Value of the load in the direction X.
   * \param[in] LoadX - Value of the load in the direction Y.
   * \param[in] LoadX - Value of the load in the direction Z.
   */
  vector<passivedouble> GetFlowLoad_Sensitivity(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Get the flow load (from the extra step - the repeated methods should be unified once the postprocessing
   * strategy is in place).
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   */
  vector<passivedouble> GetFlowLoad(unsigned short iMarker, unsigned long iVertex);

  /*!
   * \brief Set the adjoint of the flow tractions (from the extra step -
   * the repeated methods should be unified once the postprocessing strategy is in place).
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_AdjointX - Value of the adjoint in the direction X.
   * \param[in] val_AdjointY - Value of the adjoint in the direction Y.
   * \param[in] val_AdjointZ - Value of the adjoint in the direction Z.
   */
  void SetFlowLoad_Adjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                    passivedouble val_AdjointY, passivedouble val_AdjointZ);

  /*!
   * \brief Set the adjoint of the structural displacements (from an outside source)
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_AdjointX - Value of the adjoint in the direction X.
   * \param[in] val_AdjointY - Value of the adjoint in the direction Y.
   * \param[in] val_AdjointZ - Value of the adjoint in the direction Z.
   */
  void SetSourceTerm_DispAdjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                 passivedouble val_AdjointY, passivedouble val_AdjointZ);

  /*!
   * \brief Get the undeformed mesh coordinates
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Undeformed Vertex Coordinates
   */
  vector<passivedouble> GetVertex_UndeformedCoord(unsigned short iMarker, unsigned long iVertex);

};

/*!
 * \class CFluidDriver
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon, G. Gori
 */
class CFluidDriver : public CDriver {

protected:
   unsigned long Max_Iter;

public:

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

  /*!
   * \brief Destructor of the class.
   */
  ~CFluidDriver(void);

  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  void StartSolver();

  /*!
   * \brief Run a single iteration of the physics within multiple zones.
   */
  void Run();

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update();

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long InnerIter);

  /*!
   * \brief Monitor the computation.
   */
  bool Monitor(unsigned long ExtIter);

  /*!
   * \brief Perform some pre-processing before an iteration of the physics.
   */
  void Preprocess(unsigned long Iter);

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure (multiple zone).
   */
  void DynamicMeshUpdate(unsigned long TimeIter);

  /*!
   * \brief Perform a static mesh deformation, without considering grid velocity (multiple zone).
   */
  void StaticMeshUpdate();

  /*!
   * \brief Perform a mesh deformation as initial condition (multiple zone).
   */
  void SetInitialMesh();

  /*!
   * \brief Process the boundary conditions and update the multigrid structure.
   */
  void BoundaryConditionsUpdate();

  /*!
   * \brief Transfer data among different zones (multiple zone).
   */
  void Transfer_Data(unsigned short donorZone, unsigned short targetZone);

  /*!
   * \brief Set the total temperature of a vertex on a specified inlet marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_Ttotal - Value of the total (stagnation) temperature.
   */
  void SetVertexTtotal(unsigned short iMarker, unsigned long iVertex, passivedouble val_Ttotal);

  /*!
   * \brief Set the total pressure of a vertex on a specified inlet marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_Ptotal - Value of the total (stagnation) pressure.
   */
  void SetVertexPtotal(unsigned short iMarker, unsigned long iVertex, passivedouble val_Ptotal);

  /*!
   * \brief Set the flow direction of a vertex on a specified inlet marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the flow direction unit vector
   * \param[in] val_FlowDir - Component of a unit vector representing the flow direction
   */
  void SetVertexFlowDir(unsigned short iMarker, unsigned long iVertex, unsigned short iDim, passivedouble val_FlowDir);

  /*!
   * \brief Set a turbulence variable on a specified inlet marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the turbulence variable (i.e. k is 0 in SST)
   * \param[in] val_turb_var - Value of the turbulence variable to be used.
   */
  void SetVertexTurbVar(unsigned short iMarker, unsigned long iVertex, unsigned short iDim, passivedouble val_tub_var);

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
  ~CTurbomachineryDriver(void);

  /*!
   * \brief Run a single iteration of the physics within multiple zones.
   */

  void Run();

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
  bool Monitor(unsigned long TimeIter);



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
  ~CHBDriver(void);

  /*!
   * \brief Run a single iteration of a Harmonic Balance problem.
   */
  void Run();

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
  void Update();

  /*!
   * \brief Reset the convergence flag (set to false) of the solver for the Harmonic Balance.
   */
  void ResetConvergence();
};

/*!
 * \class CDiscAdjFSIDriver
 * \brief Overload: Class for driving a discrete adjoint FSI iteration.
 * \author R. Sanchez.
 * \version 7.0.0 "Blackbird"
 */
class CDiscAdjFSIDriver : public CDriver {

  COutputLegacy* output_legacy;

  CIteration** direct_iteration;
  unsigned short RecordingState;
  unsigned short CurrentRecording;          /*!< \brief Stores the current status of the recording. */
  unsigned short Kind_Objective_Function;   /*!< \brief Stores the kind of objective function of the recording. */

  su2double *init_res_flow,     /*!< \brief Stores the initial residual for the flow. */
            *init_res_struct,   /*!< \brief Stores the initial residual for the structure. */
            *residual_flow,     /*!< \brief Stores the current residual for the flow. */
            *residual_struct,   /*!< \brief Stores the current residual for the structure. */
            *residual_flow_rel,
            *residual_struct_rel;

  su2double flow_criteria,
            flow_criteria_rel,
            structure_criteria,
            structure_criteria_rel;


  enum OF_KIND{
    NO_OBJECTIVE_FUNCTION = 0,               /*!< \brief Indicates that there is no objective function. */
    FLOW_OBJECTIVE_FUNCTION = 1,            /*!< \brief Indicates that the objective function is only flow-dependent. */
    FEM_OBJECTIVE_FUNCTION = 2              /*!< \brief Indicates that the objective function is only structural-dependent. */
  };

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Total number of dimensions.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDiscAdjFSIDriver(char* confFile,
                    unsigned short val_nZone,
                    SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFSIDriver(void);

  /*!
   * \brief Launch the computation for FSI adjoint (legacy) driver
   */
  inline void StartSolver(){

      /*--- Run the solver. ---*/
      if (rank == MASTER_NODE)
        cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
      Run();
  }

  /*!
   * \brief Run a Discrete Adjoint iteration for the FSI problem.
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   */

  void Run();

  /*!
   * \brief Iterate the direct solver for recording.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */

  void Iterate_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT, unsigned short kind_recording);

  /*!
   * \brief Run a direct flow iteration.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   */
  void Fluid_Iteration_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT);

  /*!
   * \brief Run a direct structural iteration.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   */
  void Structural_Iteration_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT);

  /*!
   * \brief Run a direct mesh deformation.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   */
  void Mesh_Deformation_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT);

  /*!
   * \brief Set the recording for a Discrete Adjoint iteration for the FSI problem.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */

  void SetRecording(unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Load the restarts for fluid, structure and mesh.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void Preprocess(unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Iterate a certain block for adjoint FSI - may be the whole set of variables or independent and subiterate
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void Iterate_Block(unsigned short ZONE_FLOW,
                       unsigned short ZONE_STRUCT,
                       unsigned short kind_recording);

  /*!
   * \brief Initialize the adjoint - set the objective funcition and the output of the adjoint iteration
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void InitializeAdjoint(unsigned short ZONE_FLOW,
                            unsigned short ZONE_STRUCT,
                            unsigned short kind_recording);

  /*!
   * \brief Extract the adjoint solution variables
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void ExtractAdjoint(unsigned short ZONE_FLOW,
                         unsigned short ZONE_STRUCT,
                         unsigned short kind_recording);


  /*!
   * \brief Check the convergence of the problem
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  bool CheckConvergence(unsigned long IntIter,
                          unsigned short ZONE_FLOW,
                          unsigned short ZONE_STRUCT,
                          unsigned short kind_recording);

  /*!
   * \brief Check the convergence of BGS subiteration process
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  bool BGSConvergence(unsigned long IntIter,
                         unsigned short ZONE_FLOW,
                         unsigned short ZONE_STRUCT);


  /*!
   * \brief Output the convergence history
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void ConvergenceHistory(unsigned long IntIter,
                             unsigned long nIntIter,
                             unsigned short ZONE_FLOW,
                             unsigned short ZONE_STRUCT,
                             unsigned short kind_recording);

  /*!
   * \brief Load the restarts for fluid, structure and mesh.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void PrintDirect_Residuals(unsigned short ZONE_FLOW,
                                unsigned short ZONE_STRUCT,
                                unsigned short kind_recording);

  /*!
   * \brief Restart the variables to the converged solution.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void PrepareRecording(unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Register the input variables for adjoint FSI problems: flow conservative, fluid mesh position and structural displacements.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void RegisterInput(unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Register the input variables for adjoint FSI problems: flow conservative, fluid mesh position and structural displacements.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void SetDependencies(unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Restart the output variables for adjoint FSI problems: flow conservative, fluid mesh position and structural displacements.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void RegisterOutput(unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Overload, does nothing but avoids dynamic mesh updates in adjoint FSI problems before the iteration
   */
  void DynamicMeshUpdate(unsigned long TimeIter);

  /*!
   * \brief Transfer the displacements computed on the structural solver into the fluid solver.
   * \param[in] donorZone - zone in which the displacements will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  void Transfer_Displacements(unsigned short donorZone, unsigned short targetZone);

  /*!
   * \brief Transfer the tractions computed on the fluid solver into the structural solver.
   * \param[in] donorZone - zone from which the tractions will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  void Transfer_Tractions(unsigned short donorZone, unsigned short targetZone);

};
