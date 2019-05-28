/*!
 * \file driver_structure.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../Common/include/mpi_structure.hpp"
#include "iteration_structure.hpp"
#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "transfer_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/interpolation_structure.hpp"

using namespace std;

/*! 
 * \class CDriver
 * \brief Parent class for driving an iteration of a single or multi-zone problem.
 * \author T. Economon
 */
class CDriver {
protected:
  int rank, 	/*!< \brief MPI Rank. */
  size;       	/*!< \brief MPI Size. */
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
  su2double MDOFs;                              /*!< \brief Total number of DOFs in millions in the calculation (including ghost points).*/
  su2double MDOFsDomain;                        /*!< \brief Total number of DOFs in millions in the calculation (excluding ghost points).*/
  unsigned long ExtIter;                        /*!< \brief External iteration.*/
  ofstream **ConvHist_file;                       /*!< \brief Convergence history file.*/
  ofstream FSIHist_file;                        /*!< \brief FSI convergence history file.*/
  unsigned short iMesh,                         /*!< \brief Iterator on mesh levels.*/
                iZone,                          /*!< \brief Iterator on zones.*/
                nZone,                          /*!< \brief Total number of zones in the problem. */
                nDim,                           /*!< \brief Number of dimensions.*/
                iInst,                          /*!< \brief Iterator on instance levels.*/
                *nInst,                         /*!< \brief Total number of instances in the problem (per zone). */
                **transfer_types;               /*!< \brief Type of coupling between the distinct (physical) zones.*/
  bool StopCalc,                                /*!< \brief Stop computation flag.*/
       mixingplane,                             /*!< \brief mixing-plane simulation flag.*/
       fsi,                                     /*!< \brief FSI simulation flag.*/
       fem_solver;                              /*!< \brief FEM fluid solver simulation flag. */
  CIteration ***iteration_container;             /*!< \brief Container vector with all the iteration methods. */
  COutput *output;                              /*!< \brief Pointer to the COutput class. */
  CIntegration ****integration_container;        /*!< \brief Container vector with all the integration methods. */
  CGeometry ****geometry_container;              /*!< \brief Geometrical definition of the problem. */
  CSolver *****solver_container;                 /*!< \brief Container vector with all the solutions. */
  CNumerics ******numerics_container;            /*!< \brief Description of the numerical method (the way in which the equations are solved). */
  CConfig **config_container;                   /*!< \brief Definition of the particular problem. */
  CConfig *driver_config;                       /*!< \brief Definition of the driver configuration. */
  CSurfaceMovement **surface_movement;          /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement ***grid_movement;          /*!< \brief Volume grid movement classes of the problem. */
  CFreeFormDefBox*** FFDBox;                    /*!< \brief FFD FFDBoxes of the problem. */
  CInterpolator ***interpolator_container;      /*!< \brief Definition of the interpolation method between non-matching discretizations of the interface. */
  CTransfer ***transfer_container;              /*!< \brief Definition of the transfer of information and the physics involved in the interface. */
  su2double PyWrapVarCoord[3],                  /*!< \brief This is used to store the VarCoord of each vertex. */
            PyWrapNodalForce[3],                /*!< \brief This is used to store the force at each vertex. */
            PyWrapNodalForceDensity[3],         /*!< \brief This is used to store the force density at each vertex. */
            PyWrapNodalHeatFlux[3];             /*!< \brief This is used to store the heat flux at each vertex. */

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
          unsigned short val_nDim,
          SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CDriver(void);

  /*!
   * \brief A virtual member.
   */  
  virtual void Run() { };

  /*!
   * \brief Read in the config and mesh files.
   */
  void Input_Preprocessing(SU2_Comm MPICommunicator);

  /*!
   * \brief Construction of the edge-based data structure and the multigrid structure.
   */
  void Geometrical_Preprocessing();
  
  /*!
   * \brief Do the geometrical preprocessing for the DG FEM solver.
   */
  void Geometrical_Preprocessing_DGFEM();
  
  /*!
   * \brief Definition of the physics iteration class or within a single zone.
   * \param[in] iteration_container - Pointer to the iteration container to be instantiated.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void Iteration_Preprocessing();

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Preprocessing(CSolver ****solver_container, CGeometry ***geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Restart of the solvers from the restart files.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Restart(CSolver ****solver_container, CGeometry ***geometry, CConfig *config, bool update_geo, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Postprocessing(CSolver ****solver_container, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Preprocessing(CIntegration ***integration_container, CGeometry ***geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Postprocessing(CIntegration ***integration_container, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all interface classes.
   */
  void Interface_Preprocessing();

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Preprocessing(CNumerics *****numerics_container, CSolver ****solver_container, CGeometry ***geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Postprocessing(CNumerics *****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Initialize Python interface functionalities
   */
  void PythonInterface_Preprocessing();

  /*!
   * \brief Deallocation routine
   */
  void Postprocessing();

  /*!
   * \brief Initiate value for static mesh movement such as the gridVel for the ROTATING frame.
   */
  void InitStaticMeshMovement();

  /*!
   * \brief Initiate value for static mesh movement such as the gridVel for the ROTATING frame.
   */
  void TurbomachineryPreprocessing(void);

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the displacements will be predicted.
   * \param[in] targetZone - zone which receives the predicted displacements.
   */
  virtual void Predict_Displacements(unsigned short donorZone, unsigned short targetZone) {};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the tractions will be predicted.
   * \param[in] targetZone - zone which receives the predicted traction.
   */
  virtual void Predict_Tractions(unsigned short donorZone, unsigned short targetZone) {};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the displacements will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void Transfer_Displacements(unsigned short donorZone, unsigned short targetZone) {};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone from which the tractions will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void Transfer_Tractions(unsigned short donorZone, unsigned short targetZone) {};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Displacements(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter) {};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Tractions(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter) {};

  /*!
   * \brief A virtual member.
   */
  virtual void Update() {};

  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  virtual void StartSolver();

  /*!
   * \brief A virtual member.
   */
  virtual void ResetConvergence() { };

  /*!
   * \brief Perform some pre-processing before an iteration of the physics.
   */
  void PreprocessExtIter(unsigned long ExtIter);

  /*!
   * \brief Monitor the computation.
   */
  virtual bool Monitor(unsigned long ExtIter);

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long ExtIter);

  /*!
   * \brief Perform a dynamic mesh deformation, including grid velocity computation and update of the multigrid structure.
   */
  virtual void DynamicMeshUpdate(unsigned long ExtIter) { };

  /*!
   * \brief Perform a dynamic mesh deformation, including grid velocity computation and update of the multigrid structure.
   */
  virtual void DynamicMeshUpdate(unsigned short val_iZone, unsigned long ExtIter) { };

  /*!
   * \brief Perform a static mesh deformation, without considering grid velocity.
   */
  virtual void StaticMeshUpdate() { };

  /*!
   * \brief Perform a mesh deformation as initial condition.
   */
  virtual void SetInitialMesh() { };

  /*!
   * \brief Process the boundary conditions and update the multigrid structure.
   */
  virtual void BoundaryConditionsUpdate() { };

  /*!
   * \brief Get the total drag.
   * \return Total drag.
   */
  su2double Get_Drag();

  /*!
   * \brief Get the total lift.
   * \return Total lift.
   */
  su2double Get_Lift();

  /*!
   * \brief Get the total x moment.
   * \return Total x moment.
   */
  su2double Get_Mx();

  /*!
   * \brief Get the total y moment.
   * \return Total y moment.
   */
  su2double Get_My();

  /*!
   * \brief Get the total z moment.
   * \return Total z moment.
   */
  su2double Get_Mz();

  /*!
   * \brief Get the total drag coefficient.
   * \return Total drag coefficient.
   */
  su2double Get_DragCoeff();

  /*!
   * \brief Get the total lift coefficient.
   * \return Total lift coefficient.
   */
  su2double Get_LiftCoeff();

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
  bool IsAHaloNode(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the number of external iterations.
   * \return Number of external iterations.
   */
  unsigned long GetnExtIter();

  /*!
   * \brief Get the current external iteration.
   * \return Current external iteration.
   */
  unsigned long GetExtIter();

  /*!
   * \brief Get the unsteady time step.
   * \return Unsteady time step.
   */
  su2double GetUnsteady_TimeStep();

  /*!
   * \brief Get the global index of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vertex global index.
   */
  unsigned long GetVertexGlobalIndex(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the x coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x coordinate of the vertex.
   */
  su2double GetVertexCoordX(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the y coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y coordinate of the vertex.
   */
  su2double GetVertexCoordY(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the z coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z coordinate of the vertex.
   */
  su2double GetVertexCoordZ(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Compute the total force (pressure and shear stress) at a vertex on a specified marker (3 components).
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return True if the vertex is a halo node (non physical force).
   */
  bool ComputeVertexForces(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the x component of the force at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x component of the force at the vertex.
   */
  su2double GetVertexForceX(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the y component of the force at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y component of the force at the vertex.
   */
  su2double GetVertexForceY(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the z component of the force at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z component of the force at the vertex.
   */
  su2double GetVertexForceZ(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the x component of the force density at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x component of the force density at the vertex.
   */
  su2double GetVertexForceDensityX(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the y component of the force density at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y component of the force density at the vertex.
   */
  su2double GetVertexForceDensityY(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the z component of the force density at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z component of the force density at the vertex.
   */
  su2double GetVertexForceDensityZ(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Set the x coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] newPosX - New x coordinate of the vertex.
   */
  void SetVertexCoordX(unsigned short iMarker, unsigned short iVertex, su2double newPosX);

  /*!
   * \brief Set the y coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] newPosY - New y coordinate of the vertex.
   */
  void SetVertexCoordY(unsigned short iMarker, unsigned short iVertex, su2double newPosY);

  /*!
   * \brief Set the z coordinate of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] newPosZ - New z coordinate of the vertex.
   */
  void SetVertexCoordZ(unsigned short iMarker, unsigned short iVertex, su2double newPosZ);

  /*!
   * \brief Set the VarCoord of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Norm of the VarCoord.
   */
  su2double SetVertexVarCoord(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the temperature at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Temperature of the vertex.
   */
  su2double GetVertexTemperature(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Set the temperature of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_WallTemp - Value of the temperature.
   */
  void SetVertexTemperature(unsigned short iMarker, unsigned short iVertex, su2double val_WallTemp);

  /*!
   * \brief Compute the heat flux at a vertex on a specified marker (3 components).
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return True if the vertex is a halo node.
   */
  bool ComputeVertexHeatFluxes(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the x component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return x component of the heat flux at the vertex.
   */
  su2double GetVertexHeatFluxX(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the y component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return y component of the heat flux at the vertex.
   */
  su2double GetVertexHeatFluxY(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the z component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return z component of the heat flux at the vertex.
   */
  su2double GetVertexHeatFluxZ(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Get the wall normal component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Wall normal component of the heat flux at the vertex.
   */
  su2double GetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Set the wall normal component of the heat flux at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_WallHeatFlux - Value of the normal heat flux.
   */
  void SetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex, su2double val_WallHeatFlux);

  /*!
   * \brief Get the thermal conductivity at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Thermal conductivity at the vertex.
   */
  su2double GetThermalConductivity(unsigned short iMarker, unsigned short iVertex);

  /*!
   * \brief Preprocess the inlets via file input for all solvers.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Inlet_Preprocessing(CSolver ***solver_container, CGeometry **geometry,
                                    CConfig *config);

  /*!
   * \brief Get the unit normal (vector) at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Unit normal (vector) at the vertex.
   */
  vector<su2double> GetVertexUnitNormal(unsigned short iMarker, unsigned short iVertex);

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
   * \brief A virtual member to run a Block Gauss-Seidel iteration in multizone problems.
   */
  virtual void Run_GaussSeidel(){};

  /*!
   * \brief A virtual member to run a Block-Jacobi iteration in multizone problems.
   */
  virtual void Run_Jacobi(){};

};

/*!
 * \class CFluidDriver
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon, G. Gori
 */
class CFluidDriver : public CDriver {
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
               unsigned short val_nDim,
               SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CFluidDriver(void);

  /*!
   * \brief Run a single iteration of the physics within multiple zones.
   */
  void Run();

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update();

  /*!
   * \brief Reset the convergence flag (set to false) of the multizone solver.
   */
  void ResetConvergence();

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure (multiple zone).
   */
  void DynamicMeshUpdate(unsigned long ExtIter);

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
  void SetVertexTtotal(unsigned short iMarker, unsigned short iVertex, su2double val_Ttotal);

  /*!
   * \brief Set the total pressure of a vertex on a specified inlet marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_Ptotal - Value of the total (stagnation) pressure.
   */
  void SetVertexPtotal(unsigned short iMarker, unsigned short iVertex, su2double val_Ptotal);

  /*!
   * \brief Set the flow direction of a vertex on a specified inlet marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the flow direction unit vector
   * \param[in] val_FlowDir - Component of a unit vector representing the flow direction
   */
  void SetVertexFlowDir(unsigned short iMarker, unsigned short iVertex, unsigned short iDim, su2double val_FlowDir);

  /*!
   * \brief Set a turbulence variable on a specified inlet marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] iDim - Index of the turbulence variable (i.e. k is 0 in SST)
   * \param[in] val_turb_var - Value of the turbulence variable to be used.
   */
  void SetVertexTurbVar(unsigned short iMarker, unsigned short iVertex, unsigned short iDim, su2double val_tub_var);

};


/*!
 * \class CTurbomachineryDriver
 * \brief Class for driving an iteration for turbomachinery flow analysis.
 * \author S. Vitale
 */
class CTurbomachineryDriver : public CFluidDriver {
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
                        unsigned short val_nDim,
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
  bool Monitor(unsigned long ExtIter);



};

/*!
 * \class CHBDriver
 * \brief Class for driving an iteration of Harmonic Balance (HB) method problem using multiple time zones.
 * \author T. Economon
 */
class CHBDriver : public CDriver {

private:

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
            unsigned short val_nDim,
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
 * \class CFSIDriver
 * \brief Class for driving a BGS iteration for a fluid-structure interaction problem in multiple zones.
 * \author R. Sanchez.
 */
class CFSIDriver : public CDriver {

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

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CFSIDriver(char* confFile,
             unsigned short val_nZone,
             unsigned short val_nDim,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CFSIDriver(void);

  /*!
   * \brief Run a Block Gauss-Seidel iteration of the FSI problem.
   */
  void Run();

  /*!
   * \brief Predict the structural displacements to pass them into the fluid solver on a BGS implementation.
   * \param[in] donorZone - zone in which the displacements will be predicted.
   * \param[in] targetZone - zone which receives the predicted displacements.
   */
  void Predict_Displacements(unsigned short donorZone, unsigned short targetZone);

  /*!
   * \brief Predict the fluid tractions to pass them into the structural solver on a BGS implementation.
   * \param[in] donorZone - zone in which the tractions will be predicted.
   * \param[in] targetZone - zone which receives the predicted traction.
   */
  void Predict_Tractions(unsigned short donorZone, unsigned short targetZone);

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

  /*!
   * \brief Apply a relaxation method into the computed displacements.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
   */
  void Relaxation_Displacements(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter);

  /*!
   * \brief Apply a relaxation method into the computed tractions.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iOuterIter - Fluid-Structure Interaction subiteration.
   */
  void Relaxation_Tractions(unsigned short donorZone, unsigned short targetZone, unsigned long iOuterIter);

  /*!
   * \brief Check the convergence of BGS subiteration process
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  bool BGSConvergence(unsigned long IntIter, unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT);

  /*!
   * \brief Enforce the coupling condition at the end of the time step
   */
  void Update(void);

  /*!
   * \brief Overload, does nothing but avoids dynamic mesh updates in FSI problems before the iteration
   */
  void DynamicMeshUpdate(unsigned long ExtIter);

};

/*!
 * \class CDiscAdjFSIDriver
 * \brief Overload: Class for driving a discrete adjoint FSI iteration.
 * \author R. Sanchez.
 * \version 6.2.0 "Falcon"
 */
class CDiscAdjFSIDriver : public CDriver {

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
                    unsigned short val_nDim,
                    SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFSIDriver(void);

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
   * \brief Run the post-processing routines.
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   */
  void Postprocess(unsigned short ZONE_FLOW,
                     unsigned short ZONE_STRUCT);

  /*!
   * \brief Overload, does nothing but avoids updates in adjoint FSI problems before the iteration
   */
  void Update(void);

  /*!
   * \brief Overload, does nothing but avoids dynamic mesh updates in adjoint FSI problems before the iteration
   */
  void DynamicMeshUpdate(unsigned long ExtIter);

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

/*!
 * \class CMultiphysicsZonalDriver
 * \brief Class for driving zone-specific iterations.
 * \author O. Burghardt
 * \version 6.2.0 "Falcon"
 */
class CMultiphysicsZonalDriver : public CDriver {
protected:

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CMultiphysicsZonalDriver(char* confFile,
                           unsigned short val_nZone,
                           unsigned short val_nDim,
                           SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CMultiphysicsZonalDriver(void);

  /*!
   * \brief Run one iteration in all physical zones.
   */
  void Run();

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update();

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure (multiple zone).
   */
  void DynamicMeshUpdate(unsigned long ExtIter);

  /*!
   * \brief Routine to provide all the desired physical transfers between the different zones during one iteration.
   */
  void Transfer_Data(unsigned short donorZone, unsigned short targetZone);

};

/*!
 * \class CSinglezoneDriver
 * \brief Class for driving single-zone solvers.
 * \author R. Sanchez
 * \version 6.0.1 "Falcon"
 */
class CSinglezoneDriver : public CDriver {
protected:

  unsigned long TimeIter;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CSinglezoneDriver(char* confFile,
             unsigned short val_nZone,
             unsigned short val_nDim,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CSinglezoneDriver(void);

  /*!
   * \brief [Overload] Launch the computation for single-zone problems.
   */
  void StartSolver();

  /*!
   * \brief Preprocess the single-zone iteration
   */
  virtual void Preprocess(unsigned long TimeIter);

  /*!
   * \brief Run the iteration for ZONE_0.
   */
  virtual void Run();

  /*!
   * \brief Postprocess the iteration for ZONE_0.
   */
  virtual void Postprocess();

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update();

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long ExtIter);

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure.
   */
  void DynamicMeshUpdate(unsigned long ExtIter);


};

/*!
 * \class CDiscAdjSinglezoneDriver
 * \brief Class for driving single-zone adjoint solvers.
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 */
class CDiscAdjSinglezoneDriver : public CSinglezoneDriver {
protected:

  unsigned long nAdjoint_Iter;                  /*!< \brief The number of adjoint iterations that are run on the fixed-point solver.*/
  unsigned short RecordingState;                /*!< \brief The kind of recording the tape currently holds.*/
  unsigned short MainVariables,                 /*!< \brief The kind of recording linked to the main variables of the problem.*/
                 SecondaryVariables;            /*!< \brief The kind of recording linked to the secondary variables of the problem.*/
  su2double ObjFunc;                            /*!< \brief The value of the objective function.*/
  CIteration* direct_iteration;                 /*!< \brief A pointer to the direct iteration.*/

  CConfig *config;                              /*!< \brief Definition of the particular problem. */
  CIteration *iteration;                        /*!< \brief Container vector with all the iteration methods. */
  CIntegration **integration;                   /*!< \brief Container vector with all the integration methods. */
  CGeometry *geometry;                          /*!< \brief Geometrical definition of the problem. */
  CSolver **solver;                             /*!< \brief Container vector with all the solutions. */


public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Total number of dimensions.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDiscAdjSinglezoneDriver(char* confFile,
             unsigned short val_nZone,
             unsigned short val_nDim,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjSinglezoneDriver(void);

  /*!
   * \brief Preprocess the single-zone iteration
   * \param[in] TimeIter - index of the current time-step.
   */
  void Preprocess(unsigned long TimeIter);

  /*!
   * \brief Run a single iteration of the discrete adjoint solver with a single zone.
   */
  void Run(void);

  /*!
   * \brief Postprocess the adjoint iteration for ZONE_0.
   */
  void Postprocess(void);

  /*!
   * \brief Record one iteration of a flow iteration in within multiple zones.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void SetRecording(unsigned short kind_recording);

  /*!
   * \brief Run one iteration of the solver.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void DirectRun(unsigned short kind_recording);

  /*!
   * \brief Set the objective function.
   */
  void SetObjFunction(void);

  /*!
   * \brief Initialize the adjoint value of the objective function.
   */
  void SetAdj_ObjFunction(void);

  /*!
   * \brief Print out the direct residuals.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void Print_DirectResidual(unsigned short kind_recording);

  /*!
   * \brief Record the main computational path.
   */
  void MainRecording(void);

  /*!
   * \brief Record the secondary computational path.
   */
  void SecondaryRecording(void);

};


/*!
 * \class CMultizoneDriver
 * \brief Class for driving zone-specific iterations.
 * \author R. Sanchez, O. Burghardt
 * \version 6.0.1 "Falcon"
 */
class CMultizoneDriver : public CDriver {
protected:

  bool fsi;
  bool cht;

  unsigned long TimeIter;

  unsigned short *nVarZone;
  su2double **init_res,      /*!< \brief Stores the initial residual. */
            **residual,      /*!< \brief Stores the current residual. */
            **residual_rel;  /*!< \brief Stores the residual relative to the initial. */

  su2double flow_criteria,
            flow_criteria_rel,
            structure_criteria,
            structure_criteria_rel;

  bool *prefixed_motion;     /*!< \brief Determines if a fixed motion is imposed in the config file. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CMultizoneDriver(char* confFile,
             unsigned short val_nZone,
             unsigned short val_nDim,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CMultizoneDriver(void);

  /*!
   * \brief [Overload] Launch the computation for multizone problems.
   */
  void StartSolver();

  /*!
   * \brief Preprocess the multizone iteration
   */
  void Preprocess(unsigned long TimeIter);

  /*!
   * \brief Use a corrector step to prevent convergence issues.
   */
  void Corrector(unsigned short val_iZone);

  /*!
   * \brief Run a Block Gauss-Seidel iteration in all physical zones.
   */
  void Run_GaussSeidel();

  /*!
   * \brief Run a Block-Jacobi iteration in all physical zones.
   */
  void Run_Jacobi();

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update();

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long TimeIter);

  /*!
   * \brief Check the convergence at the outer level.
   */
  bool OuterConvergence(unsigned long OuterIter);

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure (multiple zone).
   */
  void DynamicMeshUpdate(unsigned long ExtIter);

  /*!
   * \brief Perform a dynamic mesh deformation, including grid velocity computation and update of the multigrid structure.
   */
  void DynamicMeshUpdate(unsigned short val_iZone, unsigned long ExtIter);

  /*!
   * \brief Routine to provide all the desired physical transfers between the different zones during one iteration.
   * \return Boolean that determines whether the mesh needs to be updated for this particular transfer
   */
  bool Transfer_Data(unsigned short donorZone, unsigned short targetZone);

};

