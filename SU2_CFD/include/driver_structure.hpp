/*!
 * \file driver_structure.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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
 * \version 4.3.0 "Cardinal"
 */
class CDriver {
protected:
  char* config_file_name;                       /*!< \brief Configuration file name of the problem.*/
  char runtime_file_name[MAX_STRING_SIZE];
  su2double StartTime,                          /*!< \brief Start point of the timer for performance benchmarking.*/
            StopTime,                           /*!< \brief Stop point of the timer for performance benchmarking.*/
            UsedTime;                           /*!< \brief Elapsed time between Start and Stop point of the timer.*/
  unsigned long ExtIter;                        /*!< \brief External iteration.*/
  ofstream ConvHist_file;                       /*!< \brief Convergence history file.*/
  unsigned short iMesh,                         /*!< \brief Iterator on mesh levels.*/
                iZone,                          /*!< \brief Iterator on zones.*/
                nZone,                          /*!< \brief Total number of zones in the problem. */
                nDim;                           /*!< \brief Number of dimensions.*/
  bool StopCalc,                                /*!< \brief Stop computation flag.*/
       fsi;                                     /*!< \brief FSI simulation flag.*/
  CIteration **iteration_container;             /*!< \brief Container vector with all the iteration methods. */
  COutput *output;                              /*!< \brief Pointer to the COutput class. */
  CIntegration ***integration_container;        /*!< \brief Container vector with all the integration methods. */
  CGeometry ***geometry_container;              /*!< \brief Geometrical definition of the problem. */
  CSolver ****solver_container;                 /*!< \brief Container vector with all the solutions. */
  CNumerics *****numerics_container;            /*!< \brief Description of the numerical method (the way in which the equations are solved). */
  CConfig **config_container;                   /*!< \brief Definition of the particular problem. */
  CSurfaceMovement **surface_movement;          /*!< \brief Surface movement classes of the problem. */
  CVolumetricMovement **grid_movement;          /*!< \brief Volume grid movement classes of the problem. */
  CFreeFormDefBox*** FFDBox;                    /*!< \brief FFD FFDBoxes of the problem. */
  CInterpolator ***interpolator_container;      /*!< \brief Definition of the interpolation method between non-matching discretizations of the interface. */
  CTransfer ***transfer_container;              /*!< \brief Definition of the transfer of information and the physics involved in the interface. */
  //Those are used to store the VarCoord of each node during FSI communications
  su2double APIVarCoord[3];
  su2double APINodalForce[3];
  su2double APINodalForceDensity[3];

public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] confFile - Configuration file name.
   	 * \param[in] val_nZone - Total number of zones.
	 * \param[in] val_nDim - Number of dimensions.
	 */
  CDriver(char* confFile,
          unsigned short val_nZone,
          unsigned short val_nDim);
	
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CDriver(void);

	/*!
	 * \brief A virtual member.
	 */  
  virtual void Run() { };

    /*!
     * \brief Construction of the edge-based data structure and the multigrid structure.
     */
  void Geometrical_Preprocessing();

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
  void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Postprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Preprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Postprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config);

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
  void Numerics_Preprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config);


  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Postprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config);


  /*!
   * \brief Deallocation routine
   */
  void Postprocessing();

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
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Displacements(unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter) {};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Tractions(unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter) {};

  /*!
   * \brief A virtual member.
   */
  virtual void Update() {};

  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  void StartSolver();

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
  bool Monitor(unsigned long ExtIter);

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long ExtIter);

  /*!
   * \brief Perform a dynamic mesh deformation, included grid velocity computation and update of the multigrid structure.
   */
  virtual void DynamicMeshUpdate(unsigned long ExtIter) { };

  /*!
   * \brief Perform a static mesh deformation, without considering grid velocity.
   */
  virtual void StaticMeshUpdate() { };

  /*!
   * \brief Perform a mesh deformation as initial condition.
   */
  virtual void SetInitialMesh() { };

  /*--- External communication layer ---*/
  
  su2double Get_Drag();
  su2double Get_Lift();
  su2double Get_Mz();
  unsigned short GetMovingMarker();
  unsigned long GetNumberVertices(unsigned short iMarker);
  unsigned long GetVertexGlobalIndex(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexCoordX(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexCoordY(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexCoordZ(unsigned short iMarker, unsigned short iVertex);
  bool ComputeVertexForces(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexForceX(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexForceY(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexForceZ(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexForceDensityX(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexForceDensityY(unsigned short iMarker, unsigned short iVertex);
  su2double GetVertexForceDensityZ(unsigned short iMarker, unsigned short iVertex);
  void SetVertexCoordX(unsigned short iMarker, unsigned short iVertex, su2double newPosX);
  void SetVertexCoordY(unsigned short iMarker, unsigned short iVertex, su2double newPosY);
  void SetVertexCoordZ(unsigned short iMarker, unsigned short iVertex, su2double newPosZ);
  su2double SetVertexVarCoord(unsigned short iMarker, unsigned short iVertex);

};
/*!
 * \class CSingleZoneDriver
 * \brief Class for driving an iteration of the physics within a single zone.
 * \author T. Economon
 * \version 4.3.0 "Cardinal"
 */
class CSingleZoneDriver : public CDriver {
public:
	
	/*! 
	 * \brief Constructor of the class.
	 * \param[in] confFile - Configuration file name.
	 * \param[in] val_nZone - Total number of zones.
	 * \param[in] val_nDim - Number of dimensions.
	 */
  CSingleZoneDriver(char* confFile,
                    unsigned short val_nZone,
                    unsigned short val_nDim);
	
	/*!
	 * \brief Destructor of the class.
	 */
	~CSingleZoneDriver(void);
	
	/*! 
	 * \brief Run a single iteration of the physics within a single zone.
	 */
  
  void Run();

    /*!
     * \brief Update the dual-time solution for a single zone.
     */
  void Update();

    /*!
     * \brief Reset the convergence flag (set to false) of the single zone solver.
     */
  void ResetConvergence();

    /*!
     * \brief Perform a dynamic mesh deformation, included grid velocity computation and the update of the multigrid structure (single zone).
     */
  void DynamicMeshUpdate(unsigned long ExtIter);

    /*!
     * \brief Perform a static mesh deformation, without considering grid velocity (single zone).
     */
  void StaticMeshUpdate();

    /*!
     * \brief Perform a mesh deformation as initial condition (single zone).
     */
  void SetInitialMesh();
};


/*!
 * \class CMultiZoneDriver
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon
 * \version 4.3.0 "Cardinal"
 */
class CMultiZoneDriver : public CDriver {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Number of dimensions.
   */
  CMultiZoneDriver(char* confFile,
                   unsigned short val_nZone,
                   unsigned short val_nDim);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMultiZoneDriver(void);
  
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
};


/*!
 * \class CHBDriver
 * \brief Class for driving an iteration of Harmonic Balance (HB) method problem using multiple time zones.
 * \author T. Economon
 * \version 4.3.0 "Cardinal"
 */
class CHBDriver : public CDriver {

private:

	su2double **D; /*!< \brief Harmonic Balance operator. */

public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] confFile - Configuration file name.
	 * \param[in] val_nZone - Total number of zones.
	 * \param[in] val_nDim - Number of dimensions.
	 */
	CHBDriver(char* confFile,
			unsigned short val_nZone,
			unsigned short val_nDim);

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
 * \version 4.3.0 "Cardinal"
 */
class CFSIDriver : public CDriver {
public:

	/*!
	 * \brief Constructor of the class.
	 * \param[in] confFile - Configuration file name.
	 * \param[in] val_nZone - Total number of zones.
	 */
	CFSIDriver(char* confFile,
			unsigned short val_nZone,
			unsigned short val_nDim);

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
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  void Relaxation_Displacements(unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter);

  /*!
   * \brief Apply a relaxation method into the computed tractions.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  void Relaxation_Tractions(unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter);

  /*!
   * \brief Enforce the coupling condition at the end of the time step
   * \param[in] zoneFlow - zone of the flow equations.
   * \param[in] zoneStruct - zone of the structural equations.
   */
  void Update(unsigned short zoneFlow, unsigned short zoneStruct);
  using CDriver::Update;
  
};

/*!
 * \class CFSIDriver
 * \brief Overload: class for driving a steady-state BGS iteration for a fluid-structure interaction problem in multiple zones.
 * \author R. Sanchez.
 * \version 4.2.0 "Cardinal"
 */
class CFSIStatDriver : public CFSIDriver {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of zones.
   */
  CFSIStatDriver(CIteration **iteration_container,
             CSolver ****solver_container,
             CGeometry ***geometry_container,
             CIntegration ***integration_container,
             CNumerics *****numerics_container,
             CInterpolator ***interpolator_container,
             CTransfer ***transfer_container,
             CConfig **config,
             unsigned short val_nZone,
             unsigned short val_nDim);

  /*!
   * \brief Destructor of the class.
   */
  ~CFSIStatDriver(void);

  /*!
   * \brief Run a Block Gauss-Seidel iteration of the FSI problem.
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

  void Run(CIteration **iteration_container,
           COutput *output,
           CIntegration ***integration_container,
           CGeometry ***geometry_container,
           CSolver ****solver_container,
           CNumerics *****numerics_container,
           CConfig **config_container,
           CSurfaceMovement **surface_movement,
           CVolumetricMovement **grid_movement,
           CFreeFormDefBox*** FFDBox,
           CInterpolator ***interpolator_container,
           CTransfer ***transfer_container);

};

/*!
 * \class CDiscAdjBlockFSIDriver
 * \brief Class for driving a BLOCK discrete adjoint FSI iteration.
 * \author R. Sanchez.
 * \version 4.2.0 "Cardinal"
 */
class CDiscAdjFSIStatDriver : public CFSIStatDriver {

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

  /*--- The first three are maintained so the values are the same as in the single-physics zone ---*/
  enum RECORDING{
    NONE = 0,               /*!< \brief Indicates that nothing is recorded. */
    FLOW_VARIABLES = 1,     /*!< \brief Indicates that the current recording can
                                        be used to compute the gradients of the flow problem with
                                        respect to the conservative flow variables. */
    GEOMETRY_VARIABLES = 2, /*!< \brief Indicates that the current recording can
                                        be used to compute the gradients of the flow problem with respect
                                        to the mesh geometry variables. */
    FEM_VARIABLES = 3,      /*!< \brief Indicates that the current recording can
                                        be used to compute the gradients of the structural problem with respect
                                        to the structural displacements. */
    ALL_VARIABLES = 4,      /*!< \brief All variables (monolithic solution) */
    FLOW_CROSS_TERM = 5,    /*!< \brief Indicates that the current recording can
                                        be used to compute the gradients of the structural problem
                                        with respect to the flow variables. */
    FEM_CROSS_TERM_GEOMETRY = 6,      /*!< \brief Indicates that the current recording can
                                        be used to compute the gradients of the mesh problem
                                        with respect to the structural displacements. */
    GEOMETRY_CROSS_TERM = 7,   /*!< \brief Indicates that the current recording can
                                        be used to compute the gradients of the structural problem
                                        with respect to the geometry variables. */
    FEM_CROSS_TERM_FLOW = 8   /*!< \brief Indicates that the current recording can
                                        be used to compute the gradients of the structural problem
                                        with respect to the flow variables. */
  };

  enum OF_KIND{
    NO_OBJECTIVE_FUNCTION = 0,               /*!< \brief Indicates that there is no objective function. */
    FLOW_OBJECTIVE_FUNCTION = 1,            /*!< \brief Indicates that the objective function is only flow-dependent. */
    FEM_OBJECTIVE_FUNCTION = 2              /*!< \brief Indicates that the objective function is only structural-dependent. */
  };

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_nZone - Total number of zones.
   */
  CDiscAdjFSIStatDriver(CIteration **iteration_container,
                            CSolver ****solver_container,
                            CGeometry ***geometry_container,
                            CIntegration ***integration_container,
                            CNumerics *****numerics_container,
                            CInterpolator ***interpolator_container,
                            CTransfer ***transfer_container,
                            CConfig **config,
                            unsigned short val_nZone,
                            unsigned short val_nDim);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFSIStatDriver(void);

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

  void Run(CIteration **iteration_container,
           COutput *output,
           CIntegration ***integration_container,
           CGeometry ***geometry_container,
           CSolver ****solver_container,
           CNumerics *****numerics_container,
           CConfig **config_container,
           CSurfaceMovement **surface_movement,
           CVolumetricMovement **grid_movement,
           CFreeFormDefBox*** FFDBox,
           CInterpolator ***interpolator_container,
           CTransfer ***transfer_container);

  /*!
   * \brief Iterate the direct solver for recording.
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

  void Iterate_Direct(CIteration **iteration_container, COutput *output, CIntegration ***integration_container,
                                CGeometry ***geometry_container, CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                                CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox, CInterpolator ***interpolator_container,
                                CTransfer ***transfer_container, unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT, unsigned short kind_recording);

  /*!
   * \brief Run a direct flow iteration.
   * \param[in] fluidZone - zone of the fluid solver.
   * \param[in] structuralZone - zone of the structural solver.
   */
  void Fluid_Iteration_Direct(CIteration **iteration_container, CTransfer ***transfer_container, COutput *output,
      CIntegration ***integration_container, CGeometry ***geometry_container, CSolver ****solver_container,
      CNumerics *****numerics_container, CConfig **config_container, CInterpolator ***interpolator_container,
      CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
      unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT);

  /*!
   * \brief Run a direct structural iteration.
   * \param[in] fluidZone - zone of the fluid solver.
   * \param[in] structuralZone - zone of the structural solver.
   */
  void Structural_Iteration_Direct(CIteration **iteration_container, CTransfer ***transfer_container, COutput *output,
      CIntegration ***integration_container, CGeometry ***geometry_container, CSolver ****solver_container,
      CNumerics *****numerics_container, CConfig **config_container, CInterpolator ***interpolator_container,
      CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
      unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT);

  /*!
   * \brief Run a direct mesh deformation.
   * \param[in] fluidZone - zone of the fluid solver.
   * \param[in] structuralZone - zone of the structural solver.
   */
  void Mesh_Deformation_Direct(CIteration **iteration_container, CTransfer ***transfer_container, COutput *output,
      CIntegration ***integration_container, CGeometry ***geometry_container, CSolver ****solver_container,
      CNumerics *****numerics_container, CConfig **config_container, CInterpolator ***interpolator_container,
      CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
      unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT);

  /*!
   * \brief Set the recording for a Discrete Adjoint iteration for the FSI problem.
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

  void SetRecording(CIteration **iteration_container,
                    COutput *output,
                    CIntegration ***integration_container,
                    CGeometry ***geometry_container,
                    CSolver ****solver_container,
                    CNumerics *****numerics_container,
                    CConfig **config_container,
                    CSurfaceMovement **surface_movement,
                    CVolumetricMovement **grid_movement,
                    CFreeFormDefBox*** FFDBox,
                    CInterpolator ***interpolator_container,
                    CTransfer ***transfer_container,
                    unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Load the restarts for fluid, structure and mesh.
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
  void Preprocess(CIteration **iteration_container,
                    COutput *output,
                    CIntegration ***integration_container,
                    CGeometry ***geometry_container,
                    CSolver ****solver_container,
                    CNumerics *****numerics_container,
                    CConfig **config_container,
                    CSurfaceMovement **surface_movement,
                    CVolumetricMovement **grid_movement,
                    CFreeFormDefBox*** FFDBox,
                    CInterpolator ***interpolator_container,
                    CTransfer ***transfer_container,
                    unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Iterate a certain block for adjoint FSI - may be the whole set of variables or independent and subiterate
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
  void Iterate_Block(CIteration **iteration_container,
                       COutput *output,
                       CIntegration ***integration_container,
                       CGeometry ***geometry_container,
                       CSolver ****solver_container,
                       CNumerics *****numerics_container,
                       CConfig **config_container,
                       CSurfaceMovement **surface_movement,
                       CVolumetricMovement **grid_movement,
                       CFreeFormDefBox*** FFDBox,
                       CInterpolator ***interpolator_container,
                       CTransfer ***transfer_container,
                       unsigned short ZONE_FLOW,
                       unsigned short ZONE_STRUCT,
                       unsigned short kind_recording);

  /*!
   * \brief Iterate a block for adjoint FSI with flow Objective Function
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
  void Iterate_Block_FlowOF(CIteration **iteration_container,
                               COutput *output,
                               CIntegration ***integration_container,
                               CGeometry ***geometry_container,
                               CSolver ****solver_container,
                               CNumerics *****numerics_container,
                               CConfig **config_container,
                               CSurfaceMovement **surface_movement,
                               CVolumetricMovement **grid_movement,
                               CFreeFormDefBox*** FFDBox,
                               CInterpolator ***interpolator_container,
                               CTransfer ***transfer_container,
                               unsigned short ZONE_FLOW,
                               unsigned short ZONE_STRUCT,
                               unsigned short kind_recording);

  /*!
   * \brief Iterate a block for adjoint FSI with structural Objective Function
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
  void Iterate_Block_StructuralOF(CIteration **iteration_container,
                                      COutput *output,
                                      CIntegration ***integration_container,
                                      CGeometry ***geometry_container,
                                      CSolver ****solver_container,
                                      CNumerics *****numerics_container,
                                      CConfig **config_container,
                                      CSurfaceMovement **surface_movement,
                                      CVolumetricMovement **grid_movement,
                                      CFreeFormDefBox*** FFDBox,
                                      CInterpolator ***interpolator_container,
                                      CTransfer ***transfer_container,
                                      unsigned short ZONE_FLOW,
                                      unsigned short ZONE_STRUCT,
                                      unsigned short kind_recording);

  /*!
   * \brief Initialize the adjoint - set the objective funcition and the output of the adjoint iteration
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   */
  void InitializeAdjoint(CIteration **iteration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CConfig **config_container,
                            unsigned short ZONE_FLOW,
                            unsigned short ZONE_STRUCT,
                            unsigned short kind_recording);

  /*!
   * \brief Extract the adjoint solution variables
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   */
  void ExtractAdjoint(CIteration **iteration_container,
                         CGeometry ***geometry_container,
                         CSolver ****solver_container,
                         CConfig **config_container,
                         unsigned short ZONE_FLOW,
                         unsigned short ZONE_STRUCT,
                         unsigned short kind_recording);


  /*!
   * \brief Check the convergence of the problem
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   */
  bool CheckConvergence(CIntegration ***integration_container,
                          CGeometry ***geometry_container,
                          CSolver ****solver_container,
                          CConfig **config_container,
                          unsigned long IntIter,
                          unsigned short ZONE_FLOW,
                          unsigned short ZONE_STRUCT,
                          unsigned short kind_recording);

  /*!
   * \brief Check the convergence of BGS subiteration process
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] numerics_container -  Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   */
  bool BGSConvergence(CIntegration ***integration_container,
                         CGeometry ***geometry_container,
                         CSolver ****solver_container,
                         CNumerics *****numerics_container,
                         CConfig **config_container,
                         unsigned long IntIter,
                         unsigned short ZONE_FLOW,
                         unsigned short ZONE_STRUCT);


  /*!
   * \brief Output the convergence history
   * \param[in] iteration_container - Container vector with all the iteration methods.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   */
  void ConvergenceHistory(CIntegration ***integration_container,
                             CGeometry ***geometry_container,
                             CSolver ****solver_container,
                             CConfig **config_container,
                             COutput *output,
                             unsigned long IntIter,
                             unsigned long nIntIter,
                             unsigned short ZONE_FLOW,
                             unsigned short ZONE_STRUCT,
                             unsigned short kind_recording);

  /*!
   * \brief Load the restarts for fluid, structure and mesh.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problem.
   * \param[in] ZONE_FLOW - Zone of the flow solver
   * \param[in] ZONE_STRUCT - Zone of the structural solver
   * \param[in] kind_recording - Kind of recording we are doing
   */
  void PrintDirect_Residuals(CSolver ****solver_container,
                                CConfig **config_container,
                                unsigned short ZONE_FLOW,
                                unsigned short ZONE_STRUCT,
                                unsigned short kind_recording);

  /*!
   * \brief Restart the variables to the converged solution.
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
  void PrepareRecording(CIteration **iteration_container,
                    COutput *output,
                    CIntegration ***integration_container,
                    CGeometry ***geometry_container,
                    CSolver ****solver_container,
                    CNumerics *****numerics_container,
                    CConfig **config_container,
                    CSurfaceMovement **surface_movement,
                    CVolumetricMovement **grid_movement,
                    CFreeFormDefBox*** FFDBox,
                    CInterpolator ***interpolator_container,
                    CTransfer ***transfer_container,
                    unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Register the input variables for adjoint FSI problems: flow conservative, fluid mesh position and structural displacements.
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
  void RegisterInput(CIteration **iteration_container,
                    COutput *output,
                    CIntegration ***integration_container,
                    CGeometry ***geometry_container,
                    CSolver ****solver_container,
                    CNumerics *****numerics_container,
                    CConfig **config_container,
                    CSurfaceMovement **surface_movement,
                    CVolumetricMovement **grid_movement,
                    CFreeFormDefBox*** FFDBox,
                    CInterpolator ***interpolator_container,
                    CTransfer ***transfer_container,
                    unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Register the input variables for adjoint FSI problems: flow conservative, fluid mesh position and structural displacements.
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
  void SetDependencies(CIteration **iteration_container,
                    COutput *output,
                    CIntegration ***integration_container,
                    CGeometry ***geometry_container,
                    CSolver ****solver_container,
                    CNumerics *****numerics_container,
                    CConfig **config_container,
                    CSurfaceMovement **surface_movement,
                    CVolumetricMovement **grid_movement,
                    CFreeFormDefBox*** FFDBox,
                    CInterpolator ***interpolator_container,
                    CTransfer ***transfer_container,
                    unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);

  /*!
   * \brief Restart the output variables for adjoint FSI problems: flow conservative, fluid mesh position and structural displacements.
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
  void RegisterOutput(CIteration **iteration_container,
                    COutput *output,
                    CIntegration ***integration_container,
                    CGeometry ***geometry_container,
                    CSolver ****solver_container,
                    CNumerics *****numerics_container,
                    CConfig **config_container,
                    CSurfaceMovement **surface_movement,
                    CVolumetricMovement **grid_movement,
                    CFreeFormDefBox*** FFDBox,
                    CInterpolator ***interpolator_container,
                    CTransfer ***transfer_container,
                    unsigned short ZONE_FLOW,
                    unsigned short ZONE_STRUCT,
                    unsigned short kind_recording);


};

