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
  virtual void Run(){ };

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
  virtual void Predict_Displacements(unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the tractions will be predicted.
   * \param[in] targetZone - zone which receives the predicted traction.
   */
  virtual void Predict_Tractions(unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone in which the displacements will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void Transfer_Displacements(unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - zone from which the tractions will be transferred.
   * \param[in] targetZone - zone which receives the tractions transferred.
   */
  virtual void Transfer_Tractions(unsigned short donorZone, unsigned short targetZone){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Displacements(unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter){};

  /*!
   * \brief A virtual member.
   * \param[in] donorZone - origin of the information.
   * \param[in] targetZone - destination of the information.
   * \param[in] iFSIIter - Fluid-Structure Interaction subiteration.
   */
  virtual void Relaxation_Tractions(unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter){};

  /*!
   * \brief A virtual member.
   */
  virtual void Update(){};

  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  void StartSolver();

  /*!
   * \brief A virtual member.
   */
  virtual void ResetConvergence(){ };

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
  virtual void DynamicMeshUpdate(unsigned long ExtIter){ };

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
  unsigned int GetVertexGlobalIndex(unsigned short iMarker, unsigned short iVertex);
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
 * \class CSpectralDriver
 * \brief Class for driving an iteration of a spectral method problem using multiple zones.
 * \author T. Economon
 * \version 4.3.0 "Cardinal"
 */
class CSpectralDriver : public CDriver {
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Number of dimensions.
   */
  CSpectralDriver(char* confFile,
                  unsigned short val_nZone,
                  unsigned short val_nDim);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSpectralDriver(void);
  
  /*!
   * \brief Run a single iteration of a spectral method problem.
   */
  
  void Run();
  
  /*!
   * \brief Computation and storage of the time spectral source terms.
   * \author T. Economon, K. Naik
   * \param[in] iZone - Current zone number.
   */
  void SetTimeSpectral(unsigned short iZone);
  
  /*!
   * \brief Computation of the Time-Spectral operator matrix.
   * \author K. Naik
   * \param[in] D - su2double pointer to the operator matrix.
   */
  void ComputeTimeSpectral_Operator(su2double **D, su2double period);
  
  /*!
   * \brief Computation and storage of the time-spectral mesh velocities.
   * \author K. Naik, T. Economon
   */
  void SetTimeSpectral_Velocities();

  /*!
   * \brief Update the solution for the spectral method.
   */
  void Update();

  /*!
   * \brief Reset the convergence flag (set to false) of the solver for the spectral method.
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
