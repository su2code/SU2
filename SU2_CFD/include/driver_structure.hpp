/*!
 * \file driver_structure.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
 * \version 5.0.0 "Raven"
 */
class CDriver {
protected:
  char* config_file_name;                       /*!< \brief Configuration file name of the problem.*/
  char runtime_file_name[MAX_STRING_SIZE];
  su2double StartTime,                          /*!< \brief Start point of the timer for performance benchmarking.*/
            StopTime,                           /*!< \brief Stop point of the timer for performance benchmarking.*/
            UsedTime;                           /*!< \brief Elapsed time between Start and Stop point of the timer.*/
  unsigned long ExtIter;                        /*!< \brief External iteration.*/
  ofstream *ConvHist_file;                       /*!< \brief Convergence history file.*/
  unsigned short iMesh,                         /*!< \brief Iterator on mesh levels.*/
                iZone,                          /*!< \brief Iterator on zones.*/
                nZone,                          /*!< \brief Total number of zones in the problem. */
                nDim;                           /*!< \brief Number of dimensions.*/
  bool StopCalc,                                /*!< \brief Stop computation flag.*/
       mixingplane,                             /*!< \brief mixing-plane simulation flag.*/
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
  su2double APIVarCoord[3];                     /*!< \brief This is used to store the VarCoord of each node. */
  su2double APINodalForce[3];                   /*!< \brief This is used to store the force at each node. */
  su2double APINodalForceDensity[3];            /*!< \brief This is used to store the force density at each node. */

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
   * \brief Restart of the solvers from the restart files.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Restart(CSolver ***solver_container, CGeometry **geometry, CConfig *config, bool update_geo);

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
   * \brief Perform a static mesh deformation, without considering grid velocity.
   */
  virtual void StaticMeshUpdate() { };

  /*!
   * \brief Perform a mesh deformation as initial condition.
   */
  virtual void SetInitialMesh() { };

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

};

/*!
 * \class CGeneralDriver
 * \brief Class for driving a structural iteration of the physics within multiple zones.
 * \author T. Economon
 * \version 5.0.0 "Raven"
 */
class CGeneralDriver : public CDriver {
public:
  
  /*! 
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Number of dimensions.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CGeneralDriver(char* confFile,
                 unsigned short val_nZone,
                 unsigned short val_nDim,
                 SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralDriver(void);

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
 * \class CFluidDriver
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon, G. Gori
 * \version 5.0.0 "Raven"
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
   * \brief Transfer data among different zones (multiple zone).
   */
  void Transfer_Data(unsigned short donorZone, unsigned short targetZone);
};


/*!
 * \class CTurbomachineryDriver
 * \brief Class for driving an iteration for turbomachinery flow analysis.
 * \author S. Vitale
 * \version 5.0.0 "Raven"
 */
class CTurbomachineryDriver : public CFluidDriver {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Number of dimensions.
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
 * \class CDiscAdjMultiZoneDriver
 * \brief Class for driving an iteration of the discrete adjoint within multiple zones.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 */
class CDiscAdjFluidDriver : public CFluidDriver {

protected:
  unsigned short RecordingState; /*!< \brief The kind of recording the tape currently holds.*/
  su2double ObjFunc;             /*!< \brief The value of the objective function.*/
  CIteration** direct_iteration; /*!< \brief A pointer to the direct iteration.*/

public:

  /*!
    * \brief Constructor of the class.
    * \param[in] confFile - Configuration file name.
    * \param[in] val_nZone - Total number of zones.
    * \param[in] val_nDim - Number of dimensions.
    */
  CDiscAdjFluidDriver(char* confFile,
                   unsigned short val_nZone,
                   unsigned short val_nDim,
                   SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFluidDriver(void);

  /*!
   * \brief Run a single iteration of the discrete adjoint solver within multiple zones.
   */

  void Run();

  /*!
   * \brief Record one iteration of a flow iteration in within multiple zones.
   * \param[in] kind_recording - Type of recording (either FLOW_CONS_VARS, MESH_COORDS, COMBINED or NONE)
   */

  void SetRecording(unsigned short kind_recording);

  /*!
   * \brief Run one iteration of the solver. It is virtual because it depends on the kind of physics.
   */
  virtual void DirectRun();

  /*!
   * \brief Set the objective function. It is virtual because it depends on the kind of physics.
   */
  virtual void SetObjFunction();

  /*!
   * \brief Initialize the adjoint value of the objective function.
   */
  void SetAdj_ObjFunction();
};

/*!
 * \class CDiscAdjTurbomachineryDriver
 * \brief Class for driving an iteration of the discrete adjoint within multiple zones.
 * \author S. Vitale, T. Albring
 * \version 5.0.0 "Raven"
 */
class CDiscAdjTurbomachineryDriver : public  CDiscAdjFluidDriver {

public:

	 /*!
	   * \brief Constructor of the class.
	   * \param[in] confFile - Configuration file name.
	   * \param[in] val_nZone - Total number of zones.
	   * \param[in] val_nDim - Number of dimensions.
	   */
  CDiscAdjTurbomachineryDriver(char* confFile,
                   unsigned short val_nZone,
                   unsigned short val_nDim, SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjTurbomachineryDriver(void);

  /*!
   * \brief Run a single iteration of the direct solver.
   */
  void DirectRun();

  /*!
   * \brief Set Obj.Function for turbomachinery design.
   */
  void SetObjFunction();

  /*!
   * \brief Set Mixing Plane interface within multiple zones.
   */
  void SetMixingPlane(unsigned short iZone);

  /*!
   * \brief Set Mixing Plane interface within multiple zones.
   */
  void SetTurboPerformance(unsigned short targetZone);


};
/*!
 * \class CHBDriver
 * \brief Class for driving an iteration of Harmonic Balance (HB) method problem using multiple time zones.
 * \author T. Economon
 * \version 5.0.0 "Raven"
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
 * \version 5.0.0 "Raven"
 */
class CFSIDriver : public CDriver {
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
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   */
  CFSIStatDriver(char* confFile,
             unsigned short val_nZone,
             unsigned short val_nDim,
             SU2_Comm MPICommunicator);

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

  void Run();

};

/*!
 * \class CDiscAdjFSIStatDriver
 * \brief Class for driving a discrete adjoint FSI iteration.
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
  CDiscAdjFSIStatDriver(char* confFile,
                           unsigned short val_nZone,
                           unsigned short val_nDim,
                           SU2_Comm MPICommunicator);

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
   * \brief Iterate a block for adjoint FSI with flow Objective Function
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void Iterate_Block_FlowOF(unsigned short ZONE_FLOW,
                               unsigned short ZONE_STRUCT,
                               unsigned short kind_recording);

  /*!
   * \brief Iterate a block for adjoint FSI with structural Objective Function
   * \param[in] ZONE_FLOW - zone of the fluid solver.
   * \param[in] ZONE_STRUCT - zone of the structural solver.
   * \param[in] kind_recording - kind of recording (flow, structure, mesh, cross terms)
   */
  void Iterate_Block_StructuralOF(unsigned short ZONE_FLOW,
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


};

