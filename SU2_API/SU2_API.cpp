#include "../SU2_CFD/include/SU2_CFD.hpp"
#include "SU2_API.h"
#include <stdlib.h>

using namespace std;

SU2Solver::SU2Solver(string str):confFile(str){

  output                       = NULL;
  config = NULL;
  integration_container = NULL;
  geometry_container       = NULL;
  solver_container          = NULL;
  numerics_container     = NULL;
  config_container            = NULL;
  surface_movement   = NULL;
  grid_movement   = NULL;
  FFDBox             = NULL;
  driver            =NULL;
  iteration_container   =NULL;

}

SU2Solver::~SU2Solver(){}

void SU2Solver::initialize(bool FSIComp){

  unsigned short iMesh, iZone, iSol;
  unsigned long ExtIter;

  int rank = MASTER_NODE;
  int size = SINGLE_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if(rank == MASTER_NODE){
    if(FSIComp)cout << endl << "***************************** Setting SU2 for FSI simulation *****************************" << endl;
    else cout << endl <<"***************************** Setting SU2 for CFD simulation *****************************" << endl;
  }

  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/

  char config_file_name[200];
  strcpy(config_file_name, confFile.c_str());

  /*--- Read the name and format of the input mesh file ---*/
  config = new CConfig(config_file_name, SU2_CFD);

  /*--- Get the number of zones and dimensions from the numerical grid
   (required for variables allocation) ---*/

  nZone = GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
  nDim  = GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());

  /*--- Definition and of the containers for all possible zones. ---*/

  iteration_container   = new CIteration*[nZone];
  solver_container      = new CSolver***[nZone];
  integration_container = new CIntegration**[nZone];
  numerics_container    = new CNumerics****[nZone];
  config_container      = new CConfig*[nZone];
  geometry_container    = new CGeometry **[nZone];
  surface_movement      = new CSurfaceMovement *[nZone];
  grid_movement         = new CVolumetricMovement *[nZone];
  FFDBox                = new CFreeFormDefBox**[nZone];

  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]       = NULL;
    integration_container[iZone]  = NULL;
    numerics_container[iZone]     = NULL;
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
    surface_movement[iZone]       = NULL;
    grid_movement[iZone]          = NULL;
    FFDBox[iZone]                 = NULL;
  }

  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {

      /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    config_container[iZone] = new CConfig(config_file_name, SU2_CFD, iZone, nZone, nDim, VERB_HIGH);

    /*--- Definition of the geometry class to store the primal grid in the
     partitioning process. ---*/

    CGeometry *geometry_aux = NULL;

    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);

    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);

    /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

    geometry_container[iZone] = new CGeometry *[config_container[iZone]->GetnMGLevels()+1];
    geometry_container[iZone][MESH_0] = new CPhysicalGeometry(geometry_aux, config_container[iZone], 1);

    /*--- Deallocate the memory of geometry_aux ---*/

    delete geometry_aux;

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone][MESH_0]->SetSendReceive(config_container[iZone]);

    /*--- Add the Send/Receive boundaries ---*/

    geometry_container[iZone][MESH_0]->SetBoundaries(config_container[iZone]);

  }

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;

  /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
   based data structure is constructed, i.e. node and cell neighbors are
   identified and linked, face areas and volumes of the dual mesh cells are
   computed, and the multigrid levels are created using an agglomeration procedure. ---*/

  Geometrical_Preprocessing(geometry_container, config_container, nZone);

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Computation of wall distances for turbulence modeling ---*/

    if (rank == MASTER_NODE)
      cout << "Computing wall distances." << endl;

    if ((config_container[iZone]->GetKind_Solver() == RANS) ||
        (config_container[iZone]->GetKind_Solver() == ADJ_RANS) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS))
      geometry_container[iZone][MESH_0]->ComputeWall_Distance(config_container[iZone]);

    /*--- Computation of positive surface area in the z-plane which is used for
     the calculation of force coefficient (non-dimensionalization). ---*/

    geometry_container[iZone][MESH_0]->SetPositive_ZArea(config_container[iZone]);

    /*--- Set the near-field, interface and actuator disk boundary conditions, if necessary. ---*/

    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
      geometry_container[iZone][iMesh]->MatchNearField(config_container[iZone]);
      geometry_container[iZone][iMesh]->MatchInterface(config_container[iZone]);
      geometry_container[iZone][iMesh]->MatchActuator_Disk(config_container[iZone]);
    }

  }

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Driver Preprocessing --------------------------" << endl;

  /*--- First, given the basic information about the number of zones and the
   solver types from the config, instantiate the appropriate driver for the problem. ---*/

  Driver_Preprocessing(&driver, iteration_container, solver_container,
                       geometry_container, integration_container, numerics_container, config_container, nZone);


  /*--- Instantiate the geometry movement classes for the solution of unsteady
   flows on dynamic meshes, including rigid mesh transformations, dynamically
   deforming meshes, and time-spectral preprocessing. ---*/


  for (iZone=0; iZone < nZone; iZone++){

    if (FSIComp || config_container[iZone]->GetGrid_Movement() ||
        (config_container[iZone]->GetDirectDiff() == D_DESIGN)) {
      if (rank == MASTER_NODE)
        cout << "Setting dynamic mesh structure." << endl;
      grid_movement[iZone] = new CVolumetricMovement(geometry_container[iZone][MESH_0], config_container[iZone]);
      FFDBox[iZone] = new CFreeFormDefBox*[MAX_NUMBER_FFD];
      surface_movement[iZone] = new CSurfaceMovement();
      surface_movement[iZone]->CopyBoundary(geometry_container[iZone][MESH_0], config_container[iZone]);
      if (config_container[iZone]->GetUnsteady_Simulation() == TIME_SPECTRAL)
        SetGrid_Movement(geometry_container[iZone], surface_movement[iZone], grid_movement[iZone],
                         FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, 0, 0);
    }

    if (config_container[iZone]->GetDirectDiff() == D_DESIGN){
      if (rank == MASTER_NODE)
        cout << "Setting surface/volume derivatives." << endl;

      /*--- Set the surface derivatives, i.e. the derivative of the surface mesh nodes with respect to the design variables ---*/

      surface_movement[iZone]->SetSurface_Derivative(geometry_container[iZone][MESH_0],config_container[iZone]);

      /*--- Call the volume deformation routine with derivative mode enabled.
       This computes the derivative of the volume mesh with respect to the surface nodes ---*/

      grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone][MESH_0],config_container[iZone], true, true);

      /*--- Update the multi-grid structure to propagate the derivative information to the coarser levels ---*/

      geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);

      /*--- Set the derivative of the wall-distance with respect to the surface nodes ---*/

      if ( (config_container[iZone]->GetKind_Solver() == RANS) ||
          (config_container[iZone]->GetKind_Solver() == ADJ_RANS) ||
          (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS))
        geometry_container[iZone][MESH_0]->ComputeWall_Distance(config_container[iZone]);
    }
  }

  /*--- Coupling between zones (limited to two zones at the moment) ---*/

  if (nZone == 2) {
    if (rank == MASTER_NODE)
      cout << endl <<"--------------------- Setting Coupling Between Zones --------------------" << endl;
    geometry_container[ZONE_0][MESH_0]->MatchZone(config_container[ZONE_0], geometry_container[ZONE_1][MESH_0],
                                                  config_container[ZONE_1], ZONE_0, nZone);
    geometry_container[ZONE_1][MESH_0]->MatchZone(config_container[ZONE_1], geometry_container[ZONE_0][MESH_0],
                                                  config_container[ZONE_0], ZONE_1, nZone);
  }


  if(FSIComp){
  /*--- Initialisation of the FSI interface (currently only for rigid body displacements) ---*/
  interfRigidDispArray = new su2double[6];
  interfRigidDispArray[0] = 0.0;
  interfRigidDispArray[1] = 0.0;
  interfRigidDispArray[2] = 0.0;
  interfRigidDispArray[3] = 0.0;
  interfRigidDispArray[4] = 0.0;
  interfRigidDispArray[5] = 0.0;

  attitude_i1 = new su2double[3];
  attitude_i1[0] = 0.0;
  attitude_i1[1] = 0.0;
  attitude_i1[2] = 0.0;

  }

  /*--- Definition of the output class (one for all zones). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/

  output = new COutput();

  /*--- Open the convergence history file ---*/

  if(rank == MASTER_NODE)
    output->SetConvHistory_Header(&ConvHist_file, config_container[ZONE_0]);

  /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
    ExtIter = config_container[ZONE_0]->GetUnst_RestartIter();

  /*--- Initiate value at each interface for the mixing plane ---*/
  if(config_container[ZONE_0]->GetBoolMixingPlane())
  	for (iZone = 0; iZone < nZone; iZone++)
  	  iteration_container[iZone]->Preprocess(output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone);

  ofstream historyFile_FSI;
  bool writeHistFSI = config_container[ZONE_0]->GetWrite_Conv_FSI();
  if (writeHistFSI){
	  char cstrFSI[200];
	  string filenameHistFSI = config_container[ZONE_0]->GetConv_FileName_FSI();
	  strcpy (cstrFSI, filenameHistFSI.data());
	  historyFile_FSI.open (cstrFSI);
	  historyFile_FSI << "Time,Iteration,Aitken,URes,logResidual,orderMagnResidual" << endl;
	  historyFile_FSI.close();
  }

  /*--------------------------------------------------------------------------------------------------------------------------------------*/

  if(rank == MASTER_NODE){
    if(FSIComp)cout << endl << "***************************** SU2 is set for FSI simulation *****************************" << endl;
    else cout << endl <<"***************************** SU2 is set for CFD simulation *****************************" << endl;
  }

#ifdef HAVE_MPI
/*--- Synchronization point after SU2 initialization ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void SU2Solver::exit(){

int rank = MASTER_NODE;
int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif


  unsigned short iMesh, iZone, iSol;

  if (rank == MASTER_NODE) {
    cout << endl;

  /*--- Print out the number of non-physical points and reconstructions ---*/
  if (config_container[ZONE_0]->GetNonphysical_Points() > 0)
    cout << "Warning: there are " << config_container[ZONE_0]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
  if (config_container[ZONE_0]->GetNonphysical_Reconstr() > 0)
    cout << "Warning: " << config_container[ZONE_0]->GetNonphysical_Reconstr() << " reconstructed states for upwinding are non-physical." << endl;

  /*--- Close the convergence history file. ---*/
  ConvHist_file.close();
  cout << "Fluid history file is closed." << endl;


  /*--- Exit the solver cleanly ---*/
  cout << endl << "***************************** Exit SU2 *****************************" << endl;


  }

  /* --- Free the memory ---*/  // Need to fix problems when deleting some container !!!
  for(iZone = 0; iZone < nZone; iZone ++){
    for(iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
      for(iSol = 0; iSol < MAX_SOLS; iSol++){
        //delete [] solver_container[iZone][iMesh][iSol]; no need !!
      }
      delete [] solver_container[iZone][iMesh];
    }
    //cout << config_container[iZone];
    //delete config_container[iZone];
    //delete geometry_container[iZone][MESH_0];
    delete [] geometry_container[iZone];
    delete [] solver_container[iZone];
    delete [] integration_container[iZone];
    delete [] numerics_container[iZone];
    //delete grid_movement[iZone];
    delete surface_movement[iZone];
    delete [] FFDBox[iZone];
  }

  delete [] iteration_container;
  delete [] config_container;
  delete [] geometry_container;
  delete [] solver_container;
  delete [] integration_container;
  delete [] numerics_container;
  delete [] surface_movement;
  delete [] grid_movement;
  delete [] FFDBox;

  //delete config;

  delete output;

  delete [] interfRigidDispArray;
  delete [] attitude_i1;

  if(solidInterface != NULL){
    for(int iVertex=0; iVertex<nSolidInterfaceVertex; iVertex++){
      delete [] solidInterface[iVertex];
    }

    delete [] solidInterface;
    solidInterface = NULL;
  }

  if(solidInterfaceBuffer != NULL){
    delete [] solidInterfaceBuffer;
    solidInterfaceBuffer = NULL;
  }

  if(fluidSurfaceloads != NULL){
    for(int iVertex=0; iVertex<nAllFluidInterfaceVertex; iVertex++){
      delete [] fluidSurfaceloads[iVertex];
    }
    delete [] fluidSurfaceloads;
    fluidSurfaceloads = NULL;
  }

  if(partFluidSurfaceLoads !=NULL){
    for(int iVertex=0; iVertex<nLocalFluidInterfaceVertex; iVertex++){
      delete [] partFluidSurfaceLoads[iVertex];
    }
    delete [] partFluidSurfaceLoads;
  }

  //delete [] attitude_n1;


#ifdef HAVE_MPI
/*--- Synchronization point after exiting the solver ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void SU2Solver::connectToSolidSolver(unsigned long nVertex){

  int rank = MASTER_NODE;
  int size = 1;

  //bool matchingMeshes(true);

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /*--- Import the number of solid interface vertices. ---*/
  nSolidInterfaceVertex = nVertex;
  if(rank == MASTER_NODE) cout << "Solid solver has communicated a number of nodes on the moving solid interface equal to : " << nSolidInterfaceVertex << "." << endl;

  /*if(matchingMeshes){
    if(nSolidInterfaceVertex == nInterfaceVertex){
      //Check if mesh are matching (use of tolerance)
    }
    else{
      //if(rank == MASTER_NODE) cout << "For matching meshes computation, the number of vertices on each interface must be equal ! Exiting..." << endl;
      //throw(-1);
    }
  }*/


}

void SU2Solver::setFSIInterface(){

  int rank = MASTER_NODE;
  int size = 1;

  unsigned short iMarker, jMarker, Moving, iVertex;
  string Marker_Tag, Moving_Tag;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /*--- Broadcast the number of solid interface vertices. ---*/
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&nSolidInterfaceVertex, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#endif // HAVE_MPI

  /*--- Initialize the solid interface container and the associated buffer ---*/
  solidInterfaceBuffer = NULL;
  solidInterface = NULL;

  solidInterface = new su2double*[nSolidInterfaceVertex];

  for(int iVertex=0; iVertex < nSolidInterfaceVertex; iVertex++){
    solidInterface[iVertex] = new su2double[4];
  }

  if(solidInterfaceBuffer == NULL) solidInterfaceBuffer = new su2double[nSolidInterfaceVertex*4];


  /*--- Get the number of vertices on the fluid interface on each process ---*/
  nAllFluidInterfaceVertex = 0;
  nLocalFluidInterfaceVertex = 0;

  for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
    Moving = config_container[ZONE_0]->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config_container[ZONE_0]->GetnMarker_Moving(); jMarker++) {
        Moving_Tag = config_container[ZONE_0]->GetMarker_Moving(jMarker);
        Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Moving_Tag) {
          nLocalFluidInterfaceVertex = geometry_container[ZONE_0][MESH_0]->nVertex[iMarker];
          cout << "Message from process " << rank << " : there are " << nLocalFluidInterfaceVertex << " nodes on the moving fluid interface." << endl;
        }
      }
    }
  }

  /*--- Initialize the fluid interface container ---*/
  fluidSurfaceloads = NULL;
  partFluidSurfaceLoads = NULL;

  partFluidSurfaceLoads = new su2double*[nLocalFluidInterfaceVertex];
  for(iVertex = 0; iVertex<nLocalFluidInterfaceVertex; iVertex++){
    partFluidSurfaceLoads[iVertex] = new su2double[4];
  }

  /*--- Compute the total number of vertices by MPI_Reduce, including the halo nodes ---*/
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&nLocalFluidInterfaceVertex, &nAllFluidInterfaceVertex, 1, MPI_INT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#endif // HAVE_MPI

if (rank==MASTER_NODE){
  cout << "There are " << nAllFluidInterfaceVertex << " nodes on the moving fluid interface (including halo nodes)" << endl;
  fluidSurfaceloads = new su2double*[nAllFluidInterfaceVertex];
   for(int iVertex=0; iVertex < nAllFluidInterfaceVertex; iVertex++){
     fluidSurfaceloads[iVertex] = new su2double[4];
   }
}


#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif // HAVE_MPI

}

double* SU2Solver::getInterRigidDispArray() const{ //Will disappear

  return interfRigidDispArray;
}

void SU2Solver::setTemporalIteration(unsigned long ExtIter){

  unsigned short iZone;

  for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetExtIter(ExtIter);
      //config_container[iZone]->UpdateCFL(ExtIter);
  }

  /*--- Read the target pressure ---*/
  if (config_container[ZONE_0]->GetInvDesign_Cp() == YES)
    output->SetCp_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);

  /*--- Read the target heat flux ---*/
  if (config_container[ZONE_0]->GetInvDesign_HeatFlux() == YES)
    output->SetHeat_InverseDesign(solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);

  /*--- Set the initial condition ---*/
  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
  }

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void SU2Solver::dualTimeInnerLoop(unsigned long ExtIter){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short iZone,iMesh;
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);

  /*--- Reset the convergence flags before fluid iteration ---*/
  for(iZone = 0; iZone < nZone; iZone++){
    integration_container[iZone][FLOW_SOL]->SetConvergence(false);
    if(config_container[iZone]->GetKind_Solver() == RANS)
      integration_container[iZone][TURB_SOL]->SetConvergence(false);
    if(config_container[iZone]->GetKind_Trans_Model() == LM)
      integration_container[iZone][TRANS_SOL]->SetConvergence(false);
  }

  /*--- Inner loop for dual-time stepping ---*/
  for(IntIter = 0; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {

    if(IntIter > 0){
    /*--- Write the convergence history (only screen output) ---*/
    output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);
    }

    /*--- Set the value of the internal iteration ---*/
    config_container[ZONE_0]->SetIntIter(IntIter);

    /*--- All zones must be advanced and coupled with each pseudo timestep. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {

      /*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/
      if (config_container[iZone]->GetKind_Solver() == EULER)
        config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
	  if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES)
        config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
      if (config_container[iZone]->GetKind_Solver() == RANS)
        config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);

      /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
      integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_FLOW_SYS, IntIter, iZone);

      /*--- Pseudo-timestepping the turbulence model ---*/
      if (config_container[iZone]->GetKind_Solver() == RANS) {
        /*--- Solve the turbulence model ---*/
        config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
        integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_TURB_SYS, IntIter, iZone);
        /*--- Solve transition model ---*/
        if (config_container[iZone]->GetKind_Trans_Model() == LM) {
		  config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
		  integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
        }

      }

      /*--- Call Dynamic mesh update if AEROELASTIC motion was specified ---*/
      if ((config_container[ZONE_0]->GetGrid_Movement()) && (config_container[ZONE_0]->GetAeroelastic_Simulation()))
        SetGrid_Movement(geometry_container[iZone], surface_movement[iZone], grid_movement[iZone], FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, IntIter, ExtIter);
    }//ENF OF IZONE LOOP

    if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()){
      break;
    }

  } //END OF INNER LOOP

  output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

  if(rank == MASTER_NODE){
    cout << endl << "Fluid solution is converged (dual-time inner loop) !" << endl;
  }

#ifdef HAVE_MPI
/*--- Synchronization point after the dual-time inner loop ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void SU2Solver::steadyFluidIteration(unsigned long ExtIter){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short iZone, iMesh;
  unsigned long IntIter = 0;

  /*--- Set the initial condition (may be ignore if steady computation...) ---*/
  for (iZone = 0; iZone < nZone; iZone++)
    solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Set the value of the internal iteration ---*/
	IntIter = ExtIter;
    config_container[ZONE_0]->SetIntIter(ExtIter);

    /*--- Update global parameters ---*/
	if (config_container[iZone]->GetKind_Solver() == EULER){
      config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
    }
    if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES){
      config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
    }
    if (config_container[iZone]->GetKind_Solver() == RANS){
      config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
    }

    /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
    integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_FLOW_SYS, IntIter, iZone);

    if (config_container[iZone]->GetKind_Solver() == RANS) {
      /*--- Solve the turbulence model ---*/
      config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
      integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_TURB_SYS, IntIter, iZone);

      /*--- Solve transition model ---*/
      if (config_container[iZone]->GetKind_Trans_Model() == LM) {
		config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
		integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
      }
    }

  }

#ifdef HAVE_MPI
/*--- Synchronization point after the solver definition subroutine ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void SU2Solver::resetConvergence(){

  unsigned short iZone;

  for(iZone = 0; iZone < nZone; iZone++){
    integration_container[iZone][FLOW_SOL]->SetConvergence(false);
    if(config_container[iZone]->GetKind_Solver() == RANS)
      integration_container[iZone][TURB_SOL]->SetConvergence(false);
    if(config_container[iZone]->GetKind_Trans_Model() == LM)
      integration_container[iZone][TRANS_SOL]->SetConvergence(false);
  }

}

void SU2Solver::updateDualTime(){

  unsigned short iZone, iMesh;

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Update dual time solver on all mesh levels ---*/
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
	  integration_container[iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][FLOW_SOL], config_container[iZone], iMesh);
	  integration_container[iZone][FLOW_SOL]->SetConvergence(false);
    }

    /*--- Update dual time solver for the turbulence model ---*/
    if (config_container[iZone]->GetKind_Solver() == RANS) {
      integration_container[iZone][TURB_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][TURB_SOL], config_container[iZone], MESH_0);
      integration_container[iZone][TURB_SOL]->SetConvergence(false);
    }

    /*--- Update dual time solver for the transition model ---*/
	if (config_container[iZone]->GetKind_Trans_Model() == LM) {
	  integration_container[iZone][TRANS_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][TRANS_SOL], config_container[iZone], MESH_0);
	  integration_container[iZone][TRANS_SOL]->SetConvergence(false);
	}
  }

  /*--- Reset the rigid body position, only relevant for steady computations using purely rigid body method ---*/
  attitude_i1[0] = 0.0;
  attitude_i1[1] = 0.0;
  attitude_i1[2] = 0.0;

  /*--- Update the center of reference for the next time step, only relevant for unsteady computations using purely rigid body method  ---*/
  Center_n[0] += interfRigidDispArray[0];
  Center_n[1] += interfRigidDispArray[1];
  Center_n[2] += interfRigidDispArray[2];

#ifdef HAVE_MPI
/*--- Synchronization point after updating the dual-time solution ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

bool SU2Solver::writeSolution(unsigned long ExtIter){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  bool StopCalc = false;
  su2double UsedTime(0.0);

  /*--- For specific applications, evaluate and plot the equivalent area. ---*/
  if (config_container[ZONE_0]->GetEquivArea() == YES) {
    output->SetEquivalentArea(solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);
  }

  /*--- Check if there is any change in the runtime parameters ---*/
  //CConfig *runtime = NULL;
  //strcpy(runtime_file_name, "runtime.dat");
  //runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);

  /*--- Update the convergence history file (serial and parallel computations). ---*/
  output->SetConvHistory_Body(&ConvHist_file, geometry_container, solver_container, config_container, integration_container, false, UsedTime, ZONE_0);

  /*--- Evaluate the new CFL number (adaptive). ---*/
  if (config_container[ZONE_0]->GetCFL_Adapt() == YES) {
    output->SetCFL_Number(solver_container, config_container, ZONE_0);
  }

  /*--- Check whether the current simulation has reached the specified
     convergence criteria, and set StopCalc to true, if so. ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
      case EULER: case NAVIER_STOKES: case RANS:
        StopCalc = integration_container[ZONE_0][FLOW_SOL]->GetConvergence(); break;
      case WAVE_EQUATION:
        StopCalc = integration_container[ZONE_0][WAVE_SOL]->GetConvergence(); break;
      case HEAT_EQUATION:
        StopCalc = integration_container[ZONE_0][HEAT_SOL]->GetConvergence(); break;
      case LINEAR_ELASTICITY:
        // This is a temporal fix, while we code the non-linear solver
        //	        StopCalc = integration_container[ZONE_0][FEA_SOL]->GetConvergence(); break;
        StopCalc = false; break;
      case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
        StopCalc = integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence(); break;
    }

    /*--- Solution output. Determine whether a solution needs to be written
     after the current iteration, and if so, execute the output file writing
     routines. ---*/

  if ((ExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
        ((ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (ExtIter != 0) &&
         !((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
           (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
        (StopCalc) ||
        (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
          (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
         ((ExtIter == 0) || (ExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {

          /*--- Low-fidelity simulations (using a coarser multigrid level
           approximation to the solution) require an interpolation back to the
           finest grid. ---*/

          if (config_container[ZONE_0]->GetLowFidelitySim()) {
            integration_container[ZONE_0][FLOW_SOL]->SetProlongated_Solution(RUNTIME_FLOW_SYS, solver_container[ZONE_0][MESH_0][FLOW_SOL], solver_container[ZONE_0][MESH_1][FLOW_SOL], geometry_container[ZONE_0][MESH_0], geometry_container[ZONE_0][MESH_1], config_container[ZONE_0]);
            integration_container[ZONE_0][FLOW_SOL]->Smooth_Solution(RUNTIME_FLOW_SYS, solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], 3, 1.25, config_container[ZONE_0]);
            solver_container[ZONE_0][MESH_0][config_container[ZONE_0]->GetContainerPosition(RUNTIME_FLOW_SYS)]->Set_MPI_Solution(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
            solver_container[ZONE_0][MESH_0][config_container[ZONE_0]->GetContainerPosition(RUNTIME_FLOW_SYS)]->Preprocessing(geometry_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_0], config_container[ZONE_0], MESH_0, 0, RUNTIME_FLOW_SYS, false);
          }

          if (rank == MASTER_NODE) cout << endl << "-------------------------- File Output Summary --------------------------";

          /*--- Execute the routine for writing restart, volume solution,
           surface solution, and surface comma-separated value files. ---*/
          output->SetResult_Files(solver_container, geometry_container, config_container, ExtIter, nZone);

          /*--- Output a file with the forces breakdown. ---*/
          output->SetForces_Breakdown(geometry_container, solver_container,
                                      config_container, integration_container, ZONE_0);

          /*--- Compute the forces at different sections. ---*/
          if (config_container[ZONE_0]->GetPlot_Section_Forces())
            output->SetForceSections(solver_container[ZONE_0][MESH_0][FLOW_SOL], geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ExtIter);

          if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;

  }

  return StopCalc;

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void SU2Solver::mapRigidDisplacementOnFluidMesh_Old(){

  int rank = MASTER_NODE;
  int size = 1;

  unsigned short iZone;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if(rank == MASTER_NODE){
    cout << "Displacment to map" << endl;
    cout << interfRigidDispArray[0] << endl;
    cout << interfRigidDispArray[1] << endl;
    cout << interfRigidDispArray[2] << endl;
    cout << interfRigidDispArray[3] << endl;
    cout << interfRigidDispArray[4] << endl;
    cout << interfRigidDispArray[5] << endl;
    cout << " " << endl;
  }

  for(iZone = 0; iZone < nZone; iZone++){

    if(rank == MASTER_NODE) cout << "Updating node locations with rigid body motion." << endl;
    surface_movement[iZone]->SetFSIRigidBodyMotion(geometry_container[iZone][MESH_0], config_container[iZone], interfRigidDispArray, Center_n);
    if(rank == MASTER_NODE) cout << " Deforming the volume grid." << endl;
    grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone][MESH_0], config_container[iZone], true);

  }

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void SU2Solver::getSolidDisplacement(double** solidCoordinate,const double* CenterCoordinate){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /*--- Convert the matrix form of the solid interface into a buffer (vector) for MPI communication ---*/
  MatrixToVec(ROW_MAJ,solidCoordinate,solidInterfaceBuffer,nSolidInterfaceVertex,4,nSolidInterfaceVertex*4);

  Center[0] = CenterCoordinate[0];
  Center[1] = CenterCoordinate[1];
  Center[2] = CenterCoordinate[2];

}

void SU2Solver::mapRigidDisplacementOnFluidMesh(){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short iMarker, Moving, equivIndex, jMarker, iDim;
  unsigned long iVertex, iPoint, GlobalIndex;
  string Marker_Tag, Moving_Tag;
  su2double* Coord;
  su2double newCoord[3], VarCoord[3];

  /*--- Broadcast the solid interface position (buffer form). ---*/
#ifdef HAVE_MPI
  MPI_Bcast(solidInterfaceBuffer, nSolidInterfaceVertex*4, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Bcast(Center,3,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
#endif // HAVE_MPI

  /*--- Convert the buffer form of the solid interface into a matrix after MPI communication. ---*/
  VecToMatrix(ROW_MAJ,solidInterface, solidInterfaceBuffer,nSolidInterfaceVertex,4,nSolidInterfaceVertex*4);


  /*--- Loop over markers and find the particular marker(s) (surface) to deform or move. ---*/
  for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
    Moving = config_container[ZONE_0]->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config_container[ZONE_0]->GetnMarker_Moving(); jMarker++) {
        Moving_Tag = config_container[ZONE_0]->GetMarker_Moving(jMarker);
        Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Moving_Tag) {

          /*--- Loop over vertices and compute the matching between fluid and solid interface meshes. ---*/
          //NB : Currently only a one-to-one relation is enabled (= matching mesh)
          for(iVertex = 0; iVertex < geometry_container[ZONE_0][MESH_0]->nVertex[iMarker]; iVertex++) {
            iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
            GlobalIndex = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetGlobalIndex();
            Coord = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();
            equivIndex = searchEquivIndex(solidInterface, nSolidInterfaceVertex, GlobalIndex); //This currently represents the mesh interpolation procedure

            /*--- Calculate delta change in the x, y, & z directions ---*/
            for(iDim=0; iDim < nDim; iDim++){
              newCoord[iDim] = solidInterface[equivIndex][iDim+1];
              VarCoord[iDim] = newCoord[iDim] - Coord[iDim];
            }

            /*--- Set node displacement for volume deformation ---*/
            geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
          }
        }
      }
    }
  }

  for (jMarker=0; jMarker<config_container[ZONE_0]->GetnMarker_Moving(); jMarker++) {

    /*-- Check if we want to update the motion origin for the given marker ---*/

    if (config_container[ZONE_0]->GetMoveMotion_Origin(jMarker) == YES) {
      if(rank == MASTER_NODE){
        cout << "Updating the motion center" << endl;
      }
      config_container[ZONE_0]->SetMotion_Origin_X(jMarker, Center[0]);
      config_container[ZONE_0]->SetMotion_Origin_Y(jMarker, Center[1]);
      config_container[ZONE_0]->SetMotion_Origin_Z(jMarker, Center[2]);
    }
  }

  /*--- Set the moment computation center to the new location after
   incrementing the position with the plunging. ---*/

  for (jMarker=0; jMarker<config_container[ZONE_0]->GetnMarker_Monitoring(); jMarker++) {
    if (rank == MASTER_NODE){
      cout << "Updating the moment computation center" << endl;
      cout << "New location of the center : " << Center[0] << " " << Center[1] << " " << Center[2] << endl;
    }
    config_container[ZONE_0]->SetRefOriginMoment_X(jMarker, Center[0]);
    config_container[ZONE_0]->SetRefOriginMoment_Y(jMarker, Center[1]);
    config_container[ZONE_0]->SetRefOriginMoment_Z(jMarker, Center[2]);
  }

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if(rank == MASTER_NODE) cout << " Deforming the volume grid." << endl;
  grid_movement[ZONE_0]->SetVolume_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], true);

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

unsigned short SU2Solver::searchEquivIndex(double** target, unsigned long size_target, unsigned long iPoint){

  bool indexFound(false);
  int ii(0);

  while(indexFound == false && ii<(size_target)){
    if(target[ii][0] == iPoint){
      indexFound = true;
      break;
    }
    else ii++;
  }

  if(indexFound == false){
    cout << "INDEX IS NOT FOUND" << endl;
  }

  return ii;

}

int SU2Solver::searchNumberOfOccurence(double** target, unsigned long size_target, unsigned long iPoint){

  int nboccur(0);

  for(int ii=0; ii<size_target; ii++){
    if(target[ii][0] == iPoint)
      nboccur++;
  }

  return nboccur;
}

void SU2Solver::dynamicMeshUpdate(unsigned long ExtIter){

  int rank = MASTER_NODE;
  int size = 1;

  unsigned short iZone;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  for(iZone = 0; iZone < nZone; iZone++){

    if(rank == MASTER_NODE) cout << " Computing grid velocities by finite differencing." << endl;
    geometry_container[iZone][MESH_0]->SetGridVelocity(config_container[iZone], ExtIter);

    if(rank == MASTER_NODE) cout << " Updating multigrid structure." << endl;
    grid_movement[iZone]->UpdateMultiGrid(geometry_container[iZone], config_container[iZone]);
  }

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void SU2Solver::staticMeshUpdate(){

  int rank = MASTER_NODE;
  int size = 1;

  unsigned short iZone;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  for(iZone = 0; iZone < nZone; iZone++){

    if(rank == MASTER_NODE) cout << " Updating multigrid structure." << endl;
    grid_movement[iZone]->UpdateMultiGrid(geometry_container[iZone], config_container[iZone]);
  }

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void SU2Solver::setInitialMesh_Old(bool restart){

  unsigned short iZone, iMesh;
  unsigned long iPoint;

  /*--- Get the center of rotation (rigid body motion) as it is located in the input mesh ---*/
  Center_n[0] = config->GetMotion_Origin_X(0);
  Center_n[1] = config->GetMotion_Origin_Y(0);
  Center_n[2] = config->GetMotion_Origin_Z(0);

  /*--- First static mesh deformation according to the initial position of the solid interface ---*/
  mapRigidDisplacementOnFluidMesh_Old();
  staticMeshUpdate();

  /*--- Propagate the initial deformation to the past ---*/
  //if (!restart){
    for(iZone = 0; iZone < nZone; iZone++){
      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
        for(iPoint = 0; iPoint < geometry_container[iZone][iMesh]->GetnPoint(); iPoint++){
          //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
          geometry_container[iZone][iMesh]->node[iPoint]->SetVolume_n();
          geometry_container[iZone][iMesh]->node[iPoint]->SetVolume_nM1();
          geometry_container[iZone][iMesh]->node[iPoint]->SetCoord_n();
          geometry_container[iZone][iMesh]->node[iPoint]->SetCoord_n1();
        }
      }
    }
  //}

  /*--- Update the center of rotation with the initial displacement of the solide body ---*/
  Center_n[0] += interfRigidDispArray[0];
  Center_n[1] += interfRigidDispArray[1];
  Center_n[2] += interfRigidDispArray[2];

  attitude_i1[0] = 0.0;
  attitude_i1[1]  = 0.0;
  attitude_i1[2]  = 0.0;

  /*--- Reset rigid displacement values of the solid body ---*/
  interfRigidDispArray[0] = 0.0;
  interfRigidDispArray[1] = 0.0;
  interfRigidDispArray[2] = 0.0;
  interfRigidDispArray[3] = 0.0;
  interfRigidDispArray[4] = 0.0;
  interfRigidDispArray[5] = 0.0;

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void SU2Solver::setInitialMesh(bool restart){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short iZone, iMesh,iDim;
  unsigned long iPoint;

  /*--- First static mesh deformation according to the initial position of the solid interface ---*/
  mapRigidDisplacementOnFluidMesh();
  staticMeshUpdate();

  /*--- Propagate the initial deformation to the past ---*/
  //if (!restart){
    for(iZone = 0; iZone < nZone; iZone++){
      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
        for(iPoint = 0; iPoint < geometry_container[iZone][iMesh]->GetnPoint(); iPoint++){
          //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
          geometry_container[iZone][iMesh]->node[iPoint]->SetVolume_n();
          geometry_container[iZone][iMesh]->node[iPoint]->SetVolume_nM1();
          geometry_container[iZone][iMesh]->node[iPoint]->SetCoord_n();
          geometry_container[iZone][iMesh]->node[iPoint]->SetCoord_n1();
        }
      }
    }
  //}

#ifdef HAVE_MPI
/*--- Synchronization point ---*/
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

double SU2Solver::outputFluidLoads_Drag(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag, RefDensity, RefAreaCoeff, RefVel2(0.0), factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefAreaCoeff = config_container[val_iZone]->GetRefAreaCoeff();

  /*--- Calculate free-stream velocity (squared) ---*/
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = 0.5*RefDensity*RefAreaCoeff*RefVel2;
  CDrag = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();

  return CDrag*factor;
}

double SU2Solver::outputFluidLoads_Lift(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, RefDensity, RefAreaCoeff, RefVel2(0.0), factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefAreaCoeff = config_container[val_iZone]->GetRefAreaCoeff();

  /*--- Calculate free-stream velocity (squared) ---*/
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = 0.5*RefDensity*RefAreaCoeff*RefVel2;
  CLift = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();

  return CLift*factor;
}

void SU2Solver::outputGlobalFluidLoads(double* globalFluidLoad){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, CDrag, CMz, RefDensity, RefAreaCoeff, RefLengthCoeff, RefVel2(0.0), factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefAreaCoeff = config_container[val_iZone]->GetRefAreaCoeff();
  RefLengthCoeff = config_container[val_iZone]->GetRefLengthMoment();

  /*--- Calculate free-stream velocity (squared) ---*/
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate drag and lift force based on coefficients ---*/
  factor = 0.5*RefDensity*RefAreaCoeff*RefVel2;
  CLift = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
  CDrag = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
  CMz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();

  /*--- Fill the global fluid loads array coming from external solid solver ---*/
  globalFluidLoad[0] = CLift*factor;
  globalFluidLoad[1] = CDrag*factor;
  globalFluidLoad[2] = CMz*factor*RefLengthCoeff;

}

void SU2Solver::mergeSurfaceLoads(){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    unsigned short iMarker, Moving, jMarker, iDim, jDim;
	unsigned long iPoint, iVertex, GlobalIndex;
	su2double *Normal, RefDensity, RefAreaCoeff, RefVel2, factor;
	string Marker_Tag, Moving_Tag;
	int *rcvcounts, *displ, somDispl(0);

    unsigned short FinestMesh = config_container[ZONE_0]->GetFinestMesh();

	/*--- Check the kind of fluid problem ---*/
	bool compressible       = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
	bool incompressible     = (config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
	bool viscous_flow       = ((config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES) ||
							   (config_container[ZONE_0]->GetKind_Solver() == RANS) );

	/*--- Parameters for the calculations ---*/
	// Pn: Pressure
	// Pinf: Pressure_infinite
	// div_vel: Velocity divergence
	// Dij: Dirac delta
	su2double Pn = 0.0, div_vel = 0.0, Dij = 0.0;
	su2double Viscosity = 0.0;
	su2double Grad_Vel[3][3] = { {0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} } ;
	su2double Tau[3][3] = { {0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} ,
							{0.0, 0.0, 0.0} } ;

	su2double Pinf = solver_container[ZONE_0][FinestMesh][FLOW_SOL]->GetPressure_Inf();

    /*--- Loop over markers and find the particular marker(s) (surface) to deform or move. ---*/
    for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
    Moving = config_container[ZONE_0]->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config_container[ZONE_0]->GetnMarker_Moving(); jMarker++) {
        Moving_Tag = config_container[ZONE_0]->GetMarker_Moving(jMarker);
        Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Moving_Tag) {
          /*--- Loop over vertices and compute the traction vector in the cartesian frame. ---*/
          for(iVertex = 0; iVertex < nLocalFluidInterfaceVertex; iVertex++) {
	        iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
	        GlobalIndex = geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetGlobalIndex();
	        /*--- It necessary to distinguish the halo nodes from the others, since they introduice non physical forces. The halo nodes are marked with a global index equal to -1 ---*/
	        if(geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetDomain()) partFluidSurfaceLoads[iVertex][0] = GlobalIndex;
	        else partFluidSurfaceLoads[iVertex][0] = -1.0;

		    /*--- Get the normal at the vertex: this normal goes inside the fluid domain. ---*/
	        Normal = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

	        /*--- Get the values of pressure and viscosity ---*/
	        if (incompressible){
		      Pn = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetPressureInc();
		      if (viscous_flow){
                for(iDim=0; iDim<nDim; iDim++){
                  for(jDim=0; jDim<nDim; jDim++){
                    Grad_Vel[iDim][jDim] = solver_container[ZONE_0][FinestMesh][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
                  }
                }
			    Viscosity = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
		      }
	        }
	        else if (compressible){
		      Pn = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetPressure();
		      if (viscous_flow){
			    for(iDim=0; iDim<nDim; iDim++){
			      for(jDim=0; jDim<nDim; jDim++){
			        Grad_Vel[iDim][jDim] = solver_container[ZONE_0][FinestMesh][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
			      }
			    }
			    Viscosity = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
		      }
	        }

	       /*--- Calculate the inviscid (pressure) part of tn in the fluid nodes (force units) ---*/
           for (iDim = 0; iDim < nDim; iDim++) {
		     partFluidSurfaceLoads[iVertex][iDim+1] = -(Pn-Pinf)*Normal[iDim];   //NB : norm(Normal) = Area
	       }

	       /*--- Calculate the viscous (shear stress) part of tn in the fluid nodes (force units ---*/
	       if ((incompressible || compressible) && viscous_flow){
		     div_vel = 0.0;
		     for (iDim = 0; iDim < nDim; iDim++)
		       div_vel += Grad_Vel[iDim][iDim];
		     if (incompressible) div_vel = 0.0;

		     for (iDim = 0; iDim < nDim; iDim++) {
               for (jDim = 0 ; jDim < nDim; jDim++) {
				 Dij = 0.0; if (iDim == jDim) Dij = 1.0;
				 Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*Dij;
				 partFluidSurfaceLoads[iVertex][iDim+1] += Tau[iDim][jDim]*Normal[jDim];
			   }
		     }
	       }
	       /*// Redimensionalize and take into account ramp transfer of the loads
	       for (iDim = 0; iDim < nDim; iDim++){
		     //Donor_Variable[iVar] = Donor_Variable[iVar] * Physical_Constants[0] * Physical_Constants[1];
           }*/
          }
        }
       }
    }
  }

  /*--- Be sure that all z components are zero in case of 2D computation ---*/
  if(nDim == 2){
    for(iVertex = 0; iVertex<nLocalFluidInterfaceVertex; iVertex++)
      partFluidSurfaceLoads[iVertex][3] = 0.0;
  }

  /*--- Create the sending and receiving buffers for MPI communications ---*/
  su2double* SendBuff;
  unsigned int sizeSendBuff = nLocalFluidInterfaceVertex*4;
  su2double* RecBuff;
  unsigned int sizeRecBuff = nAllFluidInterfaceVertex*4;

  SendBuff = new su2double [sizeSendBuff];

  /*--- Convert the local surface loads array into the sending buffer before MPI communication. ---*/
  MatrixToVec(ROW_MAJ, partFluidSurfaceLoads, SendBuff, nLocalFluidInterfaceVertex,4,sizeSendBuff);

  /*--- Pre-processing of the MPI_Gatherv ---*/
  if(rank==MASTER_NODE){
    displ = new int [size];
    rcvcounts = new int [size];
    RecBuff = new su2double [sizeRecBuff];
  }

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(&sizeSendBuff, 1, MPI_INT, rcvcounts, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
#endif // HAVE_MPI

if(rank == MASTER_NODE){
    displ[0] = 0;
    for(int i=0; i<size; i++){
      somDispl += rcvcounts[i-1];
      displ[i] = somDispl;
    }
  }

/*--- Gather the partitioned surface loads to the MASTER_NODE ---*/
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gatherv(SendBuff, sizeSendBuff, MPI_DOUBLE, RecBuff, rcvcounts, displ, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#endif // HAVE_MPI

if(rank == MASTER_NODE){
    VecToMatrix(ROW_MAJ, fluidSurfaceloads, RecBuff, nAllFluidInterfaceVertex, 4, sizeRecBuff);
    delete [] RecBuff;
    delete [] displ;
    delete [] rcvcounts;
  }

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif // HAVE_MPI

  delete [] SendBuff;

}

void SU2Solver::outputSurfaceLoads(double** solidSurfaceLoads){

  unsigned long iVertex, GlobalIndex, equivIndex;
  //int NbOccur;

  for(iVertex=0; iVertex<nSolidInterfaceVertex; iVertex++){
    //NbOccur = 0;
    GlobalIndex = solidSurfaceLoads[iVertex][0];
    //NbOccur = searchNumberOfOccurence(fluidSurfaceloads, nAllFluidInterfaceVertex, GlobalIndex);
    //if(NbOccur > 1) cout << "There are " << NbOccur << " occurences of the node " << GlobalIndex << endl;
    equivIndex = searchEquivIndex(fluidSurfaceloads, nAllFluidInterfaceVertex, GlobalIndex);
    //cout << GlobalIndex << "<---->" << fluidSurfaceloads[equivIndex][0] << endl;
    solidSurfaceLoads[iVertex][1] = fluidSurfaceloads[equivIndex][1];
    solidSurfaceLoads[iVertex][2] = fluidSurfaceloads[equivIndex][2];
    solidSurfaceLoads[iVertex][3] = fluidSurfaceloads[equivIndex][3];
  }
}

/*double SU2Solver::computeRigidDisplacementNorm() const {

  double norm(0.0);
  int i;

  for(int i=0; i<5; i++){
    norm += interfRigidDispArray[i]*interfRigidDispArray[i];
  }

  return sqrt(norm);
}*/

void SU2Solver::broadcastInterfDisp(){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /*--- Broadcast an array of 6 doubles (for rigid body displacement) ---*/
  MPI_Bcast(interfRigidDispArray, 6, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);


}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void MatrixToVec(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec){

    /*--- Convert a matrix into a vector with row or column order ---*/
    if(order == COL_MAJ){
      for(int j=0; j<sizeVec; j++)
        vecteur[j] = matrix[j%Nrow][j/Nrow];
    }
    else if(order == ROW_MAJ){
      for(int j=0; j<sizeVec; j++)
        vecteur[j] = matrix[j/Ncol][j%Ncol];
    }
    else{
      cerr << "Wrong storage order" << endl;
      throw(-1);
    }
}

void VecToMatrix(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec){

    /*--- Convert a vector into a matrix with row or column order ---*/
    if(order == COL_MAJ){
      for(int i=0; i<Nrow;i++)
        for(int j=0; j<Ncol; j++)
          matrix[i][j] = vecteur[j*Nrow+i];
    }
    else if(order == ROW_MAJ){
      for(int i=0; i<Nrow;i++)
        for(int j=0; j<Ncol; j++)
          matrix[i][j] = vecteur[i*Ncol+j];
    }
    else{
      cerr << "Wrong storage order" << endl;
      throw(-1);
    }
}
