/*!
 * \file SU2_ERR.cpp
 * \brief Main file of the Gradient-Norm Error Estimation Code (SU2_ERR).
 * \author B. Munguia
 * \version 5.0.0 "Raven"
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

#include "../include/SU2_ERR.hpp"
using namespace std;


int main(int argc, char *argv[]) {
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  
  char config_file_name[MAX_STRING_SIZE], *cstr;
  ofstream Gradient_file;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  /*--- MPI initialization, and buffer setting ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Init(&argc,&argv);
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
  MPI_Comm_rank(MPICommunicator,&rank);
  MPI_Comm_size(MPICommunicator,&size);
#else
  SU2_Comm MPICommunicator(0);
#endif
  
  /*--- Pointer to different structures that will be used throughout the entire code ---*/
  
  CConfig **config_container          = NULL;
  CGeometry **geometry_container      = NULL;
  CSurfaceMovement *surface_movement  = NULL;
  CVolumetricMovement *mesh_movement  = NULL;
  
  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/
  
  if (argc == 2) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }
  
  /*--- Definition of the containers per zones ---*/
  
  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];
  
  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
  }
  
  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/
    
    config_container[iZone] = new CConfig(config_file_name, SU2_ERR, iZone, nZone, 0, VERB_HIGH);

    /*--- Set the MPI communicator ---*/
    config_container[iZone]->SetMPICommunicator(MPICommunicator);
        
    /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/
    
    CGeometry *geometry_aux = NULL;
    
    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/
    
    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);
    
    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/
    
    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);
    
    /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/
    
    geometry_container[iZone] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);
    
    /*--- Deallocate the memory of geometry_aux ---*/
    
    delete geometry_aux;
    
    /*--- Add the Send/Receive boundaries ---*/
    
    geometry_container[iZone]->SetSendReceive(config_container[iZone]);
    
    /*--- Add the Send/Receive boundaries ---*/
    
    geometry_container[iZone]->SetBoundaries(config_container[iZone]);
    
  }
  
  /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/
  
#ifdef HAVE_MPI
  StartTime = MPI_Wtime();
#else
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
  
  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting local point connectivity." <<endl;
  geometry_container[ZONE_0]->SetPoint_Connectivity();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." <<endl;
  geometry_container[ZONE_0]->SetBoundVolume();
  geometry_container[ZONE_0]->Check_IntElem_Orientation(config_container[ZONE_0]);
  geometry_container[ZONE_0]->Check_BoundElem_Orientation(config_container[ZONE_0]);
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
  geometry_container[ZONE_0]->SetEdges(); geometry_container[ZONE_0]->SetVertex(config_container[ZONE_0]);
  
  /*--- Compute center of gravity ---*/
  
  if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
  geometry_container[ZONE_0]->SetCoord_CG();
  
  /*--- Create the dual control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);

  /*--- Store the global to local mapping after preprocessing. ---*/
 
  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
  geometry_container[ZONE_0]->SetGlobal_to_Local_Point();
 
  /*--- Load the surface sensitivities from file. This is done only
   once: if this is an unsteady problem, a time-average of the surface
   sensitivities at each node is taken within this routine. ---*/
  
  if (!config_container[ZONE_0]->GetDiscrete_Adjoint()){
    if (rank == MASTER_NODE) cout << "Reading surface sensitivities at each node from file." << endl;
    geometry_container[ZONE_0]->SetBoundSensitivity(config_container[ZONE_0]);
  }
  else {
    if (rank == MASTER_NODE) cout << "Reading volume sensitivities at each node from file." << endl;
    mesh_movement = new CVolumetricMovement(geometry_container[ZONE_0], config_container[ZONE_0]);
    geometry_container[ZONE_0]->SetSensitivity(config_container[ZONE_0]);
    
    if (rank == MASTER_NODE)
      cout << endl <<"---------------------- Mesh sensitivity computation ---------------------" << endl;
    mesh_movement->SetVolume_Deformation(geometry_container[ZONE_0], config_container[ZONE_0], false, true);
    
    COutput *output = new COutput();
    output->SetSensitivity_Files(geometry_container, config_container, nZone);
  }
  
  /*--- Definition of the Class for surface deformation ---*/
  
  surface_movement = new CSurfaceMovement();
  
  /*--- Copy coordinates to the surface structure ---*/
  
  surface_movement->CopyBoundary(geometry_container[ZONE_0], config_container[ZONE_0]);
  
  if (config_container[ZONE_0]->GetDesign_Variable(0) != NONE){
    if (rank == MASTER_NODE)
      cout << endl <<"---------- Start gradient evaluation using sensitivity information ----------" << endl;
    
    /*--- Write the gradient in a external file ---*/
    
    if (rank == MASTER_NODE) {
      cstr = new char [config_container[ZONE_0]->GetObjFunc_Grad_FileName().size()+1];
      strcpy (cstr, config_container[ZONE_0]->GetObjFunc_Grad_FileName().c_str());
      Gradient_file.open(cstr, ios::out);
    }
    
    /*--- If AD mode is enabled we can use it to compute the projection,
     * otherwise we use finite differences. ---*/
    
    if (config_container[ZONE_0]->GetAD_Mode()){
      SetProjection_AD(geometry_container[ZONE_0], config_container[ZONE_0], surface_movement, Gradient_file);
    }/*else{
      SetProjection_FD(geometry_container[ZONE_0], config_container[ZONE_0], surface_movement, Gradient_file);
    }*/
    
    if (rank == MASTER_NODE)
      Gradient_file.close();
  }
  
  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/
  
#ifdef HAVE_MPI
  StopTime = MPI_Wtime();
#else
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
  
  /*--- Compute/print the total time for performance benchmarking. ---*/
  
  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }
  
  /*--- Exit the solver cleanly ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_ERR) ------------------------" << endl << endl;
  
  /*--- Finalize MPI parallelization ---*/
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}

void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, ofstream& Gradient_file){
  
  su2double DV_Value, *VarCoord, Sensitivity, **Gradient, my_Gradient, *Normal, Area = 0.0;
  unsigned short iDV_Value = 0, iMarker, nMarker, iDim, nDim, iDV, nDV, nDV_Value;
  unsigned long iVertex, nVertex, iPoint;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();
  
  VarCoord = NULL;
  
  /*--- Structure to store the gradient ---*/
  
  Gradient = new su2double*[nDV];
  
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    Gradient[iDV] = new su2double[nDV_Value];
  }
  
  /*--- Discrete adjoint gradient computation ---*/
  
  if (rank == MASTER_NODE)
    cout << "Evaluate functional gradient using Algorithmic Differentiation." << endl;
  
  /*--- Start recording of operations ---*/
  
  AD::StartRecording();
  
  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/
  
  
  
  for (iDV = 0; iDV < nDV; iDV++){
    
    nDV_Value =  config->GetnDV_Value(iDV);
    
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      
      /*--- Initilization with su2double resets the index ---*/
      
      DV_Value = 0.0;
      
      AD::RegisterInput(DV_Value);
      
      config->SetDV_Value(iDV, iDV_Value, DV_Value);
    }
  }
  
  /*--- Call the surface deformation routine ---*/
  
  surface_movement->SetSurface_Deformation(geometry, config);
  
  /*--- Stop the recording --- */
  
  AD::StopRecording();
  
  /*--- Initialize the derivatives of the output of the surface deformation routine
   * with the discrete adjoints from the CFD solution ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint      = geometry->vertex[iMarker][iVertex]->GetNode();
        VarCoord    = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        Normal      = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++){
          Area += Normal[iDim]*Normal[iDim];
        }
        Area = sqrt(Area);
        
        for (iDim = 0; iDim < nDim; iDim++){
          if (config->GetDiscrete_Adjoint()){
            Sensitivity = geometry->GetSensitivity(iPoint, iDim);
          } else {
            Sensitivity = -Normal[iDim]*geometry->vertex[iMarker][iVertex]->GetAuxVar()/Area;
          }
          SU2_TYPE::SetDerivative(VarCoord[iDim], SU2_TYPE::GetValue(Sensitivity));
        }
      }
    }
  }
  
  /*--- Compute derivatives and extract gradient ---*/
  
  AD::ComputeAdjoint();
  
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      DV_Value = config->GetDV_Value(iDV, iDV_Value);
      my_Gradient = SU2_TYPE::GetDerivative(DV_Value);
#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&my_Gradient, &Gradient[iDV][iDV_Value], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      Gradient[iDV][iDV_Value] = my_Gradient;
#endif
      
      /*--- Angle of Attack design variable (this is different,
       the value comes form the input file) ---*/
      
      if ((config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) ||
          (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK))  {
        Gradient[iDV][iDV_Value] = config->GetAoA_Sens();
      }
      
    }
  }
  
  /*--- Print gradients to screen and file ---*/
  
  OutputGradient(Gradient, config, Gradient_file);
  
  for (iDV = 0; iDV  < nDV; iDV++){
    delete [] Gradient[iDV];
  }
  delete [] Gradient;
}
