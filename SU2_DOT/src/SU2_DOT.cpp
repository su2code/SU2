/*!
 * \file SU2_DOT.cpp
 * \brief Main file of the Gradient Projection Code (SU2_DOT).
 * \author F. Palacios, T. Economon
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

#include "../include/SU2_DOT.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  
  char config_file_name[MAX_STRING_SIZE], *cstr;
  ofstream Gradient_file;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;

  su2double** Gradient;
  unsigned short iDV, iDV_Value;

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
  CSurfaceMovement **surface_movement = NULL;
  CVolumetricMovement **grid_movement = NULL;
  
  /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
   file is specified, default.cfg is used) ---*/
  
  if (argc == 2) { strcpy(config_file_name,argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig *config = NULL;
  config = new CConfig(config_file_name, SU2_DEF);

  nZone = CConfig::GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);


  /*--- Definition of the containers per zones ---*/
  
  config_container = new CConfig*[nZone];
  geometry_container = new CGeometry*[nZone];
  surface_movement   = new CSurfaceMovement*[nZone];
  grid_movement      = new CVolumetricMovement*[nZone];
  
  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]       = NULL;
    geometry_container[iZone]     = NULL;
    grid_movement [iZone]     = NULL;
    surface_movement[iZone]   = NULL;
  }
  
  /*--- Loop over all zones to initialize the various classes. In most
   cases, nZone is equal to one. This represents the solution of a partial
   differential equation on a single block, unstructured mesh. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/
    
    config_container[iZone] = new CConfig(config_file_name, SU2_DOT, iZone, nZone, 0, VERB_HIGH);

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
  
  for (iZone = 0; iZone < nZone; iZone++){
  
  if (rank == MASTER_NODE)
    cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
  
  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting local point connectivity." <<endl;
    geometry_container[iZone]->SetPoint_Connectivity();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation of the elements." <<endl;
    geometry_container[iZone]->SetBoundVolume();
    geometry_container[iZone]->Check_IntElem_Orientation(config_container[iZone]);
    geometry_container[iZone]->Check_BoundElem_Orientation(config_container[iZone]);
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
    geometry_container[iZone]->SetEdges(); geometry_container[iZone]->SetVertex(config_container[iZone]);
  
  /*--- Compute center of gravity ---*/
  
  if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
  geometry_container[iZone]->SetCoord_CG();
  
  /*--- Create the dual control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry_container[iZone]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);

  /*--- Store the global to local mapping after preprocessing. ---*/
 
  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
  geometry_container[iZone]->SetGlobal_to_Local_Point();
 
  /*--- Load the surface sensitivities from file. This is done only
   once: if this is an unsteady problem, a time-average of the surface
   sensitivities at each node is taken within this routine. ---*/
    if (!config_container[iZone]->GetDiscrete_Adjoint()){
      if (rank == MASTER_NODE) cout << "Reading surface sensitivities at each node from file." << endl;
      geometry_container[iZone]->SetBoundSensitivity(config_container[iZone]);
    } else {
      if (rank == MASTER_NODE) cout << "Reading volume sensitivities at each node from file." << endl;
      grid_movement[iZone] = new CVolumetricMovement(geometry_container[iZone], config_container[iZone]);
      geometry_container[iZone]->SetSensitivity(config_container[iZone]);

      if (rank == MASTER_NODE)
        cout << endl <<"---------------------- Mesh sensitivity computation ---------------------" << endl;
      grid_movement[iZone]->SetVolume_Deformation(geometry_container[iZone], config_container[iZone], false, true);

    }
  }

   if (config_container[ZONE_0]->GetDiscrete_Adjoint()){
     if (rank == MASTER_NODE)
       cout << endl <<"------------------------ Mesh sensitivity Output ------------------------" << endl;
     COutput *output = new COutput(config_container[ZONE_0]);
     output->SetSensitivity_Files(geometry_container, config_container, nZone);
   }

   if (config_container[ZONE_0]->GetDesign_Variable(0) != NONE){

     /*--- Initialize structure to store the gradient ---*/

     Gradient = new su2double*[config_container[ZONE_0]->GetnDV()];

     for (iDV = 0; iDV  < config_container[ZONE_0]->GetnDV(); iDV++){
       Gradient[iDV] = new su2double[config_container[ZONE_0]->GetnDV_Value(iDV)];
       for (iDV_Value = 0; iDV_Value < config_container[ZONE_0]->GetnDV_Value(iDV); iDV_Value++){
         Gradient[iDV][iDV_Value] = 0.0;
       }
     }

     if (rank == MASTER_NODE)
       cout << endl <<"---------- Start gradient evaluation using sensitivity information ----------" << endl;

     /*--- Write the gradient in a external file ---*/

     if (rank == MASTER_NODE) {
       cstr = new char [config_container[ZONE_0]->GetObjFunc_Grad_FileName().size()+1];
       strcpy (cstr, config_container[ZONE_0]->GetObjFunc_Grad_FileName().c_str());
       Gradient_file.open(cstr, ios::out);
     }

     /*--- Loop through each zone and add it's contribution to the gradient array ---*/

     for (iZone = 0; iZone < nZone; iZone++){

       /*--- Definition of the Class for surface deformation ---*/

       surface_movement[iZone] = new CSurfaceMovement();

       /*--- Copy coordinates to the surface structure ---*/

       surface_movement[iZone]->CopyBoundary(geometry_container[iZone], config_container[iZone]);

       /*--- If AD mode is enabled we can use it to compute the projection,
        *    otherwise we use finite differences. ---*/

       if (config_container[iZone]->GetAD_Mode()){
         SetProjection_AD(geometry_container[iZone], config_container[iZone], surface_movement[iZone] , Gradient);
       }else{
         SetProjection_FD(geometry_container[iZone], config_container[iZone], surface_movement[iZone] , Gradient);
       }
     }

     /*--- Compute sensitivities of transpiration boundary velocity ---*/
     if(config_container[ZONE_0]->GetnMarker_Transpiration() > 0 && config_container[ZONE_0]->GetAD_Mode()){
      if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Transpiration boundary sensitivitiy -----------------------" << endl;
      /*--- Initialize solver class ---*/
      CSolver **solver_container = new CSolver*[nZone];
       for (iZone = 0; iZone < nZone; iZone++) {
         solver_container[iZone] = NULL;
       }

       SolutionPostprocessing(geometry_container, config_container, solver_container, nZone);
       AD::StartRecording();

       /*--- Compute pressure at boundary using transpiration velocity ---*/
       ComputeTranspirationPressure(geometry_container[ZONE_0], config_container[ZONE_0], solver_container[ZONE_0]);
       solver_container[ZONE_0]->Pressure_Forces(geometry_container[ZONE_0], config_container[ZONE_0]);
       su2double Objective_Function = Compute_TotalObjFunc(config_container[ZONE_0], solver_container[ZONE_0]);

       if (rank==MASTER_NODE){
         SU2_TYPE::SetDerivative(Objective_Function,1.0);
       }else{
         SU2_TYPE::SetDerivative(Objective_Function,0.0);
       }
       AD::StopRecording();
       AD::ComputeAdjoint();

       su2double extracted_derivative;
       su2double dJdalpha = SU2_TYPE::GetDerivative(extracted_derivative);

       cout << "partial_J/partial_alpha = " << dJdalpha << "on rank " << rank << endl;
       
     }

     /*--- Print gradients to screen and file ---*/

     OutputGradient(Gradient, config_container[ZONE_0], Gradient_file);

     if (rank == MASTER_NODE)
       Gradient_file.close();

     for (iDV = 0; iDV  < config_container[ZONE_0]->GetnDV(); iDV++){
       delete [] Gradient[iDV];
     }
     delete [] Gradient;

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
    cout << endl <<"------------------------- Exit Success (SU2_DOT) ------------------------" << endl << endl;
  
  /*--- Finalize MPI parallelization ---*/
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}

void SetProjection_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double** Gradient){
  
  unsigned short iDV, nDV, iFFDBox, nDV_Value, iMarker, iDim;
  unsigned long iVertex, iPoint;
  su2double delta_eps, my_Gradient, localGradient, *Normal, dS, *VarCoord, Sensitivity,
  dalpha[3], deps[3], dalpha_deps;
  bool *UpdatePoint, MoveSurface, Local_MoveSurface;
  CFreeFormDefBox **FFDBox;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  
  nDV = config->GetnDV();
  
  /*--- Boolean controlling points to be updated ---*/
  
  UpdatePoint = new bool[geometry->GetnPoint()];
  
  /*--- Definition of the FFD deformation class ---*/
  
  unsigned short nFFDBox = MAX_NUMBER_FFD;
  FFDBox = new CFreeFormDefBox*[nFFDBox];

  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (nDV_Value != 1){
      cout << "The projection using finite differences currently only supports a fixed direction of movement for FFD points." << endl;
      exit(EXIT_FAILURE);
    }
  }

  /*--- Continuous adjoint gradient computation ---*/
  
  if (rank == MASTER_NODE)
    cout << "Evaluate functional gradient using Finite Differences." << endl;
  
  for (iDV = 0; iDV < nDV; iDV++) {
    
    MoveSurface = true;
    Local_MoveSurface = true;
    
    /*--- Free Form deformation based ---*/
    
    if ((config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
        (config->GetDesign_Variable(iDV) == FFD_NACELLE) ||
        (config->GetDesign_Variable(iDV) == FFD_GULL) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST) ||
        (config->GetDesign_Variable(iDV) == FFD_ROTATION) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS) ||
        (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK)) {
      
      /*--- Read the FFD information in the first iteration ---*/
      
      if (iDV == 0) {
        
        if (rank == MASTER_NODE)
          cout << "Read the FFD information from mesh file." << endl;
        
        /*--- Read the FFD information from the grid file ---*/
        
        surface_movement->ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());
        
        /*--- If the FFDBox was not defined in the input file ---*/
        if (!surface_movement->GetFFDBoxDefinition() && (rank == MASTER_NODE)) {
          cout << "The input grid doesn't have the entire FFD information!" << endl;
          cout << "Press any key to exit..." << endl;
          cin.get();
        }
        
        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
          
          if (rank == MASTER_NODE) cout << "Checking FFD box dimension." << endl;
          surface_movement->CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);
          
          if (rank == MASTER_NODE) cout << "Check the FFD box intersections with the solid surfaces." << endl;
          surface_movement->CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);
          
        }
        
        if (rank == MASTER_NODE)
          cout <<"-------------------------------------------------------------------------" << endl;
        
      }
      
      if (rank == MASTER_NODE) {
        cout << endl << "Design variable number "<< iDV <<"." << endl;
        cout << "Performing 3D deformation of the surface." << endl;
      }
      
      /*--- Apply the control point change ---*/
      
      MoveSurface = false;
      
      for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
        
        /*--- Reset FFD box ---*/
        
        switch (config->GetDesign_Variable(iDV) ) {
          case FFD_CONTROL_POINT_2D : Local_MoveSurface = surface_movement->SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER_2D :        Local_MoveSurface = surface_movement->SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS_2D :     Local_MoveSurface = surface_movement->SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST_2D :         Local_MoveSurface = surface_movement->SetFFDTwist_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_POINT :    Local_MoveSurface = surface_movement->SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_NACELLE :          Local_MoveSurface = surface_movement->SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_GULL :             Local_MoveSurface = surface_movement->SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST :            Local_MoveSurface = surface_movement->SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ROTATION :         Local_MoveSurface = surface_movement->SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER :           Local_MoveSurface = surface_movement->SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS :        Local_MoveSurface = surface_movement->SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_SURFACE :  Local_MoveSurface = surface_movement->SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ANGLE_OF_ATTACK :  Gradient[iDV][0] = config->GetAoA_Sens(); break;
        }
        
        /*--- Recompute cartesian coordinates using the new control points position ---*/
        
        if (Local_MoveSurface) {
          MoveSurface = true;
          surface_movement->SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, true);
        }
        
      }
      
    }
    
    /*--- Hicks Henne design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == HICKS_HENNE) {
      surface_movement->SetHicksHenne(geometry, config, iDV, true);
    }
    
    /*--- Surface bump design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SURFACE_BUMP) {
      surface_movement->SetSurface_Bump(geometry, config, iDV, true);
    }

    /*--- Kulfan (CST) design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == CST) {
      surface_movement->SetCST(geometry, config, iDV, true);
    }
    
    /*--- Displacement design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == TRANSLATION) {
      surface_movement->SetTranslation(geometry, config, iDV, true);
    }
    
    /*--- Angle of Attack design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) {
      Gradient[iDV][0] = config->GetAoA_Sens();
    }
    
    /*--- Scale design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == SCALE) {
      surface_movement->SetScale(geometry, config, iDV, true);
    }
    
    /*--- Rotation design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == ROTATION) {
      surface_movement->SetRotation(geometry, config, iDV, true);
    }
    
    /*--- NACA_4Digits design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == NACA_4DIGITS) {
      surface_movement->SetNACA_4Digits(geometry, config);
    }
    
    /*--- Parabolic design variable ---*/
    
    else if (config->GetDesign_Variable(iDV) == PARABOLIC) {
      surface_movement->SetParabolic(geometry, config);
    }
    
    /*--- Design variable not implement ---*/
    
    else {
      if (rank == MASTER_NODE)
        cout << "Design Variable not implement yet" << endl;
    }
    
    /*--- Load the delta change in the design variable (finite difference step). ---*/
    
    if ((config->GetDesign_Variable(iDV) != ANGLE_OF_ATTACK) &&
        (config->GetDesign_Variable(iDV) != FFD_ANGLE_OF_ATTACK)) {
      
      /*--- If the Angle of attack is not involved, reset the value of the gradient ---*/
      
      my_Gradient = 0.0; Gradient[iDV][0] = 0.0;
      
      if (MoveSurface) {
        
        delta_eps = config->GetDV_Value(iDV);
        
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
          UpdatePoint[iPoint] = true;
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if ((iPoint < geometry->GetnPointDomain()) && UpdatePoint[iPoint]) {
                
                Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
                Sensitivity = geometry->vertex[iMarker][iVertex]->GetAuxVar();
                
                dS = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dS += Normal[iDim]*Normal[iDim];
                  deps[iDim] = VarCoord[iDim] / delta_eps;
                }
                dS = sqrt(dS);
                
                dalpha_deps = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dalpha[iDim] = Normal[iDim] / dS;
                  dalpha_deps -= dalpha[iDim]*deps[iDim];
                }
                
                my_Gradient += Sensitivity*dalpha_deps;
                UpdatePoint[iPoint] = false;
              }
            }
          }
        }
        
      }
      
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    localGradient = my_Gradient;
#endif
    Gradient[iDV][0] += localGradient;
    }
  }
}
  

void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double** Gradient){

  su2double DV_Value, *VarCoord, Sensitivity, my_Gradient, localGradient, *Normal, Area = 0.0;
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

  /*--- Discrete adjoint gradient computation ---*/
  
  if (rank == MASTER_NODE)
    cout  << endl << "Evaluate functional gradient using Algorithmic Differentiation (ZONE " << config->GetiZone() << ")." << endl;

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
    SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      localGradient = my_Gradient;
#endif
      /*--- Angle of Attack design variable (this is different,
       the value comes form the input file) ---*/
      
      if ((config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) ||
          (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK))  {
        Gradient[iDV][iDV_Value] = config->GetAoA_Sens();
      }

      Gradient[iDV][iDV_Value] += localGradient;
    }
  }

  AD::Reset();

}

void OutputGradient(su2double** Gradient, CConfig* config, ofstream& Gradient_file){
  
  unsigned short nDV, iDV, iDV_Value, nDV_Value;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  
  nDV = config->GetnDV();
  
  /*--- Loop through all design variables and their gradients ---*/
  
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (rank == MASTER_NODE){
      
      /*--- Print the kind of design variable on screen ---*/
      
      cout << endl << "Design variable (";
      for (std::map<string, ENUM_PARAM>::const_iterator it = Param_Map.begin(); it != Param_Map.end(); ++it ){
        if (it->second == config->GetDesign_Variable(iDV)){
          cout << it->first << ") number "<< iDV << "." << endl;
        }
      }
      
      /*--- Print the kind of objective function to screen ---*/
      
      for (std::map<string, ENUM_OBJECTIVE>::const_iterator it = Objective_Map.begin(); it != Objective_Map.end(); ++it ){
        if (it->second == config->GetKind_ObjFunc()){
          cout << it->first << " gradient : ";
          if (iDV == 0) Gradient_file << it->first << " gradient " << endl;
        }
      }
      
      /*--- Print the gradient to file and screen ---*/
      
      for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
        cout << Gradient[iDV][iDV_Value];
        if (iDV_Value != nDV_Value-1 ){
          cout << ", ";
        }
        Gradient_file << Gradient[iDV][iDV_Value] << endl;
      }
      cout << endl;
      cout <<"-------------------------------------------------------------------------" << endl;
    }
  }
}

void SolutionPostprocessing(CGeometry **geometry_container, CConfig **config_container, CSolver **solver_container, unsigned short nZone){

  unsigned short iZone;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;

#ifdef HAVE_MPI
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
  MPI_Comm_rank(MPICommunicator,&rank);
  MPI_Comm_size(MPICommunicator,&size);
#else
  SU2_Comm MPICommunicator(0);
#endif

  /*--- Currently no FSI support ---*/
  if (config_container[ZONE_0]->GetWrt_Unsteady()) {

      /*--- Unsteady simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long iExtIter = 0;
      bool StopCalc = false;
      bool *SolutionInstantiated = new bool[nZone];

      for (iZone = 0; iZone < nZone; iZone++)
        SolutionInstantiated[iZone] = false;

      /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetWrt_Unsteady() && config_container[ZONE_0]->GetRestart())
        iExtIter = config_container[ZONE_0]->GetUnst_RestartIter();

      while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

        /*--- Check several conditions in order to merge the correct time step files. ---*/
        Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
        Physical_t  = (iExtIter+1)*Physical_dt;
        if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
          StopCalc = true;

        if ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
            ((iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq() == 0) && (iExtIter != 0) &&
             !((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
               (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND))) ||
            (StopCalc) ||
            (((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
              (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
             ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {

              /*--- Read in the restart file for this time step ---*/
              for (iZone = 0; iZone < nZone; iZone++) {

                /*--- Set the current iteration number in the config class. ---*/
                config_container[iZone]->SetExtIter(iExtIter);

                /*--- Either instantiate the solution class or load a restart file. ---*/
                if (SolutionInstantiated[iZone] == false &&
                    (iExtIter == 0 || (config_container[ZONE_0]->GetRestart() && ((long)iExtIter == config_container[ZONE_0]->GetUnst_RestartIter() ||
                                                                                  iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
                                                                                  iExtIter+1 == config_container[ZONE_0]->GetnExtIter())))) {
                  solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
                  SolutionInstantiated[iZone] = true;
                }
                  solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);
              }

            }

        iExtIter++;
        if (StopCalc) break;
      }

    }

    else if (config_container[ZONE_0]->GetUnsteady_Simulation() == HARMONIC_BALANCE) {

      /*--- Read in the restart file for this time step ---*/
      for (iZone = 0; iZone < nZone; iZone++) {

        /*--- Either instantiate the solution class or load a restart file. ---*/
        solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
        solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);

        /*--- Print progress in solution writing to the screen. ---*/
        if (rank == MASTER_NODE) {
          cout << "Storing the volume solution for time instance " << iZone << "." << endl;
        }

      }

    }

    else if (config_container[ZONE_0]->GetWrt_Dynamic()){

      /*--- Dynamic simulation: merge all unsteady time steps. First,
       find the frequency and total number of files to write. ---*/

      su2double Physical_dt, Physical_t;
      unsigned long iExtIter = 0;
      bool StopCalc = false;
      bool SolutionInstantiated = false;



      /*--- Check for an dynamic restart (structural analysis). Update ExtIter if necessary. ---*/
      if (config_container[ZONE_0]->GetKind_Solver() == FEM_ELASTICITY &&
          config_container[ZONE_0]->GetWrt_Dynamic() && config_container[ZONE_0]->GetRestart())
        iExtIter = config_container[ZONE_0]->GetDyn_RestartIter();

      while (iExtIter < config_container[ZONE_0]->GetnExtIter()) {

        /*--- Check several conditions in order to merge the correct time step files. ---*/
        /*--- If the solver is structural, the total and delta_t are obtained from different functions. ---*/

        Physical_dt = config_container[ZONE_0]->GetDelta_DynTime();
        Physical_t  = (iExtIter+1)*Physical_dt;
        if (Physical_t >=  config_container[ZONE_0]->GetTotal_DynTime())
          StopCalc = true;

        if ((iExtIter+1 == config_container[ZONE_0]->GetnExtIter()) ||
            (StopCalc) ||
            ((config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC) &&
             ((iExtIter == 0) || (iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0)))) {

              /*--- Set the current iteration number in the config class. ---*/
              config_container[ZONE_0]->SetExtIter(iExtIter);

              /*--- Read in the restart file for this time step ---*/
              for (iZone = 0; iZone < nZone; iZone++) {

                /*--- Either instantiate the solution class or load a restart file. ---*/
                if (SolutionInstantiated == false &&
                    (iExtIter == 0 || ((config_container[ZONE_0]->GetRestart() && (SU2_TYPE::Int(iExtIter) == config_container[ZONE_0]->GetDyn_RestartIter())) ||
                                       iExtIter % config_container[ZONE_0]->GetWrt_Sol_Freq_DualTime() == 0 ||
                                       iExtIter+1 == config_container[ZONE_0]->GetnExtIter()))) {
                  solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
                  SolutionInstantiated = true;
                }
                  solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);
              }

            }
        
        iExtIter++;
        if (StopCalc) break;
      }
      
      }
    
    else {

      /*--- Steady simulation: merge the single solution file. ---*/

      for (iZone = 0; iZone < nZone; iZone++) {
        /*--- Definition of the solution class ---*/
        solver_container[iZone] = new CBaselineSolver(geometry_container[iZone], config_container[iZone]);
        solver_container[iZone]->LoadRestart(geometry_container, &solver_container, config_container[iZone], SU2_TYPE::Int(MESH_0), true);
      }

      }
}

void ComputeTranspirationPressure(CGeometry *geometry, CConfig *config, CSolver *solver){
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned short iMarker, nMarker = config->GetnMarker_All();
  unsigned long iPoint, iVertex;
  su2double eps, *Vel_b, Vel_i, Vel_2, Density_b, Energy_b, Pressure_b;
  su2double Area, *Normal, *NormalArea, *UnitNormal;

  Vel_b = new su2double[nDim];
  Normal = new su2double[nDim];
  NormalArea = new su2double[nDim];
  UnitNormal = new su2double[nDim];

  string Marker_Tag;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if(config->GetMarker_All_KindBC(iMarker) == TRANSPIRATION){
     Marker_Tag  = config->GetMarker_All_TagBound(iMarker);
     for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
       iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
         if (geometry->node[iPoint]->GetDomain()) {

           geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
      
           Area = 0.0;
           for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
           Area = sqrt (Area);
      
           for (iDim = 0; iDim < nDim; iDim++) {
             NormalArea[iDim] = -Normal[iDim];
             UnitNormal[iDim] = -Normal[iDim]/Area;
           }

           for(iDim = 0; iDim < nDim; iDim++){
             Vel_b[iDim] = solver->node[iPoint]->GetVelocity(iDim);
           }
           eps = config->GetTranspiration(Marker_Tag);
           Density_b = solver->node[iPoint]->GetDensity();
           Energy_b = solver->node[iPoint]->GetEnergy();
           Vel_2  = 0.0;
           AD::RegisterInput(eps);
           for(iDim = 0; iDim < nDim; iDim++){
             Vel_i = Vel_b[iDim] - eps*UnitNormal[iDim];
             Vel_2 += Vel_i*Vel_i;
           }
           Pressure_b = (Energy_b - 0.5*Vel_2)*Density_b*(config->GetGamma()-1);

           solver->node[iPoint]->SetPressure(Pressure_b);
         }
       }
    }
  }

  delete [] Vel_b;
  delete [] Normal;
  delete [] NormalArea;
  delete [] UnitNormal;

}

su2double Compute_TotalObjFunc(CConfig *config, CSolver *solver) {
  
  unsigned short iMarker_Monitoring;
  su2double Weight_ObjFunc, Total_ComboObj = 0.0;
  
  /*--- Loop over all monitored markers, add to the 'combo' objective ---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    
    Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
    
    switch(config->GetKind_ObjFunc(iMarker_Monitoring)) {
      case DRAG_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CD(iMarker_Monitoring);
        break;
      case LIFT_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CL(iMarker_Monitoring);
        break;
      case SIDEFORCE_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CSF(iMarker_Monitoring);
        break;
      case EFFICIENCY:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CEff(iMarker_Monitoring);
        break;
      case MOMENT_X_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CMx(iMarker_Monitoring);
        break;
      case MOMENT_Y_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CMy(iMarker_Monitoring);
        break;
      case MOMENT_Z_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CMz(iMarker_Monitoring);
        break;
      case FORCE_X_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CFx(iMarker_Monitoring);
        break;
      case FORCE_Y_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CFy(iMarker_Monitoring);
        break;
      case FORCE_Z_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_CFz(iMarker_Monitoring);
        break;
      case TOTAL_HEATFLUX:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_HF_Visc(iMarker_Monitoring);
        break;
      case MAXIMUM_HEATFLUX:
        Total_ComboObj+=Weight_ObjFunc*solver->GetSurface_MaxHF_Visc(iMarker_Monitoring);
        break;

    }
  }

  return Total_ComboObj;
  
}
