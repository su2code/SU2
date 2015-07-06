/*!
 * \file SU2_DOT.cpp
 * \brief Main file of the Gradient Projection Code (SU2_DOT).
 * \author F. Palacios
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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
	unsigned short iMarker, iDim, iDV, iFFDBox;
	unsigned long iVertex, iPoint;
	su2double delta_eps, my_Gradient, Gradient, *Normal, dS, *VarCoord, Sensitivity,
  dalpha[3], deps[3], dalpha_deps;
	char config_file_name[MAX_STRING_SIZE], *cstr;
	ofstream Gradient_file, Jacobian_file;
	bool *UpdatePoint;
	int rank = MASTER_NODE;
	int size = SINGLE_NODE;

  /*--- MPI initialization, and buffer setting ---*/

#ifdef HAVE_MPI
	SU2_MPI::Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
	
	/*--- Pointer to different structures that will be used throughout the entire code ---*/
  
	CConfig **config_container          = NULL;
	CGeometry **geometry_container      = NULL;
	CSurfaceMovement *surface_movement  = NULL;
  CFreeFormDefBox** FFDBox            = NULL;
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
    
    config_container[iZone] = new CConfig(config_file_name, SU2_DOT, iZone, nZone, 0, VERB_HIGH);
        
    /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/
    
    CGeometry *geometry_aux = NULL;
    
    /*--- All ranks process the grid and call ParMETIS for partitioning ---*/
    
    geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);
    
    /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/
    
    geometry_aux->SetColorGrid_Parallel(config_container[iZone]);
    
    /*--- Allocate the memory of the current domain, and
     divide the grid between the nodes ---*/
    
    geometry_container[iZone] = new CPhysicalGeometry(geometry_aux, config_container[iZone], 1);
    
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
  
  if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation of the interior elements." <<endl;
  geometry_container[ZONE_0]->Check_IntElem_Orientation(config_container[ZONE_0]);
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identify edges and vertices." <<endl;
  geometry_container[ZONE_0]->SetEdges(); geometry_container[ZONE_0]->SetVertex(config_container[ZONE_0]);
  
  /*--- Compute center of gravity ---*/
  
  if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
  geometry_container[ZONE_0]->SetCG();
  
  /*--- Create the dual control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);
  
  /*--- Load the surface sensitivities from file. This is done only
   once: if this is an unsteady problem, a time-average of the surface
   sensitivities at each node is taken within this routine. ---*/
  if (!config_container[ZONE_0]->GetDiscrete_Adjoint()){
    if (rank == MASTER_NODE) cout << "Reading surface sensitivities at each node from file." << endl;
    geometry_container[ZONE_0]->SetBoundSensitivity(config_container[ZONE_0]);
  } else {
    mesh_movement = new CVolumetricMovement(geometry_container[ZONE_0]);
    geometry_container[ZONE_0]->SetSensitivity(config_container[ZONE_0]);

    if (rank == MASTER_NODE) cout << "Setting mesh sensitivity." << endl;
    mesh_movement->SetVolume_Deformation(geometry_container[ZONE_0], config_container[ZONE_0], false, true);
  }
  
  /*--- Boolean controlling points to be updated ---*/
  
	UpdatePoint = new bool[geometry_container[ZONE_0]->GetnPoint()];
	
	/*--- Definition of the Class for surface deformation ---*/
  
	surface_movement = new CSurfaceMovement();
  
  /*--- Copy coordinates to the surface structure ---*/
  
  surface_movement->CopyBoundary(geometry_container[ZONE_0], config_container[ZONE_0]);

	/*--- Definition of the FFD deformation class ---*/
  
	unsigned short nFFDBox = MAX_NUMBER_FFD;
	FFDBox = new CFreeFormDefBox*[nFFDBox];
	
	if (rank == MASTER_NODE) 
		cout << endl <<"---------- Start gradient evaluation using surface sensitivity ----------" << endl;
	
	/*--- Write the gradient in a external file ---*/

	if (rank == MASTER_NODE) {
		cstr = new char [config_container[ZONE_0]->GetObjFunc_Grad_FileName().size()+1];
		strcpy (cstr, config_container[ZONE_0]->GetObjFunc_Grad_FileName().c_str());
		Gradient_file.open(cstr, ios::out);
    
//    /*--- Write an additional file with the geometric Jacobian ---*/
//    /*--- WARNING: This is only for serial calculations!!! ---*/
//    /*--- WARNING: We should generalize this... ---*/
//    if (size == SINGLE_NODE) {
//      Jacobian_file.open("geo_jacobian.csv", ios::out);
//      Jacobian_file.precision(15);
//      
//      /*--- Write the CSV file header ---*/
//      Comma = false;
//      for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
//        if (config_container[ZONE_0]->GetMarker_All_DV(iMarker) == YES) {
//          for (iVertex = 0; iVertex < geometry_container[ZONE_0]->nVertex[iMarker]; iVertex++) {
//            iPoint = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNode();
//            if (!Comma) { Jacobian_file << "\t\"DesignVariable\""; Comma = true;}
//            Jacobian_file  << ", " << "\t\"" << iPoint << "\"";
//          }
//        }
//      }
//      Jacobian_file << endl;
//    }

	}

  /*--- For the discrete adjoint method we use AD to compute the derivatives---*/
  
  if (config_container[ZONE_0]->GetDiscrete_Adjoint()){

    su2double DV_Value = 0.0;

    /*--- Start recording of operations ---*/

    AD::StartRecording();

    /*--- Register design variables as input and set them to zero ---*/

    for (iDV = 0; iDV < config_container[ZONE_0]->GetnDV(); iDV++){
      AD::RegisterInput(DV_Value);

      config_container[ZONE_0]->SetDV_Value(iDV, DV_Value);
    }

    /*--- Call the surface deformation routine ---*/

    surface_movement->SetSurface_Deformation(geometry_container[ZONE_0], config_container[ZONE_0]);

    /*--- Stop the recording --- */

    AD::StopRecording();

    /*--- Initialize the derivatives of the output of the surface deformation routine
     * with the discrete adjoints from the CFD solution ---*/

    for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
      if (config_container[ZONE_0]->GetMarker_All_DV(iMarker) == YES) {
        for (iVertex = 0; iVertex < geometry_container[ZONE_0]->nVertex[iMarker]; iVertex++) {
          iPoint = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNode();
          for (iDim = 0; iDim < geometry_container[ZONE_0]->GetnDim(); iDim++){
            SU2_TYPE::SetDerivative(geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetVarCoord()[iDim],
                                  SU2_TYPE::GetPrimary(geometry_container[ZONE_0]->GetSensitivity(iPoint, iDim)));

          }
        }
      }
    }

    /*--- Compute derivatives ---*/

    AD::ComputeAdjoint();

    for (iDV = 0; iDV  < config_container[ZONE_0]->GetnDV(); iDV++){
      my_Gradient = SU2_TYPE::GetDerivative(DV_Value);
#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&my_Gradient, &Gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      Gradient = my_Gradient;
#endif
      /*--- Print result to screen and file ---*/

      if (rank == MASTER_NODE) {
        switch (config_container[ZONE_0]->GetKind_ObjFunc()) {
          case LIFT_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Lift coeff. grad. using disc. adj." << endl;
            cout << "Lift coefficient gradient: "<< Gradient << "." << endl; break;
          case DRAG_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Drag coeff. grad. using disc. adj." << endl;
            cout << "Drag coefficient gradient: "<< Gradient << "." << endl; break;
          case SIDEFORCE_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Sideforce coeff. grad. using disc. adj." << endl;
            cout << "Sideforce coefficient gradient: "<< Gradient << "." << endl; break;
          case INVERSE_DESIGN_PRESSURE :
            if (iDV == 0) Gradient_file << "Pressure inverse design using disc. adj."<< endl;
            cout << "Pressure inverse design gradient: "<< Gradient << "." << endl; break;
          case INVERSE_DESIGN_HEATFLUX :
            if (iDV == 0) Gradient_file << "Heat inverse design using disc. adj."<< endl;
            cout << "Heat flux inverse design gradient: "<< Gradient << "." << endl; break;
          case TOTAL_HEATFLUX :
            if (iDV == 0) Gradient_file << "Integrated surface heat flux. using disc. adj."<< endl;
            cout << "Total heat flux gradient: "<< Gradient << "." << endl; break;
          case MAXIMUM_HEATFLUX :
            if (iDV == 0) Gradient_file << "Integrated surface heat flux. using disc. adj."<< endl;
            cout << "Maximum heat flux gradient: "<< Gradient << "." << endl; break;
          case MOMENT_X_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Moment x coeff. grad. using disc. adj." << endl;
            cout << "Moment x coefficient gradient: "<< Gradient << "." << endl; break;
          case MOMENT_Y_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Moment y coeff. grad. using disc. adj." << endl;
            cout << "Moment y coefficient gradient: "<< Gradient << "." << endl; break;
          case MOMENT_Z_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Moment z coeff. grad. using disc. adj." << endl;
            cout << "Moment z coefficient gradient: "<< Gradient << "." << endl; break;
          case EFFICIENCY :
            if (iDV == 0) Gradient_file << "Efficiency coeff. grad. using disc. adj." << endl;
            cout << "Efficiency coefficient gradient: "<< Gradient << "." << endl; break;
          case EQUIVALENT_AREA :
            if (iDV == 0) Gradient_file << "Equivalent area coeff. grad. using disc. adj." << endl;
            cout << "Equivalent Area coefficient gradient: "<< Gradient << "." << endl; break;
          case NEARFIELD_PRESSURE :
            if (iDV == 0) Gradient_file << "Near-field pressure coeff. grad. using disc. adj." << endl;
            cout << "Near-field pressure coefficient gradient: "<< Gradient << "." << endl; break;
          case FORCE_X_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Force x coeff. grad. using disc. adj." << endl;
            cout << "Force x coefficient gradient: "<< Gradient << "." << endl; break;
          case FORCE_Y_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Force y coeff. grad. using disc. adj." << endl;
            cout << "Force y coefficient gradient: "<< Gradient << "." << endl; break;
          case FORCE_Z_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Force z coeff. grad. using disc. adj." << endl;
            cout << "Force z coefficient gradient: "<< Gradient << "." << endl; break;
          case THRUST_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Thrust coeff. grad. using disc. adj."<< endl;
            cout << "Thrust coefficient gradient: "<< Gradient << "." << endl; break;
          case TORQUE_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Torque coeff. grad. using disc. adj."<< endl;
            cout << "Torque coefficient gradient: "<< Gradient << "." << endl; break;
          case FIGURE_OF_MERIT :
            if (iDV == 0) Gradient_file << "Rotor Figure of Merit grad. using disc. adj."<< endl;
            cout << "Rotor Figure of Merit gradient: "<< Gradient << "." << endl; break;
          case FREE_SURFACE :
            if (iDV == 0) Gradient_file << "Free-Surface grad. using disc. adj."<< endl;
            cout << "Free-surface gradient: "<< Gradient << "." << endl; break;
          case MASS_FLOW_RATE :
            if (iDV == 0) Gradient_file << "Mass flow rate grad. using disc. adj."<< endl;
            cout << "Mass flow rate gradient: "<< Gradient << "." << endl; break;
          case AVG_OUTLET_PRESSURE :
            if (iDV == 0) Gradient_file << "Average outlet presure grad. using disc. adj."<< endl;
            cout << "Average outlet pressure gradient: "<< Gradient << "." << endl; break;
          case AVG_TOTAL_PRESSURE :
            if (iDV == 0) Gradient_file << "Average total presure grad. using disc. adj."<< endl;
            cout << "Average total pressure gradient: "<< Gradient << "." << endl; break;

        }

        Gradient_file << Gradient << endl;

        cout <<"-------------------------------------------------------------------------" << endl;

      }

    }

    /*--- For continuous adjoint we use Finite differences ---*/

  }else{


    for (iDV = 0; iDV < config_container[ZONE_0]->GetnDV(); iDV++) {

      if (size == SINGLE_NODE) { Jacobian_file << iDV; }

      /*--- Free Form deformation based ---*/

      if ((config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CAMBER_2D) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_THICKNESS_2D) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_DIHEDRAL_ANGLE) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_TWIST_ANGLE) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_ROTATION) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CAMBER) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_THICKNESS) ) {

        /*--- Read the FFD information in the first iteration ---*/

        if (iDV == 0) {

          if (rank == MASTER_NODE)
            cout << "Read the FFD information from mesh file." << endl;

          /*--- Read the FFD information from the grid file ---*/

          surface_movement->ReadFFDInfo(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox, config_container[ZONE_0]->GetMesh_FileName());

          /*--- If the FFDBox was not defined in the input file ---*/
          if (!surface_movement->GetFFDBoxDefinition() && (rank == MASTER_NODE)) {
            cout << "The input grid doesn't have the entire FFD information!" << endl;
            cout << "Press any key to exit..." << endl;
            cin.get();
          }

          for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

            if (rank == MASTER_NODE)
              cout << "Check the FFD box intersections with the solid surfaces." << endl;

            surface_movement->CheckFFDIntersections(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iFFDBox);

          }

          if (rank == MASTER_NODE)
            cout <<"-------------------------------------------------------------------------" << endl;

        }

        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 3D deformation of the surface." << endl;
        }

        /*--- Apply the control point change ---*/

        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

          /*--- Reset FFD box ---*/

          switch (config_container[ZONE_0]->GetDesign_Variable(iDV) ) {
            case FFD_CONTROL_POINT_2D : surface_movement->SetFFDCPChange_2D(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_CAMBER_2D :        surface_movement->SetFFDCamber_2D(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_THICKNESS_2D :     surface_movement->SetFFDThickness_2D(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_CONTROL_POINT :    surface_movement->SetFFDCPChange(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_DIHEDRAL_ANGLE :   surface_movement->SetFFDDihedralAngle(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_TWIST_ANGLE :      surface_movement->SetFFDTwistAngle(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_ROTATION :         surface_movement->SetFFDRotation(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_CAMBER :           surface_movement->SetFFDCamber(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_THICKNESS :        surface_movement->SetFFDThickness(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
            case FFD_CONTROL_SURFACE :  surface_movement->SetFFDControl_Surface(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iDV, true); break;
          }

          /*--- Recompute cartesian coordinates using the new control points position ---*/

          surface_movement->SetCartesianCoord(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iFFDBox);

        }

      }

      /*--- Hicks Henne design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == HICKS_HENNE) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        surface_movement->SetHicksHenne(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Displacement design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == TRANSLATION) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        surface_movement->SetTranslation(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Scale design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == SCALE) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        surface_movement->SetScale(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Rotation design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == ROTATION) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        surface_movement->SetRotation(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- NACA_4Digits design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == NACA_4DIGITS) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        surface_movement->SetNACA_4Digits(geometry_container[ZONE_0], config_container[ZONE_0]);
      }

      /*--- Parabolic design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == PARABOLIC) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        surface_movement->SetParabolic(geometry_container[ZONE_0], config_container[ZONE_0]);
      }

      /*--- Design variable not implement ---*/

      else { cout << "Design Variable not implement yet" << endl; }


      /*--- Continuous adjoint gradient computation ---*/
      if (rank == MASTER_NODE)
        cout << "Evaluate functional gradient using the continuous adjoint strategy." << endl;

      /*--- Load the delta change in the design variable (finite difference step). ---*/

      delta_eps = config_container[ZONE_0]->GetDV_Value(iDV);
      my_Gradient = 0.0; Gradient = 0.0;
      
      /*--- Reset update points ---*/

      for (iPoint = 0; iPoint < geometry_container[ZONE_0]->GetnPoint(); iPoint++)
        UpdatePoint[iPoint] = true;
      
      for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
        if (config_container[ZONE_0]->GetMarker_All_DV(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry_container[ZONE_0]->nVertex[iMarker]; iVertex++) {

            iPoint = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNode();
            if ((iPoint < geometry_container[ZONE_0]->GetnPointDomain()) && UpdatePoint[iPoint]) {

              Normal = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetNormal();
              VarCoord = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetVarCoord();
              Sensitivity = geometry_container[ZONE_0]->vertex[iMarker][iVertex]->GetAuxVar();

              dS = 0.0;
              for (iDim = 0; iDim < geometry_container[ZONE_0]->GetnDim(); iDim++) {
                dS += Normal[iDim]*Normal[iDim];
                deps[iDim] = VarCoord[iDim] / delta_eps;
              }
              dS = sqrt(dS);

              dalpha_deps = 0.0;
              for (iDim = 0; iDim < geometry_container[ZONE_0]->GetnDim(); iDim++) {
                dalpha[iDim] = Normal[iDim] / dS;
                dalpha_deps -= dalpha[iDim]*deps[iDim];
              }

              /*--- Store the geometric sensitivity for this DV (rows) & this node (column) ---*/
              
              if (size == SINGLE_NODE) {
                Jacobian_file  << ", " << dalpha_deps;
              }
              
              my_Gradient += Sensitivity*dalpha_deps;
              UpdatePoint[iPoint] = false;
            }
          }
        }
      }

#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&my_Gradient, &Gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      Gradient = my_Gradient;
#endif

      if (rank == MASTER_NODE) {
        switch (config_container[ZONE_0]->GetKind_ObjFunc()) {
          case LIFT_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Lift coeff. grad. using cont. adj." << endl;
            cout << "Lift coefficient gradient: "<< Gradient << "." << endl; break;
          case DRAG_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Drag coeff. grad. using cont. adj." << endl;
            cout << "Drag coefficient gradient: "<< Gradient << "." << endl; break;
          case SIDEFORCE_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Sideforce coeff. grad. using cont. adj." << endl;
            cout << "Sideforce coefficient gradient: "<< Gradient << "." << endl; break;
          case INVERSE_DESIGN_PRESSURE :
            if (iDV == 0) Gradient_file << "Pressure inverse design using cont. adj."<< endl;
            cout << "Pressure inverse design gradient: "<< Gradient << "." << endl; break;
          case INVERSE_DESIGN_HEATFLUX :
            if (iDV == 0) Gradient_file << "Heat inverse design using cont. adj."<< endl;
            cout << "Heat flux inverse design gradient: "<< Gradient << "." << endl; break;
          case TOTAL_HEATFLUX :
            if (iDV == 0) Gradient_file << "Integrated surface heat flux. using cont. adj."<< endl;
            cout << "Total heat flux gradient: "<< Gradient << "." << endl; break;
          case MAXIMUM_HEATFLUX :
            if (iDV == 0) Gradient_file << "Integrated surface heat flux. using cont. adj."<< endl;
            cout << "Maximum heat flux gradient: "<< Gradient << "." << endl; break;
          case MOMENT_X_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Moment x coeff. grad. using cont. adj." << endl;
            cout << "Moment x coefficient gradient: "<< Gradient << "." << endl; break;
          case MOMENT_Y_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Moment y coeff. grad. using cont. adj." << endl;
            cout << "Moment y coefficient gradient: "<< Gradient << "." << endl; break;
          case MOMENT_Z_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Moment z coeff. grad. using cont. adj." << endl;
            cout << "Moment z coefficient gradient: "<< Gradient << "." << endl; break;
          case EFFICIENCY :
            if (iDV == 0) Gradient_file << "Efficiency coeff. grad. using cont. adj." << endl;
            cout << "Efficiency coefficient gradient: "<< Gradient << "." << endl; break;
          case EQUIVALENT_AREA :
            if (iDV == 0) Gradient_file << "Equivalent area coeff. grad. using cont. adj." << endl;
            cout << "Equivalent Area coefficient gradient: "<< Gradient << "." << endl; break;
          case NEARFIELD_PRESSURE :
            if (iDV == 0) Gradient_file << "Near-field pressure coeff. grad. using cont. adj." << endl;
            cout << "Near-field pressure coefficient gradient: "<< Gradient << "." << endl; break;
          case FORCE_X_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Force x coeff. grad. using cont. adj." << endl;
            cout << "Force x coefficient gradient: "<< Gradient << "." << endl; break;
          case FORCE_Y_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Force y coeff. grad. using cont. adj." << endl;
            cout << "Force y coefficient gradient: "<< Gradient << "." << endl; break;
          case FORCE_Z_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Force z coeff. grad. using cont. adj." << endl;
            cout << "Force z coefficient gradient: "<< Gradient << "." << endl; break;
          case THRUST_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Thrust coeff. grad. using cont. adj."<< endl;
            cout << "Thrust coefficient gradient: "<< Gradient << "." << endl; break;
          case TORQUE_COEFFICIENT :
            if (iDV == 0) Gradient_file << "Torque coeff. grad. using cont. adj."<< endl;
            cout << "Torque coefficient gradient: "<< Gradient << "." << endl; break;
          case FIGURE_OF_MERIT :
            if (iDV == 0) Gradient_file << "Rotor Figure of Merit grad. using cont. adj."<< endl;
            cout << "Rotor Figure of Merit gradient: "<< Gradient << "." << endl; break;
          case FREE_SURFACE :
            if (iDV == 0) Gradient_file << "Free-Surface grad. using cont. adj."<< endl;
            cout << "Free-surface gradient: "<< Gradient << "." << endl; break;
          case MASS_FLOW_RATE :
            if (iDV == 0) Gradient_file << "Mass flow rate grad. using cont. adj."<< endl;
            cout << "Mass flow rate gradient: "<< Gradient << "." << endl; break;
          case AVG_OUTLET_PRESSURE :
            if (iDV == 0) Gradient_file << "Average outlet presure grad. using cont. adj."<< endl;
            cout << "Average outlet pressure gradient: "<< Gradient << "." << endl; break;
          case AVG_TOTAL_PRESSURE :
            if (iDV == 0) Gradient_file << "Average total presure grad. using cont. adj."<< endl;
            cout << "Average total pressure gradient: "<< Gradient << "." << endl; break;

        }

        Gradient_file << Gradient << endl;

        cout <<"-------------------------------------------------------------------------" << endl;

      }

      /*--- End the line for the current DV in the geometric Jacobian file ---*/

      if (size == SINGLE_NODE) Jacobian_file << endl;
    }
  }
	if (rank == MASTER_NODE)
		Gradient_file.close();
  
  if (size == SINGLE_NODE)
    Jacobian_file.close();
	
	delete [] UpdatePoint;
	
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
