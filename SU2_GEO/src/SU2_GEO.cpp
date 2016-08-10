/*!
 * \file SU2_GEO.cpp
 * \brief Main file of the Geometry Definition Code (SU2_GEO).
 * \author F. Palacios, T. Economon
 * \version 4.2.0 "Cardinal"
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

#include "../include/SU2_GEO.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
	unsigned short iDV, iFFDBox, iPlane, nPlane, iVar;
	su2double *ObjectiveFunc, *ObjectiveFunc_New, *Gradient, delta_eps, MinPlane, MaxPlane, MinXCoord, MaxXCoord,
  **Plane_P0, **Plane_Normal, Volume = 0.0, Volume_New = 0.0, Volume_Grad = 0.0;
  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  char config_file_name[MAX_STRING_SIZE];
 	char *cstr;
	ofstream Gradient_file, ObjFunc_file;
	int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  /*--- MPI initialization ---*/

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
    
    config_container[iZone] = new CConfig(config_file_name, SU2_GEO, iZone, nZone, 0, VERB_HIGH);
        
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
  
  /*--- Evaluation of the objective function ---*/
  
  if (rank == MASTER_NODE)
		cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;

  /*--- Set the number of sections, and allocate the memory ---*/
  
  if (geometry_container[ZONE_0]->GetnDim() == 2) nPlane = 1;
  else nPlane = config_container[ZONE_0]->GetnSections();

  Xcoord_Airfoil = new vector<su2double>[nPlane];
  Ycoord_Airfoil = new vector<su2double>[nPlane];
  Zcoord_Airfoil = new vector<su2double>[nPlane];
  Variable_Airfoil = new vector<su2double>[nPlane];

  Plane_P0 = new su2double*[nPlane];
  Plane_Normal = new su2double*[nPlane];
  for(iPlane = 0; iPlane < nPlane; iPlane++ ) {
    Plane_P0[iPlane] = new su2double[3];
    Plane_Normal[iPlane] = new su2double[3];
  }
  
  ObjectiveFunc = new su2double[nPlane*20];
  ObjectiveFunc_New = new su2double[nPlane*20];
  Gradient = new su2double[nPlane*20];

  for (iVar = 0; iVar < nPlane*20; iVar++) {
    ObjectiveFunc[iVar] = 0.0;
    ObjectiveFunc_New[iVar] = 0.0;
    Gradient[iVar] = 0.0;
  }
  
  
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
  geometry_container[ZONE_0]->SetCoord_CG();
  
  /*--- Create the dual control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the bound control volume structure." << endl;
  geometry_container[ZONE_0]->SetBoundControlVolume(config_container[ZONE_0], ALLOCATE);
  
  /*--- Compute the surface curvature ---*/
  
  if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
  geometry_container[ZONE_0]->ComputeSurf_Curvature(config_container[ZONE_0]);
  
//  if (rank == MASTER_NODE) cout << "Writing a Tecplot file of the surface curvature." << endl;
//  if (size > 1) SPRINTF (buffer_char, "_%d.plt", rank+1); else SPRINTF (buffer_char, ".plt");
//  strcpy (out_file, "Surface_Curvature"); strcat(out_file, buffer_char); geometry_container[ZONE_0]->SetBoundTecPlot(out_file, true, config);
  
	/*--- Create plane structure ---*/
  
  if (rank == MASTER_NODE) cout << "Set plane structure." << endl;
  if (geometry_container[ZONE_0]->GetnDim() == 2) {
    MinXCoord = -1E6; MaxXCoord = 1E6;
    Plane_Normal[0][0] = 0.0;   Plane_P0[0][0] = 0.0;
    Plane_Normal[0][1] = 1.0;   Plane_P0[0][1] = 0.0;
    Plane_Normal[0][2] = 0.0;   Plane_P0[0][2] = 0.0;
  }
  else if (geometry_container[ZONE_0]->GetnDim() == 3) {
    MinPlane = config_container[ZONE_0]->GetSection_Location(0); MaxPlane = config_container[ZONE_0]->GetSection_Location(1);
    MinXCoord = -1E6; MaxXCoord = 1E6;
    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
      Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
      Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;
      Plane_Normal[iPlane][config_container[ZONE_0]->GetAxis_Orientation()] = 1.0;
      Plane_P0[iPlane][config_container[ZONE_0]->GetAxis_Orientation()] = MinPlane + iPlane*(MaxPlane - MinPlane)/su2double(nPlane-1);
    }
  }

  /*--- Create airfoil section structure ---*/
  
  if (rank == MASTER_NODE) cout << "Set airfoil section structure." << endl;
  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    geometry_container[ZONE_0]->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], MinXCoord, MaxXCoord, NULL,
                                     Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], Variable_Airfoil[iPlane], true, config_container[ZONE_0]);
  }
  
  /*--- Compute the internal volume of a 3D body. ---*/
  
  if (rank == MASTER_NODE)  cout << "Computing the internal volume." << endl;
  if (geometry_container[ZONE_0]->GetnDim() == 3) Volume = geometry_container[ZONE_0]->Compute_Volume(config_container[ZONE_0], true);
  
  if (rank == MASTER_NODE)
    cout << endl <<"-------------------- Objective function evaluation ----------------------" << endl;

  if (rank == MASTER_NODE) {
    
    /*--- Evaluate objective function ---*/
    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      if (Xcoord_Airfoil[iPlane].size() != 0) {
        
        cout << "\nSection " << (iPlane+1) << ". Plane (yCoord): " << Plane_P0[iPlane][1] << "." << endl;
        
        ObjectiveFunc[iPlane]           = geometry_container[ZONE_0]->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[1*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.250000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[2*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.333333, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[3*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.500000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[4*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.666666, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[5*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.750000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[6*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[7*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_AoA(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        ObjectiveFunc[8*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
        
        cout << "Maximum thickness: "   << ObjectiveFunc[iPlane] << "." << endl;
        cout << "1/4 chord thickness: " << ObjectiveFunc[1*nPlane+iPlane] << "." << endl;
        cout << "1/3 chord thickness: " << ObjectiveFunc[2*nPlane+iPlane] << "." << endl;
        cout << "1/2 chord thickness: " << ObjectiveFunc[3*nPlane+iPlane] << "." << endl;
        cout << "2/3 chord thickness: " << ObjectiveFunc[4*nPlane+iPlane] << "." << endl;
        cout << "3/4 chord thickness: " << ObjectiveFunc[5*nPlane+iPlane] << "." << endl;
        cout << "Area: "                << ObjectiveFunc[6*nPlane+iPlane] << "." << endl;
        cout << "Twist angle: "     << ObjectiveFunc[7*nPlane+iPlane] << "." << endl;
        cout << "Chord: "               << ObjectiveFunc[8*nPlane+iPlane] << "." << endl;
        
      }
      
    }
    
    /*--- The value for Volume to the area of the first section
     in case this is a 2D calculation ---*/
    
    if (geometry_container[ZONE_0]->GetnDim() == 2) Volume = ObjectiveFunc[6];
    cout << "\nInternal volume: " << Volume << "." << endl;
    
    /*--- Write the objective function in a external file ---*/
		cstr = new char [config_container[ZONE_0]->GetObjFunc_Value_FileName().size()+1];
		strcpy (cstr, config_container[ZONE_0]->GetObjFunc_Value_FileName().c_str());
		ObjFunc_file.open(cstr, ios::out);
    ObjFunc_file << "TITLE = \"SU2_GEO Simulation\"" << endl;
    
    if (geometry_container[ZONE_0]->GetnDim() == 2) {
      ObjFunc_file << "VARIABLES = \"MAX_THICKNESS\",\"1/4_THICKNESS\",\"1/3_THICKNESS\",\"1/2_THICKNESS\",\"2/3_THICKNESS\",\"3/4_THICKNESS\",\"AREA\",\"AOA\",\"CHORD\",\"VOLUME\"" << endl;
    }
    else if (geometry_container[ZONE_0]->GetnDim() == 3) {
      ObjFunc_file << "VARIABLES = ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"MAX_THICKNESS_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"1/4_THICKNESS_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"1/3_THICKNESS_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"1/2_THICKNESS_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"2/3_THICKNESS_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"3/4_THICKNESS_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"AREA_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"AOA_SEC"<< (iPlane+1) << "\", ";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"CHORD_SEC"<< (iPlane+1) << "\", ";
      ObjFunc_file << "\"VOLUME\"" << endl;
    }
    
    ObjFunc_file << "ZONE T= \"Geometrical variables (value)\"" << endl;
    
    for (iPlane = 0; iPlane < nPlane*9; iPlane++)
      ObjFunc_file << ObjectiveFunc[iPlane] <<", ";
    ObjFunc_file << Volume << endl;

    ObjFunc_file.close();
    
	}
	
	if (config_container[ZONE_0]->GetGeometryMode() == GRADIENT) {
		
		/*--- Definition of the Class for surface deformation ---*/
		surface_movement = new CSurfaceMovement();
    
    /*--- Copy coordinates to the surface structure ---*/
    surface_movement->CopyBoundary(geometry_container[ZONE_0], config_container[ZONE_0]);
		
		/*--- Definition of the FFD deformation class ---*/
		FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
		
		if (rank == MASTER_NODE)
			cout << endl <<"------------- Gradient evaluation using finite differences --------------" << endl;

		/*--- Write the gradient in a external file ---*/
		if (rank == MASTER_NODE) {
			cstr = new char [config_container[ZONE_0]->GetObjFunc_Grad_FileName().size()+1];
			strcpy (cstr, config_container[ZONE_0]->GetObjFunc_Grad_FileName().c_str());
			Gradient_file.open(cstr, ios::out);
		}
		
		for (iDV = 0; iDV < config_container[ZONE_0]->GetnDV(); iDV++) {
			   
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
          
          if (rank == MASTER_NODE) cout << "Read the FFD information from mesh file." << endl;
          
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
          
          switch ( config_container[ZONE_0]->GetDesign_Variable(iDV) ) {
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
      
      /*--- Translation design variable ---*/
      
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
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == CUSTOM and rank==MASTER_NODE)
        cout <<"Custom design variable will be used in external script" << endl;

      /*--- Design variable not implement ---*/
      
      else {
        if (rank==MASTER_NODE)
          cout << "Design Variable not implemented yet" << endl;
      }
      
      /*--- Create airfoil structure ---*/
      for (iPlane = 0; iPlane < nPlane; iPlane++) {
        geometry_container[ZONE_0]->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], MinXCoord, MaxXCoord, NULL,
                                         Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], Variable_Airfoil[iPlane], false, config_container[ZONE_0]);
      }
      
      /*--- Compute the gradient for the volume. In 2D this is just
       the gradient of the area. ---*/
      
      if (geometry_container[ZONE_0]->GetnDim() == 3) Volume_New = geometry_container[ZONE_0]->Compute_Volume(config_container[ZONE_0], false);
      
			/*--- Compute gradient ---*/
			if (rank == MASTER_NODE) {
        
        delta_eps = config_container[ZONE_0]->GetDV_Value(iDV);
        
        if (delta_eps == 0) {
          cout << "The finite difference steps is zero!!" << endl;
          cout << "Press any key to exit..." << endl;
          cin.get();
#ifdef HAVE_MPI
          MPI_Abort(MPI_COMM_WORLD,1);
          MPI_Finalize();
#else
          exit(1);
#endif
        }

        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          
          if (Xcoord_Airfoil[iPlane].size() != 0) {
            
            cout << "\nSection " << (iPlane+1) << ". Plane (yCoord): " << Plane_P0[iPlane][1] << "." << endl;
            
            ObjectiveFunc_New[iPlane] = geometry_container[ZONE_0]->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[iPlane] = (ObjectiveFunc_New[iPlane] - ObjectiveFunc[iPlane]) / delta_eps;
            
            ObjectiveFunc_New[1*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.250000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[1*nPlane + iPlane] = (ObjectiveFunc_New[1*nPlane + iPlane] - ObjectiveFunc[1*nPlane + iPlane]) / delta_eps;
            
            ObjectiveFunc_New[2*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.333333, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[2*nPlane + iPlane] = (ObjectiveFunc_New[2*nPlane + iPlane] - ObjectiveFunc[2*nPlane + iPlane]) / delta_eps;
            
            ObjectiveFunc_New[3*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.500000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[3*nPlane + iPlane] = (ObjectiveFunc_New[3*nPlane + iPlane] - ObjectiveFunc[3*nPlane + iPlane]) / delta_eps;
            
            ObjectiveFunc_New[4*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.666666, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[4*nPlane + iPlane] = (ObjectiveFunc_New[4*nPlane + iPlane] - ObjectiveFunc[4*nPlane + iPlane]) / delta_eps;
            
            ObjectiveFunc_New[5*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, 0.750000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[5*nPlane + iPlane] = (ObjectiveFunc_New[5*nPlane + iPlane] - ObjectiveFunc[5*nPlane + iPlane]) / delta_eps;
            
            ObjectiveFunc_New[6*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[6*nPlane + iPlane] = (ObjectiveFunc_New[6*nPlane + iPlane] - ObjectiveFunc[6*nPlane + iPlane]) / delta_eps;
            
            ObjectiveFunc_New[7*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_AoA(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[7*nPlane + iPlane] = (ObjectiveFunc_New[7*nPlane + iPlane] - ObjectiveFunc[7*nPlane + iPlane]) / delta_eps;
            
            ObjectiveFunc_New[8*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
            Gradient[8*nPlane + iPlane] = (ObjectiveFunc_New[8*nPlane + iPlane] - ObjectiveFunc[8*nPlane + iPlane]) / delta_eps;
            
            cout << "Maximum thickness gradient: "    << Gradient[iPlane] << "." << endl;
            cout << "1/4 chord thickness gradient: "  << Gradient[1*nPlane + iPlane] << "." << endl;
            cout << "1/3 chord thickness gradient: "  << Gradient[2*nPlane + iPlane] << "." << endl;
            cout << "1/2 chord thickness gradient: "  << Gradient[3*nPlane + iPlane] << "." << endl;
            cout << "2/3 chord thickness gradient: "  << Gradient[4*nPlane + iPlane] << "." << endl;
            cout << "3/4 chord thickness gradient: "  << Gradient[5*nPlane + iPlane] << "." << endl;
            cout << "Area gradient: "                 << Gradient[6*nPlane + iPlane] << "." << endl;
            cout << "Twist angle gradient: "      << Gradient[7*nPlane + iPlane] << "." << endl;
            cout << "Chord gradient: "                << Gradient[8*nPlane + iPlane] << "." << endl;
            
          }
          
        }
        
        /*--- Compute the gradient for the volume. In 2D this is just
         the gradient of the area. ---*/
        
        if (geometry_container[ZONE_0]->GetnDim() == 2) Volume_New = ObjectiveFunc_New[6];
        Volume_Grad = (Volume_New - Volume)/ delta_eps;
        cout << "\nInternal volume gradient: " << Volume_Grad << "." << endl;
 				
        if (iDV == 0) {
          Gradient_file << "TITLE = \"SU2_GEO Simulation\"" << endl;

          if (geometry_container[ZONE_0]->GetnDim() == 2) {
            Gradient_file << "VARIABLES = \"DESIGN_VARIABLE\",\"MAX_THICKNESS\",\"1/4_THICKNESS\",\"1/3_THICKNESS\",\"1/2_THICKNESS\",\"2/3_THICKNESS\",\"3/4_THICKNESS\",\"AREA\",\"AOA\",\"CHORD\",\"VOLUME\"" << endl;
          }
          else if (geometry_container[ZONE_0]->GetnDim() == 3) {
            Gradient_file << "VARIABLES = \"DESIGN_VARIABLE\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"MAX_THICKNESS_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"1/4_THICKNESS_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"1/3_THICKNESS_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"1/2_THICKNESS_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"2/3_THICKNESS_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"3/4_THICKNESS_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"AREA_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"AOA_SEC"<< (iPlane+1) << "\", ";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"CHORD_SEC"<< (iPlane+1) << "\", ";
            Gradient_file << "\"VOLUME\"" << endl;
          }
          
          Gradient_file << "ZONE T= \"Geometrical variables (gradient)\"" << endl;
          
        }
        
        Gradient_file << (iDV) <<", ";
        for (iPlane = 0; iPlane < nPlane*9; iPlane++)
          Gradient_file << Gradient[iPlane] <<", ";
        Gradient_file << Volume_Grad << endl;
        
				if (iDV != (config_container[ZONE_0]->GetnDV()-1)) cout <<"-------------------------------------------------------------------------" << endl;
				
			}

		}
		
		if (rank == MASTER_NODE)
			Gradient_file.close();
    
	}
		
  /*--- Deallocate memory ---*/
  
  delete [] Xcoord_Airfoil;
  delete [] Ycoord_Airfoil;
  delete [] Zcoord_Airfoil;
  
  delete [] ObjectiveFunc;
  delete [] ObjectiveFunc_New;
  delete [] Gradient;
  
  for(iPlane = 0; iPlane < nPlane; iPlane++ ) {
    delete Plane_P0[iPlane];
    delete Plane_Normal[iPlane];
  }
  delete [] Plane_P0;
  delete [] Plane_Normal;
  
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
		cout << endl <<"------------------------- Exit Success (SU2_GEO) ------------------------" << endl << endl;

  
  /*--- Finalize MPI parallelization ---*/
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
	return EXIT_SUCCESS;
	
}
