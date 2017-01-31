/*!
 * \file SU2_GEO.cpp
 * \brief Main file of the Geometry Definition Code (SU2_GEO).
 * \author F. Palacios, T. Economon
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

#include "../include/SU2_GEO.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  
  unsigned short iZone, nZone = SINGLE_ZONE;
  su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  unsigned short iDV, iFFDBox, iPlane, nPlane, iVar;
  su2double FanRadius_Diff_Grad, *ObjectiveFunc, *ObjectiveFunc_New, *Gradient, delta_eps, MinPlane, MaxPlane, MinXCoord, MaxXCoord,
  **Plane_P0, **Plane_Normal,
  
  Wing_Volume = 0.0, Wing_MinMaxThickness = 0.0, Wing_MaxChord = 0.0, Wing_MinToC = 0.0, Wing_MaxTwist = 0.0, Wing_MaxCurvature = 0.0, Wing_MaxDihedral = 0.0,
  Wing_Volume_New = 0.0,Wing_MinMaxThickness_New = 0.0, Wing_MaxChord_New = 0.0, Wing_MinToC_New = 0.0, Wing_MaxTwist_New = 0.0, Wing_MaxCurvature_New = 0.0, Wing_MaxDihedral_New = 0.0,
  Wing_Volume_Grad = 0.0, Wing_MinMaxThickness_Grad = 0.0, Wing_MaxChord_Grad = 0.0, Wing_MinToC_Grad = 0.0, Wing_MaxTwist_Grad = 0.0, Wing_MaxCurvature_Grad = 0.0, Wing_MaxDihedral_Grad = 0.0;
  
  
  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  vector<su2double> Xcoord_Fan, Ycoord_Fan, Zcoord_Fan;
  char config_file_name[MAX_STRING_SIZE];
 	char *cstr;
  bool Local_MoveSurface, MoveSurface = false;
	ofstream Gradient_file, ObjFunc_file;
	int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  /*--- MPI initialization ---*/

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
  
  /*--- Evaluation of the objective function ---*/
  
  if (rank == MASTER_NODE)
		cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;

  /*--- Set the number of sections, and allocate the memory ---*/
  
  if (geometry_container[ZONE_0]->GetnDim() == 2) nPlane = 1;
  else nPlane = config_container[ZONE_0]->GetnLocationStations();

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
  
  /*--- Create plane structure ---*/
  
  if (rank == MASTER_NODE) cout << "Set plane structure." << endl;
  if (geometry_container[ZONE_0]->GetnDim() == 2) {
    MinXCoord = -1E6; MaxXCoord = 1E6;
    Plane_Normal[0][0] = 0.0;   Plane_P0[0][0] = 0.0;
    Plane_Normal[0][1] = 1.0;   Plane_P0[0][1] = 0.0;
    Plane_Normal[0][2] = 0.0;   Plane_P0[0][2] = 0.0;
  }
  else if (geometry_container[ZONE_0]->GetnDim() == 3) {
    
    MinPlane = config_container[ZONE_0]->GetSection_WingBounds(0); MaxPlane = config_container[ZONE_0]->GetSection_WingBounds(1);
    MinXCoord = -1E6; MaxXCoord = 1E6;
        
      for (iPlane = 0; iPlane < nPlane; iPlane++) {
        Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
        Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
        Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;
        Plane_Normal[iPlane][config_container[ZONE_0]->GetAxis_Stations()] = 1.0;
        Plane_P0[iPlane][config_container[ZONE_0]->GetAxis_Stations()] = config_container[ZONE_0]->GetLocationStations(iPlane);
      }
  }
  
  /*--- Compute the wing and fan description (only 3D). ---*/
  
  if (geometry_container[ZONE_0]->GetnDim() == 3) {
    
    if (rank == MASTER_NODE)  cout << "Computing the wing continuous description." << endl << endl;
    
    geometry_container[ZONE_0]->Compute_Wing(config_container[ZONE_0], true,
                                             Wing_Volume, Wing_MinMaxThickness, Wing_MaxChord, Wing_MinToC,
                                             Wing_MaxTwist, Wing_MaxCurvature, Wing_MaxDihedral);
    
    /*--- Screen output for the wing definition ---*/
    
    if (rank == MASTER_NODE) {
      if (config_container[ZONE_0]->GetSystemMeasurements() == US) cout << "Wing volume: "    << Wing_Volume << " in^3. ";
      else cout << "Wing volume: "    << Wing_Volume << " m^3. ";
      if (config_container[ZONE_0]->GetSystemMeasurements() == US) cout << "Wing min. max. thickness: "  << Wing_MinMaxThickness << " in. ";
      else cout << "Wing min. max. thickness: "  << Wing_MinMaxThickness << " m. ";
      if (config_container[ZONE_0]->GetSystemMeasurements() == US) cout << "Wing max. chord: "  << Wing_MaxChord << " in." << endl;
      else cout << "Wing max. chord: "  << Wing_MaxChord << " m." << endl;
      cout << "Wing min. ToC: "  << Wing_MinToC*100 << "%. ";
      cout << "Wing max. twist: "  << Wing_MaxTwist << " deg. ";
      if (config_container[ZONE_0]->GetSystemMeasurements() == US) cout << "Wing max. curvature: "  << Wing_MaxCurvature << " 1/in. " << endl;
      else cout << "Wing max. curvature: "  << Wing_MaxCurvature << " 1/m. " << endl;
      cout << "Wing max. dihedral: "  << Wing_MaxDihedral << " deg." << endl;
    }
    
  }

  /*--- Create airfoil section structure ---*/
  
  if (rank == MASTER_NODE) cout << "Set airfoil section structure." << endl;
  
  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    geometry_container[ZONE_0]->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], MinXCoord, MaxXCoord, NULL,
                                     Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], Variable_Airfoil[iPlane], true, config_container[ZONE_0]);
  }
  
  if (rank == MASTER_NODE)
    cout << endl <<"-------------------- Objective function evaluation ----------------------" << endl;

  if (rank == MASTER_NODE) {
    
    /*--- Evaluate objective function ---*/
    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      if (Xcoord_Airfoil[iPlane].size() != 0) {
        
        cout << "\nStation " << (iPlane+1) << ". Plane (yCoord): " << Plane_P0[iPlane][1] << "." << endl;
        
        ObjectiveFunc[iPlane]           = geometry_container[ZONE_0]->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                           config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[1*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                        0.250000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[2*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                        0.333333, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[3*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                        0.500000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[4*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                        0.666666, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[5*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                        0.750000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[6*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                   config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[7*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Twist(Plane_P0[iPlane], Plane_Normal[iPlane],
                                                                                    iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        ObjectiveFunc[8*nPlane+iPlane]  = geometry_container[ZONE_0]->Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane],
                                                                                    iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
        
        cout << "Max. thickness: "   << ObjectiveFunc[iPlane] << ", ";
        cout << "1/3C thickness: " << ObjectiveFunc[2*nPlane+iPlane] << ", ";
        cout << "2/3C thickness: " << ObjectiveFunc[4*nPlane+iPlane] << endl;
        cout << "1/4C thickness: " << ObjectiveFunc[1*nPlane+iPlane] << ", ";
        cout << "1/2C thickness: " << ObjectiveFunc[3*nPlane+iPlane] << ", ";
        cout << "3/4C thickness: " << ObjectiveFunc[5*nPlane+iPlane] << endl;
        cout << "Area: "                << ObjectiveFunc[6*nPlane+iPlane] << ", ";
        cout << "Twist angle: "     << ObjectiveFunc[7*nPlane+iPlane] << ", ";
        cout << "Chord: "               << ObjectiveFunc[8*nPlane+iPlane] << endl;
        
      }
      
    }
    
    /*--- Write the objective function in a external file ---*/
    
    cstr = new char [config_container[ZONE_0]->GetObjFunc_Value_FileName().size()+1];
    strcpy (cstr, config_container[ZONE_0]->GetObjFunc_Value_FileName().c_str());
    ObjFunc_file.open(cstr, ios::out);
    ObjFunc_file << "TITLE = \"SU2_GEO Evaluation\"" << endl;
    
    if (geometry_container[ZONE_0]->GetnDim() == 2) {
      ObjFunc_file << "VARIABLES = \"MAX_THICKNESS\",\"1/4_THICKNESS\",\"1/3_THICKNESS\",\"1/2_THICKNESS\",\"2/3_THICKNESS\",\"3/4_THICKNESS\",\"AREA\",\"AOA\",\"CHORD\"";
    }
    else if (geometry_container[ZONE_0]->GetnDim() == 3) {
      ObjFunc_file << "VARIABLES = ";
      ObjFunc_file << "\"WING_VOLUME\",\"WING_MIN_MAXTHICKNESS\",\"WING_MAX_CHORD\",\"WING_MIN_TOC\",\"WING_MAX_TWIST\",\"WING_MAX_CURVATURE\",\"WING_MAX_DIHEDRAL\",";
      
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"MAX_THICKNESS_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"1/4_THICKNESS_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"1/3_THICKNESS_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"1/2_THICKNESS_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"2/3_THICKNESS_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"3/4_THICKNESS_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"AREA_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) ObjFunc_file << "\"AOA_SEC"<< (iPlane+1) << "\",";
      for (iPlane = 0; iPlane < nPlane; iPlane++) {
        ObjFunc_file << "\"CHORD_SEC"<< (iPlane+1) << "\"";
        if (iPlane != nPlane-1) ObjFunc_file << ",";
      }
      
    }
    
    ObjFunc_file << "\nZONE T= \"Geometrical variables (value)\"" << endl;
    
    if (geometry_container[ZONE_0]->GetnDim() == 3) {
      ObjFunc_file << Wing_Volume <<", "<< Wing_MinMaxThickness <<", "<< Wing_MaxChord
      <<", "<< Wing_MinToC <<", "<< Wing_MaxTwist <<", "<< Wing_MaxCurvature
      <<", "<< Wing_MaxDihedral <<", ";
    }
    for (iPlane = 0; iPlane < nPlane*9; iPlane++) {
      ObjFunc_file << ObjectiveFunc[iPlane];
      if (iPlane != (nPlane*9)-1) ObjFunc_file <<", ";
    }
    
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
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_TWIST_2D) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_NACELLE) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_GULL) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_TWIST) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_ROTATION) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_CAMBER) ||
          (config_container[ZONE_0]->GetDesign_Variable(iDV) == FFD_THICKNESS) ) {
        
        /*--- Read the FFD information in the first iteration ---*/
        
        if (iDV == 0) {
          
          if (rank == MASTER_NODE) cout << "Read the FFD information from mesh file." << endl;
          
          /*--- Read the FFD information from the grid file ---*/
          
          surface_movement->ReadFFDInfo(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox, config_container[ZONE_0]->GetMesh_FileName());
          
          /*--- Modify the control points for polar based computations ---*/
          
          if (config_container[ZONE_0]->GetFFD_CoordSystem() == CYLINDRICAL) {
            for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCart2Cyl_ControlPoints(config_container[ZONE_0]);
            }
          }
          else if (config_container[ZONE_0]->GetFFD_CoordSystem() == SPHERICAL) {
            for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCart2Sphe_ControlPoints(config_container[ZONE_0]);
            }
          }

          /*--- If the FFDBox was not defined in the input file ---*/
          
          if (!surface_movement->GetFFDBoxDefinition() && (rank == MASTER_NODE)) {
            cout << "The input grid doesn't have the entire FFD information!" << endl;
            cout << "Press any key to exit..." << endl;
            cin.get();
          }
          
          for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
            
            if (rank == MASTER_NODE) cout << "Checking FFD box dimension." << endl;
            surface_movement->CheckFFDDimension(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iFFDBox);

            
            if (rank == MASTER_NODE) cout << "Check the FFD box intersections with the solid surfaces." << endl;
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
        
        MoveSurface = false;

        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {
          
          switch ( config_container[ZONE_0]->GetDesign_Variable(iDV) ) {
            case FFD_CONTROL_POINT_2D : Local_MoveSurface = surface_movement->SetFFDCPChange_2D(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_CAMBER_2D :        Local_MoveSurface = surface_movement->SetFFDCamber_2D(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_THICKNESS_2D :     Local_MoveSurface = surface_movement->SetFFDThickness_2D(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_TWIST_2D :         Local_MoveSurface = surface_movement->SetFFDTwist_2D(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_CONTROL_POINT :    Local_MoveSurface = surface_movement->SetFFDCPChange(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_NACELLE :    Local_MoveSurface = surface_movement->SetFFDNacelle(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_GULL :    Local_MoveSurface = surface_movement->SetFFDGull(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_TWIST :            Local_MoveSurface = surface_movement->SetFFDTwist(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_ROTATION :         Local_MoveSurface = surface_movement->SetFFDRotation(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_CAMBER :           Local_MoveSurface = surface_movement->SetFFDCamber(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_THICKNESS :        Local_MoveSurface = surface_movement->SetFFDThickness(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
            case FFD_CONTROL_SURFACE :  Local_MoveSurface = surface_movement->SetFFDControl_Surface(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], FFDBox, iDV, true); break;
          }
          
          /*--- Recompute cartesian coordinates using the new control points position ---*/
          
          if (Local_MoveSurface) {
            MoveSurface = true;
            surface_movement->SetCartesianCoord(geometry_container[ZONE_0], config_container[ZONE_0], FFDBox[iFFDBox], iFFDBox, true);
          }
          
        }
        
      }
      
      /*--- Hicks Henne design variable ---*/
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == HICKS_HENNE) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetHicksHenne(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- Surface bump design variable ---*/

      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == SURFACE_BUMP) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetSurface_Bump(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }

      /*--- CST design variable ---*/
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == CST) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetCST(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }
      
      /*--- Translation design variable ---*/
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == TRANSLATION) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetTranslation(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }
      
      /*--- Scale design variable ---*/
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == SCALE) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetScale(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }
      
      /*--- Rotation design variable ---*/
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == ROTATION) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetRotation(geometry_container[ZONE_0], config_container[ZONE_0], iDV, true);
      }
      
      /*--- NACA_4Digits design variable ---*/
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == NACA_4DIGITS) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetNACA_4Digits(geometry_container[ZONE_0], config_container[ZONE_0]);
      }
      
      /*--- Parabolic design variable ---*/
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == PARABOLIC) {
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 2D deformation of the surface." << endl;
        }
        MoveSurface = true;
        surface_movement->SetParabolic(geometry_container[ZONE_0], config_container[ZONE_0]);
      }
      
      else if (config_container[ZONE_0]->GetDesign_Variable(iDV) == CUSTOM && rank == MASTER_NODE) {
        cout <<"Custom design variable will be used in external script" << endl;
      }

      /*--- Design variable not implement ---*/
      
      else {
        if (rank==MASTER_NODE)
          cout << "Design Variable not implemented yet" << endl;
      }
      
      if (MoveSurface) {
        
        /*--- Compute the gradient for the volume. In 2D this is just
         the gradient of the area. ---*/
        
        if (geometry_container[ZONE_0]->GetnDim() == 3) {
          
          geometry_container[ZONE_0]->Compute_Wing(config_container[ZONE_0], false,
                                                   Wing_Volume_New, Wing_MinMaxThickness_New, Wing_MaxChord_New, Wing_MinToC_New,
                                                   Wing_MaxTwist_New, Wing_MaxCurvature_New, Wing_MaxDihedral_New);
          
        }
        
        /*--- Create airfoil structure ---*/
        
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          geometry_container[ZONE_0]->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], MinXCoord, MaxXCoord, NULL,
                                                             Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], Variable_Airfoil[iPlane], false, config_container[ZONE_0]);
        }
        
      }
      
      /*--- Compute gradient ---*/
      
      if (rank == MASTER_NODE) {
        
        delta_eps = config_container[ZONE_0]->GetDV_Value(iDV);
        
        if (delta_eps == 0) {
          cout << "The finite difference steps is zero!!" << endl;
          cout << "Press any key to exit..." << endl;
          cin.get();
#ifdef HAVE_MPI
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Abort(MPI_COMM_WORLD,1);
          MPI_Finalize();
#else
          exit(EXIT_FAILURE);
#endif
        }

        if (MoveSurface) {
          
          Wing_Volume_Grad = (Wing_Volume_New - Wing_Volume) / delta_eps;
          Wing_MinMaxThickness_Grad = (Wing_MinMaxThickness_New - Wing_MinMaxThickness) / delta_eps;
          Wing_MaxChord_Grad = (Wing_MaxChord_New - Wing_MaxChord) / delta_eps;
          Wing_MinToC_Grad = (Wing_MinToC_New - Wing_MinToC) / delta_eps;
          Wing_MaxTwist_Grad = (Wing_MaxTwist_New - Wing_MaxTwist) / delta_eps;
          Wing_MaxCurvature_Grad = (Wing_MaxCurvature_New - Wing_MaxCurvature) / delta_eps;
          Wing_MaxDihedral_Grad = (Wing_MaxDihedral_New - Wing_MaxDihedral) / delta_eps;
          
          for (iPlane = 0; iPlane < nPlane; iPlane++) {
            if (Xcoord_Airfoil[iPlane].size() != 0) {
              ObjectiveFunc_New[iPlane] = geometry_container[ZONE_0]->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                           config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[iPlane] = (ObjectiveFunc_New[iPlane] - ObjectiveFunc[iPlane]) / delta_eps;
              ObjectiveFunc_New[1*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                                   0.250000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[1*nPlane + iPlane] = (ObjectiveFunc_New[1*nPlane + iPlane] - ObjectiveFunc[1*nPlane + iPlane]) / delta_eps;
              ObjectiveFunc_New[2*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                                   0.333333, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[2*nPlane + iPlane] = (ObjectiveFunc_New[2*nPlane + iPlane] - ObjectiveFunc[2*nPlane + iPlane]) / delta_eps;
              ObjectiveFunc_New[3*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                                   0.500000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[3*nPlane + iPlane] = (ObjectiveFunc_New[3*nPlane + iPlane] - ObjectiveFunc[3*nPlane + iPlane]) / delta_eps;
              ObjectiveFunc_New[4*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                                   0.666666, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[4*nPlane + iPlane] = (ObjectiveFunc_New[4*nPlane + iPlane] - ObjectiveFunc[4*nPlane + iPlane]) / delta_eps;
              ObjectiveFunc_New[5*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Thickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                                   0.750000, config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[5*nPlane + iPlane] = (ObjectiveFunc_New[5*nPlane + iPlane] - ObjectiveFunc[5*nPlane + iPlane]) / delta_eps;
              ObjectiveFunc_New[6*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                              config_container[ZONE_0], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[6*nPlane + iPlane] = (ObjectiveFunc_New[6*nPlane + iPlane] - ObjectiveFunc[6*nPlane + iPlane]) / delta_eps;
              ObjectiveFunc_New[7*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Twist(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                               Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[7*nPlane + iPlane] = (ObjectiveFunc_New[7*nPlane + iPlane] - ObjectiveFunc[7*nPlane + iPlane]) / delta_eps;
              ObjectiveFunc_New[8*nPlane + iPlane] = geometry_container[ZONE_0]->Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                                                                                               Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);
              Gradient[8*nPlane + iPlane] = (ObjectiveFunc_New[8*nPlane + iPlane] - ObjectiveFunc[8*nPlane + iPlane]) / delta_eps;
            }
          }
          
        }
        
        else {
          
          Wing_Volume_Grad = 0.0;
          Wing_MinMaxThickness_Grad = 0.0;
          Wing_MaxChord_Grad = 0.0;
          Wing_MinToC_Grad = 0.0;
          Wing_MaxTwist_Grad = 0.0;
          Wing_MaxCurvature_Grad = 0.0;
          Wing_MaxDihedral_Grad = 0.0;
          FanRadius_Diff_Grad = 0.0;
          
          for (iPlane = 0; iPlane < nPlane; iPlane++) {
            Gradient[iPlane] = 0.0;
            Gradient[1*nPlane + iPlane] = 0.0;
            Gradient[2*nPlane + iPlane] = 0.0;
            Gradient[3*nPlane + iPlane] = 0.0;
            Gradient[4*nPlane + iPlane] = 0.0;
            Gradient[5*nPlane + iPlane] = 0.0;
            Gradient[6*nPlane + iPlane] = 0.0;
            Gradient[7*nPlane + iPlane] = 0.0;
            Gradient[8*nPlane + iPlane] = 0.0;
          }
          
        }
        
        /*--- Screen output ---*/
        
        if (geometry_container[ZONE_0]->GetnDim() == 3) {
          cout << "\nWing volume grad.: "    << Wing_Volume_Grad << ". ";
          cout << "Wing min. max. thickness grad.: "  << Wing_MinMaxThickness_Grad << ". ";
          cout << "Wing max. chord grad.: "  << Wing_MaxChord_Grad << "." << endl;
          cout << "Wing min. ToC grad.: "  << Wing_MinToC_Grad << ". ";
          cout << "Wing max. twist grad.: "  << Wing_MaxTwist_Grad << ". ";
          cout << "Wing max. curvature grad.: "  << Wing_MaxCurvature_Grad << "." << endl;
          cout << "Wing max. dihedral grad.: "  << Wing_MaxDihedral_Grad << "." << endl;
        }
        
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          if (Xcoord_Airfoil[iPlane].size() != 0) {
            cout << "\nStation " << (iPlane+1) << ". Plane (yCoord): " << Plane_P0[iPlane][1] << "." << endl;
            cout << "Max. thick. grad.: "    << Gradient[iPlane] << ". ";
            cout << "1/3c thick. grad.: "  << Gradient[2*nPlane + iPlane] << ". ";
            cout << "2/3c thick. grad.: "  << Gradient[4*nPlane + iPlane] << "." << endl;
            cout << "1/4c thick. grad.: "  << Gradient[1*nPlane + iPlane] << ". ";
            cout << "1/2c thick. grad.: "  << Gradient[3*nPlane + iPlane] << ". ";
            cout << "3/4c thick. grad.: "  << Gradient[5*nPlane + iPlane] << "." << endl;
            cout << "Area grad.: "                 << Gradient[6*nPlane + iPlane] << ". ";
            cout << "Twist angle grad.: "      << Gradient[7*nPlane + iPlane] << ". ";
            cout << "Chord grad.: "                << Gradient[8*nPlane + iPlane] << "." << endl;
          }
        }
        cout << endl;
        
        
        if (iDV == 0) {
          Gradient_file << "TITLE = \"SU2_GEO Gradient\"" << endl;
          
          if (geometry_container[ZONE_0]->GetnDim() == 2) {
            Gradient_file << "VARIABLES = \"DESIGN_VARIABLE\",\"MAX_THICKNESS\",\"1/4_THICKNESS\",\"1/3_THICKNESS\",\"1/2_THICKNESS\",\"2/3_THICKNESS\",\"3/4_THICKNESS\",\"AREA\",\"AOA\",\"CHORD\"";
          }
          else if (geometry_container[ZONE_0]->GetnDim() == 3) {
            Gradient_file << "VARIABLES = \"DESIGN_VARIABLE\",";
            Gradient_file << "\"WING_VOLUME\",\"WING_MIN_MAXTHICKNESS\",\"WING_MAX_CHORD\",\"WING_MIN_TOC\",\"WING_MAX_TWIST\",\"WING_MAX_CURVATURE\",\"WING_MAX_DIHEDRAL\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"MAX_THICKNESS_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"1/4_THICKNESS_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"1/3_THICKNESS_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"1/2_THICKNESS_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"2/3_THICKNESS_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"3/4_THICKNESS_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"AREA_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) Gradient_file << "\"AOA_SEC"<< (iPlane+1) << "\",";
            for (iPlane = 0; iPlane < nPlane; iPlane++) {
              Gradient_file << "\"CHORD_SEC"<< (iPlane+1) << "\"";
              if (iPlane != nPlane-1) Gradient_file << ",";
            }
          }
          
          Gradient_file << "\nZONE T= \"Geometrical variables (gradient)\"" << endl;
          
        }
        
        Gradient_file << (iDV) <<",";
        
        if (geometry_container[ZONE_0]->GetnDim() == 3) {
          Gradient_file << Wing_Volume_Grad <<","<< Wing_MinMaxThickness_Grad <<","<< Wing_MaxChord_Grad <<","<< Wing_MinToC_Grad
          <<","<< Wing_MaxTwist_Grad <<","<< Wing_MaxCurvature_Grad <<","<< Wing_MaxDihedral_Grad <<",";
        }
        
        for (iPlane = 0; iPlane < nPlane*9; iPlane++) {
          Gradient_file << Gradient[iPlane];
          if (iPlane != (nPlane*9)-1) Gradient_file <<",";
        }
        
        Gradient_file << endl;
        
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
