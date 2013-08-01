/*!
 * \file SU2_GDC.cpp
 * \brief Main file of the Geometry Definition Code (SU2_GDC).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/SU2_GDC.hpp"
using namespace std;

int main(int argc, char *argv[]) {
  
  /*--- Local variables ---*/
	unsigned short iDV, nZone = 1, iFFDBox, iPlane, nPlane = 5;
	double ObjectiveFunc, ObjectiveFunc_New, Gradient, delta_eps, Plane_P0[3] = {0.0, 0.0, 0.0}, Plane_Normal[3] = {0.0, 1.0, 0.0}, Max_Thickness[100];
  vector<double> Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil;
  double MinPlane = 1.5, MaxPlane = 3.5;

	char *cstr;
	ofstream Gradient_file, ObjFunc_file;
	int rank = MASTER_NODE;
  int size = 1;
  
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
	MPI::Init(argc,argv);
	rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#endif
	
	/*--- Pointer to different structures that will be used throughout the entire code ---*/
	CFreeFormDefBox** FFDBox = NULL;
	CConfig *config = NULL;
	CGeometry *boundary = NULL;
	CSurfaceMovement *surface_mov = NULL;
	
	/*--- Definition of the Class for the definition of the problem ---*/
	if (argc == 2) config = new CConfig(argv[1], SU2_GDC, ZONE_0, nZone, VERB_HIGH);
	else {
		char grid_file[200];
		strcpy (grid_file, "default.cfg");
		config = new CConfig(grid_file, SU2_GPC, ZONE_0, nZone, VERB_HIGH);
	}
	
#ifndef NO_MPI
	/*--- Change the name of the input-output files for the
	 parallel computation ---*/
	config->SetFileNameDomain(rank+1);
#endif
	
	/*--- Definition of the Class for the boundary of the geometry ---*/
	boundary = new CBoundaryGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());
  
  /*--- Perform the non-dimensionalization, in case any values are needed ---*/
  config->SetNondimensionalization(boundary->GetnDim(), ZONE_0);
  
	if (rank == MASTER_NODE)
		cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
  /*--- Boundary geometry preprocessing ---*/
	if (rank == MASTER_NODE) cout << "Identify vertices." <<endl;
	boundary->SetVertex();
	
	/*--- Create the control volume structures ---*/
	if (rank == MASTER_NODE) cout << "Set boundary control volume structure." << endl;
	boundary->SetBoundControlVolume(config, ALLOCATE);
	

  /*--- Evaluate objective function ---*/
  if (boundary->GetnDim() == 2) {
    
    /*--- The 2D geometrical constraints are not defined in parallel ---*/
    if ((rank == MASTER_NODE) && (size == 1)) {
      switch (config->GetKind_GeoObjFunc()) {
        case MAX_THICKNESS :
          ObjectiveFunc = boundary->GetMaxThickness(config, true);
          cout << "Maximum thickness: "<< ObjectiveFunc << "." << endl;
          break;
        case MIN_THICKNESS :
          ObjectiveFunc = boundary->GetMinThickness(config, true);
          cout << "Minimum thickness: "<< ObjectiveFunc << "." << endl;
          break;
        case TOTAL_VOLUME :
          ObjectiveFunc = boundary->GetTotalVolume(config, true);
          cout << "Total volume: "<< ObjectiveFunc << "." << endl;
          break;
        case CLEARANCE :
          ObjectiveFunc = boundary->GetClearance(config, true);
          cout << "Clearance: "<< ObjectiveFunc << "." << endl;
          break;
      }
    }
    else {
      cout << "The 2D geometrical constraints are not defined in parallel!!!!!" << endl;
      exit(1);
    }
    
  }
  else if (boundary->GetnDim() == 3) {
    
    /*--- Create airfoil structure and compute the thickness ---*/
    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      Plane_P0[1] = MinPlane + iPlane*(MaxPlane - MinPlane)/double(nPlane-1);
      boundary->ComputeAirfoil_Section(Plane_P0, Plane_Normal, iPlane, config, Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil, true);
      Max_Thickness[iPlane] = boundary->Compute_MaxThickness(Plane_P0, Plane_Normal, iPlane, Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil);
    }
    
    switch (config->GetKind_GeoObjFunc()) {
      case MAX_THICK_SEC1 : ObjectiveFunc = Max_Thickness[0]; break;
      case MAX_THICK_SEC2 : ObjectiveFunc = Max_Thickness[1]; break;
      case MAX_THICK_SEC3 : ObjectiveFunc = Max_Thickness[2]; break;
      case MAX_THICK_SEC4 : ObjectiveFunc = Max_Thickness[3]; break;
      case MAX_THICK_SEC5 : ObjectiveFunc = Max_Thickness[4]; break;
    }
    
  }
	
	/*--- Write the objective function in a external file ---*/
	if (rank == MASTER_NODE) {
		cstr = new char [config->GetObjFunc_Eval_FileName().size()+1];
		strcpy (cstr, config->GetObjFunc_Eval_FileName().c_str());
		ObjFunc_file.open(cstr, ios::out);
		ObjFunc_file << ObjectiveFunc << endl;
		ObjFunc_file.close();
	}
	
	if (config->GetGeometryMode() == GRADIENT) {
		
		/*--- Definition of the Class for surface deformation ---*/
		surface_mov = new CSurfaceMovement();
		
		/*--- Definition of the FFD deformation class ---*/
		FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
		
		if (rank == MASTER_NODE)
			cout << endl <<"---------- Start gradient evaluation using finite differences -----------" << endl;
		
		/*--- Write the gradient in a external file ---*/
		if (rank == MASTER_NODE) {
			cstr = new char [config->GetObjFunc_Grad_FileName().size()+1];
			strcpy (cstr, config->GetObjFunc_Grad_FileName().c_str());
			Gradient_file.open(cstr, ios::out);
		}
		
		for (iDV = 0; iDV < config->GetnDV(); iDV++) {
			
			/*--- Bump deformation for 2D problems ---*/
			if (boundary->GetnDim() == 2) {
				
				if (rank == MASTER_NODE)
					cout << "Perform 2D deformation of the surface." << endl;
				
				switch ( config->GetDesign_Variable(iDV) ) {
					case HICKS_HENNE : surface_mov->SetHicksHenne(boundary, config, iDV, true); break;
					case DISPLACEMENT : surface_mov->SetDisplacement(boundary, config, iDV, true); break;
					case ROTATION : surface_mov->SetRotation(boundary, config, iDV, true); break;
					case NACA_4DIGITS : surface_mov->SetNACA_4Digits(boundary, config); break;
					case PARABOLIC : surface_mov->SetParabolic(boundary, config); break;
				}
				
			}
			
			/*--- Free Form deformation for 3D problems ---*/
			else if (boundary->GetnDim() == 3) {
        
        /*--- Read the FFD information in the first iteration ---*/
        if (iDV == 0) {
          
          if (rank == MASTER_NODE) cout << "Read the FFD information from mesh file." << endl;
          
          /*--- Read the FFD information from the grid file ---*/
          surface_mov->ReadFFDInfo(boundary, config, FFDBox, config->GetMesh_FileName(), true);
          
          /*--- If the FFDBox was not defined in the input file ---*/
          if (!surface_mov->GetFFDBoxDefinition() && (rank == MASTER_NODE)) {
            cout << "The input grid doesn't have the entire FFD information!" << endl;
            cout << "Press any key to exit..." << endl;
            cin.get();
          }
          
          if (rank == MASTER_NODE)
            cout <<"-------------------------------------------------------------------------" << endl;
          
        }
        
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 3D deformation of the surface." << endl;
        }
        
        /*--- Apply the control point change ---*/
        for (iFFDBox = 0; iFFDBox < surface_mov->GetnFFDBox(); iFFDBox++) {
          
          switch ( config->GetDesign_Variable(iDV) ) {
            case FFD_CONTROL_POINT : surface_mov->SetFFDCPChange(boundary, config, FFDBox[iFFDBox], iFFDBox, iDV, true); break;
            case FFD_DIHEDRAL_ANGLE : surface_mov->SetFFDDihedralAngle(boundary, config, FFDBox[iFFDBox], iFFDBox, iDV, true); break;
            case FFD_TWIST_ANGLE : surface_mov->SetFFDTwistAngle(boundary, config, FFDBox[iFFDBox], iFFDBox, iDV, true); break;
            case FFD_ROTATION : surface_mov->SetFFDRotation(boundary, config, FFDBox[iFFDBox], iFFDBox, iDV, true); break;
            case FFD_CAMBER : surface_mov->SetFFDCamber(boundary, config, FFDBox[iFFDBox], iFFDBox, iDV, true); break;
            case FFD_THICKNESS : surface_mov->SetFFDThickness(boundary, config, FFDBox[iFFDBox], iFFDBox, iDV, true); break;
            case FFD_VOLUME : surface_mov->SetFFDVolume(boundary, config, FFDBox[iFFDBox], iFFDBox, iDV, true); break;
          }
          
          /*--- Recompute cartesian coordinates using the new control points position ---*/
          surface_mov->SetCartesianCoord(boundary, config, FFDBox[iFFDBox], iFFDBox);
          
        }
        
        /*--- Create airfoil structure and compute the thickness ---*/
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          Plane_P0[1] = MinPlane + iPlane*(MaxPlane - MinPlane)/double(nPlane-1);
          boundary->ComputeAirfoil_Section(Plane_P0, Plane_Normal, iPlane, config, Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil, false);
          Max_Thickness[iPlane] = boundary->Compute_MaxThickness(Plane_P0, Plane_Normal, iPlane, Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil);
        }
        
 			}
			
			/*--- Compute gradient ---*/
			if (rank == MASTER_NODE) {
        
				switch (config->GetKind_GeoObjFunc()) {
					case MAX_THICKNESS :
						ObjectiveFunc_New = boundary->GetMaxThickness(config, false);
						delta_eps = config->GetDV_Value_New(iDV);
						Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
						if (iDV == 0) Gradient_file << "Max thickness grad. using fin. dif." << endl;
						cout << "Max thickness gradient: "<< Gradient << "." << endl;
						break;
          case MIN_THICKNESS :
						ObjectiveFunc_New = boundary->GetMinThickness(config, false);
						delta_eps = config->GetDV_Value_New(iDV);
						Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
						if (iDV == 0) Gradient_file << "Min thickness grad. using fin. dif." << endl;
						cout << "Min thickness gradient: "<< Gradient << "." << endl;
						break;
					case TOTAL_VOLUME :
						ObjectiveFunc_New = boundary->GetTotalVolume(config, false);
						delta_eps = config->GetDV_Value_New(iDV);
						Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
						if (iDV == 0) Gradient_file << "Total volume grad. using fin. dif." << endl;
						cout << "Total volume gradient: "<< Gradient << "." << endl;
						break;
          case CLEARANCE :
						ObjectiveFunc_New = boundary->GetClearance(config, false);
						delta_eps = config->GetDV_Value_New(iDV);
						Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
						if (iDV == 0) Gradient_file << "Clearance grad. using fin. dif." << endl;
						cout << "Clearance gradient: "<< Gradient << "." << endl;
						break;
          case MAX_THICK_SEC1 :
            ObjectiveFunc_New = Max_Thickness[0];
            delta_eps = config->GetDV_Value_New(iDV);
            Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
            if (iDV == 0) Gradient_file << "Max thickness grad. using fin. dif." << endl;
            cout << "Max thickness gradient (section 1): "<< Gradient << "." << endl;
            break;
          case MAX_THICK_SEC2 :
            ObjectiveFunc_New = Max_Thickness[1];
            delta_eps = config->GetDV_Value_New(iDV);
            Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
            if (iDV == 0) Gradient_file << "Max thickness grad. using fin. dif." << endl;
            cout << "Max thickness gradient (section 2): "<< Gradient << "." << endl;
            break;
          case MAX_THICK_SEC3 :
            ObjectiveFunc_New = Max_Thickness[2];
            delta_eps = config->GetDV_Value_New(iDV);
            Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
            if (iDV == 0) Gradient_file << "Max thickness grad. using fin. dif." << endl;
            cout << "Max thickness gradient (section 3): "<< Gradient << "." << endl;
            break;
          case MAX_THICK_SEC4 :
            ObjectiveFunc_New = Max_Thickness[3];
            delta_eps = config->GetDV_Value_New(iDV);
            Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
            if (iDV == 0) Gradient_file << "Max thickness grad. using fin. dif." << endl;
            cout << "Max thickness gradient (section 4): "<< Gradient << "." << endl;
            break;
          case MAX_THICK_SEC5 :
            ObjectiveFunc_New = Max_Thickness[4];
            delta_eps = config->GetDV_Value_New(iDV);
            Gradient = (ObjectiveFunc_New - ObjectiveFunc) / (delta_eps + EPS);
            if (iDV == 0) Gradient_file << "Max thickness grad. using fin. dif." << endl;
            cout << "Max thickness gradient (section 5): "<< Gradient << "." << endl;
            break;            
				}
				
				Gradient_file << Gradient << endl;
				
				cout <<"-------------------------------------------------------------------------" << endl;
				
			}
		}
		
		if (rank == MASTER_NODE)
			Gradient_file.close();
    
	}
	
	
#ifndef NO_MPI
	/*--- Finalize MPI parallelization ---*/
	MPI::Finalize();
#endif
	
	/*--- End solver ---*/
	if (rank == MASTER_NODE)
		cout << endl <<"------------------------- Exit Success (SU2_GDC) ------------------------" << endl << endl;
  
	return EXIT_SUCCESS;
	
}
