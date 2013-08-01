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
	double ObjectiveFunc[5], ObjectiveFunc_New[5], Gradient[5], delta_eps, Max_Thickness[100], MinPlane, MaxPlane;
  vector<double> Xcoord_Airfoil[5], Ycoord_Airfoil[5], Zcoord_Airfoil[5];
  double Plane_P0[5][3] = {0.0, 0.0, 0.0}, Plane_Normal[5][3];
  
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
  
	if (rank == MASTER_NODE)
		cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
  /*--- Boundary geometry preprocessing ---*/
	if (rank == MASTER_NODE) cout << "Identify vertices." <<endl;
	boundary->SetVertex();
	
	/*--- Create the control volume structures ---*/
	if (rank == MASTER_NODE) cout << "Set boundary control volume structure." << endl;
	boundary->SetBoundControlVolume(config, ALLOCATE);
	
	/*--- Create airfoil section structure ---*/
  if (boundary->GetnDim() == 3) {
    
    if (rank == MASTER_NODE) cout << "Set airfoil section structure." << endl;
    MinPlane = config->GetSection_Limit(0); MaxPlane = config->GetSection_Limit(1);
    
    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      
      Plane_P0[iPlane][0] = 0.0;
      Plane_P0[iPlane][1] = MinPlane + iPlane*(MaxPlane - MinPlane)/double(nPlane-1);
      Plane_P0[iPlane][2] = 0.0;
      
      Plane_Normal[iPlane][0] = 0.0;
      Plane_Normal[iPlane][1] = 1.0;
      Plane_Normal[iPlane][2] = 0.0;
      
      boundary->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config,
                                       Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
    }
    
  }
  
  if (rank == MASTER_NODE) {
    
    /*--- Evaluate objective function ---*/
    if (boundary->GetnDim() == 2) {
      
      /*--- The 2D geometrical constraints are not defined in parallel ---*/
      if (size == 1) {
        
        ObjectiveFunc[0] = boundary->GetMaxThickness(config, true);
        ObjectiveFunc[1] = boundary->GetTotalVolume(config, true);
        
        cout << "Maximum thickness: "<< ObjectiveFunc[0] << "." << endl;
        cout << "Total volume: "<< ObjectiveFunc[1] << "." << endl;
        
      }
      else {
        cout << "The 2D geometrical constraints are not defined in parallel!!!!!" << endl;
        exit(1);
      }
      
    }
    
    else if (boundary->GetnDim() == 3) {
      
      /*--- Evaluate the thickness ---*/
      for (iPlane = 0; iPlane < nPlane; iPlane++) {
        Max_Thickness[iPlane] = boundary->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], true);
      }
      
      ObjectiveFunc[0] = Max_Thickness[0];
      ObjectiveFunc[1] = Max_Thickness[1];
      ObjectiveFunc[2] = Max_Thickness[2];
      ObjectiveFunc[3] = Max_Thickness[3];
      ObjectiveFunc[4] = Max_Thickness[4];
      
    }
    
    /*--- Write the objective function in a external file ---*/
    
		cstr = new char [config->GetObjFunc_Value_FileName().size()+1];
		strcpy (cstr, config->GetObjFunc_Value_FileName().c_str());
		ObjFunc_file.open(cstr, ios::out);
    ObjFunc_file << "TITLE = \"SU2_GDC Simulation\"" << endl;
    if (boundary->GetnDim() == 2) {
      ObjFunc_file << "VARIABLES = \"MAX_THICKNESS\",\"TOTAL_VOLUME\"" << endl;
    }
    else if (boundary->GetnDim() == 3) {
      ObjFunc_file << "VARIABLES = \"MAX_THICK_SEC1\",\"MAX_THICK_SEC2\",\"MAX_THICK_SEC3\",\"MAX_THICK_SEC4\",\"MAX_THICK_SEC5\"" << endl;
    }
    ObjFunc_file << "ZONE T= \"Geometrical variables (value)\"" << endl;
    if (boundary->GetnDim() == 2) {
      ObjFunc_file << ObjectiveFunc[0] <<", "<< ObjectiveFunc[1] << endl;
    }
    else if (boundary->GetnDim() == 3) {
      ObjFunc_file << ObjectiveFunc[0] <<", "<< ObjectiveFunc[1] <<", "<< ObjectiveFunc[2] <<", "<< ObjectiveFunc[3] <<", "<< ObjectiveFunc[4] << endl;
		}
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
        
        /*--- Create airfoil structure ---*/
        for (iPlane = 0; iPlane < nPlane; iPlane++) {
          boundary->ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
        }
        
 			}
			
      
			/*--- Compute gradient ---*/
			if (rank == MASTER_NODE) {
        
        delta_eps = config->GetDV_Value_New(iDV);

        if (boundary->GetnDim() == 2) {
          
          ObjectiveFunc_New[0] = boundary->GetMaxThickness(config, false);
          Gradient[0] = (ObjectiveFunc_New[0] - ObjectiveFunc[0]) / delta_eps;
          
          ObjectiveFunc_New[2] = boundary->GetTotalVolume(config, false);
          Gradient[1] = (ObjectiveFunc_New[1] - ObjectiveFunc[1]) / delta_eps;
          
          cout << "Max thickness gradient: "<< Gradient[0] << "." << endl;
          cout << "Total volume gradient: "<< Gradient[1] << "." << endl;
          
          if (iDV == 0) {
            Gradient_file << "TITLE = \"SU2_GDC Simulation\"" << endl;
            Gradient_file << "VARIABLES = \"MAX_THICKNESS\",\"TOTAL_VOLUME\"" << endl;
            Gradient_file << "ZONE T= \"Geometrical variables (gradient)\"" << endl;
          }
          
          Gradient_file << Gradient[0] << ", " << Gradient[1] << endl;
          
        }
        else if (boundary->GetnDim() == 3) {
                    
          for (iPlane = 0; iPlane < nPlane; iPlane++) {
            Max_Thickness[iPlane] = boundary->Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], false);
          }
          
          ObjectiveFunc_New[0] = Max_Thickness[0];
          Gradient[0] = (ObjectiveFunc_New[0] - ObjectiveFunc[0]) / delta_eps;
          ObjectiveFunc_New[1] = Max_Thickness[1];
          Gradient[1] = (ObjectiveFunc_New[1] - ObjectiveFunc[1]) / delta_eps;
          ObjectiveFunc_New[2] = Max_Thickness[2];
          Gradient[2] = (ObjectiveFunc_New[2] - ObjectiveFunc[2]) / delta_eps;
          ObjectiveFunc_New[3] = Max_Thickness[3];
          Gradient[3] = (ObjectiveFunc_New[3] - ObjectiveFunc[3]) / delta_eps;
          ObjectiveFunc_New[4] = Max_Thickness[4];
          Gradient[4] = (ObjectiveFunc_New[4] - ObjectiveFunc[4]) / delta_eps;
          
          cout << "Section 1. MaxThickness gradient: " << Gradient[0] << "." << endl;
          cout << "Section 2. MaxThickness gradient: " << Gradient[1] << "." << endl;
          cout << "Section 3. MaxThickness gradient: " << Gradient[2] << "." << endl;
          cout << "Section 4. MaxThickness gradient: " << Gradient[3] << "." << endl;
          cout << "Section 5. MaxThickness gradient: " << Gradient[4] << "." << endl;
          
          if (iDV == 0) {
            Gradient_file << "TITLE = \"SU2_GDC Simulation\"" << endl;
            Gradient_file << "VARIABLES = \"MAX_THICK_SEC1\",\"MAX_THICK_SEC2\",\"MAX_THICK_SEC3\",\"MAX_THICK_SEC4\",\"MAX_THICK_SEC5\"" << endl;
            Gradient_file << "ZONE T= \"Geometrical variables (gradient)\"" << endl;
          }
          
          Gradient_file << Gradient[0] << ", " << Gradient[1] << ", " << Gradient[2] << ", " << Gradient[3] << ", " << Gradient[4] << endl;
          
        }
 				
				
				if (iDV != (config->GetnDV()-1)) cout <<"-------------------------------------------------------------------------" << endl;
				
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
