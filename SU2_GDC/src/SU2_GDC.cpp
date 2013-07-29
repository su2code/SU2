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
	unsigned short iDV, nZone = 1;
	double ObjectiveFunc, ObjectiveFunc_New, Gradient, delta_eps;
	char *cstr;
	ofstream Gradient_file, ObjFunc_file;
	int rank = MASTER_NODE;
	
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
	void *buffer, *old_buffer;
	int bufsize;
	bufsize = MAX_MPI_BUFFER;
	buffer = new char[bufsize];
	MPI::Init(argc,argv);
	MPI::Attach_buffer(buffer, bufsize);
	rank = MPI::COMM_WORLD.Get_rank();
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
		unsigned short nFFDBox = MAX_NUMBER_FFD;
		FFDBox = new CFreeFormDefBox*[nFFDBox];
		
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
			if (boundary->GetnDim() == 3) {
				if (rank == MASTER_NODE)
					cout << "No geometrical constraints defined for 3D problems." << endl;
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
	old_buffer = buffer;
	MPI::Detach_buffer(old_buffer);
	//	delete [] buffer;
	MPI::Finalize();
#endif
	
	/*--- End solver ---*/
	if (rank == MASTER_NODE) 
		cout << endl <<"------------------------- Exit Success (SU2_GDC) ------------------------" << endl << endl;
		
	return EXIT_SUCCESS;
	
}
