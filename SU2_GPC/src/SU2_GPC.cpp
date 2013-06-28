/*!
 * \file SU2_GPC.cpp
 * \brief Main file of the Gradient Projection Code (SU2_GPC).
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009). 
 * \version 1.1.
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

#include "../include/SU2_GPC.hpp"
using namespace std;

int main(int argc, char *argv[])
{	
	unsigned short iMarker, iDim, iDV, nZone = 1;
	unsigned long iVertex, iPoint;
	double delta_eps, my_Gradient, Gradient, *Normal, dS, *VarCoord, Sensitivity, dalpha_dx, dalpha_dy, 
	dalpha_dz, dx_deps, dy_deps, dz_deps, dalpha_deps;
	char *cstr;
	ofstream Gradient_file;
	bool *UpdatePoint;
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
	CFreeFormChunk** chunk = NULL;
	CConfig *config = NULL;
	CGeometry *boundary = NULL;
	CSurfaceMovement *surface_mov = NULL;
	
	/*--- Definition of the Class for the definition of the problem ---*/
	if (argc == 2) config = new CConfig(argv[1], SU2_GPC);
	else {
		char grid_file[200];
		strcpy (grid_file, "default.cfg");
		config = new CConfig(grid_file, SU2_GPC);
	}
	
#ifndef NO_MPI
	/*--- Change the name of the input-output files for the 
	 parallel computation ---*/
	config->SetFileNameDomain(rank+1);
#endif
	
	/*--- Definition of the Class for the boundary of the geometry ---*/
	boundary = new CBoundaryGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());
  
  /*--- Perform the non-dimensionalization, in case any values are needed ---*/
  config->SetNondimensionalization(boundary->GetnDim(), MASTER_NODE, nZone);
  
	if (rank == MASTER_NODE)
		cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
  /*--- Boundary geometry preprocessing ---*/
	if (rank == MASTER_NODE) cout << "Identify vertices." <<endl; 
	boundary->SetVertex();
	
	/*--- Create the control volume structures ---*/
	if (rank == MASTER_NODE) cout << "Set boundary control volume structure." << endl; 
	boundary->SetBoundControlVolume(config, ALLOCATE);
	
	/*--- Create the control volume structures ---*/
	if (rank == MASTER_NODE) cout << "Set boundary sensitivity." << endl; 
	boundary->SetBoundSensitivity(config);
	
	UpdatePoint = new bool[boundary->GetnPoint()];
	
	/*--- Definition of the Class for surface deformation ---*/
	surface_mov = new CSurfaceMovement();
	
	/*--- Definition of the FFD deformation class ---*/
	unsigned short nChunk = MAX_NUMBER_CHUNCK;
	chunk = new CFreeFormChunk*[nChunk];
	
	if (rank == MASTER_NODE) 
		cout << endl <<"---------- Start gradient evaluation using surface sensitivity ----------" << endl;
	
	/*--- Write the gradient in a external file ---*/
	if (rank == MASTER_NODE) {
		cstr = new char [config->GetObjFunc_Grad_FileName().size()+1];
		strcpy (cstr, config->GetObjFunc_Grad_FileName().c_str());
		Gradient_file.open(cstr, ios::out);
	}
	
	for (iDV = 0; iDV < config->GetnDV(); iDV++) {
		
		if (rank == MASTER_NODE)
			cout << "Design variable number "<< iDV <<"." << endl;
		
		/*--- Bump deformation for 2D problems ---*/
		if (boundary->GetnDim() == 2) {
			
			if (rank == MASTER_NODE)
				cout << "Perform 2D deformation of the surface." << endl;
			
			switch ( config->GetDesign_Variable(iDV) ) {
				case HICKS_HENNE :
					surface_mov->SetHicksHenne(boundary, config, iDV, true); break;
				case DISPLACEMENT :
					surface_mov->SetDisplacement(boundary, config, iDV, true); break;
				case ROTATION :
					surface_mov->SetRotation(boundary, config, iDV, true); break;
				case NACA_4DIGITS :
					surface_mov->SetNACA_4Digits(boundary, config); break;
				case PARABOLIC :
					surface_mov->SetParabolic(boundary, config); break;
			}
			
		}
		
		/*--- Free Form deformation for 3D problems ---*/
		if (boundary->GetnDim() == 3) {
			
			/*--- Read the FFD information in hte first iteration ---*/
			if (iDV == 0) {
				
				if (rank == MASTER_NODE)
					cout << "Read the FFD information from mesh file." << endl;
				
				/*--- Read the FFD information from the grid file ---*/
				surface_mov->ReadFFDInfo(config, boundary, chunk, config->GetMesh_FileName());
				
				/*--- If the chunk was not defined in the input file ---*/
				if (!surface_mov->GetChunkDefinition() && (rank == MASTER_NODE)) {
					cout << "The input grid doesn't have the entire FFD information!" << endl;
					cout << "Press any key to exit..." << endl;
					cin.get();
				}
				
			}
			
			if (rank == MASTER_NODE)
				cout << "Perform 3D deformation of the surface." << endl;
			
			/*--- Apply the control point change ---*/
			switch ( config->GetDesign_Variable(iDV) ) {
				case FFD_CONTROL_POINT : 
					surface_mov->SetFFDCPChange(boundary, config, chunk, iDV, true); break;
				case FFD_DIHEDRAL_ANGLE : 
					surface_mov->SetFFDDihedralAngle(boundary, config, chunk, iDV, true); break;
				case FFD_TWIST_ANGLE : 
					surface_mov->SetFFDTwistAngle(boundary, config, chunk, iDV, true); break;
				case FFD_ROTATION : 
					surface_mov->SetFFDRotation(boundary, config, chunk, iDV, true); break;
				case FFD_CAMBER : 
					surface_mov->SetFFDCamber(boundary, config, chunk, iDV, true); break;
				case FFD_THICKNESS : 
					surface_mov->SetFFDThickness(boundary, config, chunk, iDV, true); break;
				case FFD_VOLUME : 
					surface_mov->SetFFDVolume(boundary, config, chunk, iDV, true); break;
			}
			
			/*--- Recompute cartesian coordinates using the new control points position ---*/
			surface_mov->SetCartesianCoord(boundary, config, chunk);
			
		}
		
		/*--- Continuos adjoint gradient computation ---*/
		if (rank == MASTER_NODE)
			cout << "Evaluate functional gradient using the continuous adjoint strategy." << endl;
		
		for (iPoint = 0; iPoint < boundary->GetnPoint(); iPoint++)
			UpdatePoint[iPoint] = true;
		
		delta_eps = config->GetDV_Value_New(iDV);
		
		my_Gradient = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
			for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				iPoint = boundary->vertex[iMarker][iVertex]->GetNode();
				if ((iPoint < boundary->GetnPointDomain()) && 
						(config->GetMarker_All_Moving(iMarker) == YES) && UpdatePoint[iPoint]) {
					
					Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
					dS = 0.0; for (iDim = 0; iDim < boundary->GetnDim(); iDim++) dS += Normal[iDim]*Normal[iDim];
					dS = sqrt(dS);
					VarCoord = boundary->vertex[iMarker][iVertex]->GetVarCoord();
					Sensitivity = boundary->vertex[iMarker][iVertex]->GetAuxVar();
					
					dalpha_dx = Normal[0] / dS; 
					dalpha_dy = Normal[1] / dS;
					
					dx_deps = VarCoord[0] / delta_eps;
					dy_deps = VarCoord[1] / delta_eps;
					
					dalpha_deps = -dalpha_dx*dx_deps - dalpha_dy*dy_deps;
					
					if (boundary->GetnDim() == 3) {
						dalpha_dz = Normal[2] / dS;
						dz_deps = VarCoord[2] / delta_eps;
						dalpha_deps -= dalpha_dz*dz_deps;
					}
					
					my_Gradient += Sensitivity * dalpha_deps / config->GetScale_GradOF();
					UpdatePoint[iPoint] = false;
				}
			}
		}
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Allreduce(&my_Gradient, &Gradient, 1, MPI::DOUBLE, MPI::SUM); 
#else
		Gradient = my_Gradient;
#endif
		
		if (rank == MASTER_NODE) {
			switch (config->GetKind_ObjFunc()) {
				case LIFT_COEFFICIENT : 
					if (iDV == 0) Gradient_file << "Lift coeff. grad. using cont. adj." << endl;
					cout << "Lift coefficient gradient: "<< Gradient << "." << endl; break;
				case DRAG_COEFFICIENT : 
					if (iDV == 0) Gradient_file << "Drag coeff. grad. using cont. adj." << endl;
					cout << "Drag coefficient gradient: "<< Gradient << "." << endl; break;
				case SIDEFORCE_COEFFICIENT :
					if (iDV == 0) Gradient_file << "Sideforce coeff. grad. using cont. adj." << endl;
					cout << "Sideforce coefficient gradient: "<< Gradient << "." << endl; 
					break;
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
					cout << "Equivalent area coefficient gradient: "<< Gradient << "." << endl; break;
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
				case FREESURFACE :
					if (iDV == 0) Gradient_file << "Free-Surface grad. using cont. adj."<< endl;
					cout << "Free-surface gradient: "<< Gradient << "." << endl; break;
			}
			
			
			Gradient_file << Gradient << endl;
			
			cout << endl <<"-------------------------------------------------------------------------" << endl;
			
		}
		
	}
	
	if (rank == MASTER_NODE)
		Gradient_file.close();
	
	delete [] UpdatePoint;
	
#ifndef NO_MPI
	/*--- Finalize MPI parallelization ---*/	
	old_buffer = buffer;
	MPI::Detach_buffer(old_buffer);
	//	delete [] buffer;
	MPI::Finalize();
#endif
	
	return 1;
	
}
