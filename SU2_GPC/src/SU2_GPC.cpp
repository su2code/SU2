/*!
 * \file SU2_GPC.cpp
 * \brief Main file of the Gradient Projection Code (SU2_GPC).
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

#include "../include/SU2_GPC.hpp"
using namespace std;

int main(int argc, char *argv[]) {	
  
  /*--- Local variables ---*/
	unsigned short iMarker, iDim, iDV, iFFDBox, nZone = 1;
	unsigned long iVertex, iPoint;
	double delta_eps, my_Gradient, Gradient, *Normal, dS;
  double *VarCoord, Sensitivity;
	double dalpha[3], deps[3], dalpha_deps;
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
	CFreeFormDefBox** FFDBox = NULL;
	CConfig *config = NULL;
	CGeometry *boundary = NULL;
	CSurfaceMovement *surface_mov = NULL;
	
	/*--- Definition of the Class for the definition of the problem ---*/
	if (argc == 2) config = new CConfig(argv[1], SU2_GPC, ZONE_0, nZone, VERB_HIGH);
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
	
	/*--- Definition of the Class for the boundary of the geometry,
   note that the orientation of the elements is not checked ---*/
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
  
  /*--- Load the surface sensitivities from file. This is done only
   once: if this is an unsteady problem, a time-average of the surface
   sensitivities at each node is taken within this routine. ---*/
  if (rank == MASTER_NODE) cout << "Reading surface sensitivities at each node from file." << endl; 
  boundary->SetBoundSensitivity(config);
  
  /*--- Boolean controlling points to be updated ---*/
	UpdatePoint = new bool[boundary->GetnPoint()];
	
	/*--- Definition of the Class for surface deformation ---*/
	surface_mov = new CSurfaceMovement();
	
	/*--- Definition of the FFD deformation class ---*/
	unsigned short nFFDBox = MAX_NUMBER_FFD;
	FFDBox = new CFreeFormDefBox*[nFFDBox];
	
	if (rank == MASTER_NODE) 
		cout << endl <<"---------- Start gradient evaluation using surface sensitivity ----------" << endl;
	
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
				case COSINE_BUMP : surface_mov->SetCosBump(boundary, config, iDV, true); break;
				case FOURIER : surface_mov->SetFourier(boundary, config, iDV, true); break;
				case HICKS_HENNE : surface_mov->SetHicksHenne(boundary, config, iDV, true); break;
				case DISPLACEMENT : surface_mov->SetDisplacement(boundary, config, iDV, true); break;
				case ROTATION : surface_mov->SetRotation(boundary, config, iDV, true); break;
				case NACA_4DIGITS : surface_mov->SetNACA_4Digits(boundary, config); break;
				case PARABOLIC : surface_mov->SetParabolic(boundary, config); break;
			}
			
		}
		
		/*--- Free Form deformation for 3D problems ---*/
		if (boundary->GetnDim() == 3) {
			
      if (config->GetDesign_Variable(0) == SPHERICAL)  {
        
        if (rank == MASTER_NODE) {
          cout << endl << "Design variable number "<< iDV <<"." << endl;
          cout << "Perform 3D deformation of the surface." << endl;
        }
        surface_mov->SetSpherical(boundary, config, iDV, true);
        
      }
      else {
        
        /*--- Read the FFD information in the first iteration ---*/
        if (iDV == 0) {
          
          if (rank == MASTER_NODE)
            cout << "Read the FFD information from mesh file." << endl;
          
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
      }
			
		}
		
		/*--- Continuous adjoint gradient computation ---*/
		if (rank == MASTER_NODE)
			cout << "Evaluate functional gradient using the continuous adjoint strategy." << endl;
		
    /*--- Load the delta change in the design variable (finite difference step). 
     Note that this assumes DV_Value_New = finite_diff_step and DV_Value_Old = 0.0
     in the config file. ---*/
		delta_eps = config->GetDV_Value_New(iDV);
    my_Gradient = 0.0; Gradient = 0.0;
      
      /*--- Reset update points ---*/
      for (iPoint = 0; iPoint < boundary->GetnPoint(); iPoint++)
        UpdatePoint[iPoint] = true;
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
				if (config->GetMarker_All_Moving(iMarker) == YES) {
					for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
						
						iPoint = boundary->vertex[iMarker][iVertex]->GetNode();
						if ((iPoint < boundary->GetnPointDomain()) && UpdatePoint[iPoint]) {
							
							Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
							VarCoord = boundary->vertex[iMarker][iVertex]->GetVarCoord();
							Sensitivity = boundary->vertex[iMarker][iVertex]->GetAuxVar();
							
							dS = 0.0; 
							for (iDim = 0; iDim < boundary->GetnDim(); iDim++) {
								dS += Normal[iDim]*Normal[iDim];
								deps[iDim] = VarCoord[iDim] / delta_eps;
							}
							dS = sqrt(dS);
							
							dalpha_deps = 0.0;
							for (iDim = 0; iDim < boundary->GetnDim(); iDim++) {
								dalpha[iDim] = Normal[iDim] / dS;
								dalpha_deps -= dalpha[iDim]*deps[iDim];
							}
							
							my_Gradient += Sensitivity*dalpha_deps;
							UpdatePoint[iPoint] = false;
						}
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
					cout << "Sideforce coefficient gradient: "<< Gradient << "." << endl; break;
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
        case NOISE :
					if (iDV == 0) Gradient_file << "Noise grad. using cont. adj."<< endl;
					cout << "Noise gradient: "<< Gradient << "." << endl; break;
        case HEAT_LOAD :
					if (iDV == 0) Gradient_file << "Integrated surface heat flux. using cont. adj."<< endl;
					cout << "Heat load gradient: "<< Gradient << "." << endl; break;
          
			}
			
			Gradient_file << Gradient << endl;
			
			cout <<"-------------------------------------------------------------------------" << endl;
			      
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
	
	/*--- End solver ---*/
	if (rank == MASTER_NODE) 
	  cout << endl <<"------------------------- Exit Success (SU2_GPC) ------------------------" << endl << endl;
	
	return EXIT_SUCCESS;
	
}
