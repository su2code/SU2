/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

#include "../include/output_structure.hpp"

COutput::COutput(CConfig *config) {
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

COutput::~COutput(void) { }

void COutput::SetDomain_Flow(CConfig *config, CGeometry *geometry, CSolution **solution_container, 
																string val_filename, unsigned long iExtIter) {
	char cstr[200], buffer[50];
	ofstream ConsVar_file;
	unsigned short iSpecies, iVar, loc, iDim;
	unsigned long iPoint, iElem;
	double coord[3] = {0.0,0.0,0.0}, velocity[3] = {0.0,0.0,0.0}, vorticity[3] = {0.0,0.0,0.0}, mach_number;
	
	unsigned short solver = config->GetKind_Solver();
	bool incompressible = config->GetIncompressible();
	bool scalar = true;
	
	
	/*--- Paraview Output ---*/
	if (config->GetOutput_FileFormat() == PARAVIEW) {
		
		/*--- Write file name with extension ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.vtk", int(iExtIter));
		else sprintf (buffer, ".vtk");
		strcat(cstr,buffer);
		
		/*--- Write header and open the file ---*/
		geometry->SetParaView(cstr);
		ConsVar_file.precision(15);
		ConsVar_file.open(cstr, ios::out | ios::app);
		ConsVar_file << "POINT_DATA " << geometry->GetnPoint() << endl;

		if (solver == ELECTRIC_POTENTIAL) {
			WriteInOutputFile(geometry, solution_container[ELEC_SOL], ConsVar_file, "Electric_Potential", scalar, 0, "GetSolution", config);
			return;
		}
		
		else if (solver == MULTI_SPECIES_NAVIER_STOKES) {
			unsigned short nDim = geometry->GetnDim();
			unsigned short nSpecies = solution_container[PLASMA_SOL]->GetnSpecies();
			unsigned short nDiatomics = solution_container[PLASMA_SOL]->GetnDiatomics();
			
			
			if (config->GetKind_GasModel() == ARGON)
				WriteReactingOutputFile(config, geometry,solution_container, ConsVar_file);
			else if (config->GetKind_GasModel()==AIR7) {
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density", scalar,  loc+0, "GetSolution", config);
					//					WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "", scalar,  5, "GetSolution", config);
				}
			}
		
			return;
		}

		else if (solver == POTENTIAL_FLOW) {
			WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Flow_Potential", scalar,  0, "GetSolution", config);
			WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density_x_Velocity", !scalar,  0, "GetGradient", config);
			return;
		}
		
		else if  ((solver == EULER) || (solver == NAVIER_STOKES) || (solver == RANS) ||
							(solver == TWO_PHASE_FLOW) || (solver == COMBUSTION)) {
			
			if (incompressible) {
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Pressure", scalar, 0, "GetSolution", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Velocity", !scalar, 0, "GetSolution", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density", scalar,  0, "GetDensityInc", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Viscosity", scalar, 0, "GetLaminarViscosityInc", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "UndLaplacian", scalar,  0, "GetUndLaplacian", config);
			}
			
			else {
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density", scalar,  0, "GetSolution", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density_x_Velocity", !scalar, 0, "GetSolution", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density_x_Energy", scalar, geometry->GetnDim()+1, "GetSolution", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Velocity", !scalar,  0, "GetVelocity", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Pressure", scalar,  0, "GetPressure", config);
				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Mach_Number", scalar,  0, "GetMach", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "UndLaplacian", scalar,  0, "GetUndLaplacian", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Limiter", scalar,  0, "GetLimiter", config);
				
				if (config->GetGrid_Movement())
					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Grid_Velocity", !scalar,  0, "GetGridVel", config);
				
				if (config->GetRotating_Frame())
					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Rotational_Velocity", !scalar,  0, "GetRotVel", config);
				
				if ((solver == NAVIER_STOKES) || (solver == RANS) || (solver == TWO_PHASE_FLOW) || (solver == COMBUSTION)) {
					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Temperature", scalar,  0, "GetTemperature", config);
					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Laminar_Viscosity", scalar,  0, "GetLaminarViscosity", config);
				}
				
				if (solver == RANS) {
					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Eddy_Viscosity", scalar,  0, "GetEddyViscosity", config);
					WriteInOutputFile(geometry, solution_container[TURB_SOL], ConsVar_file,  "Nu_Tilde", scalar, 0, "GetSolution", config);
				}
			}
			
			if (solver == TWO_PHASE_FLOW)
				WriteInOutputFile(geometry, solution_container[LEVELSET_SOL], ConsVar_file,   "LevelSet", scalar, 0, "GetSolution", config);
			
			if (solver == COMBUSTION)
				WriteInOutputFile(geometry, solution_container[COMBUSTION_SOL], ConsVar_file,  "Lambda", scalar, 0, "GetSolution", config);
		
			return;
		}
	}
	
	/*--- Tecplot Output ---*/
	if (config->GetOutput_FileFormat() == TECPLOT) {
		
		/*--- Write file name with extension ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.plt", int(iExtIter));
		else sprintf (buffer, ".plt");
		strcat(cstr,buffer);
		
		/*--- Write header and open the file ---*/
		ConsVar_file.precision(15);
		ConsVar_file.open(cstr, ios::out);
		
		ConsVar_file << "TITLE = \"Visualization of the volumetric grid\"" << endl;
		ConsVar_file << "VARIABLES = \"x\", \"y\", \"z\"";
		
		if ((solver == POTENTIAL_FLOW) || (solver == EULER)
				|| (solver == NAVIER_STOKES) || (solver == RANS))
			ConsVar_file << ", \"Density\", \"Velocity_x\", \"Velocity_y\", \"Velocity_z\", \"Pressure\", \"Mach_Number\"";
		
		if (solver == TWO_PHASE_FLOW)
			ConsVar_file << ", \"Density\", \"Velocity_x\", \"Velocity_y\", \"Velocity_z\", \"Pressure\", \"Mach_Number\", \"Theta\"";

		if ((solver == NAVIER_STOKES) || (solver == RANS))
			ConsVar_file << ", \"Temperature\", \"Vorticity_x\", \"Vorticity_y\", \"Vorticity_z\", \"Distance\", \"Laminar_Viscosity\"";
		
		if (solver == RANS)
			ConsVar_file << ", \"Eddy_Viscosity\", \"Nu_Tilde\"";
		
		if (solver == MULTI_SPECIES_NAVIER_STOKES) {
			for ( iSpecies = 0; iSpecies < solution_container[PLASMA_SOL]->GetnSpecies(); iSpecies++ ) {
				ConsVar_file << "\"Density("<< iSpecies << ")\", \"DensityVelocity_x("<< iSpecies << ")\", \"DensityVelocity_y("<< iSpecies << ")\", \"DensityEnergy("<< iSpecies << ")\"";
				if (iSpecies < solution_container[PLASMA_SOL]->GetnSpecies() - 1)
					ConsVar_file << ", ";
			}
		}

		ConsVar_file << endl;
		
		ConsVar_file << "ZONE NODES="<< geometry->GetnPoint() <<" , ELEMENTS="<< geometry->GetnElem() <<", DATAPACKING=POINT";
		if (geometry->GetnDim() == 2) ConsVar_file << ", ZONETYPE=FEQUADRILATERAL" << endl;
		if (geometry->GetnDim() == 3) ConsVar_file << ", ZONETYPE=FEBRICK"<< endl;
		
    /*--- Cycle through all points and write the solution. ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
				
			if ((solver == POTENTIAL_FLOW) || (solver == EULER) || (solver == TWO_PHASE_FLOW) ||
					(solver == NAVIER_STOKES) || (solver == RANS)) {
				
				/*--- Compute vectorial quatities for 2D and 3D ---*/
				for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
					velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim);
					if ((solver == NAVIER_STOKES) || (solver == RANS))
						vorticity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(iDim);
				}
				
				/*--- Compute the mach number ---*/
				mach_number = sqrt(solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2())/
				solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
			}

			/*--- Write the coordinates ---*/
			ConsVar_file << coord[0] <<" "<< coord[1] <<" "<< coord[2];
			
			/*--- Write the Euler variables ---*/
			if ((solver == POTENTIAL_FLOW) || (solver == EULER) || (solver == TWO_PHASE_FLOW) ||
					(solver == NAVIER_STOKES) || (solver == RANS)) {
				if (!incompressible) {
					mach_number = sqrt(solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2())/
					solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
					for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
						velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim);
					ConsVar_file << " " << solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0) <<" "<<
					velocity[0] <<" "<< velocity[1] <<" "<< velocity[2]<<" "<<
					solution_container[FLOW_SOL]->node[iPoint]->GetPressure() <<" "<< mach_number; 					
				}
				else {
					double Density = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref();
					double SoundSpeed = sqrt(config->GetBulk_Modulus()/Density);
					double SoundSpeedND = SoundSpeed / config->GetVelocity_Ref();
					mach_number = sqrt(solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2())/SoundSpeedND;
					for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
						velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iDim+1);
					if (solver == TWO_PHASE_FLOW)
						ConsVar_file << " " << solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc() <<" "<<
						velocity[0] <<" "<< velocity[1] <<" "<< velocity[2]<<" "<<
						solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0) * solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc() 
						<<" "<< mach_number <<" "<< solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
					else
						ConsVar_file << " " << solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc() <<" "<<
						velocity[0] <<" "<< velocity[1] <<" "<< velocity[2]<<" "<<
						solution_container[FLOW_SOL]->node[iPoint]->GetSolution(0) * solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc() 
						<<" "<< mach_number;
				}
			}
			
			/*--- Write the Navier-Stokes variables ---*/
			if ((solver == NAVIER_STOKES) || (solver == RANS))
				ConsVar_file << " "<< solution_container[FLOW_SOL]->node[iPoint]->GetTemperature() <<" "<<
				vorticity[0] <<" "<< vorticity[1] <<" "<< vorticity[2]<<" "<<
				geometry->node[iPoint]->GetWallDistance()<<" "<<solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
			
			/*--- Write the Navier-Stokes variables ---*/
			if (solver == RANS)
				ConsVar_file << " "<< solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity()<< " " <<
        solution_container[TURB_SOL]->node[iPoint]->GetSolution(0);
			
			/*--- Write the multi-species flow variables ---*/
			if (solver == MULTI_SPECIES_NAVIER_STOKES )
				for ( iVar = 0; iVar < solution_container[PLASMA_SOL]->GetnVar(); iVar++ ) {
					ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(iVar);
				}
			/*--- End line ---*/
			ConsVar_file << endl;
			
		}
		
		/*--- Write the grid connectivity. ---*/
		for(iElem = 0; iElem < geometry->GetnElem(); iElem++) {
			if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
				ConsVar_file <<
				geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
				geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) {
				ConsVar_file <<
				geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
				geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
				ConsVar_file <<
				geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
				geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
				geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
				geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
				ConsVar_file <<
				geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
				geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
				geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 <<" "<<
				geometry->elem[iElem]->GetNode(6)+1 <<" "<< geometry->elem[iElem]->GetNode(7)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID) {
				ConsVar_file <<
				geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
				geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
				geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
				geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == WEDGE) {
				ConsVar_file <<
				geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
				geometry->elem[iElem]->GetNode(1)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
				geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
				geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 << endl;
			}
		}
		
		ConsVar_file.close();
	}
}


void COutput::SetDomain_Adjoint(CConfig *config, CGeometry *geometry, CSolution ***solution_container,
																	 string val_filename, unsigned long iExtIter) {
  char cstr[200], buffer [50];
  ofstream AdjVar_file;
	unsigned long iPoint, iElem;
	unsigned short iDim;
	double coord[3] = {0.0,0.0,0.0}, phi[3] = {0.0,0.0,0.0};
	
  unsigned short solver = config->GetKind_Solver();
	bool scalar = true;
	
	/*--- Paraview Output ---*/
	if (config->GetOutput_FileFormat() == PARAVIEW) {
		
		/*--- Write file name with extension ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.vtk", int(iExtIter));
		else sprintf (buffer, ".vtk");
		strcat(cstr,buffer);
		
		/*--- Write header and open the file ---*/
		geometry->SetParaView(cstr);
		AdjVar_file.precision(15);
		AdjVar_file.open(cstr, ios::out | ios::app);
		AdjVar_file << "POINT_DATA " << geometry->GetnPoint() << endl;
		
		if ((solver == ADJ_EULER) || (solver == ADJ_NAVIER_STOKES) || (solver == ADJ_RANS)) {
			WriteInOutputFile(geometry, solution_container[MESH_0][ADJFLOW_SOL], AdjVar_file,  "PsiRho", scalar, 0, "GetSolution", config);
			WriteInOutputFile(geometry, solution_container[MESH_0][ADJFLOW_SOL], AdjVar_file, "Phi", !scalar, 0, "GetSolution", config);
			WriteInOutputFile(geometry, solution_container[MESH_0][ADJFLOW_SOL], AdjVar_file,  "PsiE",  scalar,  geometry->GetnDim()+1, "GetSolution", config);
		}
		
		if (solver == ADJ_RANS) {
			WriteInOutputFile(geometry, solution_container[MESH_0][ADJTURB_SOL], AdjVar_file,   "PsiTurb",  scalar, 0, "GetSolution", config);
		}
	}
	
  /*--- Tecplot Output ---*/
	if(config->GetOutput_FileFormat() == TECPLOT) {
    
		/*--- Write file name with extension ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.plt", int(iExtIter));
		else sprintf (buffer, ".plt");
		strcat(cstr,buffer);
		
		/*--- Write header and open the file ---*/
		AdjVar_file.precision(15);
		AdjVar_file.open(cstr, ios::out);
    
		AdjVar_file << "TITLE = \"Visualization of the volumetric grid\"" << endl;
		AdjVar_file << "VARIABLES = \"x\", \"y\", \"z\"";
    
    /*--- Add variable labels depending on adjoint solver type. ---*/
    if ((solver == ADJ_EULER) || (solver == ADJ_NAVIER_STOKES) || (solver == ADJ_RANS))
      AdjVar_file << ", \"PsiRho\", \"Phi_x\", \"Phi_y\", \"Phi_z\", \"PsiE\"";
		
		if (solver == ADJ_RANS)
			AdjVar_file << ", \"PsiTurb\"";
    
		AdjVar_file << endl;
    
		AdjVar_file << "ZONE NODES="<< geometry->GetnPoint() <<" , ELEMENTS="<< geometry->GetnElem() <<", DATAPACKING=POINT";
		if (geometry->GetnDim() == 2) AdjVar_file << ", ZONETYPE=FEQUADRILATERAL" << endl;
		if (geometry->GetnDim() == 3) AdjVar_file << ", ZONETYPE=FEBRICK"<< endl;
		
    /*--- Cycle through all points and write the solution. ---*/
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
			/*--- Compute vectorial quatities for 2D and 3D ---*/
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
				coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
				phi[iDim] = solution_container[MESH_0][ADJFLOW_SOL]->node[iPoint]->GetSolution(iDim+1);
			}
      
			/*--- Write the coordinates ---*/
			AdjVar_file << coord[0] <<" "<< coord[1] <<" "<< coord[2];
      
      /*--- Write the adjoint variables. ---*/
      if ((solver == ADJ_EULER) || (solver == ADJ_NAVIER_STOKES) || (solver == ADJ_RANS)) {
        AdjVar_file << " " << solution_container[MESH_0][ADJFLOW_SOL]->node[iPoint]->GetSolution(0) 
        <<" "<< phi[0] <<" "<< phi[1] <<" "<< phi[2]<<" "<<
        solution_container[MESH_0][ADJFLOW_SOL]->node[iPoint]->GetSolution(geometry->GetnDim()+1);
			}
			
			/*--- If turbulent, add the turbulent adjoint variable. ---*/
			if (solver == ADJ_RANS) {
				AdjVar_file << " " << solution_container[MESH_0][ADJTURB_SOL]->node[iPoint]->GetSolution(0);
      }
      
      /*--- End line ---*/
      AdjVar_file << endl;
    }
    
    /*--- Write the grid connectivity. ---*/
		for(iElem = 0; iElem < geometry->GetnElem(); iElem++) {
			if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
				AdjVar_file <<
        geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
        geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) {
				AdjVar_file <<
        geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
        geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
				AdjVar_file <<
        geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
        geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
        geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
        geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
				AdjVar_file <<
        geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
        geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
        geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 <<" "<<
        geometry->elem[iElem]->GetNode(6)+1 <<" "<< geometry->elem[iElem]->GetNode(7)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID) {
				AdjVar_file <<
        geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
        geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
        geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
        geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 << endl;
			}
			if (geometry->elem[iElem]->GetVTK_Type() == WEDGE) {
				AdjVar_file <<
        geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
        geometry->elem[iElem]->GetNode(1)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
        geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
        geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 << endl;
			}
		}
    
		AdjVar_file.close();
	}
}  

void COutput::SetLinearized_Variables(CConfig *config, CGeometry *geometry, CSolution ***solution_container, string val_filename, unsigned long iExtIter) {
	char cstr[200], buffer[50];
	ofstream LinVar_file;
	bool scalar = true;

	/*--- Write file name with extension ---*/
	strcpy (cstr, val_filename.c_str());
	sprintf (buffer, ".vtk");
	strcat(cstr,buffer);

	/*--- Write header and open the file ---*/
	geometry->SetParaView(cstr);
	LinVar_file.open(cstr, ios::out | ios::app);

	/*--- Write information ---*/
	LinVar_file << "POINT_DATA " << geometry->GetnPoint() << endl;

	if (config->GetKind_Solver() != LIN_ELECTRIC_POTENTIAL)
		WriteInOutputFile(geometry, solution_container[MESH_0][LINFLOW_SOL], LinVar_file,  "DeltaRho", scalar, 0, "GetSolution", config);

	else {
		WriteInOutputFile(geometry, solution_container[MESH_0][LINELEC_SOL], LinVar_file,  "DeltaPot", scalar, 0, "GetSolution", config);
		return;
	}
	WriteInOutputFile(geometry, solution_container[MESH_0][LINFLOW_SOL], LinVar_file, "DeltaRhoVel", !scalar, 0, "GetSolution", config);

	WriteInOutputFile(geometry, solution_container[MESH_0][LINFLOW_SOL], LinVar_file,   "DeltaRhoE",  scalar,  geometry->GetnDim()+1, "GetSolution", config);
}

void COutput::SetSurface_Flow(CConfig *config, CGeometry *geometry, CSolution *FlowSolution, string val_filename, unsigned long iExtIter) {
	unsigned long iPoint, iVertex;
	unsigned short iMarker, iDim;
	double PressCoeff = 0.0, SkinFrictionCoeff, Mach, *aux_press, *aux_friction;
	char cstr[200], buffer [50];
	ofstream SurfFlow_file;
	unsigned short solver = config->GetKind_Solver();
	bool incompressible = config->GetIncompressible();
	
	if (config->GetOutput_FileFormat() == PARAVIEW) {
		
		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.vtk", int(iExtIter));
		else sprintf (buffer, ".vtk");
		strcat(cstr, buffer);
		
		/*--- Write the geometrical information on the output file ---*/
		geometry->SetBoundParaView(config, cstr);
		
		/*--- Open the output file ---*/
		SurfFlow_file.precision(15);
		SurfFlow_file.open(cstr, ios::out | ios::app);
		
		/*--- It is necessary to go from point to vertex, and an auxiliar variable is created ---*/
		aux_press = new double [geometry->GetnPoint()];
		
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_press[iPoint] = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					aux_press[iPoint] = FlowSolution->GetCPressure(iMarker,iVertex);
				}
		
		SurfFlow_file << "SCALARS Pressure_Coefficient float 1" << endl;
		SurfFlow_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			if (geometry->node[iPoint]->GetBoundary()) {
				PressCoeff = aux_press[iPoint];
				SurfFlow_file << scientific << PressCoeff << endl;
			}
		delete [] aux_press;
		
		switch (solver) {
				
			case POTENTIAL_FLOW: case EULER: case TWO_PHASE_FLOW:
				SurfFlow_file << "SCALARS Mach_Number float 1" << endl;
				SurfFlow_file << "LOOKUP_TABLE default" << endl;
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					if (geometry->node[iPoint]->GetBoundary()) {
						Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2())/FlowSolution->node[iPoint]->GetSoundSpeed();
						SurfFlow_file << scientific << Mach << endl;
					}
				
				break;
				
			case NAVIER_STOKES: case RANS:
				aux_friction = new double [geometry->GetnPoint()];
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_friction[iPoint] = 0.0;
				for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
					if (config->GetMarker_All_Plotting(iMarker) == YES)
						for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
							iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
							aux_friction[iPoint] = FlowSolution->GetCSkinFriction(iMarker,iVertex);
						}
				SurfFlow_file << "SCALARS Skin_Friction_Coefficient float 1" << endl;
				SurfFlow_file << "LOOKUP_TABLE default" << endl;
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					if (geometry->node[iPoint]->GetBoundary()) {
						SkinFrictionCoeff = aux_friction[iPoint];
						SurfFlow_file << scientific << SkinFrictionCoeff << endl;
					}
				delete [] aux_friction;				
				break;
		}
		
		SurfFlow_file.close();
		
	}
	
	if (config->GetOutput_FileFormat() == TECPLOT) {
		
		double coord[3] = {0.0,0.0,0.0};
		
		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.plt", int(iExtIter));
		else sprintf (buffer, ".plt");
		strcat(cstr, buffer);
		
		SurfFlow_file.precision(15);
		SurfFlow_file.open(cstr, ios::out);
		
		
		SurfFlow_file << "TITLE = \"Visualization of the surface grid\"" << endl;
		SurfFlow_file << "VARIABLES = \"x\", \"y\", \"z\"";
		
		if ((solver == POTENTIAL_FLOW) || (solver == EULER) || (solver == TWO_PHASE_FLOW)
				|| (solver == NAVIER_STOKES) || (solver == RANS))
			SurfFlow_file << ", \"Pressure_Coefficient\"";
		
		if ((solver == POTENTIAL_FLOW) || (solver == EULER) || (solver == TWO_PHASE_FLOW))
			SurfFlow_file << ", \"Mach_Number\"";
		
		if ((solver == NAVIER_STOKES) || (solver == RANS))
			SurfFlow_file << ", \"Skin_Friction_Coefficient\"";
		
		SurfFlow_file << endl;
		
		
		/*--- It is important to do a renumering to don't add points that do not belong to the surfaces ---*/
		unsigned long *PointSurface = new unsigned long[geometry->GetnPoint()];
		
		unsigned long nPointSurface = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			if (geometry->node[iPoint]->GetBoundary()) {
				PointSurface[iPoint] = nPointSurface;
				nPointSurface++;
			}
		
		/*--- Compute the total number of elements ---*/
		unsigned long nElemSurface = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES) {
				nElemSurface += geometry->GetnElem_Bound(iMarker);
			}
		
		SurfFlow_file << "ZONE NODES="<< nPointSurface <<" , ELEMENTS="<< nElemSurface <<", DATAPACKING=POINT";
		if (geometry->GetnDim() == 2) SurfFlow_file << ", ZONETYPE=FELINESEG" << endl;
		if (geometry->GetnDim() == 3) SurfFlow_file << ", ZONETYPE=FEQUADRILATERAL"<< endl;
		
		/*--- It is necessary to go from point to vertex, and an auxiliar variable is created ---*/
		aux_press = new double [geometry->GetnPoint()];
		aux_friction = new double [geometry->GetnPoint()];
		
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_press[iPoint] = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					aux_press[iPoint] = FlowSolution->GetCPressure(iMarker,iVertex);
				}
		
		if ((solver == NAVIER_STOKES) || (solver == RANS)) {
			for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_friction[iPoint] = 0.0;
			for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
				if (config->GetMarker_All_Plotting(iMarker) == YES)
					for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						aux_friction[iPoint] = FlowSolution->GetCSkinFriction(iMarker,iVertex);
					}
		}
		
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			
			/*--- Compute vectorial quatities for 2D and 3D ---*/
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
				coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
			}
			
			/*--- Write the coordinates ---*/
			if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << coord[0] <<" "<< coord[1] <<" "<< coord[2];
			
			/*--- Write the Euler variables ---*/
			if ((solver == POTENTIAL_FLOW) || (solver == EULER) ||
					(solver == NAVIER_STOKES) || (solver == RANS)) {
				if (incompressible) PressCoeff = FlowSolution->node[iPoint]->GetSolution(0);
				else PressCoeff = aux_press[iPoint];
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << PressCoeff;
			}
			
			/*--- Write the Mach number ---*/
			if ((solver == POTENTIAL_FLOW) || (solver == EULER)) {
				if (incompressible) Mach = 0;
				else Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2())/FlowSolution->node[iPoint]->GetSoundSpeed();
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << Mach;
			}
			
			/*--- Write the skin friction coefficient variables ---*/
			if ((solver == NAVIER_STOKES) || (solver == RANS)) {
				SkinFrictionCoeff = aux_friction[iPoint];
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << PressCoeff;
			}
			
			/*--- End line ---*/
			if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << endl;
		}
		
		/*--- Write the cells using the new numbering ---*/
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(unsigned long iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
					if (geometry->bound[iMarker][iElem]->GetVTK_Type() == LINE) {
						SurfFlow_file <<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)]+1 << endl;
					}
					if (geometry->bound[iMarker][iElem]->GetVTK_Type() == TRIANGLE) {
						SurfFlow_file <<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)]+1 <<" "<<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)]+1 << endl;
					}
					if (geometry->bound[iMarker][iElem]->GetVTK_Type() == RECTANGLE) {
						SurfFlow_file <<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)]+1 <<" "<<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(3)]+1 << endl;
					}
				}
		
		delete [] aux_press;
		delete [] aux_friction;
		delete [] PointSurface;
		SurfFlow_file.close();
		
	}
	
}

void COutput::SetSurface_Adjoint(CConfig *config, CGeometry *geometry, CSolution *AdjSolution, CSolution *FlowSolution, string val_filename, unsigned long iExtIter) {
	
	unsigned long iPoint, iVertex, iElem;
	double *aux_sens, sens_value;
	unsigned short iMarker, iDim;
	char cstr[200], buffer[50];
	ofstream SurfAdj_file;
	
	if (config->GetOutput_FileFormat() == PARAVIEW) {
		
		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.vtk", int(iExtIter));
		else sprintf (buffer, ".vtk");
		strcat(cstr, buffer);
		
		/*--- Write the geometrical information on the output file ---*/
		geometry->SetBoundParaView (config, cstr);
		
		/*--- Open the output file ---*/
		SurfAdj_file.precision(15);
		SurfAdj_file.open(cstr, ios::out | ios::app);
		
		/*--- It is necessary to go from point to vertex, and an auxiliar variable is created ---*/
		aux_sens = new double [geometry->GetnPoint()];
		
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_sens[iPoint] = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					aux_sens[iPoint] = AdjSolution->GetCSensitivity(iMarker,iVertex);
				}
		
		SurfAdj_file << "SCALARS Shape_Sensitivity float 1" << endl;
		SurfAdj_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			if (geometry->node[iPoint]->GetBoundary()) {
				sens_value = aux_sens[iPoint];
				SurfAdj_file << scientific << sens_value << endl;
			}
		delete [] aux_sens;
		
		SurfAdj_file << "SCALARS PsiRho float 1" << endl;
		SurfAdj_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			if (geometry->node[iPoint]->GetBoundary())
				SurfAdj_file << AdjSolution->node[iPoint]->GetSolution(0) << endl;
		
		SurfAdj_file.close();
	}

	if (config->GetOutput_FileFormat() == TECPLOT) {
		
		double coord[3] = {0.0,0.0,0.0};
		
		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, val_filename.c_str());
		if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.plt", int(iExtIter));
		else sprintf (buffer, ".plt");
		strcat(cstr, buffer);
		
		/*--- Open the output file ---*/
		SurfAdj_file.precision(15);
		SurfAdj_file.open(cstr, ios::out);
		
		SurfAdj_file << "TITLE = \"Visualization of the surface grid\"" << endl;
		SurfAdj_file << "VARIABLES = \"x\", \"y\", \"z\"";
		SurfAdj_file << ", \"Shape_Sensitivity\", \"PsiRho\"";
		SurfAdj_file << endl;
		
		/*--- It is important to do a renumering to don't add points that do not belong to the surfaces ---*/
		unsigned long *PointSurface = new unsigned long[geometry->GetnPoint()];
		
		unsigned long nPointSurface = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			if (geometry->node[iPoint]->GetBoundary()) {
				PointSurface[iPoint] = nPointSurface;
				nPointSurface++;
			}
		
		/*--- Compute the total number of elements ---*/
		unsigned long nElemSurface = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES) {
				nElemSurface += geometry->GetnElem_Bound(iMarker);
			}
		
		SurfAdj_file << "ZONE NODES="<< nPointSurface <<" , ELEMENTS="<< nElemSurface <<", DATAPACKING=POINT";
		if (geometry->GetnDim() == 2) SurfAdj_file << ", ZONETYPE=FELINESEG" << endl;
		if (geometry->GetnDim() == 3) SurfAdj_file << ", ZONETYPE=FEQUADRILATERAL"<< endl;
		
		/*--- It is necessary to go from point to vertex, and an auxiliar variable is created ---*/
		aux_sens = new double [geometry->GetnPoint()];
		
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_sens[iPoint] = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					aux_sens[iPoint] = AdjSolution->GetCSensitivity(iMarker,iVertex);
				}
		
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			
			/*--- Compute vectorial quatities for 2D and 3D ---*/
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
			
			/*--- Write the coordinates ---*/
			if (geometry->node[iPoint]->GetBoundary()) SurfAdj_file << coord[0] <<" "<< coord[1] <<" "<< coord[2];
			
			/*--- Write the sensitivity variables ---*/
			if (geometry->node[iPoint]->GetBoundary()) SurfAdj_file << " " << aux_sens[iPoint];
			
			/*--- Write the adjoint density ---*/
			if (geometry->node[iPoint]->GetBoundary()) SurfAdj_file << " " << AdjSolution->node[iPoint]->GetSolution(0);
			
			/*--- End line ---*/
			if (geometry->node[iPoint]->GetBoundary()) SurfAdj_file << endl;
		}
		
		/*--- Write the cells using the new numbering ---*/
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
					if (geometry->bound[iMarker][iElem]->GetVTK_Type() == LINE) {
						SurfAdj_file <<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)]+1 << endl;
					}
					if (geometry->bound[iMarker][iElem]->GetVTK_Type() == TRIANGLE) {
						SurfAdj_file <<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)]+1 <<" "<<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)]+1 << endl;
					}
					if (geometry->bound[iMarker][iElem]->GetVTK_Type() == RECTANGLE) {
						SurfAdj_file <<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)]+1 <<" "<<
						PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)]+1 <<" "<< PointSurface[geometry->bound[iMarker][iElem]->GetNode(3)]+1 << endl;
					}
				}
		
		delete [] aux_sens;
		delete [] PointSurface;
		SurfAdj_file.close();
		
	}
	
}

void COutput::SetSurface_Linearized(CConfig *config, CGeometry *geometry, CSolution *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::SetSurfaceCSV_Flow_Serial(CConfig *config, CGeometry *geometry, CSolution *FlowSolution, string val_filename, unsigned long iExtIter) {
	unsigned long iPoint, iVertex;
	unsigned short iMarker;
	double PressCoeff = 0.0, SkinFrictionCoeff, xCoord, yCoord, zCoord, Mach;
	char cstr[200], buffer [50];
	ofstream SurfFlow_file;
	unsigned short solver = config->GetKind_Solver();

	/*--- Write the surface .csv file ---*/
	strcpy (cstr, val_filename.c_str());
	if (config->GetUnsteady_Simulation() != NO) sprintf (buffer, "_%d.csv", int(iExtIter));
	else sprintf (buffer, ".csv");
	strcat(cstr, buffer);

	SurfFlow_file.precision(15);
	SurfFlow_file.open(cstr, ios::out);

	if (geometry->GetnDim() == 2) {
		switch (solver) {
		case POTENTIAL_FLOW: case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\"" << endl; break;
		case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\"" << endl; break;
		}
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					PressCoeff = FlowSolution->GetCPressure(iMarker,iVertex);
					switch (solver) {
					case POTENTIAL_FLOW: case EULER :
						Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2()) / FlowSolution->node[iPoint]->GetSoundSpeed();
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << Mach << "," << yCoord << endl;
						break;
					case NAVIER_STOKES: case RANS:
						SkinFrictionCoeff = FlowSolution->GetCSkinFriction(iMarker,iVertex);
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << SkinFrictionCoeff << "," << yCoord << endl;
						break;
					}
				}
	}

	if (geometry->GetnDim() == 3) {
		switch (solver) {
		case POTENTIAL_FLOW: case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\"" << endl; break;
		case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"z_coord\"" << endl; break;
		}
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					zCoord = geometry->node[iPoint]->GetCoord(2);
					PressCoeff = FlowSolution->GetCPressure(iMarker,iVertex);
					switch (solver) {
					case POTENTIAL_FLOW: case EULER :
						Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2()) / FlowSolution->node[iPoint]->GetSoundSpeed();
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << Mach <<"," << yCoord << "," << zCoord << endl;
						break;
					case NAVIER_STOKES: case RANS:
						SkinFrictionCoeff = FlowSolution->GetCSkinFriction(iMarker,iVertex);
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << SkinFrictionCoeff << "," << yCoord << "," << zCoord << endl;
						break;
					}
				}
	}

	SurfFlow_file.close();
}

void COutput::SetSurfaceCSV_Flow_Parallel(CConfig *config, CGeometry *geometry, CSolution *FlowSolution,
																						string val_filename, unsigned long iExtIter) {
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();
	double PressCoeff = 0.0, SkinFrictionCoeff, xCoord, yCoord, zCoord, Mach;
	char cstr[200], buffer [50];
	unsigned short iMarker;
	unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
	MaxLocalVertex_Surface = 0, nBuffer_Scalar, *Buffer_Receive_nVertex = NULL, position;
	ofstream SurfFlow_file;
		
	/*--- Write the surface .csv file, the information of allthe vertices is send to the MASTER_NODE ---*/
	if (rank != MASTER_NODE) {
		nLocalVertex_Surface = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
				}
	}
	else {
		nLocalVertex_Surface = 0;
		Buffer_Receive_nVertex = new unsigned long [nProcessor*sizeof(unsigned long)];
	}
	
	Buffer_Send_nVertex[0] = nLocalVertex_Surface;
	
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	
	double Buffer_Send_Coord_x[MaxLocalVertex_Surface];
	double Buffer_Send_Coord_y[MaxLocalVertex_Surface];
	double Buffer_Send_Coord_z[MaxLocalVertex_Surface];
	double Buffer_Send_Press[MaxLocalVertex_Surface];
	double Buffer_Send_Mach[MaxLocalVertex_Surface];
	double Buffer_Send_SkinFriction[MaxLocalVertex_Surface];
	
	if (rank != MASTER_NODE) {
		nVertex_Surface = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					if (geometry->node[iPoint]->GetDomain()) {
						Buffer_Send_Press[nVertex_Surface] = FlowSolution->GetCPressure(iMarker,iVertex);
						Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
						Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
						if (geometry->GetnDim() == 3) Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2);
						if ((config->GetKind_Solver() == POTENTIAL_FLOW) || (config->GetKind_Solver() == EULER))
							Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolution->node[iPoint]->GetVelocity2()) / FlowSolution->node[iPoint]->GetSoundSpeed();
						if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS))
							Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolution->GetCSkinFriction(iMarker,iVertex);
						nVertex_Surface++;
					}
				}
	}
	
	double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Press = NULL,
	*Buffer_Receive_Mach = NULL, *Buffer_Receive_SkinFriction = NULL;
	
	if (rank == MASTER_NODE) {
		Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		if (geometry->GetnDim() == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_Press = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_Mach = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_SkinFriction = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
	}
	
	nBuffer_Scalar = MaxLocalVertex_Surface;
	
	/*--- Send the information to the Master node ---*/
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(&Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (geometry->GetnDim() == 3) MPI::COMM_WORLD.Gather(&Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_Press, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Press, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if ((config->GetMPI_Kind_Solver() == POTENTIAL_FLOW) || (config->GetMPI_Kind_Solver() == EULER))
		MPI::COMM_WORLD.Gather(&Buffer_Send_Mach, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Mach, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if ((config->GetMPI_Kind_Solver() == NAVIER_STOKES) || (config->GetMPI_Kind_Solver() == RANS))
		MPI::COMM_WORLD.Gather(&Buffer_Send_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	
	/*--- The master node is the one who writes the surface files ---*/
	if (rank == MASTER_NODE) {
		
		strcpy (cstr, val_filename.c_str());
		sprintf (buffer, ".csv");
		strcat(cstr,buffer);
		
		SurfFlow_file.precision(15);
		SurfFlow_file.open(cstr, ios::out);
		
		/*--- Write the 2D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 2) {
			switch (config->GetMPI_Kind_Solver()) {
				case POTENTIAL_FLOW: case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\"" << endl; break;
				case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\"" << endl; break;
			}
			for (iProcessor = 1; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					PressCoeff = Buffer_Receive_Press[position];
					switch (config->GetMPI_Kind_Solver()) {
						case POTENTIAL_FLOW: case EULER :
							Mach = Buffer_Receive_Mach[position];
							SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << Mach << ", " << yCoord << endl;
							break;
						case NAVIER_STOKES: case RANS:
							SkinFrictionCoeff = Buffer_Receive_SkinFriction[position];
							SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << SkinFrictionCoeff << ", " << yCoord << endl;
							break;
					}
				}
		}
		
		/*--- Write the 3D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 3) {
			switch (config->GetMPI_Kind_Solver()) {
				case POTENTIAL_FLOW: case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\"" << endl; break;
				case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"z_coord\"" << endl; break;
			}
			for (iProcessor = 1; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					zCoord = Buffer_Receive_Coord_z[position];
					PressCoeff = Buffer_Receive_Press[position];
					switch (config->GetMPI_Kind_Solver()) {
						case POTENTIAL_FLOW: case EULER :
							Mach = Buffer_Receive_Mach[position];
							SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << Mach << ", " << yCoord << ", " << zCoord <<endl;
							break;
						case NAVIER_STOKES: case RANS:
							SkinFrictionCoeff = Buffer_Receive_SkinFriction[position];
							SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << SkinFrictionCoeff << ", " << yCoord << ", " << zCoord <<endl;
							break;
					}
				}
		}
	}
	
	if (rank == MASTER_NODE) {
		delete [] Buffer_Receive_Coord_x;
		delete [] Buffer_Receive_Coord_y;
		if (geometry->GetnDim() == 3) delete [] Buffer_Receive_Coord_z;
		delete [] Buffer_Receive_Press;
		delete [] Buffer_Receive_Mach;
		delete [] Buffer_Receive_SkinFriction;
	}
	
	SurfFlow_file.close();
	
#endif
}

void COutput::SetSurfaceCSV_Adjoint_Serial(CConfig *config, CGeometry *geometry, CSolution *AdjSolution, CSolution *FlowSolution, string val_filename, unsigned long iExtIter) {

	unsigned long iPoint, iVertex;
	double *Solution, xCoord, yCoord, zCoord, *Normal, *d, *IntBoundary_Jump;
	unsigned short iMarker, iDim;
	char cstr[200], buffer[50];
	ofstream SurfAdj_file;

	/*--- Write the surface .csv file ---*/
	strcpy (cstr, val_filename.c_str());
	sprintf (buffer, ".csv");
	strcat (cstr, buffer);

	SurfAdj_file.precision(15);
	SurfAdj_file.open(cstr, ios::out);

	if (geometry->GetnDim() == 2) {
		SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\",\"psi_terms\",\"d_dot_n\"" << endl;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					Solution = AdjSolution->node[iPoint]->GetSolution();

					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

					d = AdjSolution->node[iPoint]->GetForceProj_Vector();
					IntBoundary_Jump = AdjSolution->node[iPoint]->GetIntBoundary_Jump();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);

					/*--- Testing surface quantities for rotating adjoint ---*/
					double LHS = 0.0; double RHS = 0.0;
					double *FlowSol = FlowSolution->node[iPoint]->GetSolution();
					for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
						LHS += Normal[iDim]*Solution[iDim+1] + Solution[geometry->GetnDim()+1]*(Normal[iDim]*FlowSol[iDim+1]);
						RHS += d[iDim]*Normal[iDim];
					}

					SurfAdj_file << scientific << iPoint << ", " << AdjSolution->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
							<< Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord <<", "<< LHS << ", " << RHS << endl;

				}
		}
	}

	if (geometry->GetnDim() == 3) {
		SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Solution = AdjSolution->node[iPoint]->GetSolution();
					Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					d = AdjSolution->node[iPoint]->GetForceProj_Vector();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					zCoord = geometry->node[iPoint]->GetCoord(2);

					SurfAdj_file << scientific << iPoint << ", " << AdjSolution->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
							<< Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord << endl;
				}
		}
	}

	SurfAdj_file.close();

}

void COutput::SetSurfaceCSV_Adjoint_Parallel(CConfig *config, CGeometry *geometry, CSolution *AdjSolution, string val_filename, unsigned long iExtIter) {
#ifndef NO_MPI

	int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();
	unsigned short nDim = geometry->GetnDim(), iMarker;
	double *Solution, *Normal, *d, *Coord;
	unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
			MaxLocalVertex_Surface = 0, nBuffer_Scalar;
	unsigned long *Buffer_Receive_nVertex = NULL;
	ofstream SurfAdj_file;

	/*--- Write the surface .csv file ---*/
	if (rank != MASTER_NODE) {
		nLocalVertex_Surface = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
				}
	}
	else {
		nLocalVertex_Surface = 0;
		Buffer_Receive_nVertex = new unsigned long [nProcessor*sizeof(unsigned long)];
	}

	Buffer_Send_nVertex[0] = nLocalVertex_Surface;

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);

	double Buffer_Send_Coord_x[MaxLocalVertex_Surface];
	double Buffer_Send_Coord_y[MaxLocalVertex_Surface];
	double Buffer_Send_Coord_z[MaxLocalVertex_Surface];
	unsigned long Buffer_Send_GlobalPoint[MaxLocalVertex_Surface];
	double Buffer_Send_Sensitivity[MaxLocalVertex_Surface];
	double Buffer_Send_PsiRho[MaxLocalVertex_Surface];
	double Buffer_Send_Phi_x[MaxLocalVertex_Surface];
	double Buffer_Send_Phi_y[MaxLocalVertex_Surface];
	double Buffer_Send_Phi_z[MaxLocalVertex_Surface];
	double Buffer_Send_PsiE[MaxLocalVertex_Surface];

	if (rank != MASTER_NODE) {
		nVertex_Surface = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					if (geometry->node[iPoint]->GetDomain()) {
						Solution = AdjSolution->node[iPoint]->GetSolution();
						Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
						Coord = geometry->node[iPoint]->GetCoord();
						d = AdjSolution->node[iPoint]->GetForceProj_Vector();
						Buffer_Send_GlobalPoint[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
						Buffer_Send_Coord_x[nVertex_Surface] = Coord[0];
						Buffer_Send_Coord_y[nVertex_Surface] = Coord[1];
						Buffer_Send_Sensitivity[nVertex_Surface] =  AdjSolution->GetCSensitivity(iMarker,iVertex);
						Buffer_Send_PsiRho[nVertex_Surface] = Solution[0];
						Buffer_Send_Phi_x[nVertex_Surface] = Solution[1];
						Buffer_Send_Phi_y[nVertex_Surface] = Solution[2];
						if (nDim == 2) Buffer_Send_PsiE[nVertex_Surface] = Solution[3];
						if (nDim == 3) {
							Buffer_Send_Coord_z[nVertex_Surface] = Coord[2];
							Buffer_Send_Phi_z[nVertex_Surface] = Solution[3];
							Buffer_Send_PsiE[nVertex_Surface] = Solution[4];
						}
						nVertex_Surface++;
					}
				}
	}

	double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
			*Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
			*Buffer_Receive_PsiE = NULL;
	unsigned long *Buffer_Receive_GlobalPoint = NULL;

	if (rank == MASTER_NODE) {
		Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		if (nDim == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Surface*sizeof(unsigned long)];
		Buffer_Receive_Sensitivity = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_PsiRho = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_Phi_x = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_Phi_y = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		if (nDim == 3) Buffer_Receive_Phi_z = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
		Buffer_Receive_PsiE = new double [nProcessor*MaxLocalVertex_Surface*sizeof(double)];
	}

	nBuffer_Scalar = MaxLocalVertex_Surface;

	/*--- Send the information to the Master node ---*/
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(&Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (nDim == 3) MPI::COMM_WORLD.Gather(&Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_PsiRho, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_Phi_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_Phi_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (nDim == 3) MPI::COMM_WORLD.Gather(&Buffer_Send_Phi_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_PsiE, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);

	/*--- The master node is the one who writes the surface files ---*/
	if (rank == MASTER_NODE) {
		unsigned long iVertex, GlobalPoint, position;
		char cstr[200], buffer[50];
		ofstream SurfAdj_file;

		strcpy (cstr, val_filename.c_str());
		sprintf (buffer, ".csv");
		strcat(cstr,buffer);
		SurfAdj_file.open(cstr, ios::out);
		SurfAdj_file.precision(15);

		/*--- Write the 2D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 2) {

			SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;

			for (iProcessor = 1; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {

					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					GlobalPoint = Buffer_Receive_GlobalPoint[position];

					SurfAdj_file << scientific << GlobalPoint <<
							", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
							", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
							", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
							", "<< Buffer_Receive_Coord_y[position]  << endl;
				}
		}

		/*--- Write the 3D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 3) {

			SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;

			for (iProcessor = 1; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					GlobalPoint = Buffer_Receive_GlobalPoint[position];

					SurfAdj_file << scientific << GlobalPoint <<
							", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
							", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
							", " << Buffer_Receive_PsiE[position] <<", "<< Buffer_Receive_Coord_x[position] <<
							", "<< Buffer_Receive_Coord_y[position] <<", "<< Buffer_Receive_Coord_z[position] << endl;
				}
		}

	}

	if (rank == MASTER_NODE) {
		delete [] Buffer_Receive_nVertex;
		delete [] Buffer_Receive_Coord_x;
		delete [] Buffer_Receive_Coord_y;
		if (nDim == 3) delete [] Buffer_Receive_Coord_z;
		delete [] Buffer_Receive_Sensitivity;
		delete [] Buffer_Receive_PsiRho;
		delete [] Buffer_Receive_Phi_x;
		delete [] Buffer_Receive_Phi_y;
		if (nDim == 3) delete [] Buffer_Receive_Phi_z;
		delete [] Buffer_Receive_PsiE;
		delete [] Buffer_Receive_GlobalPoint;
	}

	SurfAdj_file.close();

#endif
}

void COutput::SetSurfaceCSV_Linearized_Serial(CConfig *config, CGeometry *geometry, CSolution *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::SetSurfaceCSV_Linearized_Parallel(CConfig *config, CGeometry *geometry, CSolution *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::SetReStart(CConfig *config, CGeometry *geometry, CSolution **solution, string mesh_filename) {
	unsigned short iVar, nVar_First = 0, nVar_Second = 0, FirstIndex = NONE, SecondIndex = NONE;
	unsigned long iPoint;
	ofstream restart_file;
	string filename, AdjExt;
	
	unsigned short KindSolver = config->GetKind_Solver();
	unsigned short KindOF = config->GetKind_ObjFunc();
	bool AdjProblem = ((KindSolver == ADJ_EULER) || (KindSolver == ADJ_NAVIER_STOKES) || (KindSolver == ADJ_RANS) ||
										 (KindSolver == ONE_SHOT_ADJ_EULER) || (KindSolver == ONE_SHOT_ADJ_NAVIER_STOKES));
	
	switch (KindSolver) {
		case POTENTIAL_FLOW: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
		case EULER : case NAVIER_STOKES: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
		case MULTI_SPECIES_NAVIER_STOKES : FirstIndex = PLASMA_SOL; SecondIndex = ELEC_SOL; break;
		case RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
		case TWO_PHASE_FLOW : FirstIndex = FLOW_SOL; SecondIndex = LEVELSET_SOL; break;
		case COMBUSTION : FirstIndex = FLOW_SOL; SecondIndex = COMBUSTION_SOL; break;
		case ELECTRIC_POTENTIAL: FirstIndex = ELEC_SOL; SecondIndex = NONE; break;
		case ADJ_ELECTRIC_POTENTIAL: FirstIndex = ADJELEC_SOL; SecondIndex = NONE; break;
		case LIN_ELECTRIC_POTENTIAL: FirstIndex = LINELEC_SOL; SecondIndex = NONE; break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES : FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; break;
		case ADJ_RANS : FirstIndex = ADJFLOW_SOL; SecondIndex = ADJTURB_SOL; break;
		case ONE_SHOT_ADJ_EULER : case ONE_SHOT_ADJ_NAVIER_STOKES : FirstIndex = FLOW_SOL; SecondIndex = ADJFLOW_SOL; break;
		case LIN_EULER : case LIN_NAVIER_STOKES : FirstIndex = LINFLOW_SOL; SecondIndex = NONE; break;
	}
	nVar_First = solution[FirstIndex]->GetnVar();
	if (SecondIndex != NONE) nVar_Second = solution[SecondIndex]->GetnVar();

	/*--- Create a copy of the filename ---*/
	filename.assign(mesh_filename);
	
	/*--- The adjoint problem requires a particular treatment ---*/
	if (AdjProblem) {
		filename.erase (filename.end()-4, filename.end());
		switch (KindOF) {
			case DRAG_COEFFICIENT: AdjExt = "_cd.dat"; break;
			case LIFT_COEFFICIENT: AdjExt = "_cl.dat"; break;
			case SIDEFORCE_COEFFICIENT: AdjExt = "_csf.dat"; break;
			case PRESSURE_COEFFICIENT: AdjExt = "_cp.dat"; break;
			case MOMENT_X_COEFFICIENT: AdjExt = "_cmx.dat"; break;
			case MOMENT_Y_COEFFICIENT: AdjExt = "_cmy.dat"; break;
			case MOMENT_Z_COEFFICIENT: AdjExt = "_cmz.dat"; break;
			case EFFICIENCY: AdjExt = "_eff.dat"; break;
			case ELECTRIC_CHARGE: AdjExt = "_ec.dat"; break;
			case EQUIVALENT_AREA: AdjExt = "_ea.dat"; break;
			case NEARFIELD_PRESSURE: AdjExt = "_nfp.dat"; break;
		}
		filename.append(AdjExt);
	}
	
	/*--- Open the restart file and write the solution ---*/
	restart_file.open(filename.data(), ios::out);
	restart_file.precision(15);

	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		
		/*--- Index of the point ---*/
		restart_file << iPoint << "\t";
		
		/*--- Solution (first, and second system of equations) ---*/
		for (iVar = 0; iVar < nVar_First; iVar++)
			restart_file << scientific << solution[FirstIndex]->node[iPoint]->GetSolution(iVar) << "\t";

		for (iVar = 0; iVar < nVar_Second; iVar++)
			restart_file << scientific << solution[SecondIndex]->node[iPoint]->GetSolution(iVar) << "\t";
		
		/*--- Residual (first, and second system of equations) ---*/
		for (iVar = 0; iVar < nVar_First; iVar++)
			restart_file << scientific << solution[FirstIndex]->node[iPoint]->GetResidual(iVar) << "\t";
		
		for (iVar = 0; iVar < nVar_Second; iVar++)
			restart_file << scientific << solution[SecondIndex]->node[iPoint]->GetResidual(iVar) << "\t";
		
		restart_file << geometry->node[iPoint]->GetVolume() << "\t";
		restart_file << endl;
	}

	/*--- Close the restart file ---*/
	restart_file.close();

}

void COutput::SetHistory_file(ofstream *ConvHist_file, CConfig *config) {
	char cstr[200], buffer[50];

	/*--- Write file name with extension ---*/
	strcpy (cstr, config->GetConv_FileName().data());
	if (config->GetOutput_FileFormat() == PARAVIEW)  sprintf (buffer, ".csv");
	if (config->GetOutput_FileFormat() == TECPLOT)  sprintf (buffer, ".plt");
	strcat(cstr,buffer);

	ConvHist_file->open(cstr, ios::out);
	ConvHist_file->precision(15);

	char begin[]= "\"Iteration\"";

	char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CPress\",\"CMx\",\"CMy\",\"CMz\",\"CL/CD\",\"CEquivArea\",\"CNearFieldPress\"";
	char adj_coeff[]= ",\"CSens_Geo\",\"CSens_Mach\",\"CSens_AoA\",\"CSens_AoS\"";

	char pot_coeff[]= ",\"CCharge\"";
	char pot_resid[]= ",\"Res_Pot\"";
	char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
	char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";

	char turb_resid[]= ",\"Res_Turb[0]\"";
	char adj_turb_resid[]= ",\"Res_AdjTurb[0]\"";

	char end[]= ",\"Time(min)\"\n";

	if (config->GetOutput_FileFormat() == TECPLOT) {
		char title[]= "TITLE = \"SU2 Simulation\"";
		char variables[]= "VARIABLES = ";
		ConvHist_file[0] << title << endl;
		ConvHist_file[0] << variables;
	}

	switch (config->GetKind_Solver()) {
	case POTENTIAL_FLOW:
		ConvHist_file[0] << begin << flow_coeff << pot_resid << end;
		break;
	case EULER : case NAVIER_STOKES: case TWO_PHASE_FLOW:
		ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
		break;
	case RANS :
		ConvHist_file[0] << begin << flow_coeff << flow_resid << turb_resid << end;
		break;
	case ELECTRIC_POTENTIAL:
		ConvHist_file[0] << begin << pot_coeff << pot_resid << end;
		break;
	case ADJ_EULER : case ADJ_NAVIER_STOKES :
		ConvHist_file[0] << begin << adj_coeff << flow_coeff << adj_flow_resid << end;
		break;
	case ADJ_RANS :
		ConvHist_file[0] << begin << adj_coeff << flow_coeff << adj_flow_resid << adj_turb_resid << end;
		break;
	case ONE_SHOT_ADJ_EULER : case ONE_SHOT_ADJ_NAVIER_STOKES :
		ConvHist_file[0] << begin << flow_coeff << adj_coeff << flow_resid << adj_flow_resid << end;
		break;
	}

	if (config->GetOutput_FileFormat() == TECPLOT) {
		char zone[]= "ZONE T= \"Convergence history\"";
		ConvHist_file[0] << zone << endl;
	}

}

void COutput::SetHistoryDT_Serial(CGeometry **geometry, CConfig *config, CSolution ***solution_container,
		unsigned long iExtIter) {

	unsigned short FinestMesh = config->GetFinestMesh();
	bool write_heads = ((iExtIter % (config->GetWrt_Con_Freq()*1000)) == 0);
	bool write_solution = ((iExtIter % (config->GetWrt_Con_Freq()*1)) == 0);
	bool incompressible = config->GetIncompressible();

	/*--- Write the screen header---*/
	if ((iExtIter == 1) || (write_heads)) {
		cout << endl << " Min DT: " << solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
				". Max DT: " << solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
				". Dual Time step: " << config->GetUnst_TimeStep();
		switch (config->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES: case TWO_PHASE_FLOW: case COMBUSTION:
			if (incompressible) cout << endl << " DualTime Iter" << "      Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			else cout << endl << "  DualTime Iter" << "      Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			break;
		case RANS :
			cout << endl << "  DualTime Iter" << "      Res[Rho]" << "       Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			break;
		}
	}

	/*--- Write the solution on the screen and history file ---*/
	if ((iExtIter != 0) && (write_solution)) {
		switch (config->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES: case TWO_PHASE_FLOW: case COMBUSTION:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
			if (!incompressible) {
				cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(geometry[0]->GetnDim()+1));
			}
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
			cout << endl;
			break;
		case RANS :
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[FinestMesh][TURB_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
			cout << endl;
			break;
		}
	}
	cout.unsetf(ios::fixed);
}

void COutput::SetHistory_Serial(ofstream *ConvHist_file, CGeometry **geometry, CConfig *config, CSolution ***solution_container,
		unsigned long iExtIter, unsigned long timeused) {
	char begin[200], flow_coeff[200], adj_coeff[200], lin_coeff[200], flow_resid[200], adj_flow_resid[200], lin_flow_resid[200],
	turb_resid[200], adj_turb_resid[200], end[200];

	double timeiter = double(timeused)/double(iExtIter+1), dummy = 0.0;
	bool write_heads = ((iExtIter % (config->GetWrt_Con_Freq()*20)) == 0);
	bool incompressible = config->GetIncompressible();
	unsigned short FinestMesh = config->GetFinestMesh();

	/*--- Write the begining of the history file ---*/
	sprintf (begin, "%12d", int(iExtIter));

	/*--- Write the end of the history file ---*/
	sprintf (end, ", %12.10f\n", double(timeused)/(CLOCKS_PER_SEC*60.0));

	/*--- Write the solution and residual of the history file ---*/
	switch (config->GetKind_Solver()) {
	case POTENTIAL_FLOW:
		sprintf (flow_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CSideForce(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CPress(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMx(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMy(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMz(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CEff());
		sprintf (flow_resid, ", %12.10f, 0.0, 0.0, 0.0, 0.0",
				log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)));
		break;
	case EULER : case NAVIER_STOKES: case RANS: case TWO_PHASE_FLOW: case COMBUSTION:
		sprintf (flow_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CSideForce(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CPress(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMx(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMy(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMz(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CEff(),
				config->GetWeightCd()*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CEquivArea(),
				config->GetWeightCd()*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CNearFieldPress());

		if (geometry[FinestMesh]->GetnDim() == 2) {
			if (incompressible) {
				sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
						dummy, dummy );
			}
			else {
				sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(3)),
						dummy );
			}
		}
		else {
			if (incompressible) {
				sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(3)),
						dummy );
			}
			else {
				sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(3)),
						log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(4)));
			}
		}

		if (config->GetKind_Solver() == RANS)
			sprintf (turb_resid, ", %12.10f", log10 (solution_container[FinestMesh][TURB_SOL]->GetRes_Max(0)));
		break;

	case MULTI_SPECIES_NAVIER_STOKES:
		sprintf (flow_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CLift(),
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CDrag(),
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CSideForce(),
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CPress(),
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CMx(),
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CMy(),
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CMz(),
				solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CEff(),
				config->GetWeightCd()*solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CEquivArea(),
				config->GetWeightCd()*solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][PLASMA_SOL]->GetTotal_CNearFieldPress());
		if (geometry[FinestMesh]->GetnDim() == 2) {
			sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(3)),
					dummy	 );		}
		else {
			sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(3)),
					log10 (solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(4))); }

		break;

	case ADJ_EULER: case ONE_SHOT_ADJ_EULER: case ADJ_NAVIER_STOKES: case ONE_SHOT_ADJ_NAVIER_STOKES: case ADJ_RANS:
		sprintf (flow_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CSideForce(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CPress(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMx(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMy(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMz(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CEff(),
				config->GetWeightCd()*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CEquivArea(),
				config->GetWeightCd()*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CNearFieldPress());
		sprintf (adj_coeff, ", %12.10f, %12.10f, %12.10f, 0.0",
				solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Geo(),
				solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Mach(),
				solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_AoA());
		if (geometry[FinestMesh]->GetnDim() == 2) {
			sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0",
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(3)));
			sprintf (adj_flow_resid, ", %17.15f, %17.15f, %17.15f, %17.15f, 0.0",
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(3))); }
		else {
			sprintf (flow_resid, ", 12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(3)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(4)));
			sprintf (adj_flow_resid, ", 12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(3)),
					log10 (solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(4))); }
		if (config->GetKind_Solver() == ADJ_RANS) {
			sprintf (turb_resid, ", %12.10f",
					log10 (solution_container[FinestMesh][TURB_SOL]->GetRes_Max(0)));
			sprintf (adj_turb_resid, ", %12.10f",
					log10 (solution_container[FinestMesh][ADJTURB_SOL]->GetRes_Max(0)));
		}
		break;
	case ELECTRIC_POTENTIAL:
		sprintf (flow_coeff, ", %12.10f",
				solution_container[FinestMesh][ELEC_SOL]->GetTotal_CCharge());
		sprintf (flow_resid, ", %12.10f",
				log10 (solution_container[FinestMesh][ELEC_SOL]->GetRes_Max(0)));
		break;
	case ADJ_ELECTRIC_POTENTIAL:
		sprintf (flow_coeff, ", 0.0, 0.0, 0.0, 0.0, 0.0");
		sprintf (adj_coeff, ", 0.0, 0.0, 0.0, 0.0");
		sprintf (flow_resid, ", %12.10f, 0.0, 0.0, 0.0, 0.0",
				log10 (solution_container[FinestMesh][ELEC_SOL]->GetRes_Max(0)));
		sprintf (adj_flow_resid, ", %12.10f, 0.0, 0.0, 0.0, 0.0",
				log10 (solution_container[FinestMesh][ADJELEC_SOL]->GetRes_Max(0)));
		break;
	case LIN_EULER: case LIN_NAVIER_STOKES:
		sprintf (flow_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CSideForce(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CPress(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMx(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMy(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CMz(),
				solution_container[FinestMesh][FLOW_SOL]->GetTotal_CEff(),
				config->GetWeightCd()*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CEquivArea(),
				config->GetWeightCd()*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag()
				+ (1.0-config->GetWeightCd())*solution_container[FinestMesh][FLOW_SOL]->GetTotal_CNearFieldPress());
		sprintf (lin_coeff, ", %12.10f, %12.10f, 0.0, 0.0",
				solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaDrag(),
				solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaLift());
		if (geometry[FinestMesh]->GetnDim() == 2) {
			sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0",
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(3)));
			sprintf (lin_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0",
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(3))); }
		else {
			sprintf (flow_resid, ", 12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(3)),
					log10 (solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(4)));
			sprintf (lin_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(0)),
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(1)),
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(2)),
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(3)),
					log10 (solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(4))); }
		break;
	}

	/*--- Write the screen header---*/
	/*	if (write_heads)
		cout << endl << " Min DT: " << solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Physical Time: " << solution_container[MESH_0][FLOW_SOL]->GetSum_Delta_Time();
	 */
	switch (config->GetKind_Solver()) {
	case POTENTIAL_FLOW :
		if (write_heads) cout << endl << " Iter" << "  Time(s)" << "      Res[Pot]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
		break;
	case EULER : case NAVIER_STOKES: case TWO_PHASE_FLOW : case COMBUSTION :
		if (write_heads) {
			if (incompressible) cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			else cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
		}
		break;
	case MULTI_SPECIES_NAVIER_STOKES :
			if (write_heads) {
				if (config->GetKind_GasModel() == ARGON)
					cout << endl << " Iter" << "   Time(s)" << "        Res[r1]" << "         Res[m1]" << "	   Res[n1]" << "       Res[e1)]" <<  "	 Res[r2]" << "		Res[m2]" << "	Res[n2]" << "	Res[e2]" << "	   Res[r3]" << "	   Res[m3]" << "      Res[n3]" << "      Res[e3]"<< endl;
				if (config->GetKind_GasModel() == AIR7)
					cout << endl << " Iter" << "   Time(s)" << "        Res[Rho1]" << "         Res[E]" << "	   Heat_Transfer" << endl;
			}
				
				break;
	case RANS :
		if (write_heads) cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
		break;
	case ELECTRIC_POTENTIAL :
		if (write_heads) cout << endl << " Iter" << "  Time(s)" << "      Res[Pot]" << "   CCharge(Total)" << endl;
		break;
	case ADJ_EULER : case ADJ_NAVIER_STOKES :
		if (write_heads) cout << endl << " Iter" << "   Time(s)" << "   Res[Psi_Rho]" << "     Res[Psi_E]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
		break;
	case ADJ_RANS :
		if (write_heads) cout << endl << "  Iter" << "    Time(s)" << "     Res[Psi_Rho]" << "      Res[Psi_nu]" << "       CSens_Geo" << endl;
		break;
	case ADJ_ELECTRIC_POTENTIAL :
		if (write_heads) cout << endl << " Iter" << "  Time(s)" << "      Res[Psi_Pot]" << endl;
		break;
	case LIN_ELECTRIC_POTENTIAL :
		if (write_heads) cout << endl << " Iter" << "  Time(s)" << "      Res[Delta_Pot]" << endl;
		break;
	case ONE_SHOT_ADJ_EULER : case ONE_SHOT_ADJ_NAVIER_STOKES :
		if (write_heads) cout << endl << " Iter" << "  Time(s)" << "      Res[Rho]" << "  Res[Psi_Rho]" << "           CDrag" << "       CSens_Geo" << endl;
		break;
	case LIN_EULER : case LIN_NAVIER_STOKES :
		if (write_heads) cout << endl << " Iter" << "   Time(s)" << "  Res[Delta_Rho]" << "  Res[Delta_E]" << "   CDeltaLift" << "   CDeltaDrag"<< endl;
		break;
	}

	/*--- Write the solution on the screen and history file ---*/
	if (iExtIter != 0) {
		switch (config->GetKind_Solver()) {
		case POTENTIAL_FLOW :
			ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
			cout << endl;
			break;
		case EULER : case NAVIER_STOKES: case TWO_PHASE_FLOW: case COMBUSTION:
			ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
			if (!incompressible) {
				cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(geometry[0]->GetnDim()+1));
			}
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
			cout << endl;
			break;

			case MULTI_SPECIES_NAVIER_STOKES:
				ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				if (config->GetKind_GasModel() == ARGON) {
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(0));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(1));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(2));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(3));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(4));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(5));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(6));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(7));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(8));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(9));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(10));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(11));
				}
				if (config->GetKind_GasModel() == AIR7) {
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(0));
					cout.width(14); cout << log10(solution_container[FinestMesh][PLASMA_SOL]->GetRes_Max(3));					
				}
			cout << endl;
			break;
		case RANS :
			ConvHist_file[0] << begin << flow_coeff << flow_resid << turb_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[FinestMesh][TURB_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
			cout << endl;
			break;
		case ELECTRIC_POTENTIAL :
			ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(14); cout << log10(solution_container[FinestMesh][ELEC_SOL]->GetRes_Max(0));
			cout.width(17); cout << solution_container[FinestMesh][ELEC_SOL]->GetTotal_CCharge();
			cout << endl;
			break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES :
			ConvHist_file[0] << begin << adj_coeff << flow_coeff << adj_flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(15); cout << log10(solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0));
			cout.width(15); cout << log10(solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(geometry[0]->GetnDim()+1));
			cout.width(14); cout << solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Geo();
			cout.width(14); cout << solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Mach();
			cout << endl;
			break;
		case ADJ_RANS :
			ConvHist_file[0] << begin << adj_coeff << flow_coeff << adj_flow_resid << adj_turb_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(6); cout << iExtIter;
			cout.width(11); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(17); cout << log10(solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0));
			cout.width(17); cout << log10(solution_container[FinestMesh][ADJTURB_SOL]->GetRes_Max(0));
			cout.width(16); cout << solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Geo();
			cout << endl;
			break;
		case ADJ_ELECTRIC_POTENTIAL :
			ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(18); cout << log10(solution_container[FinestMesh][ADJELEC_SOL]->GetRes_Max(0));
			cout << endl;
			break;
		case LIN_ELECTRIC_POTENTIAL :
			ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(18); cout << log10(solution_container[FinestMesh][LINELEC_SOL]->GetRes_Max(0));
			cout << endl;
			break;
		case ONE_SHOT_ADJ_EULER : case ONE_SHOT_ADJ_NAVIER_STOKES :
			ConvHist_file[0] << begin << flow_coeff << adj_coeff << flow_resid << adj_flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(14); cout << log10(solution_container[FinestMesh][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[FinestMesh][ADJFLOW_SOL]->GetRes_Max(0));
			cout.width(16); cout << solution_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
			cout.width(16); cout << solution_container[FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Geo();
			cout << endl;
			break;
		case LIN_EULER : case LIN_NAVIER_STOKES :
			ConvHist_file[0] << begin << lin_coeff << flow_coeff << lin_flow_resid << end;
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(5); cout << iExtIter;
			cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
			cout.width(16); cout << log10(solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[FinestMesh][LINFLOW_SOL]->GetRes_Max(geometry[0]->GetnDim()+1));
			cout.width(13); cout << solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaLift();
			cout.width(13); cout << solution_container[FinestMesh][LINFLOW_SOL]->GetTotal_CDeltaDrag();
			cout << endl;
			break;
		}
	}
	cout.unsetf(ios::fixed);
}

void COutput::SetHistory_Parallel(ofstream *ConvHist_file, CGeometry **geometry, CSolution ***solution_container, CConfig *config,
																	CIntegration **integration, unsigned long iExtIter, int rank, int size, unsigned long timeused) {
#ifndef NO_MPI
	
	double *sbuf_residual_flow, *sbuf_residual_turbulent, *sbuf_force, *sbuf_time, *rbuf_residual_flow = NULL, *rbuf_residual_turbulent = NULL,
	*rbuf_force = NULL, *rbuf_time = NULL, monitor = 0, *sbuf_residual_adjoint = NULL, *rbuf_residual_adjoint = NULL,
	*sbuf_residual_linearized = NULL, *rbuf_residual_linearized = NULL;
	unsigned short iVar, buf_convergence = 0;
	char begin[200], flow_coeff[200], adj_coeff[200], lin_coeff[200], flow_resid[200], adj_resid[200], lin_resid[200], turb_resid[200], end[200];
	
	double timeiter = double(timeused)/double(iExtIter+1);
	bool write_heads = ((iExtIter % (config->GetWrt_Con_Freq()*20)) == 0);
	bool incompressible = config->GetIncompressible();
	
	/*--- It is important to improve this, the problem is that the master
	 node does not know the number of variables ---*/
	unsigned short nVar_Flow = geometry[MESH_0]->GetnDim()+2;
	unsigned short nVar_Adj = geometry[MESH_0]->GetnDim()+2;
	unsigned short nVar_Lin = geometry[MESH_0]->GetnDim()+2;
	unsigned short nVar_Turb = 1;
	unsigned short nVar_Force = 10;
	unsigned short nVar_Time = 2;
	
	/*--- Allocate memory for send buffer ---*/
	sbuf_residual_flow = new double[nVar_Flow];
	sbuf_residual_adjoint = new double[nVar_Adj];
	sbuf_residual_linearized = new double[nVar_Lin];
	sbuf_residual_turbulent = new double[nVar_Turb];
	sbuf_force = new double[nVar_Force];
	sbuf_time = new double[nVar_Time];
	
	/*--- Initializate master node ---*/
	if (rank == MASTER_NODE) {
		rbuf_residual_flow = new double[nVar_Flow];
		rbuf_residual_adjoint = new double[nVar_Adj];
		rbuf_residual_linearized = new double[nVar_Lin];
		rbuf_residual_turbulent = new double[nVar_Turb];
		rbuf_force = new double[nVar_Force];
		rbuf_time = new double[nVar_Time];
		for (iVar = 0; iVar < nVar_Flow; iVar++) {
			rbuf_residual_flow[iVar] = 0.0;
			sbuf_residual_flow[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Adj; iVar++) {
			rbuf_residual_adjoint[iVar] = 0.0;
			sbuf_residual_adjoint[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Lin; iVar++) {
			rbuf_residual_linearized[iVar] = 0.0;
			sbuf_residual_linearized[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Turb; iVar++) {
			rbuf_residual_turbulent[iVar] = 0.0;
			sbuf_residual_turbulent[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Force; iVar++) {
			rbuf_force[iVar] = 0.0;
			sbuf_force[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_Time; iVar++) {
			rbuf_time[iVar] = 0.0;
			sbuf_time[iVar] = 0.0;
		}
	}
	
	/*--- Write information from slaves nodes ---*/
	if (rank != MASTER_NODE) {
		
		switch (config->GetKind_Solver()) {
			case EULER : case NAVIER_STOKES:  case TWO_PHASE_FLOW: case COMBUSTION:
				for (iVar = 0; iVar < nVar_Flow; iVar++)
					sbuf_residual_flow[iVar] = log10 (solution_container[MESH_0][FLOW_SOL]->GetRes_Max(iVar));
				sbuf_force[0] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
				sbuf_force[1] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag();
				sbuf_force[2] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CSideForce();
				sbuf_force[3] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CPress();
				sbuf_force[4] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMx();
				sbuf_force[5] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMy();
				sbuf_force[6] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMz();
				sbuf_force[7] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CEff();
				sbuf_force[8] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CEquivArea(),
				sbuf_force[9] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CNearFieldPress();
				break;
			case RANS :
				for (iVar = 0; iVar < nVar_Flow; iVar++)
					sbuf_residual_flow[iVar] = log10 (solution_container[MESH_0][FLOW_SOL]->GetRes_Max(iVar));
				for (iVar = 0; iVar < nVar_Turb; iVar++)
					sbuf_residual_turbulent[iVar] = log10 (solution_container[MESH_0][TURB_SOL]->GetRes_Max(iVar));
				sbuf_force[0] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
				sbuf_force[1] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag();
				sbuf_force[2] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CSideForce();
				sbuf_force[3] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CPress();
				sbuf_force[4] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMx();
				sbuf_force[5] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMy();
				sbuf_force[6] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMz();
				sbuf_force[7] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CEff();
				sbuf_force[8] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CEquivArea();
				sbuf_force[9] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CNearFieldPress();
				break;
			case ADJ_EULER : case ADJ_NAVIER_STOKES :
				for (iVar = 0; iVar < nVar_Adj; iVar++)
					sbuf_residual_adjoint[iVar] = log10 (solution_container[MESH_0][ADJFLOW_SOL]->GetRes_Max(iVar));
				sbuf_force[0] = solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Geo();
				sbuf_force[1] = solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Mach();
				break;
			case LIN_EULER : case LIN_NAVIER_STOKES :
				for (iVar = 0; iVar < nVar_Lin; iVar++)
					sbuf_residual_linearized[iVar] = log10 (solution_container[MESH_0][LINFLOW_SOL]->GetRes_Max(iVar));
				sbuf_force[0] = solution_container[MESH_0][LINFLOW_SOL]->GetTotal_CDeltaLift();
				sbuf_force[1] = solution_container[MESH_0][LINFLOW_SOL]->GetTotal_CDeltaDrag();
				break;
		}
		sbuf_time[0] = timeiter;
		sbuf_time[1] = timeused;
	}
	
	/*--- Send/Receive information ---*/
	switch (config->GetMPI_Kind_Solver()) {
		case EULER : case NAVIER_STOKES: case TWO_PHASE_FLOW:
			MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
			break;
		case RANS :
			MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
			MPI::COMM_WORLD.Reduce(sbuf_residual_turbulent, rbuf_residual_turbulent, nVar_Turb, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
			break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES :
			MPI::COMM_WORLD.Reduce(sbuf_residual_adjoint, rbuf_residual_adjoint, nVar_Adj, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
			break;
		case LIN_EULER : case LIN_NAVIER_STOKES :
			MPI::COMM_WORLD.Reduce(sbuf_residual_linearized, rbuf_residual_linearized, nVar_Lin, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
			break;
	}
	MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, nVar_Force, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Reduce(sbuf_time, rbuf_time, nVar_Time, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
	MPI::COMM_WORLD.Barrier();
	
	/*--- Review convergence criteria ---*/
	switch (config->GetKind_Solver()) {
		case POTENTIAL_FLOW: case EULER: case NAVIER_STOKES: case RANS: case TWO_PHASE_FLOW: case COMBUSTION:
			integration[FLOW_SOL]->SetConvergence(false); break;
		case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
			integration[ADJFLOW_SOL]->SetConvergence(false); break;
		case LIN_EULER: case LIN_NAVIER_STOKES:
			integration[LINFLOW_SOL]->SetConvergence(false); break;
		case ELECTRIC_POTENTIAL:
			integration[ELEC_SOL]->SetConvergence(false); break;
		case ADJ_ELECTRIC_POTENTIAL:
			integration[ADJELEC_SOL]->SetConvergence(false); break;
		case ONE_SHOT_ADJ_EULER: case ONE_SHOT_ADJ_NAVIER_STOKES:
			integration[FLOW_SOL]->SetConvergence(false); break;
	}
	
	/*--- Write result using the master node ---*/
	if (rank == MASTER_NODE) {
		
		/*--- Write the screen header and evaluate the convergence monitor ---*/
		switch (config->GetMPI_Kind_Solver()) {
			case EULER : case NAVIER_STOKES:  case TWO_PHASE_FLOW: case COMBUSTION:
				if (write_heads) {
					if (incompressible)
						cout << endl << " Iter" << "  MaxTime(s)" << "     Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
					else
						cout << endl << " Iter" << "  MaxTime(s)" << "     Res[Rho]" << "    Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				}
				
				if (config->GetConvCriteria() == CAUCHY) {
					if (config->GetCauchy_Func_Flow() == LIFT_COEFFICIENT) monitor = rbuf_force[0];
					if (config->GetCauchy_Func_Flow() == DRAG_COEFFICIENT) monitor = rbuf_force[1];
					if (config->GetCauchy_Func_Flow() == NEARFIELD_PRESSURE) monitor = rbuf_force[9];
				}
				if (config->GetConvCriteria() == RESIDUAL) {
					monitor = rbuf_residual_flow[0];
				}
				break;
			case RANS :
				if (write_heads) cout << endl << " Iter" << "  MaxTime(s)" << "     Res[Rho]" << "      Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				if (config->GetConvCriteria() == CAUCHY) {
					if (config->GetCauchy_Func_Flow() == LIFT_COEFFICIENT) monitor = rbuf_force[0];
					if (config->GetCauchy_Func_Flow() == DRAG_COEFFICIENT) monitor = rbuf_force[1];
					if (config->GetCauchy_Func_Flow() == NEARFIELD_PRESSURE) monitor = rbuf_force[9];
				}
				if (config->GetConvCriteria() == RESIDUAL) {
					monitor = rbuf_residual_flow[0];
				}
				break;
			case ADJ_EULER : case ADJ_NAVIER_STOKES :
				if (write_heads) cout << endl << " Iter" << "  MaxTime(s)" << "  Res[Psi_Rho]" << "    Res[Psi_E]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
				if (config->GetConvCriteria() == CAUCHY) {
					if (config->GetCauchy_Func_AdjFlow() == SENS_GEOMETRY) monitor = rbuf_force[0];
					if (config->GetCauchy_Func_AdjFlow() == SENS_MACH) monitor = rbuf_force[1];
				}
				if (config->GetConvCriteria() == RESIDUAL) {
					monitor = rbuf_residual_adjoint[0];
				}
				break;
			case ADJ_RANS :
				if (write_heads) cout << endl << " Iter" << "  MaxTime(s)" << "  Res[Psi_Rho]" << "   Res[Psi_nu]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
				if (config->GetConvCriteria() == CAUCHY) {
					if (config->GetCauchy_Func_AdjFlow() == SENS_GEOMETRY) monitor = rbuf_force[0];
					if (config->GetCauchy_Func_AdjFlow() == SENS_MACH) monitor = rbuf_force[1];
				}
				if (config->GetConvCriteria() == RESIDUAL) {
					monitor = rbuf_residual_adjoint[0];
				}
				break;
			case LIN_EULER : case LIN_NAVIER_STOKES :
				if (write_heads) cout << endl << " Iter" << "  MaxTime(s)" << "  Res[Delta_Rho]" << "  Res[Delta_E]" << "   CDeltaLift" << "   CDeltaDrag"<< endl;
				if (config->GetConvCriteria() == CAUCHY) {
					if (config->GetCauchy_Func_LinFlow() == DELTA_LIFT_COEFFICIENT) monitor = rbuf_force[0];
					if (config->GetCauchy_Func_LinFlow() == DELTA_DRAG_COEFFICIENT) monitor = rbuf_force[1];
				}
				if (config->GetConvCriteria() == RESIDUAL) {
					monitor = rbuf_residual_linearized[0];
				}
				break;
		}
		
		/*-- Compute the convergence criteria --*/
		integration[NO_SOLVER]->Convergence_Monitoring(geometry[MESH_0], config, iExtIter, monitor);
		if (integration[NO_SOLVER]->GetConvergence()) buf_convergence = 1;
		else buf_convergence = 0;
		
		/*--- Write the solution on the screen and the history file ---*/
		if (iExtIter != 0) {
			switch (config->GetMPI_Kind_Solver()) {
				case EULER : case NAVIER_STOKES: case TWO_PHASE_FLOW: case COMBUSTION:
					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(12); cout << rbuf_time[0]/CLOCKS_PER_SEC;
					cout.width(13); cout << rbuf_residual_flow[0];
					if (!incompressible) { cout.width(13); cout << rbuf_residual_flow[3]; }
					cout.width(15); cout << rbuf_force[0]; cout.width(15); cout << rbuf_force[1];
					cout << endl;
					
					/*-- Writes the history file for the global computation --*/
					sprintf (begin, "%12d", int(iExtIter));
					sprintf (flow_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
									 rbuf_force[0], rbuf_force[1], rbuf_force[2], rbuf_force[3], rbuf_force[4], rbuf_force[5],
									 rbuf_force[6], rbuf_force[7],
									 config->GetWeightCd()*rbuf_force[1] + (1.0-config->GetWeightCd())*rbuf_force[8]/ double(MPI::COMM_WORLD.Get_size()-1),
									 config->GetWeightCd()*rbuf_force[1] + (1.0-config->GetWeightCd())*rbuf_force[9]);
					if (geometry[MESH_0]->GetnDim() == 2) {
						sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0",
										 rbuf_residual_flow[0], rbuf_residual_flow[1], rbuf_residual_flow[2], rbuf_residual_flow[3]);
					}
					if (geometry[MESH_0]->GetnDim() == 3) {
						sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
										 rbuf_residual_flow[0], rbuf_residual_flow[1], rbuf_residual_flow[2], rbuf_residual_flow[3], rbuf_residual_flow[4]);
					}
					sprintf (end, ", %12.10f\n", double(rbuf_time[1])/(CLOCKS_PER_SEC*60.0));
					ConvHist_file[0] << begin << flow_coeff << flow_resid << end;
					break;
					
				case RANS:
					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(12); cout << rbuf_time[0]/CLOCKS_PER_SEC;
					cout.width(13); cout << rbuf_residual_flow[0];
					cout.width(13); cout << rbuf_residual_turbulent[0];
					cout.width(15); cout << rbuf_force[0]; cout.width(15); cout << rbuf_force[1];
					cout << endl;
					
					/*-- Writes the history file for the global computation --*/
					sprintf (begin, "%12d", int(iExtIter));
					sprintf (flow_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
									 rbuf_force[0], rbuf_force[1], rbuf_force[2], rbuf_force[3], rbuf_force[4], rbuf_force[5],
									 rbuf_force[6], rbuf_force[7],
									 config->GetWeightCd()*rbuf_force[1] + (1.0-config->GetWeightCd())*rbuf_force[8]/ double(MPI::COMM_WORLD.Get_size()-1),
									 config->GetWeightCd()*rbuf_force[1] + (1.0-config->GetWeightCd())*rbuf_force[9]);
					if (geometry[MESH_0]->GetnDim() == 2) {
						sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0",
										 rbuf_residual_flow[0], rbuf_residual_flow[1], rbuf_residual_flow[2], rbuf_residual_flow[3]);
					}
					if (geometry[MESH_0]->GetnDim() == 3) {
						sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
										 rbuf_residual_flow[0], rbuf_residual_flow[1], rbuf_residual_flow[2], rbuf_residual_flow[3],
										 rbuf_residual_flow[4]);
					}
					sprintf (turb_resid, ", %12.10f", rbuf_residual_turbulent[0]);
					sprintf (end, ", %12.10f\n", double(rbuf_time[1])/(CLOCKS_PER_SEC*60.0));
					ConvHist_file[0] << begin << flow_coeff << flow_resid << turb_resid << end;
					break;
					
				case ADJ_EULER : case ADJ_NAVIER_STOKES :
					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(12); cout << rbuf_time[0]/CLOCKS_PER_SEC;
					cout.width(14); cout << rbuf_residual_adjoint[0];
					cout.width(14); cout << rbuf_residual_adjoint[3];
					cout.width(14); cout << rbuf_force[0];
					cout.width(14); cout << rbuf_force[1];
					cout << endl;
					
					/*-- writes the history file for the global computation --*/
					sprintf (begin, "%12d", int(iExtIter));
					sprintf (adj_coeff, ", %12.10f, %12.10f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0", rbuf_force[0], rbuf_force[1]);
					if (geometry[MESH_0]->GetnDim() == 2) {
						sprintf (adj_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0",
										 rbuf_residual_adjoint[0], rbuf_residual_adjoint[1], rbuf_residual_adjoint[2], rbuf_residual_adjoint[3]);
					}
					if (geometry[MESH_0]->GetnDim() == 3) {
						sprintf (adj_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", rbuf_residual_adjoint[0], rbuf_residual_adjoint[1],
										 rbuf_residual_adjoint[2], rbuf_residual_adjoint[3], rbuf_residual_adjoint[4]);
					}
					sprintf (end, ", %12.10f\n", double(rbuf_time[1])/(CLOCKS_PER_SEC*60.0));
					ConvHist_file[0] << begin << adj_coeff << adj_resid << end;
					break;
					
				case LIN_EULER : case LIN_NAVIER_STOKES :
					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(12); cout << rbuf_time[0]/CLOCKS_PER_SEC;
					cout.width(14); cout << rbuf_residual_linearized[0];
					cout.width(14); cout << rbuf_residual_linearized[3];
					cout.width(13); cout << rbuf_force[0];
					cout.width(12); cout << rbuf_force[1];
					cout << endl;
					
					/*-- writes the history file for the global computation --*/
					sprintf (begin, "%12d", int(iExtIter));
					sprintf (lin_coeff, ", %12.10f, %12.10f, 0.0, 0.0, 0.0", rbuf_force[0], rbuf_force[1]);
					if (geometry[MESH_0]->GetnDim() == 2) {
						sprintf (lin_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0",
										 rbuf_residual_linearized[0], rbuf_residual_linearized[1], rbuf_residual_linearized[2], rbuf_residual_linearized[3]);
					}
					if (geometry[MESH_0]->GetnDim() == 3) {
						sprintf (lin_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
										 rbuf_residual_linearized[0], rbuf_residual_linearized[1], rbuf_residual_linearized[2], rbuf_residual_linearized[3],
										 rbuf_residual_linearized[4]);
					}
					sprintf (end, ", %12.10f\n", double(rbuf_time[1])/(CLOCKS_PER_SEC*60.0));
					ConvHist_file[0] << begin << lin_coeff << lin_resid << end;
					break;
			}
		}
	}
	
	MPI::COMM_WORLD.Bcast(&buf_convergence, 1, MPI::SHORT, MASTER_NODE);
	MPI::COMM_WORLD.Barrier();
	
	switch (config->GetKind_Solver()) {
		case POTENTIAL_FLOW: case EULER: case NAVIER_STOKES: case RANS: case TWO_PHASE_FLOW: case COMBUSTION:
			if (buf_convergence == 1) integration[FLOW_SOL]->SetConvergence(true); break;
		case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
			if (buf_convergence == 1) integration[ADJFLOW_SOL]->SetConvergence(true); break;
		case LIN_EULER: case LIN_NAVIER_STOKES:
			if (buf_convergence == 1) integration[LINFLOW_SOL]->SetConvergence(true); break;
		case ELECTRIC_POTENTIAL:
			if (buf_convergence == 1) integration[ELEC_SOL]->SetConvergence(true); break;
		case ADJ_ELECTRIC_POTENTIAL:
			if (buf_convergence == 1) integration[ADJELEC_SOL]->SetConvergence(true); break;
		case ONE_SHOT_ADJ_EULER: case ONE_SHOT_ADJ_NAVIER_STOKES:
			if (buf_convergence == 1) integration[FLOW_SOL]->SetConvergence(true); break;
	}
	
	delete [] sbuf_residual_flow;
	delete [] sbuf_residual_adjoint;
	delete [] sbuf_residual_linearized;
	delete [] sbuf_residual_turbulent;
	delete [] sbuf_force;
	delete [] sbuf_time;
	
	if (rank == MASTER_NODE) {
		delete [] rbuf_residual_flow;
		delete [] rbuf_residual_adjoint;
		delete [] rbuf_residual_linearized;
		delete [] rbuf_residual_turbulent;
		delete [] rbuf_force;
		delete [] rbuf_time;
	}
	
#endif
}

void COutput::SetResult_Files(CSolution ***solution_container, CGeometry *geometry,
															CConfig *config, unsigned long iExtIter) {
	
	if (config->GetKind_Solver() != NO_SOLVER) {

		switch (config->GetKind_Solver()) {
				
			case POTENTIAL_FLOW:
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				break;
				
			case EULER : case NAVIER_STOKES:
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				SetSurface_Flow(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Flow_Serial(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_FlowFileName());
				break;
				
			case MULTI_SPECIES_NAVIER_STOKES :
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_FlowFileName());
				break;
				
			case RANS :
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				SetSurface_Flow(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Flow_Serial(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_FlowFileName());
				break;
				
			case TWO_PHASE_FLOW :
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				SetSurface_Flow(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Flow_Serial(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_FlowFileName());
				break;
				
			case COMBUSTION :
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				SetSurface_Flow(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Flow_Serial(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_FlowFileName());
				break;
				
			case ELECTRIC_POTENTIAL:
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_FlowFileName());
				break;
				
			case ADJ_ELECTRIC_POTENTIAL:
				SetDomain_Adjoint(config, geometry, solution_container, config->GetAdj_FileName(), iExtIter);
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_AdjFileName());
				break;
				
			case LIN_ELECTRIC_POTENTIAL:
				SetLinearized_Variables(config, geometry, solution_container, config->GetLin_FileName(), iExtIter);
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_LinFileName());
				break;
				
			case ADJ_EULER : case ADJ_NAVIER_STOKES :
				SetDomain_Adjoint(config, geometry, solution_container, config->GetAdj_FileName(), iExtIter);
				SetSurface_Adjoint(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], solution_container[MESH_0][FLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Adjoint_Serial(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], solution_container[MESH_0][FLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_AdjFileName());
				break;
				
			case ADJ_RANS :
				SetDomain_Adjoint(config, geometry, solution_container, config->GetAdj_FileName(), iExtIter);
				SetSurface_Adjoint(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], solution_container[MESH_0][FLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Adjoint_Serial(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], solution_container[MESH_0][FLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_AdjFileName());
				break;
				
			case ONE_SHOT_ADJ_EULER : case ONE_SHOT_ADJ_NAVIER_STOKES :
				SetDomain_Flow(config, geometry, solution_container[MESH_0], config->GetFlow_FileName(), iExtIter);
				SetSurface_Flow(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
				SetDomain_Adjoint(config, geometry, solution_container, config->GetAdj_FileName(), iExtIter);
				SetSurface_Adjoint(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], solution_container[MESH_0][FLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Flow_Serial(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
				SetSurfaceCSV_Adjoint_Serial(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], solution_container[MESH_0][FLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_AdjFileName());
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_FlowFileName());
				break;
				
			case LIN_EULER : case LIN_NAVIER_STOKES :
				SetLinearized_Variables(config, geometry, solution_container, config->GetLin_FileName(), iExtIter);
				SetSurface_Linearized(config, geometry, solution_container[MESH_0][LINFLOW_SOL], config->GetSurfLinCoeff_FileName(), iExtIter);
#ifdef NO_MPI
				SetSurfaceCSV_Linearized_Serial(config, geometry, solution_container[MESH_0][LINFLOW_SOL], config->GetSurfLinCoeff_FileName(), iExtIter);
#endif
				SetReStart(config, geometry, solution_container[MESH_0], config->GetReStart_LinFileName());
				break;
		}
	}
	
#ifndef NO_MPI
	switch (config->GetMPI_Kind_Solver()) {
			
		case EULER : case NAVIER_STOKES:
			SetSurfaceCSV_Flow_Parallel(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
			break;
			
		case RANS :
			SetSurfaceCSV_Flow_Parallel(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
			break;

		case ADJ_EULER : case ADJ_NAVIER_STOKES :
			SetSurfaceCSV_Adjoint_Parallel(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
			break;

		case ADJ_RANS :
			SetSurfaceCSV_Adjoint_Parallel(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
			break;

		case ONE_SHOT_ADJ_EULER : case ONE_SHOT_ADJ_NAVIER_STOKES :
			SetSurfaceCSV_Flow_Parallel(config, geometry, solution_container[MESH_0][FLOW_SOL], config->GetSurfFlowCoeff_FileName(), iExtIter);
			SetSurfaceCSV_Adjoint_Parallel(config, geometry, solution_container[MESH_0][ADJFLOW_SOL], config->GetSurfAdjCoeff_FileName(), iExtIter);
			break;
			
		case LIN_EULER : case LIN_NAVIER_STOKES :
			SetSurfaceCSV_Linearized_Parallel(config, geometry, solution_container[MESH_0][LINFLOW_SOL], config->GetSurfLinCoeff_FileName(), iExtIter);
			break;
	}
#endif
	
}

void COutput::SetEquivalentArea(CSolution *solution_container, CGeometry *geometry,
		CConfig *config, unsigned long iExtIter) {

#ifdef NO_MPI

	unsigned long iVertex, iPoint, jVertex, nVertex_NearField, *IdPoint, auxPoint;
	unsigned short iMarker, iVar;
	double *EquivArea, *TargetArea, Coord_i, jFunction,
	jp1Function, DeltaX, MeanFuntion, *Face_Normal = NULL, InverseDesign,
	*Xcoord, *Pressure, *FaceArea, auxCoord, auxPress, auxArea, *NearFieldWeight, *Weight,
	Coord_j, jp1Coord, *Coord = NULL;
	ofstream EquivArea_file, FuncGrad_file;

	double Mach  = config->GetMach_FreeStreamND();
	double Beta = sqrt(Mach*Mach-1.0);
	double Begin_Int = config->GetEA_IntLimit(0);
	double End_Int = config->GetEA_IntLimit(1);
	double R_Plane = abs(config->GetEA_IntLimit(2));
	double Pressure_Inf = solution_container->GetPressure_Inf();
	double Density_Inf = solution_container->GetDensity_Inf();
	double RefAreaCoeff = config->GetRefAreaCoeff();
	double ModVelocity_Inf = solution_container->GetModVelocity_Inf();
	unsigned short nDim = geometry->GetnDim();
	double factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Density_Inf*Mach*Mach);

	/*--- Compute the total number of points of the near-field ---*/
	nVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				/*--- Be careful, on the 3D case we must extract a line! ---*/
				if ((Begin_Int < Coord[0]) && (Coord[0] < End_Int))
					if ((Face_Normal[nDim-1] > 0.0) &&
							((geometry->GetnDim() == 2) || ((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3)))) {
						nVertex_NearField ++;
					}

			}

	/*--- Create an array with all the coordinates, points, pressures, face area,
	 equivalent area, and nearfield weight ---*/
	Xcoord = new double[nVertex_NearField];
	IdPoint = new unsigned long[nVertex_NearField];
	Pressure = new double[nVertex_NearField];
	FaceArea = new double[nVertex_NearField];
	EquivArea = new double[nVertex_NearField];
	TargetArea = new double[nVertex_NearField];
	NearFieldWeight = new double[nVertex_NearField];
	Weight = new double[nVertex_NearField];

	/*--- Copy the information to the new array (only the upper surface of the boundary) ---*/
	nVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				/*--- Be careful, on the 3D case we must extract a line! ---*/
				if ((Begin_Int < Coord[0]) && (Coord[0] < End_Int))
					if ((Face_Normal[nDim-1] > 0.0) &&
							((geometry->GetnDim() == 2) || ((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3)))) {
						IdPoint[nVertex_NearField] = iPoint;
						Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
						Pressure[nVertex_NearField] = solution_container->node[iPoint]->GetPressure();
						FaceArea[nVertex_NearField] = abs(Face_Normal[nDim-1]);
						nVertex_NearField ++;
					}
			}

	/*--- Order the arrays (x Coordinate, Pressure, Area, and Point) ---*/
	for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
		for (jVertex = 0; jVertex < nVertex_NearField - 1 - iVertex; jVertex++)
			if (Xcoord[jVertex] > Xcoord[jVertex+1]) {
				auxCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex+1]; Xcoord[jVertex+1] = auxCoord;
				auxPress = Pressure[jVertex]; Pressure[jVertex] = Pressure[jVertex+1]; Pressure[jVertex+1] = auxPress;
				auxArea = FaceArea[jVertex]; FaceArea[jVertex] = FaceArea[jVertex+1]; FaceArea[jVertex+1] = auxArea;
				auxPoint = IdPoint[jVertex]; IdPoint[jVertex] = IdPoint[jVertex+1]; IdPoint[jVertex+1] = auxPoint;
			}

	/*--- Read target area from the configuration file ---*/
	double *x, *f;
	unsigned long nPointTargetArea;
	double CoordTA_aux, TargetArea_aux;
	ifstream index_file;
	string text_line;

	/*--- Review if there is a file ---*/
	index_file.open("TargetEA.dat", ios::in);
	if (index_file.fail()) {
		if (iExtIter == 0) { cout << "There is no Target Equivalent Area file!!"<< endl;
		cout << "Using default parameters (Target Equiv Area = 0.0)" << endl;
		}
		nPointTargetArea = 2;
		x = new double[nPointTargetArea]; f = new double[nPointTargetArea];
		x[0] = -10E6;										 f[0] = 0.0;
		x[1] = 10E6;										 f[1] = 0.0;
		goto jump_default;
	}

	/*--- Dimensionalization bucle ---*/
	nPointTargetArea = 0;
	while (!index_file.eof()) {
		getline(index_file, text_line);
		istringstream value_line(text_line);
		value_line >> CoordTA_aux >> TargetArea_aux;
		nPointTargetArea++;
	}
	index_file.close();

	x = new double [nPointTargetArea];
	f = new double [nPointTargetArea];

	/*--- Store the information ---*/
	index_file.open("TargetEA.dat", ios::in);
	nPointTargetArea = 0;
	while (!index_file.eof()) {
		getline(index_file, text_line);
		istringstream value_line(text_line);
		value_line >> CoordTA_aux >> TargetArea_aux;
		x[nPointTargetArea] = CoordTA_aux;
		f[nPointTargetArea] = TargetArea_aux;
		nPointTargetArea++;
	}
	index_file.close();

	jump_default :

	for (iVertex = 0; iVertex < nVertex_NearField; iVertex++) {
		EquivArea[iVertex] = 0.0;
		TargetArea[iVertex] = 0.0;
		NearFieldWeight[iVertex] = 0.0;
		Weight[iVertex] = 1.0;

		/*		if ((iVertex > 0) && (iVertex < nVertex_NearField-1))
		 Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]) + 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
		 if (iVertex == 0) Weight[iVertex] = 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
		 if (iVertex == nVertex_NearField-1) Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]);*/

		for (iVar = 0; iVar < nPointTargetArea-1; iVar++) {
			if ((x[iVar] <= Xcoord[iVertex]) && (Xcoord[iVertex] <= x[iVar+1])) {
				TargetArea[iVertex] = f[iVar] + (Xcoord[iVertex]-x[iVar])*(f[iVar+1]-f[iVar] )/(x[iVar+1]-x[iVar]);
			}
		}
	}

	delete [] x;
	delete [] f;

	/*--- Compute equivalent area at each node of the line ---*/
	for (iVertex = 1; iVertex < nVertex_NearField; iVertex++) {
		EquivArea[iVertex] = 0.0;
		Coord_i = Xcoord[iVertex];
		for (jVertex = 0; jVertex < iVertex-1; jVertex++) {
			Coord_j = Xcoord[jVertex];
			jp1Coord = Xcoord[jVertex+1];

			jFunction = factor*(Pressure[jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
			jp1Function = factor*(Pressure[jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);

			DeltaX = (jp1Coord-Coord_j);
			MeanFuntion = 0.5*(jp1Function + jFunction);
			EquivArea[iVertex] += DeltaX * MeanFuntion;
		}
	}

	/*--- Evaluate the objective function ---*/
	InverseDesign = 0;
	for (iVertex = 1; iVertex < nVertex_NearField; iVertex++)
		InverseDesign += Weight[iVertex]*(EquivArea[iVertex]-TargetArea[iVertex])*(EquivArea[iVertex]-TargetArea[iVertex]);

	/*--- Evaluate the weight of the nearfield pressure for this problem ---*/
	for (iVertex = 0; iVertex < nVertex_NearField; iVertex++) {
		Coord_i = Xcoord[iVertex];
		NearFieldWeight[iVertex] = 0.0;
		for (jVertex = iVertex; jVertex < nVertex_NearField; jVertex++) {
			Coord_j = Xcoord[jVertex];
			NearFieldWeight[iVertex] += Weight[iVertex]*2.0*(EquivArea[jVertex]-TargetArea[jVertex])*factor*sqrt(Coord_j-Coord_i);
		}
	}

	/*--- Write the Area Equivalent distribution, the target area, and the derivative of the functional ---*/
	EquivArea_file.precision(15);
	EquivArea_file.open("equiv_area.csv", ios::out);
	EquivArea_file << "\"x_coord\",\"Equivalent_Area\",\"Target_Area\",\"NearField_Weight\",\"Cp\"" << endl;

	for (iVertex = 0; iVertex < nVertex_NearField; iVertex++) {
		Coord_i = Xcoord[iVertex];
		EquivArea_file << scientific << Xcoord[iVertex] << ", " << EquivArea[iVertex]
		                                                                     << ", " << TargetArea[iVertex] << ", " << NearFieldWeight[iVertex] << ", " <<
		                                                                     (Pressure[iVertex]-Pressure_Inf)/(0.5*Density_Inf*RefAreaCoeff*ModVelocity_Inf*ModVelocity_Inf) << endl;
	}
	EquivArea_file.close();

	FuncGrad_file.precision(15);
	FuncGrad_file.open("WeightNF.dat", ios::out);
	for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
		FuncGrad_file << scientific << Xcoord[iVertex] <<"\t"<< NearFieldWeight[iVertex] << endl;
	FuncGrad_file.close();

	/*--- Store the value of the NearField coefficient ---*/
	solution_container->SetTotal_CEquivArea(InverseDesign);

	/*--- Delete structures ---*/
	delete [] Xcoord;
	delete [] IdPoint;
	delete [] Pressure;
	delete [] FaceArea;
	delete [] EquivArea;
	delete [] TargetArea;
	delete [] NearFieldWeight;

#else

	int rank = MPI::COMM_WORLD.Get_rank();
	int nProcessor = MPI::COMM_WORLD.Get_size();

	unsigned short iMarker = 0, iDim, iVar;
	unsigned long jVertex, nLocalVertex_NearField = 0, nGlobalVertex_NearField = 0,
			iVertex, iPoint, MaxLocalVertex_NearField = 0, Buffer_Send_nVertex[1],
			Buffer_Receive_nVertex[nProcessor][1], auxDomain, auxPoint;
	double Coord_i, Coord_j, jFunction, jp1Coord, jp1Function, DeltaX, MeanFuntion, InverseDesign,
	*Face_Normal, auxArea, auxPress, auxCoord, *Coord;
	int iProcessor;
	ofstream EquivArea_file, FuncGrad_file;

	double Mach = config->GetMach_FreeStreamND();
	double AoA = (config->GetAoA()*PI_NUMBER) / 180.0;
	double AoS = (config->GetAoS()*PI_NUMBER) / 180.0;
	double Beta = sqrt(Mach*Mach-1.0);
	double Begin_Int = config->GetEA_IntLimit(0);
	double End_Int = config->GetEA_IntLimit(1);
	double R_Plane = abs(config->GetEA_IntLimit(2));
	double Pressure_Inf = 1.0 / (Gamma*Mach*Mach);	/*--- Be carefull the adimensionalization could change ---*/
	double Density_Inf = 1.0;												/*--- Be carefull the adimensionalization could change ---*/
	double RefAreaCoeff = config->GetRefAreaCoeff();
	double Velocity_Inf[3];
	Velocity_Inf[0] = cos(AoA) * cos(AoS) * Mach * sqrt(Gamma*Pressure_Inf/Density_Inf);
	Velocity_Inf[1] = sin(AoS) * Mach * sqrt(Gamma*Pressure_Inf/Density_Inf);
	Velocity_Inf[2] = sin(AoA) * cos(AoS) * Mach * sqrt(Gamma*Pressure_Inf/Density_Inf);
	double ModVelocity_Inf = 0;
	for (iDim = 0; iDim < 3; iDim++)
		ModVelocity_Inf += Velocity_Inf[iDim]* Velocity_Inf[iDim];
	ModVelocity_Inf = sqrt(ModVelocity_Inf);
	unsigned short nDim = geometry->GetnDim();
	double factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Density_Inf*Mach*Mach);

	/*--- Compute the number of vertex without ghost nodes ---*/
	if (rank != MASTER_NODE) {
		nLocalVertex_NearField = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
				for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Coord = geometry->node[iPoint]->GetCoord();

					/*--- Be careful, on the 3D case we must extract a line! ---*/
					if (geometry->node[iPoint]->GetDomain())
						if ((Begin_Int < Coord[0]) && (Coord[0] < End_Int))
							if ((Face_Normal[nDim-1] > 0.0) && ((geometry->GetnDim() == 2) ||
									((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3))))
								nLocalVertex_NearField ++;
				}
	}
	else nLocalVertex_NearField = 0;
	Buffer_Send_nVertex[0] = nLocalVertex_NearField;

	/*--- Send Near-Field vertex information --*/
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &nGlobalVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Allgather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, &Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);

	double Buffer_Send_Xcoord[MaxLocalVertex_NearField];
	double Buffer_Receive_Xcoord[nProcessor][MaxLocalVertex_NearField];
	unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;

	unsigned long Buffer_Send_IdPoint[MaxLocalVertex_NearField];
	unsigned long Buffer_Receive_IdPoint[nProcessor][MaxLocalVertex_NearField];
	unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;

	double Buffer_Send_Pressure[MaxLocalVertex_NearField];
	double Buffer_Receive_Pressure[nProcessor][MaxLocalVertex_NearField];
	unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;

	double Buffer_Send_FaceArea[MaxLocalVertex_NearField];
	double Buffer_Receive_FaceArea[nProcessor][MaxLocalVertex_NearField];
	unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;

	for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
		Buffer_Send_IdPoint[iVertex] = 0;
		Buffer_Send_Pressure[iVertex] = 0.0;
		Buffer_Send_FaceArea[iVertex] = 0.0;
		Buffer_Send_Xcoord[iVertex] = 0.0;
	}

	/*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/
	if (rank != MASTER_NODE) {
		nLocalVertex_NearField = 0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
				for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
					Coord = geometry->node[iPoint]->GetCoord();

					if (geometry->node[iPoint]->GetDomain())
						if ((Begin_Int < Coord[0]) && (Coord[0] < End_Int))
							if ((Face_Normal[nDim-1] > 0.0) && ((geometry->GetnDim() == 2) ||
									((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3)))) {
								Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
								Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
								Buffer_Send_Pressure[nLocalVertex_NearField] = solution_container->node[iPoint]->GetPressure();
								Buffer_Send_FaceArea[nLocalVertex_NearField] = abs(Face_Normal[nDim-1]);
								nLocalVertex_NearField++;
							}
				}
	}

	/*--- Send all the information --*/
	MPI::COMM_WORLD.Gather(&Buffer_Send_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, &Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, &Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_Pressure, nBuffer_Pressure, MPI::DOUBLE, &Buffer_Receive_Pressure, nBuffer_Pressure, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(&Buffer_Send_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, &Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, MASTER_NODE);

	if (rank == MASTER_NODE) {

		/*--- Create an array with all the coordinates, points, pressures, and domains ---*/
		double *Xcoord, *Pressure, *FaceArea, *EquivArea, *TargetArea, *NearFieldWeight, *Weight;
		unsigned long *IdPoint, *IdDomain;

		Xcoord = new double[nGlobalVertex_NearField];
		IdPoint = new unsigned long[nGlobalVertex_NearField];
		IdDomain = new unsigned long[nGlobalVertex_NearField];
		Pressure = new double[nGlobalVertex_NearField];
		FaceArea = new double[nGlobalVertex_NearField];
		EquivArea = new double[nGlobalVertex_NearField];
		TargetArea = new double[nGlobalVertex_NearField];
		NearFieldWeight = new double[nGlobalVertex_NearField];
		Weight = new double[nGlobalVertex_NearField];

		nGlobalVertex_NearField = 0;
		for (iProcessor = 1; iProcessor < nProcessor; iProcessor++)
			for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor][0]; iVertex++) {
				Xcoord[nGlobalVertex_NearField] = Buffer_Receive_Xcoord[iProcessor][iVertex];
				IdPoint[nGlobalVertex_NearField] = Buffer_Receive_IdPoint[iProcessor][iVertex];
				Pressure[nGlobalVertex_NearField] = Buffer_Receive_Pressure[iProcessor][iVertex];
				FaceArea[nGlobalVertex_NearField] = Buffer_Receive_FaceArea[iProcessor][iVertex];
				IdDomain[nGlobalVertex_NearField] = iProcessor;
				nGlobalVertex_NearField++;
			}

		/*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/
		for (iVertex = 0; iVertex < nGlobalVertex_NearField; iVertex++)
			for (jVertex = 0; jVertex < nGlobalVertex_NearField - 1 - iVertex; jVertex++)
				if (Xcoord[jVertex] > Xcoord[jVertex+1]) {
					auxCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex+1]; Xcoord[jVertex+1] = auxCoord;
					auxPress = Pressure[jVertex]; Pressure[jVertex] = Pressure[jVertex+1]; Pressure[jVertex+1] = auxPress;
					auxArea = FaceArea[jVertex]; FaceArea[jVertex] = FaceArea[jVertex+1]; FaceArea[jVertex+1] = auxArea;
					auxPoint = IdPoint[jVertex]; IdPoint[jVertex] = IdPoint[jVertex+1]; IdPoint[jVertex+1] = auxPoint;
					auxDomain = IdDomain[jVertex]; IdDomain[jVertex] = IdDomain[jVertex+1]; IdDomain[jVertex+1] = auxDomain;
				}

		/*--- Read target area from the configuration file ---*/
		double *x, *f;
		unsigned long nPointTargetArea;
		double CoordTA_aux, TargetArea_aux;
		ifstream index_file;
		string text_line;

		/*--- Review if there is a file ---*/
		index_file.open("TargetEA.dat", ios::in);
		if (index_file.fail()) {
			if (iExtIter == 0) { cout << "There is no Target Equivalent Area file!!"<< endl;
			cout << "Using default parameters (Target Equiv Area = 0.0)" << endl;
			}
			nPointTargetArea = 2;
			x = new double[nPointTargetArea]; f = new double[nPointTargetArea];
			x[0] = -10E6;										 f[0] = 0.0;
			x[1] = 10E6;										 f[1] = 0.0;
			goto jump_default;
		}

		/*--- Dimensionalization bucle ---*/
		nPointTargetArea = 0;
		while (!index_file.eof()) {
			getline(index_file, text_line);
			istringstream value_line(text_line);
			value_line >> CoordTA_aux >> TargetArea_aux;
			nPointTargetArea++;
		}
		index_file.close();

		x = new double [nPointTargetArea];
		f = new double [nPointTargetArea];

		/*--- Store the information ---*/
		index_file.open("TargetEA.dat", ios::in);
		nPointTargetArea = 0;
		while (!index_file.eof()) {
			getline(index_file, text_line);
			istringstream value_line(text_line);
			value_line >> CoordTA_aux >> TargetArea_aux;
			x[nPointTargetArea] = CoordTA_aux;
			f[nPointTargetArea] = TargetArea_aux;
			nPointTargetArea++;
		}
		index_file.close();

		jump_default :

		for (iVertex = 0; iVertex < nGlobalVertex_NearField; iVertex++) {
			EquivArea[iVertex] = 0.0;
			TargetArea[iVertex] = 0.0;
			NearFieldWeight[iVertex] = 0.0;
			Weight[iVertex] = 1.0;

			/*			if ((iVertex > 0) && (iVertex < nGlobalVertex_NearField-1))
				Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]) + 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
			if (iVertex == 0) Weight[iVertex] = 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
			if (iVertex == nGlobalVertex_NearField-1) Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]);*/

			for (iVar = 0; iVar < nPointTargetArea-1; iVar++) {
				if ((x[iVar] <= Xcoord[iVertex]) && (Xcoord[iVertex] <= x[iVar+1])) {
					TargetArea[iVertex] = f[iVar] + (Xcoord[iVertex]-x[iVar])*(f[iVar+1]-f[iVar] )/(x[iVar+1]-x[iVar]);
				}
			}
		}

		/*--- Compute equivalent area at each node of the line ---*/
		EquivArea[0] = 0.0;
		for (iVertex = 1; iVertex < nGlobalVertex_NearField; iVertex++) {
			EquivArea[iVertex] = 0.0;
			Coord_i = Xcoord[iVertex];
			for (jVertex = 0; jVertex < iVertex-1; jVertex++) {
				Coord_j = Xcoord[jVertex];
				jp1Coord = Xcoord[jVertex+1];

				jFunction = factor*(Pressure[jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
				jp1Function = factor*(Pressure[jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);

				DeltaX = (jp1Coord-Coord_j);
				MeanFuntion = 0.5*(jp1Function + jFunction);
				EquivArea[iVertex] += DeltaX * MeanFuntion;
			}
		}

		/*--- Evaluate the objective function ---*/
		InverseDesign = 0.0;
		for (iVertex = 0; iVertex < nGlobalVertex_NearField; iVertex++)
			InverseDesign += Weight[iVertex]*(EquivArea[iVertex]-TargetArea[iVertex])*(EquivArea[iVertex]-TargetArea[iVertex]);

		/*--- Evaluate the weight of the nearfield pressure for this problem ---*/
		for (iVertex = 0; iVertex < nGlobalVertex_NearField; iVertex++) {
			Coord_i = Xcoord[iVertex];
			NearFieldWeight[iVertex] = 0.0;
			for (jVertex = iVertex; jVertex < nGlobalVertex_NearField; jVertex++) {
				Coord_j = Xcoord[jVertex];
				NearFieldWeight[iVertex] += Weight[iVertex]*2.0*(EquivArea[jVertex]-TargetArea[jVertex])*factor*sqrt(Coord_j-Coord_i);
			}
		}

		/*--- Write the Area Equivalent distribution, the target area, and the derivative of the functional ---*/
		EquivArea_file.precision(15);
		EquivArea_file.open("equiv_area.csv", ios::out);
		EquivArea_file << "\"x_coord\",\"Equivalent_Area\",\"Target_Area\",\"NearField_Weight\",\"Cp\"" << endl;

		for (iVertex = 0; iVertex < nGlobalVertex_NearField; iVertex++) {
			Coord_i = Xcoord[iVertex];
			EquivArea_file << scientific << Xcoord[iVertex] << ", " << EquivArea[iVertex]
			                                                                     << ", " << TargetArea[iVertex] << ", " << NearFieldWeight[iVertex] << ", " <<
			                                                                     (Pressure[iVertex]-Pressure_Inf)/(0.5*Density_Inf*RefAreaCoeff*ModVelocity_Inf*ModVelocity_Inf) << endl;
		}

		FuncGrad_file.precision(15);
		FuncGrad_file.open("WeightNF.dat", ios::out);
		for (iVertex = 0; iVertex < nGlobalVertex_NearField; iVertex++) {
			FuncGrad_file << scientific << Xcoord[iVertex] <<"\t"<< NearFieldWeight[iVertex] << endl;
		}
		FuncGrad_file.close();

		/*--- Delete structures ---*/
		delete [] Xcoord;
		delete [] IdPoint;
		delete [] Pressure;
		delete [] FaceArea;
		delete [] EquivArea;
		delete [] TargetArea;
		delete [] NearFieldWeight;
		delete [] Weight;

	}

	/*--- Send the value of the NearField coefficient to all the processors ---*/
	MPI::COMM_WORLD.Bcast (&InverseDesign, 1, MPI::DOUBLE, MASTER_NODE);

	/*--- Store the value of the NearField coefficient ---*/
	if (rank != MASTER_NODE)
		solution_container->SetTotal_CEquivArea(InverseDesign);

#endif

}

void COutput::WriteInOutputFile(CGeometry *geometry, CSolution *solution_container, ofstream &ConsVar_file, 
																string Quantity_Name, bool Scalar, short iVar, string func, CConfig *config) {
	unsigned long iPoint;
	unsigned short iDim, Function;
	unsigned short Solution = 0, Pressure = 1, Residual = 2, DensityInc = 3, LaminarViscosityInc = 4, 
	Mach = 5, Temperature = 6, LaminarViscosity = 7, EddyViscosity = 8, Source = 9, GetRotVel = 0, GetResidual = 1, 
	GetVelocity = 2, GetGridVel = 3, GetGradient = 4, GetSolution = 5, Limiter = 10, UndLaplacian = 11;
	double mach_number, *Vector = new double [geometry->GetnDim()];
	
	bool incompressible = config->GetIncompressible();

	if (Scalar) {
		
		if (func == "GetSolution")						Function = Solution;
		if (func == "GetPressure")						Function = Pressure;
		if (func == "GetResidual")						Function = Residual;
		if (func == "GetDensityInc")					Function = DensityInc;
		if (func == "GetLaminarViscosityInc")	Function = LaminarViscosityInc;
		if (func == "GetMach")								Function = Mach;
		if (func == "GetTemperature")					Function = Temperature;
		if (func == "GetLaminarViscosity")		Function = LaminarViscosity;
		if (func == "GetEddyViscosity")				Function = EddyViscosity;
		if (func == "GetSource")							Function = Source;
		if (func == "GetLimiter")							Function = Limiter;
		if (func == "GetUndLaplacian")				Function = UndLaplacian;
		
		ConsVar_file << "SCALARS " << Quantity_Name << " float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		
		switch (Function) {
			case 0:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetSolution(iVar) << endl;				
				break;
			case 1:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetPressure() << endl;
				break;
			case 2:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetResidual(iVar) << endl;
				break;
			case 3:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetDensityInc() << endl;
				break;
			case 4:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetLaminarViscosityInc() << endl;
				break;
			case 5:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					if (!incompressible)
						mach_number = sqrt(solution_container->node[iPoint]->GetVelocity2())/solution_container->node[iPoint]->GetSoundSpeed();
					else {
						double Density = solution_container->node[iPoint]->GetDensityInc()*config->GetDensity_Ref();
						double SoundSpeed = sqrt(config->GetBulk_Modulus()/Density);
						double SoundSpeedND = SoundSpeed / config->GetVelocity_Ref();
						mach_number = sqrt(solution_container->node[iPoint]->GetVelocity2())/SoundSpeedND;
					}
					ConsVar_file << fixed <<  mach_number << endl;
				}
				break;
			case 6:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetTemperature() << endl;
				break;
			case 7:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetLaminarViscosity() << endl;
				break;
			case 8:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetEddyViscosity() << endl;
				break;
			case 9:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetSource(iVar) << endl;
				break;
			case 10:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetLimiter(iVar) << endl;
				break;
			case 11:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					ConsVar_file << fixed << solution_container->node[iPoint]->GetUnd_Lapl(iVar) << endl;
				break;
			default:
				cout <<"This plot map is not defined." << endl;
				break;
		}
	}
	else {
		
		if (func == "GetRotVel")			Function = GetRotVel;
		if (func == "GetResidual")		Function = GetResidual;
		if (func == "GetVelocity")		Function = GetVelocity;
		if (func == "GetGridVel")			Function = GetGridVel;
		if (func == "GetGradient")		Function = GetGradient;
		if (func == "GetSolution")		Function = GetSolution;
		
		ConsVar_file << "VECTORS " << Quantity_Name << " float " << endl;
		
		switch (Function) {
				
			case 0:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					Vector = geometry->node[iPoint]->GetRotVel();
					if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
					if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
				}
				break;
				
			case 1:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
						Vector[iDim] = solution_container->node[iPoint]->GetResidual(iDim+1);
					
					if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
					if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
				}
				break;
				
			case 2:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
						Vector[iDim] = solution_container->node[iPoint]->GetVelocity(iDim);
					if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
					if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
				}
				break;
				
			case 3:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					Vector = geometry->node[iPoint]->GetGridVel();
					if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
					if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
				}
				break;
				
			case 4:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++)
						Vector[iDim] = solution_container->node[iPoint]->GetGradient(0, iDim);
					
					if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
					if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
				}
				break;
				
			case 5:
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++)
						Vector[iDim] = solution_container->node[iPoint]->GetSolution(iDim+1);
					if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
					if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
				}
				break;
		}
	}
	
	delete [] Vector;
}

void COutput::WriteReactingOutputFile(CConfig *config, CGeometry *geometry,CSolution **solution_container, ofstream & ConsVar_file) {

	unsigned short iPoint, iDim, nDim;
	double M1, M2, M3;
	double Vector[3];
	bool scalar = true;
	ConsVar_file.precision(15);

	M1 = config->GetParticle_Mass(0);
	M3 = config->GetParticle_Mass(2);
	M2 = M1-M3;
	nDim = geometry->GetnDim();

	/*--- PRINT OUT ELECTRIC POTENTIAL ---*/
	WriteInOutputFile(geometry, solution_container[ELEC_SOL], ConsVar_file, "Electric_Potential", scalar, 0, "GetSolution", config);

	/*--- PRINT OUT ELECTRIC FIELD Vector ---*/
	ConsVar_file << "VECTORS Electric_Field float" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			Vector[iDim] = -1.0* solution_container[ELEC_SOL]->node[iPoint]->GetGradient(0,iDim);
		}
		if (geometry->GetnDim() == 2) ConsVar_file << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
		if (geometry->GetnDim() == 3) ConsVar_file << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
	}

#ifdef OutputSource

	/*--- PRINT OUT SOURCE TERM 1 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src1", scalar, 0, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 2 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src2", scalar, 1, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 3 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src3", scalar, 2, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 4 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src4", scalar, 3, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 5 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src5", scalar, 4, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 6 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src6", scalar, 5, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 7 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src7", scalar, 6, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 8 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src8", scalar, 7, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 9 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src9", scalar, 8, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 10 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src10", scalar, 9, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 11 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src11", scalar, 10, "GetSource", config);

	/*--- PRINT OUT SOURCE TERM 12 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src12", scalar, 11, "GetSource", config);


#endif


	/*--- PRINT OUT THE NET CHARGE IN DOMAIN ---*/
	ConsVar_file << "SCALARS Net_Charge float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << ELECTRON_CHARGE*(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))- solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))*M2/M3)/M2  << endl;


	/*--- PRINT OUT THE degree of ionization  IN DOMAIN ---*/
	ConsVar_file << "SCALARS alpha float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))/M3)/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))/M3 + solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(0*(nDim+2))/M1) << endl;

	/*--- PRINT OUT THE NEGATIVE CHARGE IN DOMAIN ---*/
	ConsVar_file << "SCALARS Negative_Charge float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << ELECTRON_CHARGE*solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))/M3  << endl;

	/*--- PRINT OUT THE POSITIVE CHARGE IN DOMAIN ---*/
	ConsVar_file << "SCALARS Positive_Charge float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << ELECTRON_CHARGE*solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))/M2  << endl;

	/*--- PRINT OUT THE Total Density IN DOMAIN ---*/
	ConsVar_file << "SCALARS Total_Density float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(0*(nDim+2)) + solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))+solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))  << endl;

	/*--- PRINT OUT Total Pressure  ---*/
	ConsVar_file << "SCALARS Total_Pressure float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(0) +solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(1)+solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(2) << endl;

	int loc, iSpecies;

	/* ***************************************** ---*/
	/* SPECIES 1 DATA OUTPUT ---*/
	/* ***************************************** ---*/

	iSpecies = 0;
	loc = (nDim+2)*iSpecies;

	/*--- PRINT OUT THE DENSITY OF SPECIES 1 ---*/
	//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Argon_Density", scalar, loc + 0, "GetSolution");

	ConsVar_file << "SCALARS Argon_Density float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0)) << endl;

	/*--- PRINT OUT DENSITY* Energy OF SPECIES 1 ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species1", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);

	/*--- PRINT OUT Velocity vector OF SPECIES 1 ---*/
	ConsVar_file << "VECTORS Argon_Velocity float" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iDim = 0; iDim < nDim; iDim++)
			Vector[iDim] = solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity(iDim,iSpecies);
		if (nDim== 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
		if (nDim == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
	}

	/*--- PRINT OUT Pressure OF SPECIES 1 ---*/
	ConsVar_file << "SCALARS Argon_Pressure float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies) << endl;

	/*--- PRINT OUT Mach number OF SPECIES 1 ---*/
	ConsVar_file << "SCALARS Mach_Number_Species1 float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << fixed <<  sqrt(solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))/solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies) << endl;


	/*--- PRINT OUT Temperature OF SPECIES 1 ---*/
	ConsVar_file << "SCALARS Argon_Temperature float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies))/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(0)*8314.0/(M1*AVOGAD_CONSTANT)) << endl;

	int nSpecies = config->GetnSpecies();
	if (nSpecies > 1) {
		/* ***************************************** ---*/
		/* SPECIES 2 DATA OUTPUT ---*/
		/* ***************************************** ---*/

		iSpecies = 1;
		loc = (nDim+2)*iSpecies;

		/*--- PRINT OUT THE DENSITY OF SPECIES 2 ---*/
		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Ion_Density", scalar,loc + 0, "GetSolution");
		ConsVar_file << "SCALARS Ion_Density float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0)) << endl;

		/*--- PRINT OUT DENSITY* Energy OF SPECIES 2 ---*/
		WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species2", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);

		/*--- PRINT OUT Velocity vector OF SPECIES 2 ---*/
		ConsVar_file << "VECTORS Ion_Velocity float" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				Vector[iDim] = solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity(iDim,iSpecies);
			if (nDim == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
			if (nDim == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
		}

		/*--- PRINT OUT Pressure OF SPECIES 2 ---*/
		ConsVar_file << "SCALARS Ion_Pressure float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies) << endl;

		/*--- PRINT OUT Mach number OF SPECIES 2 ---*/
		ConsVar_file << "SCALARS Mach_Number_Species2 float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed <<  sqrt(solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))/solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies) << endl;

		/*--- PRINT OUT Temperature OF SPECIES 1 ---*/
		ConsVar_file << "SCALARS Ion_Temperature float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies))/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc + 0)*8314.0/(M2*AVOGAD_CONSTANT)) << endl;

	}

	if (nSpecies > 2) {

		/* ***************************************** ---*/
		/* SPECIES 3 DATA OUTPUT ---*/
		/* ***************************************** ---*/
		iSpecies = 2;
		loc = (nDim+2)*iSpecies;

		/*--- PRINT OUT THE DENSITY OF SPECIES 3 ---*/
		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Electron_Density", scalar, loc + 0, "GetSolution");
		ConsVar_file << "SCALARS Electron_Density float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0)) << endl;

		/*--- PRINT OUT DENSITY * Energy OF SPECIES 3 ---*/
		WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species3", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);

		/*--- PRINT OUT Velocity vector OF SPECIES 3 ---*/
		ConsVar_file << "VECTORS Electron_Velocity float" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				Vector[iDim] = solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity(iDim,iSpecies);
			if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
			if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
		}

		/*--- PRINT OUT Pressure OF SPECIES 3 ---*/
		ConsVar_file << "SCALARS Electron_Pressure float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies) << endl;

		/*--- PRINT OUT Temperature OF SPECIES 3 ---*/
		ConsVar_file << "SCALARS Electron_Temperature float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies))/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc + 0)*8314.0/(M3*AVOGAD_CONSTANT)) << endl;


		/*--- PRINT OUT Mach number OF SPECIES 3 ---*/
		ConsVar_file << "SCALARS Mach_Number_Species3 float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed <<  sqrt(solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))/solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies) << endl;

	}
}

