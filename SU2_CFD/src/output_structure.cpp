/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
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

#include "../include/output_structure.hpp"

COutput::COutput(CConfig *config) {
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

COutput::~COutput(void) { }

void COutput::SetDomain_Flow(CConfig *config, CGeometry ***geometry, CSolution ****solution_container, unsigned long iExtIter, unsigned short val_nDomain) {
	char cstr[200], buffer[50];
	ofstream ConsVar_file;
	ofstream ConsVar_file_mean;
	unsigned short iSpecies, iVar, loc, iDim, iDomain;
	unsigned long iPoint, iElem;
	double coord[3] = {0.0,0.0,0.0}, velocity[3] = {0.0,0.0,0.0}, mach_number;
	double rel_velocity[3] = {0.0,0.0,0.0}, vorticity[3] = {0.0,0.0,0.0}, rel_mach_number = 0.0, vorticity_mag = 0.0;

	unsigned short solver = config->GetKind_Solver();
	bool incompressible = config->GetIncompressible();
	bool scalar = true;
	bool rotating_frame = config->GetRotating_Frame();


	/*--- Paraview Output ---*/
	if (config->GetOutput_FileFormat() == PARAVIEW) {

		/*--- Write file name with extension ---*/
		strcpy (cstr, config->GetFlow_FileName().c_str());
		if (config->GetUnsteady_Simulation())
			sprintf (buffer, "_%d.vtk", int(iExtIter));
		else 
			sprintf (buffer, ".vtk");

		strcat(cstr,buffer);

		/*--- Write header and open the file ---*/
		geometry[DOMAIN_0][MESH_0]->SetParaView(cstr);
		ConsVar_file.precision(15);
		ConsVar_file.open(cstr, ios::out | ios::app);
		ConsVar_file << "POINT_DATA " << geometry[DOMAIN_0][MESH_0]->GetnPoint() << endl;

		if (solver == ELECTRIC_POTENTIAL) {
			WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][ELEC_SOL], ConsVar_file, "Electric_Potential", scalar, 0, "GetSolution", config);
			return;
		}

		else if (solver == NS_PLASMA) {
			unsigned short nDim = geometry[DOMAIN_0][MESH_0]->GetnDim();
			unsigned short nSpecies = solution_container[DOMAIN_0][MESH_0][PLASMA_SOL]->GetnSpecies();
			unsigned short nDiatomics = solution_container[DOMAIN_0][MESH_0][PLASMA_SOL]->GetnDiatomics();


			if (config->GetKind_GasModel() == ARGON || config->GetKind_GasModel() == AIR21)
				WriteReactingOutputFile(config, geometry[DOMAIN_0][MESH_0],solution_container[DOMAIN_0][MESH_0], ConsVar_file);
			else if (config->GetKind_GasModel()==AIR7) {
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][PLASMA_SOL], ConsVar_file, "Density", scalar,  loc+0, "GetSolution", config);
					//					WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "", scalar,  5, "GetSolution", config);
				}
			}

			return;
		}

		else if  ((solver == EULER) || (solver == NAVIER_STOKES) || (solver == RANS) ||
				(solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES)) {

			if (incompressible) {
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Pressure", scalar, 0, "GetSolution", config);
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Velocity", !scalar, 0, "GetSolution", config);
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Density", scalar,  0, "GetDensityInc", config);
				if (config->GetKind_SlopeLimit_Flow() != NONE) 
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Limiter", scalar,  0, "GetLimiter", config);
				//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Viscosity", scalar, 0, "GetLaminarViscosityInc", config);
				//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "UndLaplacian", scalar,  0, "GetUndLaplacian", config);
			}

			else {
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Density", scalar,  0, "GetSolution", config);
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Density_x_Velocity", !scalar, 0, "GetSolution", config);
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Density_x_Energy", scalar, geometry[DOMAIN_0][MESH_0]->GetnDim()+1, "GetSolution", config);
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Velocity", !scalar,  0, "GetVelocity", config);
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Pressure", scalar,  0, "GetPressure", config);
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Mach_Number", scalar,  0, "GetMach", config);
				//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "UndLaplacian", scalar,  0, "GetUndLaplacian", config);
				//				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Limiter", scalar,  0, "GetLimiter", config);

				if (config->GetGrid_Movement())
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Grid_Velocity", !scalar,  0, "GetGridVel", config);

				if (config->GetRotating_Frame())
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Rotational_Velocity", !scalar,  0, "GetRotVel", config);

				if ((solver == NAVIER_STOKES) || (solver == RANS) || (solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES)) {
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Temperature", scalar,  0, "GetTemperature", config);
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Laminar_Viscosity", scalar,  0, "GetLaminarViscosity", config);
				}

				if (solver == RANS) {
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], ConsVar_file, "Eddy_Viscosity", scalar,  0, "GetEddyViscosity", config);
					WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][TURB_SOL], ConsVar_file,  "Nu_Tilde", scalar, 0, "GetSolution", config);
				}
			}

			if ((solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES))
				WriteInOutputFile(geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][LEVELSET_SOL], ConsVar_file,   "LevelSet", scalar, 0, "GetSolution", config);

			return;
		}
	}

	/*--- Tecplot Output ---*/
	if (config->GetOutput_FileFormat() == TECPLOT) {

		/*--- Write file name with extension ---*/
		strcpy (cstr, config->GetFlow_FileName().c_str());
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.plt", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.plt", int(iExtIter));
		}
		else
			sprintf (buffer, ".plt");

		strcat(cstr,buffer);

		/*--- Write header and open the file ---*/
		ConsVar_file.precision(15);
		ConsVar_file.open(cstr, ios::out);

		ConsVar_file << "TITLE = \"Visualization of the volumetric grid\"" << endl;
		ConsVar_file << "VARIABLES = \"x\", \"y\", \"z\"";

		if ((solver == EULER) || (solver == NAVIER_STOKES) 
				|| (solver == RANS) || (solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES)) {
			ConsVar_file << ", \"Density\", \"Velocity x\", \"Velocity y\", \"Velocity z\", \"Pressure\", \"Mach Number\", \"1stVar Residual\"";
			if ((config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) && 
					(config->GetKind_SlopeLimit_Flow() != NONE)) ConsVar_file << ", \"1stVar Limiter\"";
			if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED) ConsVar_file << ", \"Dissipation Switch\", \"1stVar UndLaplacian\"";
		}

		if ((solver == NAVIER_STOKES) || (solver == RANS)) {
			if (!incompressible)
				ConsVar_file << ", \"Temperature\", \"Vorticity x\", \"Vorticity y\", \"Vorticity z\", \"Distance\", \"Laminar Viscosity\"";
			else
				ConsVar_file << ", \"Laminar Viscosity\"";
		}

		if (solver == RANS)
			ConsVar_file << ", \"Eddy Viscosity\", \"Nu Tilde\"";

		if ((solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES))
			ConsVar_file << ", \"LevelSet\", \"LevelSet Residual\", \"LevelSet Diff\"";

		if (solver == WAVE)
			ConsVar_file << ", \"Wave\", \"Wave Residual\"";

		if (rotating_frame)
			ConsVar_file << ", \"Relative Vx\", \"Relative Vy\", \"Relative Vz\", \"Relative Mach\"";

		if (solver == NS_PLASMA) {
			for ( iSpecies = 0; iSpecies < solution_container[DOMAIN_0][MESH_0][PLASMA_SOL]->GetnSpecies(); iSpecies++ ) {
				//				ConsVar_file << "\"Density("<< iSpecies << ")\", \"Velocity_x("<< iSpecies << ")\", \"Velocity_y("<< iSpecies << ")\", \"DensityEnergy("<< iSpecies << ")\",";
				ConsVar_file << "\"Density("<< iSpecies << ")\", \"Velocity_x("<< iSpecies << ")\", \"Velocity_y("<< iSpecies << ")\", \"MachNumber("<< iSpecies << ")\", Temperature_TR(" << iSpecies << ")\",";
				ConsVar_file << "\"Source_Density("<< iSpecies << ")\", \"Source_Velocity_x("<< iSpecies << ")\", \"Source_Velocity_y("<< iSpecies << ")\", \"Source_DensityEnergy("<< iSpecies << ")\"";
				if (iSpecies < solution_container[DOMAIN_0][MESH_0][PLASMA_SOL]->GetnSpecies() - 1)
					ConsVar_file << ", ";
			}

			/*--- Check for writing mean solution file ---*/
			if (config->GetWrite_Mean_Solution()) {
				char cstr_mean[200], buffer_mean[50];
				strcpy (cstr_mean, config->GetFlow_FileName().c_str());
				sprintf (buffer_mean, "_mean.plt");
				strcat(cstr_mean,buffer_mean);

				/*--- Write header and open the file ---*/
				ConsVar_file_mean.precision(15);
				ConsVar_file_mean.open(cstr_mean, ios::out);

				ConsVar_file_mean << "TITLE = \"Visualization of the volumetric grid\"" << endl;
				ConsVar_file_mean << "VARIABLES = \"x\", \"y\", \"z\"";
				ConsVar_file_mean << ", \"Density\", \"Velocity x\", \"Velocity y\", \"Energy\"";
			}
		}

		ConsVar_file << endl;
		if (config->GetWrite_Mean_Solution())
			ConsVar_file_mean << endl;

		for (iDomain = 0; iDomain < val_nDomain; iDomain++) {

			ConsVar_file << "ZONE ";
			if (config->GetWrite_Mean_Solution()) 
				ConsVar_file_mean << "ZONE ";				

			if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
				ConsVar_file << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*int(iExtIter+1)<<", ";

			ConsVar_file << "NODES="<< geometry[iDomain][MESH_0]->GetnPoint() <<" , ELEMENTS="<< geometry[iDomain][MESH_0]->GetnElem() <<", DATAPACKING=POINT";
			if (config->GetWrite_Mean_Solution())
				ConsVar_file_mean << "NODES="<< geometry[iDomain][MESH_0]->GetnPoint() <<" , ELEMENTS="<< geometry[iDomain][MESH_0]->GetnElem() <<", DATAPACKING=POINT";

			if (geometry[DOMAIN_0][MESH_0]->GetnDim() == 2) {
				ConsVar_file << ", ZONETYPE=FEQUADRILATERAL" << endl;
				if (config->GetWrite_Mean_Solution())
					ConsVar_file_mean << ", ZONETYPE=FEQUADRILATERAL" << endl;
			}
			if (geometry[DOMAIN_0][MESH_0]->GetnDim() == 3) {
				ConsVar_file << ", ZONETYPE=FEBRICK"<< endl;
				if (config->GetWrite_Mean_Solution())				
					ConsVar_file_mean << ", ZONETYPE=FEBRICK"<< endl;
			}


			/*--- Cycle through all points and write the solution. ---*/
			for (iPoint = 0; iPoint < geometry[iDomain][MESH_0]->GetnPoint(); iPoint++) {

				for (iDim = 0; iDim < geometry[iDomain][MESH_0]->GetnDim(); iDim++)
					coord[iDim] = geometry[iDomain][MESH_0]->node[iPoint]->GetCoord(iDim);

				if ((solver == EULER) || (solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES)
						|| (solver == NAVIER_STOKES) || (solver == RANS)) {

					/*--- Compute vectorial quatities for 2D and 3D ---*/
					for (iDim = 0; iDim < geometry[iDomain][MESH_0]->GetnDim(); iDim++) {
						velocity[iDim] = solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetVelocity(iDim);
						if ((solver == NAVIER_STOKES) || (solver == RANS) || rotating_frame )
							vorticity[iDim] = solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetVorticity(iDim);
					}

					/*--- Compute the mach number ---*/
					mach_number = sqrt(solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetVelocity2())/
							solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetSoundSpeed();
				}

				/*--- Compute relative velocities for rotating frame ---*/
				if (rotating_frame) {
					double *rot_vel = geometry[iDomain][MESH_0]->node[iPoint]->GetRotVel();
					double rel_velocity_mag = 0.0; vorticity_mag = 0.0;
					for (iDim = 0; iDim < geometry[iDomain][MESH_0]->GetnDim(); iDim++) {
						rel_velocity[iDim] = velocity[iDim] - rot_vel[iDim];
						rel_velocity_mag  += rel_velocity[iDim]*rel_velocity[iDim];
						vorticity_mag     += vorticity[iDim]*vorticity[iDim];
					}
					rel_velocity_mag = sqrt(rel_velocity_mag);
					vorticity_mag    = sqrt(vorticity_mag);
					rel_mach_number  = rel_velocity_mag/solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetSoundSpeed();
				}

				/*--- Write the coordinates ---*/
				ConsVar_file << coord[0] <<" "<< coord[1] <<" "<< coord[2];
				if (config->GetWrite_Mean_Solution())
					ConsVar_file_mean << coord[0] << " " << coord[1] << " " << coord[2];

				/*--- Write the Euler variables ---*/
				if ((solver == EULER) || (solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES) ||
						(solver == NAVIER_STOKES) || (solver == RANS)) {

					if (!incompressible) {
						mach_number = sqrt(solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetVelocity2())/
								solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetSoundSpeed();
						for (iDim = 0; iDim < geometry[DOMAIN_0][MESH_0]->GetnDim(); iDim++)
							velocity[iDim] = solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetVelocity(iDim);
						ConsVar_file << " " << solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution((unsigned short)0) <<" "<<
								velocity[0] <<" "<< velocity[1] <<" "<< velocity[2]<<" "<<
								solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetPressure() <<" "<< mach_number ;
					}
					else {
						double Density = solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref();
						double Pressure = solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution((unsigned short)0)*config->GetPressure_Ref(); 
						double SoundSpeed = sqrt(config->GetBulk_Modulus()/Density);
						mach_number = sqrt(solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/SoundSpeed;
						for (iDim = 0; iDim < geometry[iDomain][MESH_0]->GetnDim(); iDim++)
							velocity[iDim] = (solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution(iDim+1)/solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetDensityInc())*config->GetVelocity_Ref();
						ConsVar_file << " " << Density <<" "<< velocity[0] <<" "<< velocity[1] <<" "<< velocity[2]<<" "<< Pressure <<" "<< mach_number;
					}

					ConsVar_file <<" "<< solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetResidual((unsigned short)0);

					if ((config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) &&
							(config->GetKind_SlopeLimit_Flow() != NONE))
						ConsVar_file <<" "<< solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetLimiter((unsigned short)0);
					if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTRED)
						ConsVar_file <<" "<< solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetSensor()
						<<" "<< solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetUnd_Lapl()[0];

				}

				/*--- Write the Navier-Stokes variables ---*/
				if ((solver == NAVIER_STOKES) || (solver == RANS)) {
					if (!incompressible) {
						ConsVar_file << " "<< solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetTemperature() <<" "<<
								vorticity[0] <<" "<< vorticity[1] <<" "<< vorticity[2]<<" "<<
								geometry[iDomain][MESH_0]->node[iPoint]->GetWallDistance()<<" "<<solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
					}
					else {
						double viscosity = solution_container[iDomain][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc()*config->GetViscosity_Ref();
						ConsVar_file << " "<< viscosity;
					}
				}

				/*--- Write the Navier-Stokes variables ---*/
				if (solver == RANS)
					ConsVar_file << " "<< solution_container[iDomain][MESH_0][TURB_SOL]->node[iPoint]->GetmuT()<< " " <<
					solution_container[iDomain][MESH_0][TURB_SOL]->node[iPoint]->GetSolution(0);

				/*--- Write the level set solution ---*/
				if ((solver == FREE_SURF_EULER) || (solver == FREE_SURF_NAVIER_STOKES))
					ConsVar_file << " "<< solution_container[iDomain][MESH_0][LEVELSET_SOL]->node[iPoint]->GetSolution(0)
					<< " "<< solution_container[iDomain][MESH_0][LEVELSET_SOL]->node[iPoint]->GetResSour()[0]
					                                                                                       << " "<< solution_container[iDomain][MESH_0][LEVELSET_SOL]->node[iPoint]->GetDiffLevelSet();

				/*--- Write the wave solution ---*/
				if (solver == WAVE)
					ConsVar_file << " "<< solution_container[iDomain][MESH_0][WAVE_SOL]->node[iPoint]->GetSolution(0)
					<< " "<< solution_container[iDomain][MESH_0][WAVE_SOL]->node[iPoint]->GetResVisc()[0];

				/*--- Write the multi-species flow variables ---*/
				if (solver == NS_PLASMA ) {
					unsigned short nSpecies, nDiatomics, nDim;
					double velocity2, SoundSpeed, MachNumber, Temperature_tr;
					nSpecies = solution_container[DOMAIN_0][MESH_0][PLASMA_SOL]->GetnSpecies();
					nDiatomics = solution_container[DOMAIN_0][MESH_0][PLASMA_SOL]->GetnDiatomics();
					nDim = geometry[iDomain][MESH_0]->GetnDim();
					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
						if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);	

						/*--- Primitive Variables ---*/
						// Density
						ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);

						// Velocity & Mach number
						velocity2 = 0.0;
						for (iDim = 0; iDim < nDim; iDim++) {
							ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+iDim)/solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
							velocity2 += pow(solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+iDim)/solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+0),2.0);
						}
						SoundSpeed = solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies);
						MachNumber = sqrt(velocity2)/SoundSpeed;
						ConsVar_file << " " << MachNumber;

						// Temperature (TR)
						Temperature_tr = solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetTemperature_TR(iSpecies);
						ConsVar_file << " " << Temperature_tr;
						//						ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+nDim);

						if (iSpecies < nDiatomics)
							ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+2+nDim);

						/*--- Source terms ---*/
						ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetResSour()[loc+0];
						for (iDim = 0; iDim < nDim; iDim++)
							ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetResSour()[loc+1+iDim]/solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
						ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetResSour()[loc+1+nDim];
						if (iSpecies < nDiatomics)
							ConsVar_file << " " << solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetResSour()[loc+2+nDim];
					}

					/*--- Write mean solution file ---*/
					if (config->GetWrite_Mean_Solution()) {
						double *Mass_Ave_Solution;
						double Mixture_Density, Mass_Fraction;

						Mass_Ave_Solution = new double [nDim+2];

						for (iVar = 0; iVar < nDim+2; iVar++) {
							Mass_Ave_Solution[iVar] = 0.0;
						}

						/*--- Pre-process by mass-averaging ---*/
						Mixture_Density = 0.0;
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
							if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
							else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
							Mixture_Density += solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
						}	
						for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
							if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
							else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
							Mass_Fraction = solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+0) / Mixture_Density;
							Mass_Ave_Solution[0] += solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
							for (iDim = 0; iDim < nDim; iDim++) {
								Mass_Ave_Solution[iDim+1] += solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+iDim)*Mass_Fraction;
							}
							Mass_Ave_Solution[nDim+1] += solution_container[iDomain][MESH_0][PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+nDim)*Mass_Fraction;
						}

						ConsVar_file_mean << " " << Mass_Ave_Solution[0];
						for (iDim = 0; iDim < nDim; iDim++)
							ConsVar_file_mean << " " << Mass_Ave_Solution[iDim+1] / Mass_Ave_Solution[0];
						ConsVar_file_mean << " " << Mass_Ave_Solution[nDim+1];

						delete [] Mass_Ave_Solution;
					}

				}

				/*--- Write the rotating frame variables ---*/
				if (rotating_frame) {
					ConsVar_file <<" "<< rel_velocity[0] <<" "<< rel_velocity[1] <<" "<< rel_velocity[2]<<" "<< rel_mach_number;
				}

				/*--- End line ---*/
				ConsVar_file << endl;
				if (config->GetWrite_Mean_Solution())
					ConsVar_file_mean << endl;

			}

			/*--- Write the grid connectivity. ---*/
			for(iElem = 0; iElem < geometry[iDomain][MESH_0]->GetnElem(); iElem++) {
				if (geometry[iDomain][MESH_0]->elem[iElem]->GetVTK_Type() == TRIANGLE) {
					ConsVar_file <<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 << endl;
					if (config->GetWrite_Mean_Solution()) {
						ConsVar_file_mean <<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 << endl;
					}
				}
				if (geometry[iDomain][MESH_0]->elem[iElem]->GetVTK_Type() == RECTANGLE) {
					ConsVar_file <<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 << endl;
					if (config->GetWrite_Mean_Solution()) {
						ConsVar_file_mean <<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 << endl;
					}
				}
				if (geometry[iDomain][MESH_0]->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
					ConsVar_file <<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 << endl;
					if (config->GetWrite_Mean_Solution()) {
						ConsVar_file_mean <<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 << endl;
					}
				}
				if (geometry[iDomain][MESH_0]->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
					ConsVar_file <<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(5)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(6)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(7)+1 << endl;
					if (config->GetWrite_Mean_Solution()) {
						ConsVar_file_mean <<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(5)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(6)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(7)+1 << endl;
					}
				}
				if (geometry[iDomain][MESH_0]->elem[iElem]->GetVTK_Type() == PYRAMID) {
					ConsVar_file <<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 << endl;
					if (config->GetWrite_Mean_Solution()) {
						ConsVar_file_mean <<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 << endl;
					}
				}
				if (geometry[iDomain][MESH_0]->elem[iElem]->GetVTK_Type() == WEDGE) {
					ConsVar_file <<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<<
							geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(5)+1 << endl;
					if (config->GetWrite_Mean_Solution()) {
						ConsVar_file_mean <<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(0)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(1)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(2)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(3)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<<
								geometry[iDomain][MESH_0]->elem[iElem]->GetNode(4)+1 <<" "<< geometry[iDomain][MESH_0]->elem[iElem]->GetNode(5)+1 << endl;
					}
				}
			}

		}
		ConsVar_file.close();
		if (config->GetWrite_Mean_Solution())
			ConsVar_file_mean.close();
	}
}


void COutput::SetDomain_Adjoint(CConfig *config, CGeometry *geometry, CSolution ***solution_container, unsigned long iExtIter) {
	char cstr[200], buffer [50];
	ofstream AdjVar_file;
	unsigned long iPoint, iElem;
	unsigned short iDim;
	double coord[3] = {0.0,0.0,0.0}, phi[3] = {0.0,0.0,0.0};

	unsigned short solver = config->GetKind_Solver();
	bool incompressible = config->GetIncompressible();
	bool scalar = true;

	/*--- Paraview Output ---*/
	if (config->GetOutput_FileFormat() == PARAVIEW) {

		/*--- Write file name with extension ---*/
		strcpy (cstr, config->GetAdj_FileName().c_str());
		if (config->GetUnsteady_Simulation() != NO) 
			sprintf (buffer, "_%d.vtk", int(iExtIter));
		else 
			sprintf (buffer, ".vtk");

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
		strcpy (cstr, config->GetAdj_FileName().c_str());
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.plt", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.plt", int(iExtIter));
		}
		else
			sprintf (buffer, ".plt");

		strcat(cstr,buffer);

		/*--- Write header and open the file ---*/
		AdjVar_file.precision(15);
		AdjVar_file.open(cstr, ios::out);

		AdjVar_file << "TITLE = \"Visualization of the volumetric grid\"" << endl;
		AdjVar_file << "VARIABLES = \"x\", \"y\", \"z\"";

		/*--- Add variable labels depending on adjoint solver type. ---*/
		if ((solver == ADJ_EULER) || (solver == ADJ_NAVIER_STOKES) || (solver == ADJ_RANS) ||
				(solver == ADJ_FREE_SURF_EULER) || (solver == ADJ_FREE_SURF_NAVIER_STOKES) || (solver == ADJ_FREE_SURF_RANS)) {
			if (!incompressible) AdjVar_file << ", \"PsiRho\", \"Phi_x\", \"Phi_y\", \"Phi_z\", \"PsiE\"";
			else AdjVar_file << ", \"PsiPress\", \"Phi_x\", \"Phi_y\", \"Phi_z\"";
		}

		if (solver == ADJ_RANS)
			AdjVar_file << ", \"PsiTurb\"";

		if ((solver == ADJ_FREE_SURF_EULER) || (solver == ADJ_FREE_SURF_NAVIER_STOKES))
			AdjVar_file << ", \"PsiLevelSet\", \"Press\", \"RhoVelx\", \"RhoVely\", \"RhoVelz\", \"Density\", \"LevelSet\", \"DiffLevelSet\"";

		if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTRED)
			AdjVar_file << ", \"Diss_Switch\"";

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
			if ((solver == ADJ_EULER) || (solver == ADJ_NAVIER_STOKES) || (solver == ADJ_RANS) ||
					(solver == ADJ_FREE_SURF_EULER) || (solver == ADJ_FREE_SURF_NAVIER_STOKES) || (solver == ADJ_FREE_SURF_RANS)) {
				if (!incompressible) 
					AdjVar_file << " " << solution_container[MESH_0][ADJFLOW_SOL]->node[iPoint]->GetSolution(0)
					<<" "<< phi[0] <<" "<< phi[1] <<" "<< phi[2]<<" "<<
					solution_container[MESH_0][ADJFLOW_SOL]->node[iPoint]->GetSolution(geometry->GetnDim()+1);
				else
					AdjVar_file << " " << solution_container[MESH_0][ADJFLOW_SOL]->node[iPoint]->GetSolution(0)
					<<" "<< phi[0] <<" "<< phi[1] <<" "<< phi[2];				
			}

			/*--- If turbulent, add the turbulent adjoint variable. ---*/
			if (solver == ADJ_RANS) {
				AdjVar_file << " " << solution_container[MESH_0][ADJTURB_SOL]->node[iPoint]->GetSolution(0);
			}

			/*--- Write the adjoint level set solution ---*/
			if ((solver == ADJ_FREE_SURF_EULER) || (solver == ADJ_FREE_SURF_NAVIER_STOKES)) {
				AdjVar_file << " "<< solution_container[MESH_0][ADJLEVELSET_SOL]->node[iPoint]->GetSolution(0)
						<< " "<< solution_container[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution(0)
						<< " "<< solution_container[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution(1)
						<< " "<< solution_container[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution(2)
						<< " "<< solution_container[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution(2)
						<< " "<< solution_container[MESH_0][FLOW_SOL]->node[iPoint]->GetDensityInc()
						<< " "<< solution_container[MESH_0][LEVELSET_SOL]->node[iPoint]->GetSolution(0)
						<< " "<< solution_container[MESH_0][LEVELSET_SOL]->node[iPoint]->GetDiffLevelSet();
			}

			/*--- If rotating, visualize the dissipation sensor ---*/
			if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTRED) {
				AdjVar_file << " " << solution_container[MESH_0][ADJFLOW_SOL]->node[iPoint]->GetSensor();
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

void COutput::SetLinearized_Variables(CConfig *config, CGeometry *geometry, CSolution ***solution_container, unsigned long iExtIter) {
	char cstr[200], buffer[50];
	ofstream LinVar_file;
	bool scalar = true;

	/*--- Write file name with extension ---*/
	strcpy (cstr, config->GetLin_FileName().c_str());
	sprintf (buffer, ".vtk");
	strcat(cstr,buffer);

	/*--- Write header and open the file ---*/
	geometry->SetParaView(cstr);
	LinVar_file.open(cstr, ios::out | ios::app);

	/*--- Write information ---*/
	LinVar_file << "POINT_DATA " << geometry->GetnPoint() << endl;

	WriteInOutputFile(geometry, solution_container[MESH_0][LINFLOW_SOL], LinVar_file,  "DeltaRho", scalar, 0, "GetSolution", config);

	WriteInOutputFile(geometry, solution_container[MESH_0][LINFLOW_SOL], LinVar_file, "DeltaRhoVel", !scalar, 0, "GetSolution", config);

	WriteInOutputFile(geometry, solution_container[MESH_0][LINFLOW_SOL], LinVar_file,   "DeltaRhoE",  scalar,  geometry->GetnDim()+1, "GetSolution", config);
}

void COutput::SetSurface_Flow(CConfig *config, CGeometry *geometry, CSolution *FlowSolution, unsigned long iExtIter) {
	unsigned long iPoint, iVertex;
	unsigned short iMarker, iDim;
	double PressCoeff = 0.0, SkinFrictionCoeff, HeatTransferCoeff,WallTemperature,WallNumDensity, Mach, *aux_press = NULL, *aux_friction = NULL, *aux_rel_mach = NULL, *aux_heat_transfer = NULL;
	double *aux_wall_temperature = NULL, *aux_wall_numdensity = NULL, *aux_wall_numdensity_diff = NULL;
	char cstr[200], buffer [50];
	ofstream SurfFlow_file;

	unsigned short solver = config->GetKind_Solver();
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();

	if (config->GetOutput_FileFormat() == PARAVIEW) {
		bool adiabatic = config->GetAdiabaticWall();

		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
		if (config->GetUnsteady_Simulation() != NO) 
			sprintf (buffer, "_%d.vtk", int(iExtIter));
		else 
			sprintf (buffer, ".vtk");

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

		case NAVIER_STOKES: case RANS: case FREE_SURF_NAVIER_STOKES:
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
			if (!adiabatic) {
				aux_heat_transfer = new double [geometry->GetnPoint()];
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_heat_transfer[iPoint] = 0.0;
				for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
					if (config->GetMarker_All_Plotting(iMarker) == YES)
						for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
							iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
							aux_heat_transfer[iPoint] = FlowSolution->GetHeatTransferCoeff(iMarker,iVertex);
						}
				SurfFlow_file << "SCALARS Heat_Transfer_Coefficient float 1" << endl;
				SurfFlow_file << "LOOKUP_TABLE default" << endl;
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					if (geometry->node[iPoint]->GetBoundary()) {
						HeatTransferCoeff = aux_heat_transfer[iPoint];
						SurfFlow_file << scientific << HeatTransferCoeff << endl;
					}
				delete [] aux_heat_transfer;
			}
			break;
		case NS_PLASMA:

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
			if (!adiabatic) {
				aux_heat_transfer = new double [geometry->GetnPoint()];
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_heat_transfer[iPoint] = 0.0;
				for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
					if (config->GetMarker_All_Plotting(iMarker) == YES)
						for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
							iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
							aux_heat_transfer[iPoint] = FlowSolution->GetHeatTransferCoeff(iMarker,iVertex);
						}
				SurfFlow_file << "SCALARS Heat_Transfer_Coefficient float 1" << endl;
				SurfFlow_file << "LOOKUP_TABLE default" << endl;
				for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
					if (geometry->node[iPoint]->GetBoundary()) {
						HeatTransferCoeff = aux_heat_transfer[iPoint];
						SurfFlow_file << scientific << HeatTransferCoeff << endl;
					}

				delete [] aux_heat_transfer;
			}

			aux_wall_temperature = new double [geometry->GetnPoint()];
			aux_wall_numdensity = new double [geometry->GetnPoint()];
			aux_wall_numdensity_diff = new double [geometry->GetnPoint()];
			double *Mass = new double[3];
			Mass[0] = config->GetParticle_Mass(0);
			Mass[1] = config->GetParticle_Mass(1);
			Mass[2] = config->GetParticle_Mass(2);

			for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				aux_wall_temperature[iPoint] = 0.0;
				aux_wall_numdensity[iPoint] = 0.0;
			}
			for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
				if (config->GetMarker_All_Plotting(iMarker) == YES) {
					/* Argon Temperature and Density */
					for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						aux_wall_temperature[iPoint] = FlowSolution->node[iPoint]->GetTemperature(0);
						aux_wall_numdensity[iPoint] = FlowSolution->node[iPoint]->GetDensity(0)/Mass[0];
					}
					SurfFlow_file << "SCALARS Ar_Temperature float 1" << endl;
					SurfFlow_file << "LOOKUP_TABLE default" << endl;
					for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
						if (geometry->node[iPoint]->GetBoundary()) {
							WallTemperature = aux_wall_temperature[iPoint];
							SurfFlow_file << scientific << WallTemperature << endl;
						}

					SurfFlow_file << "SCALARS Ar_NumDen float 1" << endl;
					SurfFlow_file << "LOOKUP_TABLE default" << endl;
					for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
						if (geometry->node[iPoint]->GetBoundary()) {
							WallNumDensity = aux_wall_numdensity[iPoint];
							SurfFlow_file << scientific << WallNumDensity << endl;
						}

					/* Ion Temperature and Density*/
					for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						aux_wall_temperature[iPoint] = FlowSolution->node[iPoint]->GetTemperature(1);
						aux_wall_numdensity[iPoint] = FlowSolution->node[iPoint]->GetDensity(1)/Mass[1];

					}
					SurfFlow_file << "SCALARS Ion_Temperature float 1" << endl;
					SurfFlow_file << "LOOKUP_TABLE default" << endl;
					for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
						if (geometry->node[iPoint]->GetBoundary()) {
							WallTemperature = aux_wall_temperature[iPoint];
							SurfFlow_file << scientific << WallTemperature << endl;
						}

					SurfFlow_file << "SCALARS Ion_NumDen float 1" << endl;
					SurfFlow_file << "LOOKUP_TABLE default" << endl;
					for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
						if (geometry->node[iPoint]->GetBoundary()) {
							WallNumDensity = aux_wall_numdensity[iPoint];
							SurfFlow_file << scientific << WallNumDensity << endl;
						}

					/* Electron Temperature and Density */
					for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						aux_wall_temperature[iPoint] = FlowSolution->node[iPoint]->GetTemperature(2);
						aux_wall_numdensity[iPoint] = FlowSolution->node[iPoint]->GetDensity(2)/Mass[2];
						aux_wall_numdensity_diff[iPoint] = FlowSolution->node[iPoint]->GetDensity(2)/Mass[2] - FlowSolution->node[iPoint]->GetDensity(1)/Mass[1];
					}
					SurfFlow_file << "SCALARS Ele_Temperature float 1" << endl;
					SurfFlow_file << "LOOKUP_TABLE default" << endl;
					for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
						if (geometry->node[iPoint]->GetBoundary()) {
							WallTemperature = aux_wall_temperature[iPoint];
							SurfFlow_file << scientific << WallTemperature << endl;
						}
					SurfFlow_file << "SCALARS Ele_NumDen float 1" << endl;
					SurfFlow_file << "LOOKUP_TABLE default" << endl;
					for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
						if (geometry->node[iPoint]->GetBoundary()) {
							WallNumDensity = aux_wall_numdensity[iPoint];
							SurfFlow_file << scientific << WallNumDensity << endl;
						}

					SurfFlow_file << "SCALARS NumofNegCharge float 1" << endl;
					SurfFlow_file << "LOOKUP_TABLE default" << endl;
					for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
						if (geometry->node[iPoint]->GetBoundary()) {
							WallNumDensity = aux_wall_numdensity_diff[iPoint];
							SurfFlow_file << scientific << WallNumDensity << endl;
						}




					delete [] aux_wall_temperature;
					delete [] aux_wall_numdensity;
					delete [] aux_wall_numdensity_diff;
				}
			break;
		}

		SurfFlow_file.close();

	}
	if (config->GetOutput_FileFormat() == TECPLOT) {

		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.plt", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.plt", int(iExtIter));
		}
		else
			sprintf (buffer, ".plt");

		strcat(cstr, buffer);

		SurfFlow_file.precision(15);
		SurfFlow_file.open(cstr, ios::out);


		SurfFlow_file << "TITLE = \"Visualization of the surface grid\"" << endl;
		if (geometry->GetnDim() == 2) SurfFlow_file << "VARIABLES = \"x\", \"y\"";
		if (geometry->GetnDim() == 3) SurfFlow_file << "VARIABLES = \"x\", \"y\", \"z\"";

		if ((solver == EULER) || (solver == FREE_SURF_EULER) 
				|| (solver == FREE_SURF_NAVIER_STOKES) || (solver == NAVIER_STOKES) || (solver == RANS))
			if (incompressible) SurfFlow_file << ", \"Pressure\"";
			else SurfFlow_file << ", \"Pressure Coefficient\"";

		if ((solver == EULER) || (solver == FREE_SURF_EULER))
			SurfFlow_file << ", \"Mach Number\"";

		if ((solver == NAVIER_STOKES) || (solver == RANS) || (solver == NAVIER_STOKES))
			SurfFlow_file << ", \"Skin Friction Coefficient\"";

		if (rotating_frame && (solver == EULER))
			SurfFlow_file << ", \"Relative Mach\"";

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

		SurfFlow_file << "ZONE ";
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
			SurfFlow_file << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*int(iExtIter+1)<<", ";

		SurfFlow_file << "NODES="<< nPointSurface <<" , ELEMENTS="<< nElemSurface <<", DATAPACKING=POINT";
		if (geometry->GetnDim() == 2) SurfFlow_file << ", ZONETYPE=FELINESEG" << endl;
		if (geometry->GetnDim() == 3) SurfFlow_file << ", ZONETYPE=FEQUADRILATERAL"<< endl;

		/*--- It is necessary to go from point to vertex, and an auxiliar variable is created ---*/
		aux_press = new double [geometry->GetnPoint()];
		aux_friction = new double [geometry->GetnPoint()];
		if (rotating_frame && (solver == EULER)) aux_rel_mach = new double [geometry->GetnPoint()];

		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_press[iPoint] = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					aux_press[iPoint] = FlowSolution->GetCPressure(iMarker,iVertex);
					if (incompressible)
						aux_press[iPoint] = FlowSolution->node[iPoint]->GetSolution((unsigned short)0)*config->GetPressure_Ref();
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

		/*--- Compute relative velocities for rotating frame Euler  ---*/
		if (rotating_frame && (solver == EULER)) {
			double rel_velocity_mag, rel_velocity[3] = {0.0,0.0,0.0};
			for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_rel_mach[iPoint] = 0.0;
			for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
				if (config->GetMarker_All_Plotting(iMarker) == YES)
					for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						double *rot_vel = geometry->node[iPoint]->GetRotVel();
						rel_velocity_mag = 0.0;
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
							rel_velocity[iDim] = FlowSolution->node[iPoint]->GetVelocity(iDim) - rot_vel[iDim];
							rel_velocity_mag  += rel_velocity[iDim]*rel_velocity[iDim];
						}
						aux_rel_mach[iPoint]  = sqrt(rel_velocity_mag)/FlowSolution->node[iPoint]->GetSoundSpeed();
					}
		}


		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

			/*--- Write the coordinates ---*/
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				if (geometry->node[iPoint]->GetBoundary()) 
					SurfFlow_file <<" "<< geometry->node[iPoint]->GetCoord(iDim);

			/*--- Write the Euler variables ---*/
			if ((solver == EULER) ||
					(solver == NAVIER_STOKES) || (solver == RANS) || (solver == FREE_SURF_EULER) || 
					(solver == FREE_SURF_NAVIER_STOKES)) {
				PressCoeff = aux_press[iPoint];
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << PressCoeff;
			}

			/*--- Write the Mach number ---*/
			if ((solver == EULER || (solver == FREE_SURF_EULER))) {
				if (incompressible) {
					double Density = FlowSolution->node[iPoint]->GetDensityInc()*config->GetDensity_Ref();
					double SoundSpeed = sqrt(config->GetBulk_Modulus()/Density);					
					Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/SoundSpeed;
				}
				else Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2())/FlowSolution->node[iPoint]->GetSoundSpeed();
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << Mach;
			}

			/*--- Write the skin friction coefficient variables ---*/
			if ((solver == NAVIER_STOKES) || (solver == RANS)) {
				SkinFrictionCoeff = aux_friction[iPoint];
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << SkinFrictionCoeff;
			}

			/*--- Write the relative mach number ---*/
			if (rotating_frame && (solver == EULER)) {
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << aux_rel_mach[iPoint];
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
		if (rotating_frame && (solver == EULER)) delete[] aux_rel_mach;
		delete [] PointSurface;
		SurfFlow_file.close();

	}

}

void COutput::SetSurface_Adjoint(CConfig *config, CGeometry *geometry, CSolution *AdjSolution, CSolution *FlowSolution, unsigned long iExtIter) {

	unsigned long iPoint, iVertex, iElem;
	double *aux_sens, sens_value;
	unsigned short iMarker, iDim;
	char cstr[200], buffer[50];
	ofstream SurfAdj_file;

	if (config->GetOutput_FileFormat() == PARAVIEW) {

		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());
		if (config->GetUnsteady_Simulation() != NO) 
			sprintf (buffer, "_%d.vtk", int(iExtIter));
		else 
			sprintf (buffer, ".vtk");

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

		/*--- Write the surface .vtk file ---*/
		strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.plt", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.plt", int(iExtIter));
		}
		else
			sprintf (buffer, ".plt");

		strcat(cstr, buffer);

		/*--- Open the output file ---*/
		SurfAdj_file.precision(15);
		SurfAdj_file.open(cstr, ios::out);

		SurfAdj_file << "TITLE = \"Visualization of the surface grid\"" << endl;
		if (geometry->GetnDim() == 2) SurfAdj_file << "VARIABLES = \"x\", \"y\"";
		if (geometry->GetnDim() == 3) SurfAdj_file << "VARIABLES = \"x\", \"y\", \"z\"";
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

			/*--- Write the coordinates ---*/
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				if (geometry->node[iPoint]->GetBoundary()) 
					SurfAdj_file <<" "<< geometry->node[iPoint]->GetCoord(iDim);			

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

void COutput::SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry, CSolution *FlowSolution, unsigned long iExtIter) {

#ifdef NO_MPI
	unsigned long iPoint, iVertex;
	unsigned short iMarker;
	double PressCoeff = 0.0, SkinFrictionCoeff, xCoord, yCoord, zCoord, Mach;
	char cstr[200], buffer [50];
	ofstream SurfFlow_file;
	unsigned short solver = config->GetKind_Solver();

	/*--- Write the surface .csv file ---*/
	strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
	sprintf (buffer, ".csv");
	strcat(cstr, buffer);

	SurfFlow_file.precision(15);
	SurfFlow_file.open(cstr, ios::out);

	if (geometry->GetnDim() == 2) {
		switch (solver) {
		case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\"" << endl; break;
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
					case EULER :
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
		case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\"" << endl; break;
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
					case EULER :
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

#else

	int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();
	double PressCoeff = 0.0, SkinFrictionCoeff, xCoord, yCoord, zCoord, Mach;
	char cstr[200];
	unsigned short iMarker;
	unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
			MaxLocalVertex_Surface = 0, nBuffer_Scalar, *Buffer_Receive_nVertex = NULL, position;
	ofstream SurfFlow_file;

	/*--- Write the surface .csv file, the information of allthe vertices is send to the MASTER_NODE ---*/
	nLocalVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
			}

	if (rank == MASTER_NODE)
		Buffer_Receive_nVertex = new unsigned long [nProcessor];

	Buffer_Send_nVertex[0] = nLocalVertex_Surface;

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);

	double *Buffer_Send_Coord_x = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_y = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_z = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Press = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_Mach = new double [MaxLocalVertex_Surface];
	double *Buffer_Send_SkinFriction = new double [MaxLocalVertex_Surface];

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
					if (config->GetKind_Solver() == EULER)
						Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolution->node[iPoint]->GetVelocity2()) / FlowSolution->node[iPoint]->GetSoundSpeed();
					if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS))
						Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolution->GetCSkinFriction(iMarker,iVertex);
					nVertex_Surface++;
				}
			}

	double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Press = NULL,
			*Buffer_Receive_Mach = NULL, *Buffer_Receive_SkinFriction = NULL;

	if (rank == MASTER_NODE) {
		Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
		if (geometry->GetnDim() == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Press = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Mach = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_SkinFriction = new double [nProcessor*MaxLocalVertex_Surface];
	}

	nBuffer_Scalar = MaxLocalVertex_Surface;

	/*--- Send the information to the Master node ---*/
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (geometry->GetnDim() == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Press, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Press, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (config->GetKind_Solver() == EULER)
		MPI::COMM_WORLD.Gather(Buffer_Send_Mach, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Mach, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS))
		MPI::COMM_WORLD.Gather(Buffer_Send_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);

	/*--- The master node is the one who writes the surface files ---*/
	if (rank == MASTER_NODE) {

		strcpy (cstr, config->GetSurfFlowCSV_FileName().c_str());

		SurfFlow_file.precision(15);
		SurfFlow_file.open(cstr, ios::out);

		/*--- Write the 2D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 2) {
			switch (config->GetKind_Solver()) {
			case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\"" << endl; break;
			case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\"" << endl; break;
			}
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					PressCoeff = Buffer_Receive_Press[position];
					switch (config->GetKind_Solver()) {
					case EULER :
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
			switch (config->GetKind_Solver()) {
			case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\"" << endl; break;
			case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"z_coord\"" << endl; break;
			}
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					zCoord = Buffer_Receive_Coord_z[position];
					PressCoeff = Buffer_Receive_Press[position];
					switch (config->GetKind_Solver()) {
					case EULER :
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

	delete [] Buffer_Send_Coord_x;
	delete [] Buffer_Send_Coord_y;
	delete [] Buffer_Send_Coord_z;
	delete [] Buffer_Send_Press;
	delete [] Buffer_Send_Mach;
	delete [] Buffer_Send_SkinFriction;

	SurfFlow_file.close();

#endif
}

void COutput::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolution *AdjSolution, CSolution *FlowSolution, unsigned long iExtIter) {

#ifdef NO_MPI

	unsigned long iPoint, iVertex;
	double *Solution, xCoord, yCoord, zCoord, *IntBoundary_Jump;
	unsigned short iMarker;
	char cstr[200], buffer[50];
	ofstream SurfAdj_file;

	/*--- Write the surface .csv file ---*/
	strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());
	sprintf (buffer, ".csv");
	strcat (cstr, buffer);

	SurfAdj_file.precision(15);
	SurfAdj_file.open(cstr, ios::out);

	if (geometry->GetnDim() == 2) {
		SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					Solution = AdjSolution->node[iPoint]->GetSolution();
					IntBoundary_Jump = AdjSolution->node[iPoint]->GetIntBoundary_Jump();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					SurfAdj_file << scientific << iPoint << ", " << AdjSolution->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
							<< Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord << endl;
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
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					zCoord = geometry->node[iPoint]->GetCoord(2);

					SurfAdj_file << scientific << iPoint << ", " << AdjSolution->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
							<< Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord << endl;
				}
		}
	}

	SurfAdj_file.close();

#else

	int rank = MPI::COMM_WORLD.Get_rank(), iProcessor, nProcessor = MPI::COMM_WORLD.Get_size();
	unsigned short nDim = geometry->GetnDim(), iMarker;
	double *Solution, *Normal, *d, *Coord;
	unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
			MaxLocalVertex_Surface = 0, nBuffer_Scalar;
	unsigned long *Buffer_Receive_nVertex = NULL;
	ofstream SurfAdj_file;

	/*--- Write the surface .csv file ---*/
	nLocalVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
			}

	if (rank == MASTER_NODE)
		Buffer_Receive_nVertex = new unsigned long [nProcessor];

	Buffer_Send_nVertex[0] = nLocalVertex_Surface;

	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG, MASTER_NODE);

	double *Buffer_Send_Coord_x = new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_y= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Coord_z= new double[MaxLocalVertex_Surface];
	unsigned long *Buffer_Send_GlobalPoint= new unsigned long[MaxLocalVertex_Surface];
	double *Buffer_Send_Sensitivity= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_PsiRho= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Phi_x= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Phi_y= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_Phi_z= new double[MaxLocalVertex_Surface];
	double *Buffer_Send_PsiE= new double[MaxLocalVertex_Surface];

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

	double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
			*Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
			*Buffer_Receive_PsiE = NULL;
	unsigned long *Buffer_Receive_GlobalPoint = NULL;

	if (rank == MASTER_NODE) {
		Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
		if (nDim == 3) Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Sensitivity = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_PsiRho = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Phi_x = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Phi_y = new double [nProcessor*MaxLocalVertex_Surface];
		if (nDim == 3) Buffer_Receive_Phi_z = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_PsiE = new double [nProcessor*MaxLocalVertex_Surface];
	}

	nBuffer_Scalar = MaxLocalVertex_Surface;

	/*--- Send the information to the Master node ---*/
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (nDim == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI::UNSIGNED_LONG, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_PsiRho, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Phi_x, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Phi_y, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (nDim == 3) MPI::COMM_WORLD.Gather(Buffer_Send_Phi_z, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_PsiE, nBuffer_Scalar, MPI::DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);

	/*--- The master node is the one who writes the surface files ---*/
	if (rank == MASTER_NODE) {
		unsigned long iVertex, GlobalPoint, position;
		char cstr[200];
		ofstream SurfAdj_file;

		strcpy (cstr, config->GetSurfAdjCSV_FileName().c_str());
		SurfAdj_file.open(cstr, ios::out);
		SurfAdj_file.precision(15);

		/*--- Write the 2D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 2) {

			SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;

			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
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

			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
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

	delete [] Buffer_Send_Coord_x;
	delete [] Buffer_Send_Coord_y;
	delete [] Buffer_Send_Coord_z;
	delete [] Buffer_Send_GlobalPoint;
	delete [] Buffer_Send_Sensitivity;
	delete [] Buffer_Send_PsiRho;
	delete [] Buffer_Send_Phi_x;
	delete [] Buffer_Send_Phi_y;
	delete [] Buffer_Send_Phi_z;
	delete [] Buffer_Send_PsiE;

	SurfAdj_file.close();	

#endif
}

void COutput::SetSurfaceCSV_Linearized(CConfig *config, CGeometry *geometry, CSolution *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::SetReStart(CConfig *config, CGeometry *geometry, CSolution **solution, string mesh_filename) {
	unsigned short iVar, nVar_First = 0, nVar_Second = 0, FirstIndex = NONE, SecondIndex = NONE;
	unsigned long iPoint;
	ofstream restart_file;
	string filename, AdjExt;

	unsigned short KindSolver = config->GetKind_Solver();
	unsigned short KindOF = config->GetKind_ObjFunc();
	bool AdjProblem = ((KindSolver == ADJ_EULER) || (KindSolver == ADJ_NAVIER_STOKES) || (KindSolver == ADJ_RANS) ||
			(KindSolver == ADJ_FREE_SURF_EULER) || (KindSolver == ADJ_FREE_SURF_NAVIER_STOKES) || (KindSolver == ADJ_FREE_SURF_RANS));
	bool Write_Mean_Solution = false;

	switch (KindSolver) {
	case EULER : case NAVIER_STOKES: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
	case NS_PLASMA : FirstIndex = PLASMA_SOL; SecondIndex = NONE; Write_Mean_Solution = config->GetWrite_Mean_Solution(); break;
	case RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
	case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES: FirstIndex = FLOW_SOL; SecondIndex = LEVELSET_SOL; break;
	case ELECTRIC_POTENTIAL: FirstIndex = ELEC_SOL; SecondIndex = NONE; break;
	case WAVE: FirstIndex = WAVE_SOL; SecondIndex = NONE; break;
	case ADJ_EULER : case ADJ_NAVIER_STOKES : FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; break;
	case ADJ_FREE_SURF_EULER : case ADJ_FREE_SURF_NAVIER_STOKES : FirstIndex = ADJFLOW_SOL; SecondIndex = ADJLEVELSET_SOL; break;
	case ADJ_RANS : FirstIndex = ADJFLOW_SOL; SecondIndex = ADJTURB_SOL; break;
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
		case FORCE_X_COEFFICIENT: AdjExt = "_cfx.dat"; break;
		case FORCE_Y_COEFFICIENT: AdjExt = "_cfy.dat"; break;
		case FORCE_Z_COEFFICIENT: AdjExt = "_cfz.dat"; break;
		case THRUST_COEFFICIENT: AdjExt = "_ct.dat"; break;
		case TORQUE_COEFFICIENT: AdjExt = "_cq.dat"; break;
		case FIGURE_OF_MERIT: AdjExt = "_merit.dat"; break;
		case FREESURFACE: AdjExt = "_fs.dat"; break;
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

	/*--- Write additional restart file for mass averaged properties (Multi-species only) ---*/
	if (Write_Mean_Solution) {
		string filename_mean;
		unsigned short iSpecies, iDim, loc;
		unsigned short nDiatomics, nDim, nSpecies;
		double *Mass_Ave_Solution, *Mass_Ave_Residual;
		double Mass_Fraction, Mixture_Density;

		/*--- Create a copy of the filename ---*/
		filename_mean.assign(mesh_filename);

		filename_mean.erase (filename_mean.end()-4, filename_mean.end());
		filename_mean.append("_mean.dat");
		restart_file.open(filename_mean.data(), ios::out);
		restart_file.precision(15);

		nDiatomics = solution[FirstIndex]->GetnDiatomics();
		nDim = geometry->GetnDim();
		nSpecies = solution[FirstIndex]->GetnSpecies();
		Mass_Ave_Solution = new double [nDim+2];
		Mass_Ave_Residual = new double [nDim+2];

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

			/*--- Index of the point ---*/
			restart_file << iPoint << "\t";

			for (iVar = 0; iVar < nDim+2; iVar++) {
				Mass_Ave_Solution[iVar] = 0.0;
				Mass_Ave_Residual[iVar] = 0.0;
			}

			/*--- Pre-process by mass-averaging ---*/
			Mixture_Density = 0.0;
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				Mixture_Density += solution[FirstIndex]->node[iPoint]->GetSolution(loc+0);
			}	
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				Mass_Fraction = solution[FirstIndex]->node[iPoint]->GetSolution(loc+0) / Mixture_Density;
				Mass_Ave_Solution[0] += solution[FirstIndex]->node[iPoint]->GetSolution(loc+0);
				Mass_Ave_Residual[0] += solution[FirstIndex]->node[iPoint]->GetResidual(loc+0);
				for (iDim = 0; iDim < nDim; iDim++) {
					Mass_Ave_Solution[iDim+1] += solution[FirstIndex]->node[iPoint]->GetSolution(loc+1+iDim)*Mass_Fraction;
					Mass_Ave_Residual[iDim+1] += solution[FirstIndex]->node[iPoint]->GetResidual(loc+1+iDim)*Mass_Fraction;
				}
				Mass_Ave_Solution[nDim+1] += solution[FirstIndex]->node[iPoint]->GetSolution(loc+1+nDim)*Mass_Fraction;
				Mass_Ave_Residual[nDim+1] += solution[FirstIndex]->node[iPoint]->GetResidual(loc+1+nDim)*Mass_Fraction;
			}

			/*--- Mass-averaged Solution (first, and second system of equations) ---*/
			for (iVar = 0; iVar < nDim+2; iVar++)
				restart_file << scientific << Mass_Ave_Solution[iVar] << "\t";

			for (iVar = 0; iVar < nVar_Second; iVar++)
				restart_file << scientific << solution[SecondIndex]->node[iPoint]->GetSolution(iVar) << "\t";

			/*--- Mass-averaged Residual (first, and second system of equations) ---*/
			for (iVar = 0; iVar < nDim+2; iVar++)
				restart_file << scientific << Mass_Ave_Residual[iVar] << "\t";

			for (iVar = 0; iVar < nVar_Second; iVar++)
				restart_file << scientific << solution[SecondIndex]->node[iPoint]->GetResidual(iVar) << "\t";

			restart_file << geometry->node[iPoint]->GetVolume() << "\t";
			restart_file << endl;
		}
	}

}

void COutput::SetHistory_Header(ofstream *ConvHist_file, CConfig *config) {
	char cstr[200], buffer[50];

	bool rotating_frame = config->GetRotating_Frame();
	bool equiv_area = config->GetEquivArea();
	bool turbulent = (config->GetKind_Solver() == RANS);

	/*--- Write file name with extension ---*/
	strcpy (cstr, config->GetConv_FileName().data());
	if (config->GetOutput_FileFormat() == PARAVIEW)  sprintf (buffer, ".csv");
	if (config->GetOutput_FileFormat() == TECPLOT)  sprintf (buffer, ".plt");
	strcat(cstr,buffer);

	ConvHist_file->open(cstr, ios::out);
	ConvHist_file->precision(15);

	char begin[]= "\"Iteration\"";

	char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
	char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldPress\"";
	char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
	char free_surface_coeff[]= ",\"CFreeSurface\"";
	char plasma_coeff[]= ",\"CLift\",\"CDrag\"";
	char wave_coeff[]= ",\"CWave\"";
	char adj_coeff[]= ",\"CSens_Geo\",\"CSens_Mach\",\"CSens_AoA\",\"CSens_AoS\"";

	char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
	char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";

	char plasma_resid[]= ",\"Res_Plasma[0]\",\"Res_Plasma[1]\",\"Res_Plasma[2]\",\"Res_Plasma[3]\",\"Res_Plasma[4]\"";
	char adj_plasma_resid[]= ",\"Res_AdjPlasma[0]\",\"Res_AdjPlasma[1]\",\"Res_AdjPlasma[2]\",\"Res_AdjPlasma[3]\",\"Res_AdjPlasma[4]\"";

	char turb_resid[]= ",\"Res_Turb[0]\"";
	char adj_turb_resid[]= ",\"Res_AdjTurb[0]\"";

	char levelset_resid[]= ",\"Res_LevelSet\"";
	char adj_levelset_resid[]= ",\"Res_AdjLevelSet\"";

	char wave_resid[]= ",\"Res_Wave\"";

	char end[]= ",\"Time(min)\"\n";

	if (config->GetOutput_FileFormat() == TECPLOT) {
		char title[]= "TITLE = \"SU2 Simulation\"";
		char variables[]= "VARIABLES = ";
		ConvHist_file[0] << title << endl;
		ConvHist_file[0] << variables;
	}

	switch (config->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES: case RANS :
		ConvHist_file[0] << begin << flow_coeff;
		if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
		if (rotating_frame) ConvHist_file[0] << rotating_frame_coeff;
		ConvHist_file[0] << flow_resid;
		if (turbulent) ConvHist_file[0] << turb_resid;
		ConvHist_file[0] << end;
		break;

	case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
		ConvHist_file[0] << begin << flow_coeff << free_surface_coeff;
		ConvHist_file[0] << flow_resid << levelset_resid << end;
		break;

	case NS_PLASMA:
		ConvHist_file[0] << begin << plasma_coeff;
		ConvHist_file[0] << plasma_resid << end;
		break;

	case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
		ConvHist_file[0] << begin << adj_coeff << flow_coeff << adj_flow_resid;
		if (turbulent) ConvHist_file[0] << adj_turb_resid;
		ConvHist_file[0] << end;
		break;

	case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:
		ConvHist_file[0] << begin << adj_coeff << free_surface_coeff;
		ConvHist_file[0] << adj_flow_resid << adj_levelset_resid << end;
		break;

	case ADJ_NS_PLASMA:
		ConvHist_file[0] << begin << plasma_coeff;
		ConvHist_file[0] << adj_plasma_resid << end;
		break;

	case WAVE:
		ConvHist_file[0] << begin << wave_coeff;
		ConvHist_file[0] << wave_resid << end;
		break;
	}

	if (config->GetOutput_FileFormat() == TECPLOT) {
		char zone[]= "ZONE T= \"Convergence history\"";
		ConvHist_file[0] << zone << endl;
	}

}

void COutput::SetHistory_DualTime(CGeometry **geometry, CSolution ***solution_container, CConfig *config, 
		CIntegration **integration, unsigned long iExtIter) {
#ifdef NO_MPI

	bool write_heads = ((iExtIter % (config->GetWrt_Con_Freq()*100)) == 0);
	bool write_solution = ((iExtIter % (config->GetWrt_Con_Freq()*1)) == 0);
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();

	/*--- Write the screen header---*/
	if ((iExtIter == 1) || (write_heads)) {

		switch (config->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			if (incompressible) cout << endl << " DT Iter" << "      Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			else if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
			else cout << endl << " DT Iter" << "      Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			if (incompressible) cout << endl << " DT Iter" << "      Res[AdjRho]" << "   CSens_Geo" << "   CSens_Mach" << endl;
			else if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[AdjRho]" << "     Res[AdjRhoE]" << " CSens_Geo" << " CSens_Mach" << endl;
			else cout << endl << " DT Iter" << "      Res[AdjRho]" << "     Res[AdjRhoE]" << "   CSens_Geo" << "   CSens_Mach" << endl;
			break;
		case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "    Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "      CLevelSet" << endl;
			break;
		case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "    Res[AdjPress]" << "     Res[AdjDist]" << "      CSens_Geo" << "     CSens_Mach" << endl;
			break;
		case RANS :
			cout << endl << " Min DT: " << solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]" << " CThrust(Total)" << " CTorque(Total)" << endl;
			else cout << endl << " DT Iter" << "      Res[Rho]" << "       Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			break;
		case NS_PLASMA :
			cout << endl << " Min DT: " << solution_container[MESH_0][PLASMA_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[MESH_0][PLASMA_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "      Res[Rho_0]" << "       Res[E_0]" << endl;
			break;
		case WAVE :
			cout << endl << " Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "      Res[Wave]" << "      CWave" << endl;
			break;
		}
	}

	/*--- Write the solution on the screen and history file ---*/
	if ((iExtIter != 0) && (write_solution)) {
		switch (config->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[MESH_0][FLOW_SOL]->GetRes_Max(0));
			if (!incompressible) {
				cout.width(14); cout << log10(solution_container[MESH_0][FLOW_SOL]->GetRes_Max(geometry[0]->GetnDim()+1));
			}
			if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) {
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CT();
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CQ();
			} else {
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag();
			}
			cout << endl;
			break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[MESH_0][ADJFLOW_SOL]->GetRes_Max(0));
			if (!incompressible) {
				cout.width(14); cout << log10(solution_container[MESH_0][ADJFLOW_SOL]->GetRes_Max(geometry[0]->GetnDim()+1));
			}
			cout.width(15); cout << solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Geo();
			cout.width(15); cout << solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Mach();
			cout << endl;
			break;
		case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[MESH_0][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[MESH_0][LEVELSET_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[MESH_0][LEVELSET_SOL]->GetTotal_CFreeSurface();
			cout << endl;
			break;
		case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(17); cout << log10(solution_container[MESH_0][ADJFLOW_SOL]->GetRes_Max(0));
			cout.width(17); cout << log10(solution_container[MESH_0][ADJLEVELSET_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Geo();
			cout.width(15); cout << solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Mach();
			cout << endl;
			break;
		case RANS :
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[MESH_0][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[MESH_0][TURB_SOL]->GetRes_Max(0));
			if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) {
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CT();
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CQ();
			} else {
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
				cout.width(15); cout << solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag();
			}
			cout << endl;
			break;
		case NS_PLASMA:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			if (config->GetKind_GasModel() == AIR7) {
				cout.width(14); cout << log10(solution_container[MESH_0][PLASMA_SOL]->GetRes_Max(0));
				cout.width(14); cout << log10(solution_container[MESH_0][PLASMA_SOL]->GetRes_Max(3));
			}
			cout << endl;
			break;
		case WAVE:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[MESH_0][WAVE_SOL]->GetRes_Max(0));
			cout << endl;
			break;
		}
	}
	cout.unsetf(ios::fixed);

#else

	double *sbuf_residual_flow = NULL, *sbuf_residual_turbulent = NULL, *sbuf_residual_levelset = NULL, *sbuf_force = NULL, *sbuf_time = NULL, 
			*sbuf_residual_plasma = NULL, *rbuf_residual_flow = NULL, *rbuf_residual_turbulent = NULL, *rbuf_residual_levelset = NULL, *rbuf_force = NULL, *rbuf_time = NULL,
			*sbuf_residual_adjoint = NULL, *rbuf_residual_adjoint = NULL, *sbuf_residual_linearized = NULL,
			*rbuf_residual_linearized = NULL, *rbuf_residual_plasma = NULL;
	unsigned short iVar, buf_convergence = 0, *sbuf_conv = NULL, *rbuf_conv = NULL;

	bool write_heads = ((iExtIter % (config->GetWrt_Con_Freq()*100)) == 0);
	bool write_solution = ((iExtIter % (config->GetWrt_Con_Freq()*1)) == 0);
	bool compressible = !config->GetIncompressible();
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	unsigned short nDim = geometry[MESH_0]->GetnDim();

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();


	/*--- It is important to improve this, the problem is that the master
	 node does not know the number of variables ---*/
	unsigned short nVar_Flow;
	unsigned short nVar_Adj;
	unsigned short nVar_Lin;

	if (compressible) {
		nVar_Flow = nDim+2;
		nVar_Adj = nDim+2;
		nVar_Lin = nDim+2;
	}
	else {
		nVar_Flow = nDim+1;
		nVar_Adj = nDim+1;
		nVar_Lin = nDim+1;		
	}
	unsigned short nVar_Turb = 1;
	unsigned short nVar_LevelSet = 1;
	unsigned short nVar_Force = 16;
	unsigned short nVar_Time = 2;
	unsigned short nVar_Conv = 1;
	unsigned short nVar_Plasma = config->GetnSpecies()*(nDim+2);

	/*--- Allocate memory for send buffer ---*/
	sbuf_residual_flow = new double[nVar_Flow];
	sbuf_residual_adjoint = new double[nVar_Adj];
	sbuf_residual_linearized = new double[nVar_Lin];
	sbuf_residual_turbulent = new double[nVar_Turb];
	sbuf_residual_levelset = new double[nVar_LevelSet];
	sbuf_residual_plasma = new double[nVar_Plasma];
	sbuf_force = new double[nVar_Force];
	sbuf_time = new double[nVar_Time];
	sbuf_conv = new unsigned short[nVar_Conv];

	/*--- Initialise master node ---*/
	if (rank == MASTER_NODE) {
		rbuf_residual_flow = new double[nVar_Flow];
		rbuf_residual_adjoint = new double[nVar_Adj];
		rbuf_residual_linearized = new double[nVar_Lin];
		rbuf_residual_turbulent = new double[nVar_Turb];
		rbuf_residual_levelset = new double[nVar_LevelSet];
		rbuf_residual_plasma = new double[nVar_Plasma];
		rbuf_force = new double[nVar_Force];
		rbuf_time = new double[nVar_Time];
		rbuf_conv = new unsigned short[nVar_Conv];

		for (iVar = 0; iVar < nVar_Flow; iVar++) {
			rbuf_residual_flow[iVar] = 0.0; sbuf_residual_flow[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Adj; iVar++) {
			rbuf_residual_adjoint[iVar] = 0.0; sbuf_residual_adjoint[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Lin; iVar++) {
			rbuf_residual_linearized[iVar] = 0.0; sbuf_residual_linearized[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Turb; iVar++) {
			rbuf_residual_turbulent[iVar] = 0.0; sbuf_residual_turbulent[iVar] = -20; }
		for (iVar = 0; iVar < nVar_Plasma; iVar++) {
			rbuf_residual_plasma[iVar] = 0.0; rbuf_residual_plasma[iVar] = -20; }
		for (iVar = 0; iVar < nVar_LevelSet; iVar++) {
			rbuf_residual_levelset[iVar] = 0.0; rbuf_residual_levelset[iVar] = -20; }		
		for (iVar = 0; iVar < nVar_Force; iVar++) {
			rbuf_force[iVar] = 0.0; sbuf_force[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_Time; iVar++) {
			rbuf_time[iVar] = 0.0; sbuf_time[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_Conv; iVar++) {
			rbuf_conv[iVar] = 0; sbuf_conv[iVar] = 0; }
	}

	/*--- Write information from nodes ---*/
	switch (config->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES: case RANS :

		for (iVar = 0; iVar < nVar_Flow; iVar++)
			sbuf_residual_flow[iVar] = solution_container[MESH_0][FLOW_SOL]->GetRes_Max(iVar);

		if (config->GetKind_Solver() == RANS) {
			for (iVar = 0; iVar < nVar_Turb; iVar++)
				sbuf_residual_turbulent[iVar] = log10 (solution_container[MESH_0][TURB_SOL]->GetRes_Max(iVar));
		}

		sbuf_force[0]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
		sbuf_force[1]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag();
		sbuf_force[2]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CSideForce();
		sbuf_force[3]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CPress();
		sbuf_force[4]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMx();
		sbuf_force[5]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMy();
		sbuf_force[6]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMz();
		sbuf_force[7]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CEff();
		sbuf_force[8]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CEquivArea(),
				sbuf_force[9]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CNearFieldPress();
		sbuf_force[10] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CFx();
		sbuf_force[11] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CFy();
		sbuf_force[12] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CFz();
		sbuf_force[13] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CT();
		sbuf_force[14] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CQ();
		sbuf_force[15] = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMerit();
		sbuf_conv[0] = integration[FLOW_SOL]->GetConvergence();
		break;

	case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES: case FREE_SURF_RANS:

		for (iVar = 0; iVar < nVar_Flow; iVar++)
			sbuf_residual_flow[iVar] = solution_container[MESH_0][FLOW_SOL]->GetRes_Max(iVar);

		for (iVar = 0; iVar < nVar_LevelSet; iVar++)
			sbuf_residual_levelset[iVar] = solution_container[MESH_0][LEVELSET_SOL]->GetRes_Max(iVar);

		sbuf_force[0]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
		sbuf_force[1]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag();
		sbuf_force[2]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CSideForce();
		sbuf_force[3]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CPress();
		sbuf_force[4]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMx();
		sbuf_force[5]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMy();
		sbuf_force[6]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CMz();
		sbuf_force[7]  = solution_container[MESH_0][FLOW_SOL]->GetTotal_CEff();
		sbuf_force[8]  = solution_container[MESH_0][LEVELSET_SOL]->GetTotal_CFreeSurface();
		sbuf_conv[0] = integration[FLOW_SOL]->GetConvergence();
		break;

	case NS_PLASMA:

		for (iVar = 0; iVar < nVar_Plasma; iVar++)
			sbuf_residual_plasma[iVar] = log10 (solution_container[MESH_0][PLASMA_SOL]->GetRes_Max(iVar));

		sbuf_conv[0] = integration[PLASMA_SOL]->GetConvergence();
		break;

	case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :

		for (iVar = 0; iVar < nVar_Adj; iVar++)
			sbuf_residual_adjoint[iVar] = log10 (solution_container[MESH_0][ADJFLOW_SOL]->GetRes_Max(iVar));

		sbuf_force[0] = solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Geo();
		sbuf_force[1] = solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Mach();
		sbuf_conv[0] = integration[ADJFLOW_SOL]->GetConvergence();
		break;

	case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:

		for (iVar = 0; iVar < nVar_Adj; iVar++)
			sbuf_residual_adjoint[iVar] = log10 (solution_container[MESH_0][ADJFLOW_SOL]->GetRes_Max(iVar));

		sbuf_force[0] = solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Geo();
		sbuf_force[1] = solution_container[MESH_0][ADJFLOW_SOL]->GetTotal_CSens_Mach();
		sbuf_conv[0] = integration[ADJFLOW_SOL]->GetConvergence();
		break;

	case LIN_EULER : case LIN_NAVIER_STOKES :

		for (iVar = 0; iVar < nVar_Lin; iVar++)
			sbuf_residual_linearized[iVar] = log10 (solution_container[MESH_0][LINFLOW_SOL]->GetRes_Max(iVar));

		sbuf_force[0] = solution_container[MESH_0][LINFLOW_SOL]->GetTotal_CDeltaLift();
		sbuf_force[1] = solution_container[MESH_0][LINFLOW_SOL]->GetTotal_CDeltaDrag();
		sbuf_conv[0] = integration[LINFLOW_SOL]->GetConvergence();
		break;
	}
	sbuf_time[0] = solution_container[MESH_0][FLOW_SOL]->GetMin_Delta_Time();
	sbuf_time[1] = solution_container[MESH_0][FLOW_SOL]->GetMax_Delta_Time();

	/*--- Send/Receive information ---*/
	switch (config->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES:
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		MPI::COMM_WORLD.Reduce(sbuf_residual_levelset, rbuf_residual_levelset, nVar_LevelSet, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case RANS :
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		MPI::COMM_WORLD.Reduce(sbuf_residual_turbulent, rbuf_residual_turbulent, nVar_Turb, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	case NS_PLASMA:
		MPI::COMM_WORLD.Reduce(sbuf_residual_plasma, rbuf_residual_plasma, nVar_Plasma, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	case ADJ_EULER : case ADJ_NAVIER_STOKES :
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjoint, rbuf_residual_adjoint, nVar_Adj, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjoint, rbuf_residual_adjoint, nVar_Adj, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	case LIN_EULER : case LIN_NAVIER_STOKES :
		MPI::COMM_WORLD.Reduce(sbuf_residual_linearized, rbuf_residual_linearized, nVar_Lin, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	}
	MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, nVar_Force, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Reduce(sbuf_time, rbuf_time, nVar_Time, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
	MPI::COMM_WORLD.Reduce(sbuf_conv, rbuf_conv, nVar_Conv, MPI::UNSIGNED_SHORT, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Barrier();

	/*-- Compute global convergence criteria --*/
	if (rank == MASTER_NODE) {

		/*--- Update the value of the residual ---*/
		for (iVar = 0; iVar < nVar_Flow; iVar++)
			rbuf_residual_flow[iVar] = log10 (sqrt(rbuf_residual_flow[iVar]));
		for (iVar = 0; iVar < nVar_LevelSet; iVar++)
			rbuf_residual_levelset[iVar] = log10 (sqrt(rbuf_residual_levelset[iVar]));

		if (rbuf_conv[0] == size) buf_convergence = 1;
		else buf_convergence = 0;
	}
	MPI::COMM_WORLD.Bcast(&buf_convergence, 1, MPI::UNSIGNED_SHORT, MASTER_NODE);

	switch (config->GetKind_Solver()) {
	case EULER: case NAVIER_STOKES: case RANS:
		if (buf_convergence == 1) integration[FLOW_SOL]->SetConvergence(true);
		else integration[FLOW_SOL]->SetConvergence(false); break;
	case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
		if (buf_convergence == 1) integration[FLOW_SOL]->SetConvergence(true);
		else integration[FLOW_SOL]->SetConvergence(false);
		break;
	case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
		if (buf_convergence == 1) integration[ADJFLOW_SOL]->SetConvergence(true);
		else integration[ADJFLOW_SOL]->SetConvergence(false); break;
	case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:
		if (buf_convergence == 1) integration[ADJFLOW_SOL]->SetConvergence(true);
		else integration[ADJFLOW_SOL]->SetConvergence(false); break;
	case LIN_EULER: case LIN_NAVIER_STOKES:
		if (buf_convergence == 1) integration[LINFLOW_SOL]->SetConvergence(true);
		else integration[LINFLOW_SOL]->SetConvergence(false); break;
	case ELECTRIC_POTENTIAL:
		if (buf_convergence == 1) integration[ELEC_SOL]->SetConvergence(true);
		else integration[ELEC_SOL]->SetConvergence(false); break;
	case WAVE:
		if (buf_convergence == 1) integration[WAVE_SOL]->SetConvergence(true);
		else integration[WAVE_SOL]->SetConvergence(false); break;
	}

	/*--- Write result using the master node ---*/
	if (rank == MASTER_NODE) {

		/*--- Write the screen header---*/
		if ((iExtIter == 1) || (write_heads)) {
			cout << endl << " Min DT: " << rbuf_time[0] << 
					". Max DT: " << rbuf_time[1] <<
					". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
			switch (config->GetKind_Solver()) {
			case EULER : case NAVIER_STOKES:
				if (incompressible) cout << endl << " DT Iter" << "      Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				else if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
				else cout << endl << " DT Iter" << "      Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				break;
			case RANS :
				if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]" << " CThrust(Total)" << " CTorque(Total)" << endl;
				else cout << endl << " DT Iter" << "      Res[Rho]" << "       Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				break;
			case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
				cout << endl << " DT Iter" << "    Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "      CLevelSet" << endl;
				break;
			case ADJ_EULER : case ADJ_NAVIER_STOKES :
				if (incompressible) cout << endl << " DT Iter" << "  Res[Psi_Press]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
				else cout << endl << " DT Iter" << "  Res[Psi_Rho]" << "    Res[Psi_E]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
				break;
			case ADJ_FREE_SURF_EULER : case ADJ_FREE_SURF_NAVIER_STOKES :
				if (incompressible) cout << endl << " DT Iter" << "  Res[Psi_Press]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
				else cout << endl << " DT Iter" << "  Res[Psi_Rho]" << "    Res[Psi_E]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
				break;
			case NS_PLASMA :
				cout << endl << " DT Iter" << "      Res[Rho_0]" << "       Res[E_0]" << endl;
				break;
			}
		}
		/*--- Write the solution on the screen ---*/
		if ((iExtIter != 0) && (write_solution)) {
			switch (config->GetKind_Solver()) {
			case EULER : case NAVIER_STOKES: case RANS:
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(8); cout << iExtIter;
				cout.width(14); cout << rbuf_residual_flow[0];
				if (config->GetKind_Solver() == RANS) {
					cout.width(14); cout << rbuf_residual_turbulent[0];
				}
				if (!incompressible) {
					cout.width(13); cout << rbuf_residual_flow[3];
				}
				if (rotating_frame && geometry[MESH_0]->GetnDim() == 3 ) {
					cout.width(15); cout << rbuf_force[13];
					cout.width(15); cout << rbuf_force[14];
				} else {
					cout.width(15); cout << rbuf_force[0];
					cout.width(15); cout << rbuf_force[1];
				}
				cout << endl;
				break;
			case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(8); cout << iExtIter;
				cout.width(14); cout << rbuf_residual_flow[0];
				cout.width(14); cout << rbuf_residual_levelset[0];
				cout.width(15); cout << rbuf_force[0];
				cout.width(15); cout << rbuf_force[8];
				cout << endl;
				break;
			case ADJ_EULER : case ADJ_NAVIER_STOKES :
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(14); cout << rbuf_residual_adjoint[0];
				if (!incompressible) { cout.width(14); cout << rbuf_residual_adjoint[3]; }
				cout.precision(4);
				cout.setf(ios::scientific,ios::floatfield);
				cout.width(14); cout << rbuf_force[0];
				cout.width(14); cout << rbuf_force[1];
				cout << endl;
			case ADJ_FREE_SURF_EULER : case ADJ_FREE_SURF_NAVIER_STOKES :
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(14); cout << rbuf_residual_adjoint[0];
				if (!incompressible) { cout.width(14); cout << rbuf_residual_adjoint[3]; }
				cout.precision(4);
				cout.setf(ios::scientific,ios::floatfield);
				cout.width(14); cout << rbuf_force[0];
				cout.width(14); cout << rbuf_force[1];
				cout << endl;
			case NS_PLASMA :
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(8); cout << iExtIter;
				cout.width(14); cout << rbuf_residual_flow[0];
				cout.width(13); cout << rbuf_residual_flow[3];
			}
		}
		cout.unsetf(ios::fixed);
	}


	delete [] sbuf_residual_flow;
	delete [] sbuf_residual_adjoint;
	delete [] sbuf_residual_linearized;
	delete [] sbuf_residual_turbulent;
	delete [] sbuf_residual_plasma;
	delete [] sbuf_residual_levelset;
	delete [] sbuf_force;
	delete [] sbuf_time;

	if (rank == MASTER_NODE) {
		delete [] rbuf_residual_flow;
		delete [] rbuf_residual_adjoint;
		delete [] rbuf_residual_linearized;
		delete [] rbuf_residual_turbulent;
		delete [] rbuf_residual_plasma;
		delete [] rbuf_residual_levelset;
		delete [] rbuf_force;
		delete [] rbuf_time;
	}

#endif
}

void COutput::SetHistory_MainIter(ofstream *ConvHist_file, CGeometry ***geometry, CSolution ****solution_container, CConfig *config,
		CIntegration **integration, unsigned long iExtIter, unsigned long timeused, unsigned short val_nDomain) {

#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
#else
	int rank = MASTER_NODE;
#endif

	/*--- WARNING: These buffers have hard-coded lengths. Note that you
	 may have to adjust them to be larger if adding more entries. ---*/
	char begin[200], 
	direct_coeff[300], adjoint_coeff[200],
	flow_resid[200], adj_flow_resid[200],
	turb_resid[200], adj_turb_resid[200],
	plasma_resid[200], adj_plasma_resid[200], 
	levelset_resid[200], adj_levelset_resid[200],
	wave_coeff[300], wave_resid[200], end[200];
	unsigned short iDomain = 0;
	double dummy = 0.0;

	double timeiter = double(timeused)/double(iExtIter+1);
	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || 
			(config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
	bool write_heads = (((iExtIter % (config->GetWrt_Con_Freq()*20)) == 0) || dual_time);

	unsigned short FinestMesh = config->GetFinestMesh();
	unsigned short nDim = geometry[iDomain][FinestMesh]->GetnDim();

	bool incompressible = config->GetIncompressible();
	bool compressible = !config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool free_surface = config->GetFreeSurface();
	bool equiv_area = config->GetEquivArea();
	bool turbulent = (config->GetKind_Solver() == RANS);
	bool adjoint = ((config->GetKind_Solver() == ADJ_EULER) || (config->GetKind_Solver() == ADJ_NAVIER_STOKES) || 
			(config->GetKind_Solver() == ADJ_RANS) || (config->GetKind_Solver() == ADJ_FREE_SURF_EULER) ||
			(config->GetKind_Solver() == ADJ_FREE_SURF_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_FREE_SURF_RANS) ||
			(config->GetKind_Solver() == ADJ_NS_PLASMA));
	bool wave = (config->GetKind_Solver() == WAVE);

	/*--- Initialize variables to store information from all domains ---*/
	double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0,
			Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0,
			Total_CEff = 0.0, Total_CEquivArea = 0.0, Total_CNearFieldPress = 0.0,
			Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
			Total_CT = 0.0, Total_CQ = 0.0, Total_CFreeSurface = 0.0, Total_CWave = 0.0;

	double Total_CSens_Geo = 0.0, Total_CSens_Mach = 0.0, Total_CSens_AoA = 0.0;

	double *sbuf_residual_flow = NULL, *sbuf_residual_turbulent = NULL, *sbuf_residual_levelset = NULL, *rbuf_residual_plasma = NULL;
	double *sbuf_residual_adjflow = NULL, *sbuf_residual_adjturbulent = NULL, *sbuf_residual_adjlevelset = NULL, *rbuf_residual_adjplasma = NULL;
	double *sbuf_residual_wave = NULL;

	double *rbuf_residual_flow = NULL, *rbuf_residual_turbulent = NULL, *rbuf_residual_levelset = NULL, *sbuf_residual_plasma = NULL;
	double *rbuf_residual_adjflow = NULL, *rbuf_residual_adjturbulent = NULL, *rbuf_residual_adjlevelset = NULL, *sbuf_residual_adjplasma = NULL;
	double *rbuf_residual_wave = NULL;

	double *sbuf_force = NULL,  *rbuf_force = NULL;

	unsigned long *sbuf_time = NULL, *rbuf_time = NULL;

	unsigned short iVar, *sbuf_conv = NULL, *rbuf_conv = NULL; 

	/*--- Direct problem variables ---*/
	unsigned short nVar_Flow;
	if (compressible) nVar_Flow = nDim+2;
	else nVar_Flow = nDim+1;	
	unsigned short nVar_Plasma = config->GetnSpecies()*(nDim+2);
	unsigned short nVar_LevelSet = 1;
	unsigned short nVar_Turb = 1;
	unsigned short nVar_Wave = 1;

	/*--- Adjoint problem variables ---*/
	unsigned short nVar_AdjFlow;
	if (compressible) nVar_AdjFlow = nDim+2;
	else nVar_AdjFlow = nDim+1;
	unsigned short nVar_AdjPlasma = config->GetnSpecies()*(nDim+2);
	unsigned short nVar_AdjLevelSet = 1;
	unsigned short nVar_AdjTurb = 1;

	/*--- Other vectors ---*/
	unsigned short nVar_Force = 16;
	unsigned short nVar_Time = 2;
	unsigned short nVar_Conv = 1;

	/*--- Allocate memory for send buffer ---*/
	sbuf_residual_flow = new double[nVar_Flow];
	sbuf_residual_turbulent = new double[nVar_Turb];
	sbuf_residual_plasma = new double[nVar_Plasma];
	sbuf_residual_levelset = new double[nVar_LevelSet];
	sbuf_residual_wave = new double[nVar_Wave];

	sbuf_residual_adjflow = new double[nVar_AdjFlow];
	sbuf_residual_adjturbulent = new double[nVar_AdjTurb];
	sbuf_residual_adjplasma = new double[nVar_AdjPlasma];
	sbuf_residual_adjlevelset = new double[nVar_AdjLevelSet];

	sbuf_force = new double[nVar_Force];
	sbuf_time = new unsigned long[nVar_Time];
	sbuf_conv = new unsigned short[nVar_Conv];

	for (iVar = 0; iVar < nVar_Flow; iVar++) { sbuf_residual_flow[iVar] = -20; }
	for (iVar = 0; iVar < nVar_Turb; iVar++) { sbuf_residual_turbulent[iVar] = -20; }
	for (iVar = 0; iVar < nVar_Plasma; iVar++) { sbuf_residual_plasma[iVar] = -20; }
	for (iVar = 0; iVar < nVar_LevelSet; iVar++) { sbuf_residual_levelset[iVar] = -20; }
	for (iVar = 0; iVar < nVar_Wave; iVar++) { sbuf_residual_wave[iVar] = -20; }

	for (iVar = 0; iVar < nVar_AdjFlow; iVar++) { sbuf_residual_adjflow[iVar] = -20; }
	for (iVar = 0; iVar < nVar_AdjTurb; iVar++) { sbuf_residual_adjturbulent[iVar] = -20; }
	for (iVar = 0; iVar < nVar_AdjPlasma; iVar++) { sbuf_residual_adjplasma[iVar] = -20; }
	for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++) { sbuf_residual_adjlevelset[iVar] = -20; }

	for (iVar = 0; iVar < nVar_Force; iVar++) { sbuf_force[iVar] = 0.0; }
	for (iVar = 0; iVar < nVar_Time; iVar++) { sbuf_time[iVar] = 0; }
	for (iVar = 0; iVar < nVar_Conv; iVar++) { sbuf_conv[iVar] = 0; }	

	/*--- Allocate memory for receive buffer ---*/
	if (rank == MASTER_NODE) {

		rbuf_residual_flow = new double[nVar_Flow];
		rbuf_residual_turbulent = new double[nVar_Turb];
		rbuf_residual_plasma = new double[nVar_Plasma];
		rbuf_residual_levelset = new double[nVar_LevelSet];
		rbuf_residual_wave = new double[nVar_Wave];

		rbuf_residual_adjflow = new double[nVar_AdjFlow];
		rbuf_residual_adjturbulent = new double[nVar_AdjTurb];
		rbuf_residual_adjplasma = new double[nVar_AdjPlasma];
		rbuf_residual_adjlevelset = new double[nVar_AdjLevelSet];

		rbuf_force = new double[nVar_Force];
		rbuf_time = new unsigned long[nVar_Time];
		rbuf_conv = new unsigned short[nVar_Conv];

		for (iVar = 0; iVar < nVar_Flow; iVar++) { rbuf_residual_flow[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_Turb; iVar++) { rbuf_residual_turbulent[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_Plasma; iVar++) { rbuf_residual_plasma[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_LevelSet; iVar++) { rbuf_residual_levelset[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_Wave; iVar++) { rbuf_residual_wave[iVar] = 0.0; }

		for (iVar = 0; iVar < nVar_AdjFlow; iVar++) { rbuf_residual_adjflow[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_AdjTurb; iVar++) { rbuf_residual_adjturbulent[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_AdjPlasma; iVar++) { rbuf_residual_adjplasma[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++) { rbuf_residual_adjlevelset[iVar] = 0.0; }

		for (iVar = 0; iVar < nVar_Force; iVar++) { rbuf_force[iVar] = 0.0; }
		for (iVar = 0; iVar < nVar_Time; iVar++) { rbuf_time[iVar] = 0; }
		for (iVar = 0; iVar < nVar_Conv; iVar++) { rbuf_conv[iVar] = 0; }
	}

	/*--- Write information from nodes ---*/
	switch (config->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES: case RANS:
	case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES: case FREE_SURF_RANS:
	case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
	case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES: case ADJ_FREE_SURF_RANS:

		/*--- Flow solution coefficients (serial) ---*/
		Total_CLift += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CLift();
		Total_CDrag += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
		Total_CSideForce += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
		Total_CMx += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMx();
		Total_CMy += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMy();
		Total_CMz += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMz();
		Total_CFx += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CFx();
		Total_CFy += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CFy();
		Total_CFz += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CFz();
		Total_CEff += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CEff();

		if (equiv_area) {
			Total_CEquivArea += config->GetWeightCd()*solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CDrag()
						+ (1.0-config->GetWeightCd())*solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
			Total_CNearFieldPress += config->GetWeightCd()*solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CDrag()
						+ (1.0-config->GetWeightCd())*solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldPress();
		}

		if (rotating_frame) {
			Total_CMerit += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
			Total_CT += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CT();
			Total_CQ += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CQ();
		}

		/*--- Flow solution coefficients (parallel) ---*/
		sbuf_force[0]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CLift();
		sbuf_force[1]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
		sbuf_force[2]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
		sbuf_force[3]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMx();
		sbuf_force[4]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMy();
		sbuf_force[5]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMz();
		sbuf_force[6] += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CFx();
		sbuf_force[7] += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CFy();
		sbuf_force[8] += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CFz();

		if (free_surface) {
			sbuf_force[9]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CFreeSurface();
		}

		if (equiv_area) {
			sbuf_force[9]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea(),
					sbuf_force[10]  += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldPress();
		}

		if (rotating_frame) {
			sbuf_force[9] += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
			sbuf_force[10] += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CT();
			sbuf_force[11] += solution_container[iDomain][FinestMesh][FLOW_SOL]->GetTotal_CQ();
		}

		/*--- Convergence criteria ---*/
		sbuf_conv[0] = integration[FLOW_SOL]->GetConvergence();

		/*--- Flow Residuals ---*/
		for (iVar = 0; iVar < nVar_Flow; iVar++)
			sbuf_residual_flow[iVar] = solution_container[iDomain][FinestMesh][FLOW_SOL]->GetRes_Max(iVar);

		/*--- Turbulent residual ---*/
		if (turbulent) {
			for (iVar = 0; iVar < nVar_Turb; iVar++)
				sbuf_residual_turbulent[iVar] = solution_container[iDomain][FinestMesh][TURB_SOL]->GetRes_Max(iVar);
		}

		/*--- Free Surface residual ---*/
		if (free_surface) {
			for (iVar = 0; iVar < nVar_LevelSet; iVar++)
				sbuf_residual_levelset[iVar] = solution_container[iDomain][FinestMesh][LEVELSET_SOL]->GetRes_Max(iVar);
		}

		if (adjoint) {

			/*--- Adjoint solution coefficients (serial) ---*/
			Total_CSens_Geo  = solution_container[iDomain][FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Geo();
			Total_CSens_Mach = solution_container[iDomain][FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Mach();
			Total_CSens_AoA  = solution_container[iDomain][FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_AoA();

			/*--- Adjoint solution coefficients (parallel) ---*/
			sbuf_force[0] = solution_container[iDomain][FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Geo();
			sbuf_force[1] = solution_container[iDomain][FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_Mach();
			sbuf_force[2] = solution_container[iDomain][FinestMesh][ADJFLOW_SOL]->GetTotal_CSens_AoA();

			/*--- Convergence criteria ---*/
			sbuf_conv[0] = integration[ADJFLOW_SOL]->GetConvergence();

			/*--- Adjoint flow residuals ---*/
			for (iVar = 0; iVar < nVar_AdjFlow; iVar++)
				sbuf_residual_adjflow[iVar] = solution_container[iDomain][FinestMesh][ADJFLOW_SOL]->GetRes_Max(iVar);

			/*--- Adjoint turbulent residuals ---*/
			if (turbulent) {
				for (iVar = 0; iVar < nVar_AdjTurb; iVar++)
					sbuf_residual_adjturbulent[iVar] = solution_container[iDomain][FinestMesh][ADJTURB_SOL]->GetRes_Max(iVar);
			}

			/*--- Adjoint level set residuals ---*/
			if (free_surface) {
				for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++)
					sbuf_residual_adjlevelset[iVar] = solution_container[iDomain][FinestMesh][ADJLEVELSET_SOL]->GetRes_Max(iVar);
			}
		}

		break;

	case  NS_PLASMA:

		/*--- Plasma coefficients (serial) ---*/
		Total_CLift += solution_container[iDomain][FinestMesh][PLASMA_SOL]->GetTotal_CLift();
		Total_CDrag += solution_container[iDomain][FinestMesh][PLASMA_SOL]->GetTotal_CDrag();

		/*--- Plasma coefficients (parallel) ---*/
		sbuf_force[0] += solution_container[iDomain][FinestMesh][PLASMA_SOL]->GetTotal_CLift();
		sbuf_force[1] += solution_container[iDomain][FinestMesh][PLASMA_SOL]->GetTotal_CDrag();

		/*--- Convergence criteria ---*/
		sbuf_conv[0] = integration[PLASMA_SOL]->GetConvergence();

		/*--- Plasma Residuals ---*/
		for (iVar = 0; iVar < nVar_Plasma; iVar++)
			sbuf_residual_plasma[iVar] = solution_container[iDomain][FinestMesh][PLASMA_SOL]->GetRes_Max(iVar);

		break;

	case  WAVE:

		/*--- Plasma coefficients (serial) ---*/
		Total_CWave += 0.0; //solution_container[iDomain][FinestMesh][WAVE_SOL]->GetTotal_CWave();

		/*--- Plasma coefficients (parallel) ---*/
		sbuf_force[0] += 0.0; //solution_container[iDomain][FinestMesh][WAVE_SOL]->GetTotal_CWave();

		/*--- Convergence criteria ---*/
		sbuf_conv[0] = integration[WAVE_SOL]->GetConvergence();

		/*--- Plasma Residuals ---*/
		for (iVar = 0; iVar < nVar_Wave; iVar++) {
			sbuf_residual_wave[iVar] = solution_container[iDomain][FinestMesh][WAVE_SOL]->GetRes_Max(iVar);
		}

		break;

	}
	sbuf_time[0] = (unsigned long) timeiter; 
	sbuf_time[1] = (unsigned long) timeused;

#ifndef NO_MPI

	unsigned short buf_convergence = 0;

	/*--- Send/Receive information ---*/
	switch (config->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES: case RANS : case FREE_SURF_EULER : case FREE_SURF_NAVIER_STOKES: case FREE_SURF_RANS :
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		if (turbulent)
			MPI::COMM_WORLD.Reduce(sbuf_residual_turbulent, rbuf_residual_turbulent, nVar_Turb, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		if (free_surface)
			MPI::COMM_WORLD.Reduce(sbuf_residual_levelset, rbuf_residual_levelset, nVar_LevelSet, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case NS_PLASMA:
		MPI::COMM_WORLD.Reduce(sbuf_residual_plasma, rbuf_residual_plasma, nVar_Plasma, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case WAVE:
		MPI::COMM_WORLD.Reduce(sbuf_residual_wave, rbuf_residual_wave, nVar_Wave, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case ADJ_EULER : case ADJ_NAVIER_STOKES :
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjflow, rbuf_residual_adjflow, nVar_AdjFlow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case ADJ_FREE_SURF_EULER : case ADJ_FREE_SURF_NAVIER_STOKES :
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjflow, rbuf_residual_adjflow, nVar_AdjFlow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjlevelset, rbuf_residual_adjlevelset, nVar_AdjLevelSet, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	}

	MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, nVar_Force, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Reduce(sbuf_time, rbuf_time, nVar_Time, MPI::UNSIGNED_LONG, MPI::MAX, MASTER_NODE);
	MPI::COMM_WORLD.Reduce(sbuf_conv, rbuf_conv, nVar_Conv, MPI::UNSIGNED_SHORT, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Barrier();

	/*-- Compute global convergence criteria --*/
	if (rank == MASTER_NODE) {

		/*--- Update the value of the residual ---*/
		for (iVar = 0; iVar < nVar_Flow; iVar++) rbuf_residual_flow[iVar] = sqrt(rbuf_residual_flow[iVar]);
		for (iVar = 0; iVar < nVar_Turb; iVar++) rbuf_residual_turbulent[iVar] = sqrt(rbuf_residual_turbulent[iVar]);
		for (iVar = 0; iVar < nVar_LevelSet; iVar++) rbuf_residual_levelset[iVar] = sqrt(rbuf_residual_levelset[iVar]);
		for (iVar = 0; iVar < nVar_Plasma; iVar++) rbuf_residual_plasma[iVar] = sqrt(rbuf_residual_plasma[iVar]);
		for (iVar = 0; iVar < nVar_Wave; iVar++) rbuf_residual_wave[iVar] = sqrt(rbuf_residual_wave[iVar]);
		for (iVar = 0; iVar < nVar_AdjFlow; iVar++) rbuf_residual_adjflow[iVar] = sqrt(rbuf_residual_adjflow[iVar]);
		for (iVar = 0; iVar < nVar_AdjTurb; iVar++) rbuf_residual_adjturbulent[iVar] = sqrt(rbuf_residual_adjturbulent[iVar]);
		for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++) rbuf_residual_adjlevelset[iVar] = sqrt(rbuf_residual_adjlevelset[iVar]);
		for (iVar = 0; iVar < nVar_AdjPlasma; iVar++) rbuf_residual_adjplasma[iVar] = sqrt(rbuf_residual_adjplasma[iVar]);
		if (rbuf_conv[0] == size) buf_convergence = 1;
		else buf_convergence = 0;
	}
	MPI::COMM_WORLD.Bcast(&buf_convergence, 1, MPI::UNSIGNED_SHORT, MASTER_NODE);

	switch (config->GetKind_Solver()) {
	case EULER: case NAVIER_STOKES: case RANS:
		if (buf_convergence == 1) integration[FLOW_SOL]->SetConvergence(true);
		else integration[FLOW_SOL]->SetConvergence(false); break;
	case WAVE:
		if (buf_convergence == 1) integration[WAVE_SOL]->SetConvergence(true);
		else integration[WAVE_SOL]->SetConvergence(false); break;
	case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
		if (buf_convergence == 1) integration[FLOW_SOL]->SetConvergence(true);
		else integration[FLOW_SOL]->SetConvergence(false); break;
	case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
		if (buf_convergence == 1) integration[ADJFLOW_SOL]->SetConvergence(true);
		else integration[ADJFLOW_SOL]->SetConvergence(false); break;
	case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:
		if (buf_convergence == 1) integration[ADJFLOW_SOL]->SetConvergence(true);
		else integration[ADJFLOW_SOL]->SetConvergence(false); break;
	}


	/*--- Set the value of the funcionals, and residuals ---*/
	if (rank == MASTER_NODE) {

		Total_CLift = rbuf_force[0];
		Total_CDrag = rbuf_force[1];
		Total_CSideForce = rbuf_force[2];
		Total_CMx = rbuf_force[3];
		Total_CMy = rbuf_force[4];
		Total_CMz = rbuf_force[5];
		Total_CFx = rbuf_force[6];
		Total_CFy = rbuf_force[7];
		Total_CFz = rbuf_force[8];
		Total_CEff = rbuf_force[0]/(rbuf_force[1]+EPS);

		if (free_surface) {
			Total_CFreeSurface = rbuf_force[9];
		}
		if (equiv_area) {
			Total_CEquivArea  = rbuf_force[9]/double(size);
			Total_CNearFieldPress  = rbuf_force[10];
		}
		if (rotating_frame) {
			Total_CMerit = rbuf_force[9];
			Total_CT = rbuf_force[10];
			Total_CQ = rbuf_force[11];
		}
		if (wave) {
			Total_CWave = rbuf_force[9];
		}

		/*--- Adjoint solution coefficients ---*/
		if (adjoint) {
			Total_CSens_Geo  = rbuf_force[0];
			Total_CSens_Mach = rbuf_force[1];
			Total_CSens_AoA  = rbuf_force[2];
		}

		timeused = rbuf_time[1];

	}

#else

	/*--- Update the residual without MPI stuff ---*/
	for (iVar = 0; iVar < nVar_Flow; iVar++) rbuf_residual_flow[iVar] = sbuf_residual_flow[iVar];
	if (turbulent) for (iVar = 0; iVar < nVar_Turb; iVar++) rbuf_residual_turbulent[iVar] = sbuf_residual_turbulent[iVar];
	if (free_surface) for (iVar = 0; iVar < nVar_LevelSet; iVar++) rbuf_residual_levelset[iVar] = sbuf_residual_levelset[iVar];
	if (wave) for (iVar = 0; iVar < nVar_Wave; iVar++) rbuf_residual_wave[iVar] = sbuf_residual_wave[iVar];
	for (iVar = 0; iVar < nVar_Plasma; iVar++) rbuf_residual_plasma[iVar] = sbuf_residual_plasma[iVar];

	if (adjoint) {
		for (iVar = 0; iVar < nVar_AdjFlow; iVar++) rbuf_residual_adjflow[iVar] = sbuf_residual_adjflow[iVar];
		if (turbulent) for (iVar = 0; iVar < nVar_AdjTurb; iVar++) rbuf_residual_adjturbulent[iVar] = sbuf_residual_adjturbulent[iVar];
		if (free_surface) for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++) rbuf_residual_adjlevelset[iVar] = sbuf_residual_adjlevelset[iVar];
		for (iVar = 0; iVar < nVar_AdjPlasma; iVar++) rbuf_residual_plasma[iVar] = sbuf_residual_plasma[iVar];
	}

#endif	

	if (rank == MASTER_NODE) {
		/*--- Write the begining of the history file ---*/
		sprintf (begin, "%12d", int(iExtIter));

		/*--- Write the end of the history file ---*/
		sprintf (end, ", %12.10f\n", double(timeused)/(CLOCKS_PER_SEC*60.0));

		/*--- Write the solution and residual of the history file ---*/
		switch (config->GetKind_Solver()) {

		case EULER : case NAVIER_STOKES: case RANS: case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
		case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:

			sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
					Total_CFz, Total_CEff);
			if (equiv_area)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
						Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldPress);
			if (rotating_frame)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
						Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CMerit, Total_CT, Total_CQ);
			if (free_surface)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
						Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFreeSurface);

			if (nDim == 2) {
				if (incompressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), dummy, dummy );
				else sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), log10 (rbuf_residual_flow[3]), dummy );
			}
			else {
				if (incompressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), log10 (rbuf_residual_flow[3]), dummy );
				else sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), log10 (rbuf_residual_flow[3]), log10 (rbuf_residual_flow[4]) );
			}

			if (turbulent) sprintf (turb_resid, ", %12.10f", log10 (rbuf_residual_turbulent[0]));

			if (free_surface) sprintf (levelset_resid, ", %12.10f", log10 (rbuf_residual_levelset[0]));

			if (adjoint) {

				sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, 0.0", Total_CSens_Geo, Total_CSens_Mach, Total_CSens_AoA);

				if (nDim == 2) sprintf (adj_flow_resid, ", %17.15f, %17.15f, %17.15f, %17.15f, 0.0", log10 (rbuf_residual_adjflow[0]),log10 (rbuf_residual_adjflow[1]),log10 (rbuf_residual_adjflow[2]),log10 (rbuf_residual_adjflow[3]) );
				else sprintf (adj_flow_resid, ", 12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_adjflow[0]),log10 (rbuf_residual_adjflow[1]),log10 (rbuf_residual_adjflow[2]),log10 (rbuf_residual_adjflow[3]), log10 (rbuf_residual_adjflow[4]) );

				if (turbulent) sprintf (adj_turb_resid, ", %12.10f", log10 (rbuf_residual_adjturbulent[0]));

				if (free_surface) sprintf (adj_levelset_resid, ", %12.10f", log10 (rbuf_residual_adjlevelset[0]));

			}

			break;

		case NS_PLASMA: case ADJ_NS_PLASMA:

			sprintf (direct_coeff, ", %12.10f, %12.10f", Total_CLift, Total_CDrag);

			if (nDim == 2) sprintf (plasma_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_plasma[0]),log10 (rbuf_residual_plasma[1]),log10 (rbuf_residual_plasma[2]),log10 (rbuf_residual_plasma[3]), dummy );
			else sprintf (plasma_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(rbuf_residual_plasma[0]),log10(rbuf_residual_plasma[1]),log10 (rbuf_residual_plasma[2]),log10 (rbuf_residual_plasma[3]), log10 (rbuf_residual_plasma[4]) );

			if (adjoint) {
				sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, 0.0", Total_CSens_Geo, Total_CSens_Mach, Total_CSens_AoA);

				if (nDim == 2) sprintf (adj_plasma_resid, ", %17.15f, %17.15f, %17.15f, %17.15f, %17.15f", log10 (rbuf_residual_adjplasma[0]),log10 (rbuf_residual_adjplasma[1]),log10 (rbuf_residual_adjplasma[2]),log10 (rbuf_residual_adjplasma[3]), dummy );
				else sprintf (adj_flow_resid, ", 12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_adjplasma[0]),log10 (rbuf_residual_adjplasma[1]),log10 (rbuf_residual_adjplasma[2]),log10 (rbuf_residual_adjplasma[3]), log10 (rbuf_residual_adjplasma[4]) );

			}

			break;

		case WAVE:

			sprintf (direct_coeff, ", %12.10f", Total_CWave);

			sprintf (wave_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_wave[0]), dummy, dummy, dummy, dummy );

			break;

		}

		switch (config->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES: case FREE_SURF_EULER : case FREE_SURF_NAVIER_STOKES:
			if (write_heads) {
				if (incompressible) cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				else if (rotating_frame && nDim == 3 ) cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
				else cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			}
			break;
		case RANS :
			if (write_heads) {
				if (rotating_frame && nDim == 3 ) cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]" << " CThrust(Total)" << " CTorque(Total)" << endl;
				else cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			}
			break;
		case NS_PLASMA :
			if (write_heads) {
				if (config->GetKind_GasModel() == ARGON)
					cout << endl << " Iter" << "   Time(s)" << "      Res[r1]" << "       Res[r2]" << "   Res(r3)"<<  endl;
				if (config->GetKind_GasModel() == AIR21)
					cout << endl << " Iter" << "   Time(s)" << "      Res[r1]" << "       Res[r2]" << "   Res(r3)"<<  endl;
				if (config->GetKind_GasModel() == AIR7)
					cout << endl << " Iter" << "   Time(s)" << "        Res[Rho1]" << "         Res[E]" << "	   Heat_Transfer" << endl;
			}
			break;
		case WAVE :
			if (write_heads) {
				cout << endl << " Iter" << "   Time(s)" << "      Res[Wave]" << "   CWave(Total)"<<  endl;
			}
			break;
		case ADJ_NS_PLASMA :
			if (write_heads) {
				if (config->GetKind_GasModel() == AIR7)
					cout << endl << " Iter" << "   Time(s)" << "        Res[Psi_Rho1]" << "         Res[Psi_E]" << "	   CSens_Geo" << endl;
			}
			break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_FREE_SURF_EULER : case ADJ_FREE_SURF_NAVIER_STOKES :
			if (write_heads) {
				if (incompressible) cout << endl << " Iter" << "   Time(s)" << "   Res[Psi_Rho]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
				else cout << endl << " Iter" << "   Time(s)" << "   Res[Psi_Rho]" << "     Res[Psi_E]" << "     CSens_Geo" << "    CSens_Mach"<< endl;
			}
			break;
		case ADJ_RANS :
			if (write_heads) cout << endl << "  Iter" << "    Time(s)" << "     Res[Psi_Rho]" << "      Res[Psi_nu]" << "       CSens_Geo" << endl;
			break;
		}

		/*--- Write the solution on the screen and history file ---*/
		if (iExtIter != 0) {
			switch (config->GetKind_Solver()) {

			case EULER : case NAVIER_STOKES: case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:
				ConvHist_file[0] << begin << direct_coeff << flow_resid;
				if (free_surface) ConvHist_file[0] << levelset_resid;
				ConvHist_file[0] << end;

				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				cout.width(14); cout << log10(rbuf_residual_flow[0]);
				if (!incompressible) {
					if (nDim == 2 ) { cout.width(14); cout << log10(rbuf_residual_flow[3]); }
					else { cout.width(14); cout << log10(rbuf_residual_flow[4]); }
				}
				if (rotating_frame && nDim == 3 ) { cout.width(15); cout << Total_CT; cout.width(15); cout << Total_CQ; }
				else { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; }
				cout << endl;
				break;

			case RANS :
				ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid << end;

				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				cout.width(14); cout << log10(rbuf_residual_flow[0]);
				cout.width(14); cout << log10(rbuf_residual_turbulent[0]);
				if (rotating_frame  && nDim == 3 ) { cout.width(15); cout << Total_CT; cout.width(15); cout << Total_CQ; }
				else { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; }
				cout << endl;
				break;

			case NS_PLASMA:
				ConvHist_file[0] << begin << direct_coeff << plasma_resid << end;

				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				if (config->GetKind_GasModel() == ARGON || config->GetKind_GasModel() == AIR21) {
					cout.width(14); cout << log10(rbuf_residual_plasma[0]);
					cout.width(14); cout << log10(rbuf_residual_plasma[nDim+2]);
					cout.width(14); cout << log10(rbuf_residual_plasma[2*(nDim+2)]);
				}
				if (config->GetKind_GasModel() == AIR7) {
					cout.width(14); cout << log10(rbuf_residual_plasma[0]);
					cout.width(14); cout << log10(rbuf_residual_plasma[3]);
				}
				cout << endl;
				break;

			case WAVE:
				ConvHist_file[0] << begin << wave_coeff << wave_resid << end;

				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				cout.width(14); cout << log10(rbuf_residual_wave[0]);
				cout.width(14); cout << Total_CWave;
				cout << endl;
				break;

			case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_FREE_SURF_EULER : case ADJ_FREE_SURF_NAVIER_STOKES :
				ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				cout.width(15); cout << log10(rbuf_residual_adjflow[0]);
				if (!incompressible) {
					if (geometry[DOMAIN_0][FinestMesh]->GetnDim() == 2 ) { cout.width(15); cout << log10(rbuf_residual_adjflow[3]); }
					else { cout.width(15); cout << log10(rbuf_residual_adjflow[4]); }
				}
				cout.precision(4);
				cout.setf(ios::scientific,ios::floatfield);
				cout.width(14); cout << Total_CSens_Geo;
				cout.width(14); cout << Total_CSens_Mach;
				cout << endl;
				break;

			case ADJ_RANS :
				ConvHist_file[0] << begin << adjoint_coeff << direct_coeff << adj_flow_resid << adj_turb_resid << end;
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(6); cout << iExtIter;
				cout.width(11); cout << double(timeiter)/CLOCKS_PER_SEC;
				cout.width(17); cout << log10(rbuf_residual_adjflow[0]);
				cout.width(17); cout << log10(rbuf_residual_adjturbulent[0]);
				cout.width(16); cout << Total_CSens_Geo;
				cout << endl;
				break;

			case ADJ_NS_PLASMA:
				ConvHist_file[0] << begin << adjoint_coeff << direct_coeff << adj_plasma_resid << end;

				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				if (config->GetKind_GasModel() == AIR7) {
					cout.width(14); cout << log10(rbuf_residual_adjplasma[0]);
					cout.width(14); cout << log10(rbuf_residual_adjplasma[3]);
				}
				cout << endl;
				break;

			}
			cout.unsetf(ios::fixed);
		}
	}



	delete [] sbuf_residual_flow;
	delete [] sbuf_residual_turbulent;
	delete [] sbuf_residual_levelset;
	delete [] sbuf_residual_plasma;
	delete [] sbuf_residual_wave;

	delete [] sbuf_residual_adjflow;
	delete [] sbuf_residual_adjturbulent;
	delete [] sbuf_residual_adjlevelset;
	delete [] sbuf_residual_adjplasma;

	delete [] sbuf_force;
	delete [] sbuf_time;
	delete [] sbuf_conv;

	if (rank == MASTER_NODE) {
		delete [] rbuf_residual_flow;
		delete [] rbuf_residual_levelset;
		delete [] rbuf_residual_plasma;
		delete [] rbuf_residual_turbulent;
		delete [] rbuf_residual_wave;

		delete [] rbuf_residual_adjflow;
		delete [] rbuf_residual_adjlevelset;
		delete [] rbuf_residual_adjplasma;
		delete [] rbuf_residual_adjturbulent;

		delete [] rbuf_force;
		delete [] rbuf_time;
		delete [] rbuf_conv;

	}
}

void COutput::SetResult_Files(CSolution ****solution_container, CGeometry ***geometry, CConfig *config, unsigned long iExtIter, unsigned short val_nDomain) {

	switch (config->GetKind_Solver()) {

	case EULER : case NAVIER_STOKES : case RANS : case FREE_SURF_EULER : case FREE_SURF_NAVIER_STOKES:

		SetDomain_Flow(config, geometry, solution_container, iExtIter, val_nDomain);
		SetSurface_Flow(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], iExtIter);
		SetReStart(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0], config->GetReStart_FlowFileName());
		SetSurfaceCSV_Flow(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], iExtIter);
		break;

	case NS_PLASMA :
		SetDomain_Flow(config, geometry, solution_container, iExtIter, val_nDomain);
		SetSurface_Flow(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][PLASMA_SOL], iExtIter);
		SetReStart(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0], config->GetReStart_FlowFileName());
		break;

	case WAVE:
		SetDomain_Flow(config, geometry, solution_container, iExtIter, val_nDomain);
		SetSurface_Flow(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][WAVE_SOL], iExtIter);
		SetReStart(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0], config->GetReStart_FlowFileName());
		break;

	case ELECTRIC_POTENTIAL:
		SetDomain_Flow(config, geometry, solution_container, iExtIter, val_nDomain);
		SetReStart(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0], config->GetReStart_FlowFileName());
		break;

	case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS : case ADJ_FREE_SURF_EULER : case ADJ_FREE_SURF_NAVIER_STOKES:
		SetDomain_Adjoint(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0], iExtIter);
		SetSurface_Adjoint(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][ADJFLOW_SOL], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], iExtIter);
		SetSurfaceCSV_Adjoint(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][ADJFLOW_SOL], solution_container[DOMAIN_0][MESH_0][FLOW_SOL], iExtIter);
		SetReStart(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0], config->GetReStart_AdjFileName());
		break;

	case LIN_EULER : case LIN_NAVIER_STOKES :
		SetLinearized_Variables(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0], iExtIter);
		SetSurface_Linearized(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][LINFLOW_SOL], config->GetSurfLinCoeff_FileName(), iExtIter);
		SetSurfaceCSV_Linearized(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0][LINFLOW_SOL], config->GetSurfLinCoeff_FileName(), iExtIter);
		SetReStart(config, geometry[DOMAIN_0][MESH_0], solution_container[DOMAIN_0][MESH_0], config->GetReStart_LinFileName());
		break;
	}
}

#ifdef NEW_EQUIV_AREA
void COutput::SetEquivalentArea(CSolution *solution_container, CGeometry *geometry,
		CConfig *config, unsigned long iExtIter) {

	ofstream EquivArea_file, FuncGrad_file;
	unsigned short iMarker = 0, iDim, iVar;
	double auxCoord, InverseDesign, DeltaX, Coord_i, Coord_j, jp1Coord, 
	*Coord = NULL, MeanFuntion, *Face_Normal = NULL, auxArea, auxPress, 
	Mach, Beta, Begin_Int, End_Int, R_Plane, Pressure_Inf, Density_Inf, 
	RefAreaCoeff, ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, 
	*Ycoord = NULL, *Zcoord = NULL, *AzimutalAngle = NULL,
	*Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, 
	*NearFieldWeight = NULL, *Weight = NULL, jFunction, jp1Function;
	unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint, 
			*IdPoint = NULL, *IdDomain = NULL, auxDomain;

	unsigned short nDim = geometry->GetnDim();
	int rank = MESH_0;

	Mach  = config->GetMach_FreeStreamND();
	Beta = sqrt(Mach*Mach-1.0);
	Begin_Int = config->GetEA_IntLimit(0)*config->GetConversion_Factor();
	End_Int = config->GetEA_IntLimit(1)*config->GetConversion_Factor();
	R_Plane = abs(config->GetEA_IntLimit(2)*config->GetConversion_Factor());
	Pressure_Inf = config->GetPressure_FreeStreamND();
	Density_Inf = config->GetDensity_FreeStreamND();	
	RefAreaCoeff = config->GetRefAreaCoeff();
	Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
	Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
	Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
	ModVelocity_Inf = 0;
	for (iDim = 0; iDim < 3; iDim++)
		ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];
	ModVelocity_Inf = sqrt(ModVelocity_Inf);

	factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Density_Inf*Mach*Mach);

#ifndef NO_MPI

	/*--- Compute the total number of points of the near-field ---*/
	nVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if ((Begin_Int < Coord[0]) && (Coord[0] < End_Int))
					if (Face_Normal[nDim-1] > 0.0) nVertex_NearField ++;
			}

	/*--- Create an array with all the coordinates, points, pressures, face area,
	 equivalent area, and nearfield weight ---*/
	Xcoord = new double[nVertex_NearField];
	Ycoord = new double[nVertex_NearField];
	Zcoord = new double[nVertex_NearField];
	AzimutalAngle = new double[nVertex_NearField];
	IdPoint = new unsigned long[nVertex_NearField];
	IdDomain = new unsigned long[nVertex_NearField];
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
					if (Face_Normal[nDim-1] > 0.0) {
						IdPoint[nVertex_NearField] = iPoint;
						Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
						Ycoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
						Zcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(2);

						/*--- Compute the azymutal angle ---*/
						/*--- first cylinder rotation, and then angle computation ---*/
						AzimutalAngle[nVertex_NearField] = 0.0;

						Pressure[nVertex_NearField] = solution_container->node[iPoint]->GetPressure();
						FaceArea[nVertex_NearField] = abs(Face_Normal[nDim-1]);
						nVertex_NearField ++;
					}
			}

#else

	int nProcessor = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
	int iProcessor;

	unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor*1];
	unsigned long *Buffer_Send_nVertex = new unsigned long [1];

	/*--- Compute the total number of points of the near-field ghost nodes ---*/
	nLocalVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				/*--- Be careful, on the 3D case we must extract a line! ---*/
				if (geometry->node[iPoint]->GetDomain())
					if (((Begin_Int < Coord[0]) && (Coord[0] < End_Int)) && (Coord[nDim-1] < 0.0))
						if (Face_Normal[nDim-1] > 0.0)
							nLocalVertex_NearField ++;
			}

	Buffer_Send_nVertex[0] = nLocalVertex_NearField;

	/*--- Send Near-Field vertex information --*/
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Allgather(Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);

	double *Buffer_Send_Xcoord = new double[MaxLocalVertex_NearField];
	unsigned long *Buffer_Send_IdPoint = new unsigned long [MaxLocalVertex_NearField];
	double *Buffer_Send_Pressure = new double [MaxLocalVertex_NearField];
	double *Buffer_Send_FaceArea = new double[MaxLocalVertex_NearField];

	double *Buffer_Receive_Xcoord = new double[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_Ycoord = new double[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_Zcoord = new double[nProcessor*MaxLocalVertex_NearField];
	unsigned long *Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_Pressure = new double[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_FaceArea = new double[nProcessor*MaxLocalVertex_NearField];

	unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;	
	unsigned long nBuffer_Ycoord = MaxLocalVertex_NearField;	
	unsigned long nBuffer_Zcoord = MaxLocalVertex_NearField;	
	unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;
	unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;
	unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;

	for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
		Buffer_Send_IdPoint[iVertex] = 0; Buffer_Send_Pressure[iVertex] = 0.0;
		Buffer_Send_FaceArea[iVertex] = 0.0; Buffer_Send_Xcoord[iVertex] = 0.0;
	}

	/*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/
	nLocalVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if (geometry->node[iPoint]->GetDomain())
					if (((Begin_Int < Coord[0]) && (Coord[0] < End_Int)) && (Coord[nDim-1] < 0.0))
						if ((Face_Normal[nDim-1] > 0.0) && (geometry->GetnDim() == 2)) {
							Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
							Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
							Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
							Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
							Buffer_Send_Pressure[nLocalVertex_NearField] = solution_container->node[iPoint]->GetPressure();
							Buffer_Send_FaceArea[nLocalVertex_NearField] = abs(Face_Normal[nDim-1]);
							nLocalVertex_NearField++;
						}
			}

	/*--- Send all the information --*/
	MPI::COMM_WORLD.Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Ycoord, nBuffer_Ycoord, MPI::DOUBLE, Buffer_Receive_Ycoord, nBuffer_Ycoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Zcoord, nBuffer_Zcoord, MPI::DOUBLE, Buffer_Receive_Zcoord, nBuffer_Zcoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI::DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, MASTER_NODE);

	if (rank == MASTER_NODE) {

		Xcoord = new double[nVertex_NearField];
		Ycoord = new double[nVertex_NearField];
		Zcoord = new double[nVertex_NearField];
		AzimutalAngle = new double[nVertex_NearField];
		IdPoint = new unsigned long[nVertex_NearField];
		IdDomain = new unsigned long[nVertex_NearField];
		Pressure = new double[nVertex_NearField];
		FaceArea = new double[nVertex_NearField];
		EquivArea = new double[nVertex_NearField];
		TargetArea = new double[nVertex_NearField];
		NearFieldWeight = new double[nVertex_NearField];
		Weight = new double[nVertex_NearField];

		nVertex_NearField = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
			for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
				Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
				Ycoord[nVertex_NearField] = Buffer_Receive_Ycoord[iProcessor*MaxLocalVertex_NearField+iVertex];
				Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
				/*-Compute azimutal angle -*/

				AzimutalAngle[nVertex_NearField] = 0.0;

				/*---*/
				IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
				Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
				FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
				IdDomain[nVertex_NearField] = iProcessor;
				nVertex_NearField++;
			}

	}

	delete [] Buffer_Receive_nVertex;
	delete [] Buffer_Send_nVertex;

	delete [] Buffer_Send_Xcoord;
	delete [] Buffer_Send_Ycoord;
	delete [] Buffer_Send_Zcoord;
	delete [] Buffer_Send_IdPoint;
	delete [] Buffer_Send_Pressure;
	delete [] Buffer_Send_FaceArea;

	delete [] Buffer_Receive_Xcoord;
	delete [] Buffer_Receive_IdPoint;
	delete [] Buffer_Receive_Pressure;
	delete [] Buffer_Receive_FaceArea;

#endif

	if (rank == MASTER_NODE) {

		vector<double> AngleList;
		vector<double>::iterator IterAngleList;

		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			AngleList.push_back(AzimutalAngle[iVertex]);

		sort( AngleList.begin(), AngleList.end());
		IterAngleList = unique( AngleList.begin(), AngleList.end());
		AngleList.resize( IterAngleList - AngleList.begin() );

		nAngle = AngleList.size();

		unsigned short nAngleList[nAngle];
		for (iAngle = 0; iAngle < nAngle; iAngle++)
			nAngleList[iAngle] = 0;

		double Xcoord_Angle[nAngle][nVertex_NearField];
		unsigned long IdPoint_Angle[nAngle][nVertex_NearField];
		unsigned long IdDomain_Angle[nAngle][nVertex_NearField];
		double Pressure_Angle[nAngle][nVertex_NearField];
		double FaceArea_Angle[nAngle][nVertex_NearField];
		double EquivArea_Angle[nAngle][nVertex_NearField];
		double TargetArea_Angle[nAngle][nVertex_NearField];
		double NearFieldWeight_Angle[nAngle][nVertex_NearField];
		double Weight_Angle[nAngle][nVertex_NearField];

		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			for (iAngle = 0; iAngle < nAngle; iAngle++)
				if (AzimutalAngle[iVertex] == AngleList[iAngle]) {
					double Xcoord_Angle[iAngle][nAngleList[iAngle]];
					unsigned long IdPoint_Angle[iAngle][nAngleList[iAngle]];
					unsigned long IdDomain_Angle[iAngle][nAngleList[iAngle]];
					double Pressure_Angle[iAngle][nAngleList[iAngle]];
					double FaceArea_Angle[iAngle][nAngleList[iAngle]];
					double EquivArea_Angle[iAngle][nAngleList[iAngle]];
					double TargetArea_Angle[iAngle][nAngleList[iAngle]];
					double NearFieldWeight_Angle[iAngle][nAngleList[iAngle]];
					double Weight_Angle[iAngle][nAngleList[iAngle]];
					nAngleList[iAngle]++;					
				}

		/*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/
		for (iAngle = 0; iAngle < nAngle; iAngle++)
			for (iVertex = 0; iVertex < nAngleList[iAngle]; iVertex++)
				for (jVertex = 0; jVertex < nAngleList[iAngle] - 1 - iVertex; jVertex++)
					if (Xcoord_Angle[iAngle][jVertex] > Xcoord_Angle[iAngle][jVertex+1]) {
						auxCoord = Xcoord_Angle[iAngle][jVertex]; Xcoord_Angle[iAngle][jVertex] = Xcoord_Angle[iAngle][jVertex+1]; Xcoord_Angle[iAngle][jVertex+1] = auxCoord;
						auxPress = Pressure_Angle[iAngle][jVertex]; Pressure_Angle[iAngle][jVertex] = Pressure_Angle[iAngle][jVertex+1]; Pressure_Angle[iAngle][jVertex+1] = auxPress;
						auxArea = FaceArea_Angle[iAngle][jVertex]; FaceArea_Angle[iAngle][jVertex] = FaceArea_Angle[iAngle][jVertex+1]; FaceArea_Angle[iAngle][jVertex+1] = auxArea;
						auxPoint = IdPoint_Angle[iAngle][jVertex]; IdPoint_Angle[iAngle][jVertex] = IdPoint_Angle[iAngle][jVertex+1]; IdPoint_Angle[iAngle][jVertex+1] = auxPoint;
						auxDomain = IdDomain_Angle[iAngle][jVertex]; IdDomain_Angle[iAngle][jVertex] = IdDomain_Angle[iAngle][jVertex+1]; IdDomain_Angle[iAngle][jVertex+1] = auxDomain;
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
			if (iExtIter == 0) { cout << "There is no Target Equivalent Area file (TargetEA.dat)!!"<< endl;
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

			/*			if ((iVertex > 0) && (iVertex < nVertex_NearField-1))
			 Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]) + 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
			 if (iVertex == 0) Weight[iVertex] = 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
			 if (iVertex == nVertex_NearField-1) Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]); */

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
		EquivArea_file.open("EquivArea.plt", ios::out);
		EquivArea_file << "TITLE = \"SU2 Equivalent area computation\"" << endl;
		EquivArea_file << "VARIABLES = \"x_coord\",\"Equivalent_Area\",\"Target_Area\",\"NearField_Weight\",\"Cp\"" << endl;
		EquivArea_file << "ZONE T= \"Equivalent Area\"" << endl;

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

		/*--- Delete structures ---*/
		delete [] Xcoord; delete [] IdPoint;
		delete [] Pressure; delete [] FaceArea;
		delete [] EquivArea; delete [] TargetArea;
		delete [] NearFieldWeight; delete [] Weight;

	}

#ifdef NO_MPI

	/*--- Store the value of the NearField coefficient ---*/
	solution_container->SetTotal_CEquivArea(InverseDesign);

#else

	/*--- Send the value of the NearField coefficient to all the processors ---*/
	MPI::COMM_WORLD.Bcast (&InverseDesign, 1, MPI::DOUBLE, MASTER_NODE);

	/*--- Store the value of the NearField coefficient ---*/
	solution_container->SetTotal_CEquivArea(InverseDesign);

#endif
}

#else

void COutput::SetEquivalentArea(CSolution *solution_container, CGeometry *geometry,
		CConfig *config, unsigned long iExtIter) {

	ofstream EquivArea_file, FuncGrad_file;
	unsigned short iMarker = 0, iDim, iVar;
	double auxCoord, InverseDesign, DeltaX, Coord_i, Coord_j, jp1Coord, 
	*Coord = NULL, MeanFuntion, *Face_Normal = NULL, auxArea, auxPress, 
	Mach, Beta, Begin_Int, End_Int, R_Plane, Pressure_Inf, Density_Inf, 
	RefAreaCoeff, ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, 
	*Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, 
	*NearFieldWeight = NULL, *Weight = NULL, jFunction, jp1Function;
	unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint, 
			*IdPoint = NULL, *IdDomain = NULL, auxDomain;

	unsigned short nDim = geometry->GetnDim();
	int rank = MASTER_NODE;

	Mach  = config->GetMach_FreeStreamND();
	Beta = sqrt(Mach*Mach-1.0);
	Begin_Int = config->GetEA_IntLimit(0)*config->GetConversion_Factor();
	End_Int = config->GetEA_IntLimit(1)*config->GetConversion_Factor();
	R_Plane = abs(config->GetEA_IntLimit(2)*config->GetConversion_Factor());
	Pressure_Inf = config->GetPressure_FreeStreamND();
	Density_Inf = config->GetDensity_FreeStreamND();	
	RefAreaCoeff = config->GetRefAreaCoeff();
	Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
	Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
	Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
	ModVelocity_Inf = 0;
	for (iDim = 0; iDim < 3; iDim++)
		ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];
	ModVelocity_Inf = sqrt(ModVelocity_Inf);

	factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);

#ifdef NO_MPI

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
							((geometry->GetnDim() == 2) || 
									((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3)))) {
						nVertex_NearField ++;
					}
			}

	/*--- Create an array with all the coordinates, points, pressures, face area,
	 equivalent area, and nearfield weight ---*/
	Xcoord = new double[nVertex_NearField];
	IdPoint = new unsigned long[nVertex_NearField];
	IdDomain = new unsigned long[nVertex_NearField];
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
							((geometry->GetnDim() == 2) || 
									((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3)))) {
						IdPoint[nVertex_NearField] = iPoint;
						Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
						Pressure[nVertex_NearField] = solution_container->node[iPoint]->GetPressure();
						FaceArea[nVertex_NearField] = abs(Face_Normal[nDim-1]);
						nVertex_NearField ++;
					}
			}

#else

	int nProcessor = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
	int iProcessor;

	unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor*1];
	unsigned long *Buffer_Send_nVertex = new unsigned long [1];

	/*--- Compute the total number of points of the near-field ghost nodes ---*/
	nLocalVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				/*--- Be careful, on the 3D case we must extract a line! ---*/
				if (geometry->node[iPoint]->GetDomain())
					if (((Begin_Int < Coord[0]) && (Coord[0] < End_Int)) && (Coord[nDim-1] < 0.0))
						if ((Face_Normal[nDim-1] > 0.0) && ((geometry->GetnDim() == 2) ||
								((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3))))
							nLocalVertex_NearField ++;
			}

	Buffer_Send_nVertex[0] = nLocalVertex_NearField;

	/*--- Send Near-Field vertex information --*/
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Allgather(Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);

	double *Buffer_Send_Xcoord = new double[MaxLocalVertex_NearField];
	unsigned long *Buffer_Send_IdPoint = new unsigned long [MaxLocalVertex_NearField];
	double *Buffer_Send_Pressure = new double [MaxLocalVertex_NearField];
	double *Buffer_Send_FaceArea = new double[MaxLocalVertex_NearField];

	double *Buffer_Receive_Xcoord = new double[nProcessor*MaxLocalVertex_NearField];
	unsigned long *Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_Pressure = new double[nProcessor*MaxLocalVertex_NearField];
	double *Buffer_Receive_FaceArea = new double[nProcessor*MaxLocalVertex_NearField];

	unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;	
	unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;
	unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;
	unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;

	for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
		Buffer_Send_IdPoint[iVertex] = 0; Buffer_Send_Pressure[iVertex] = 0.0;
		Buffer_Send_FaceArea[iVertex] = 0.0; Buffer_Send_Xcoord[iVertex] = 0.0;
	}

	/*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/
	nLocalVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if (geometry->node[iPoint]->GetDomain())
					if (((Begin_Int < Coord[0]) && (Coord[0] < End_Int)) && (Coord[nDim-1] < 0.0))
						if ((Face_Normal[nDim-1] > 0.0) && ((geometry->GetnDim() == 2) ||
								((geometry->GetnDim() == 3) && (abs(Coord[1]) < 1E-3)))) {
							Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
							Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
							//							cout << Buffer_Send_Xcoord[nLocalVertex_NearField] << endl;
							Buffer_Send_Pressure[nLocalVertex_NearField] = solution_container->node[iPoint]->GetPressure();
							Buffer_Send_FaceArea[nLocalVertex_NearField] = abs(Face_Normal[nDim-1]);
							nLocalVertex_NearField++;
						}
			}

	/*--- Send all the information --*/
	MPI::COMM_WORLD.Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI::UNSIGNED_LONG, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI::DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI::DOUBLE, MASTER_NODE);

	if (rank == MASTER_NODE) {

		Xcoord = new double[nVertex_NearField];
		IdPoint = new unsigned long[nVertex_NearField];
		IdDomain = new unsigned long[nVertex_NearField];
		Pressure = new double[nVertex_NearField];
		FaceArea = new double[nVertex_NearField];
		EquivArea = new double[nVertex_NearField];
		TargetArea = new double[nVertex_NearField];
		NearFieldWeight = new double[nVertex_NearField];
		Weight = new double[nVertex_NearField];

		nVertex_NearField = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
			for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
				Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
				IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
				Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
				FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
				IdDomain[nVertex_NearField] = iProcessor;
				nVertex_NearField++;
			}

	}

	delete [] Buffer_Receive_nVertex;
	delete [] Buffer_Send_nVertex;

	delete [] Buffer_Send_Xcoord;
	delete [] Buffer_Send_IdPoint;
	delete [] Buffer_Send_Pressure;
	delete [] Buffer_Send_FaceArea;

	delete [] Buffer_Receive_Xcoord;
	delete [] Buffer_Receive_IdPoint;
	delete [] Buffer_Receive_Pressure;
	delete [] Buffer_Receive_FaceArea;

#endif

	if (rank == MASTER_NODE) {

		/*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/
		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			for (jVertex = 0; jVertex < nVertex_NearField - 1 - iVertex; jVertex++)
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
			if (iExtIter == 0) { cout << "There is no Target Equivalent Area file (TargetEA.dat)!!"<< endl;
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

			/*			if ((iVertex > 0) && (iVertex < nVertex_NearField-1))
			 Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]) + 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
			 if (iVertex == 0) Weight[iVertex] = 0.5*(Xcoord[iVertex+1]-Xcoord[iVertex]);
			 if (iVertex == nVertex_NearField-1) Weight[iVertex] = 0.5*(Xcoord[iVertex]-Xcoord[iVertex-1]); */

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
		EquivArea_file.open("EquivArea.plt", ios::out);
		EquivArea_file << "TITLE = \"SU2 Equivalent area computation\"" << endl;
		EquivArea_file << "VARIABLES = \"x_coord\",\"Equivalent_Area\",\"Target_Area\",\"NearField_Weight\",\"Cp\"" << endl;
		EquivArea_file << "ZONE T= \"Equivalent Area\"" << endl;

		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++) {
			Coord_i = Xcoord[iVertex];
			EquivArea_file << scientific << Xcoord[iVertex] << ", " << EquivArea[iVertex]
			                                                                     << ", " << TargetArea[iVertex] << ", " << NearFieldWeight[iVertex] << ", " <<
			                                                                     (Pressure[iVertex]-Pressure_Inf)/Pressure_Inf << endl;
		}
		EquivArea_file.close();


		FuncGrad_file.precision(15);
		FuncGrad_file.open("WeightNF.dat", ios::out);
		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			FuncGrad_file << scientific << Xcoord[iVertex] <<"\t"<< NearFieldWeight[iVertex] << endl;
		FuncGrad_file.close();

		/*--- Delete structures ---*/
		delete [] Xcoord; delete [] IdPoint;
		delete [] Pressure; delete [] FaceArea;
		delete [] EquivArea; delete [] TargetArea;
		delete [] NearFieldWeight; delete [] Weight;

	}

#ifdef NO_MPI

	/*--- Store the value of the NearField coefficient ---*/
	solution_container->SetTotal_CEquivArea(InverseDesign);

#else

	/*--- Send the value of the NearField coefficient to all the processors ---*/
	MPI::COMM_WORLD.Bcast (&InverseDesign, 1, MPI::DOUBLE, MASTER_NODE);

	/*--- Store the value of the NearField coefficient ---*/
	solution_container->SetTotal_CEquivArea(InverseDesign);

#endif


}
#endif

void COutput::SetFreeSurface(CSolution *solution_container, CGeometry *geometry,
		CConfig *config, unsigned long iExtIter) {

	double *coord = NULL, dist2, *iCoord = NULL, *jCoord = NULL, *U_i = NULL, *U_j = NULL, 
			**Coord_TargetLevelSet = NULL, **Coord_LevelSet = NULL, CoordTA_aux, TargetLevelSet_aux,
			*xCoord = NULL, *yCoord = NULL, auxCoordx, auxCoordy;
	unsigned short iDim;
	unsigned long iPoint, jPoint, iVertex, jVertex, nVertex_LevelSet, iEdge, nPointTargetLevelSet;
	ifstream index_file;
	string text_line;
	int rank = MASTER_NODE;

	unsigned short nDim = geometry->GetnDim();

#ifdef NO_MPI

	/*--- Identification of the 0 level set points and coordinates ---*/
	nVertex_LevelSet = 0;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = solution_container->node[iPoint]->GetSolution();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = solution_container->node[jPoint]->GetSolution();
		if (U_i[0]*U_j[0] < 0.0) nVertex_LevelSet ++;
	}

	/*--- Allocate vector of boundary coordinates ---*/
	Coord_LevelSet = new double* [nVertex_LevelSet];
	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
		Coord_LevelSet[iVertex] = new double [nDim];

	/*--- Get coordinates of the points of the surface ---*/
	nVertex_LevelSet = 0;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = solution_container->node[iPoint]->GetSolution(); iCoord = geometry->node[iPoint]->GetCoord();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = solution_container->node[jPoint]->GetSolution(); jCoord = geometry->node[jPoint]->GetCoord();
		if (U_i[0]*U_j[0] < 0.0) {
			for (iDim = 0; iDim < nDim; iDim++)
				Coord_LevelSet[nVertex_LevelSet][iDim] = iCoord[iDim]-U_i[0]*(jCoord[iDim]-iCoord[iDim])/(U_j[0]-U_i[0]);
			nVertex_LevelSet++;
		}
	}

#else

	int iProcessor;
	unsigned long *Buffer_Send_nVertex = NULL, *Buffer_Receive_nVertex = NULL, 
			nLocalVertex_LevelSet = 0, nGlobalVertex_LevelSet = 0, MaxLocalVertex_LevelSet = 0, nBuffer;
	double *Buffer_Send_Coord = NULL, *Buffer_Receive_Coord = NULL;

	int nProcessor = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	Buffer_Send_nVertex = new unsigned long [1];
	Buffer_Receive_nVertex = new unsigned long [nProcessor];

	nLocalVertex_LevelSet = 0;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = solution_container->node[iPoint]->GetSolution();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = solution_container->node[jPoint]->GetSolution();
		if (U_i[0]*U_j[0] < 0.0) nLocalVertex_LevelSet ++;
	}

	Buffer_Send_nVertex[0] = nLocalVertex_LevelSet;	

	MPI::COMM_WORLD.Allreduce(&nLocalVertex_LevelSet, &nGlobalVertex_LevelSet, 1, MPI::UNSIGNED_LONG, MPI::SUM); 	
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_LevelSet, &MaxLocalVertex_LevelSet, 1, MPI::UNSIGNED_LONG, MPI::MAX); 	
	MPI::COMM_WORLD.Allgather(Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);

	nBuffer = MaxLocalVertex_LevelSet*nDim;
	Buffer_Send_Coord = new double [nBuffer];
	Buffer_Receive_Coord = new double [nProcessor*nBuffer];

	for (iVertex = 0; iVertex < MaxLocalVertex_LevelSet; iVertex++)
		for (iDim = 0; iDim < nDim; iDim++)
			Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;

	nLocalVertex_LevelSet = 0;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0); U_i = solution_container->node[iPoint]->GetSolution(); iCoord = geometry->node[iPoint]->GetCoord();
		jPoint = geometry->edge[iEdge]->GetNode(1); U_j = solution_container->node[jPoint]->GetSolution(); jCoord = geometry->node[jPoint]->GetCoord();

		if (U_i[0]*U_j[0] < 0.0) {
			for (iDim = 0; iDim < nDim; iDim++)
				Buffer_Send_Coord[nLocalVertex_LevelSet*nDim+iDim] = iCoord[iDim]-U_i[0]*(jCoord[iDim]-iCoord[iDim])/(U_j[0]-U_i[0]);
			nLocalVertex_LevelSet++;
		}
	}

	MPI::COMM_WORLD.Allgather(Buffer_Send_Coord, nBuffer, MPI::DOUBLE, Buffer_Receive_Coord, nBuffer, MPI::DOUBLE);

	/*--- Identification of the 0 level set points and coordinates ---*/
	nVertex_LevelSet = 0;
	for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
		nVertex_LevelSet += Buffer_Receive_nVertex[iProcessor];

	/*--- Allocate vector of boundary coordinates ---*/
	Coord_LevelSet = new double* [nVertex_LevelSet];
	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
		Coord_LevelSet[iVertex] = new double [nDim];

	/*--- Set the value of the coordinates at the level set zero ---*/
	nVertex_LevelSet = 0;
	for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
		for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
			for (iDim = 0; iDim < nDim; iDim++) 
				Coord_LevelSet[nVertex_LevelSet][iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_LevelSet+iVertex)*nDim+iDim];
			nVertex_LevelSet++;
		}

	delete [] Buffer_Send_Coord;
	delete [] Buffer_Receive_Coord;
	delete [] Buffer_Send_nVertex;
	delete [] Buffer_Receive_nVertex;	

#endif


	/*--- Order the arrays (x Coordinate, y Coordinate) ---*/
	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
		for (jVertex = 0; jVertex < nVertex_LevelSet - 1 - iVertex; jVertex++)
			if (Coord_LevelSet[jVertex][0] > Coord_LevelSet[jVertex+1][0]) {
				auxCoordx = Coord_LevelSet[jVertex][0]; Coord_LevelSet[jVertex][0] = Coord_LevelSet[jVertex+1][0]; Coord_LevelSet[jVertex+1][0] = auxCoordx;
				auxCoordy = Coord_LevelSet[jVertex][1]; Coord_LevelSet[jVertex][1] = Coord_LevelSet[jVertex+1][1]; Coord_LevelSet[jVertex+1][1] = auxCoordy;
			}

	index_file.open("Target_LevelSet.dat", ios::in);
	if (index_file.fail()) {
		if ((iExtIter == 0) && (rank == MASTER_NODE)) { 
			cout << "There is no Target Level Set file (Target_LevelSet.dat)!!"<< endl;
			cout << "Using default parameters (Target Level Set = 0.0)" << endl;
		}
		nPointTargetLevelSet = 2;

		xCoord = new double [nPointTargetLevelSet];
		yCoord = new double [nPointTargetLevelSet];

		xCoord[0] = -10E6; yCoord[0] = config->GetFreeSurface_Zero();
		xCoord[1] = 10E6; yCoord[1] = config->GetFreeSurface_Zero();

	}
	else {

		/*--- Dimensionalization bucle ---*/
		nPointTargetLevelSet = 0;
		while (!index_file.eof()) {
			getline(index_file, text_line);
			istringstream value_line(text_line);
			value_line >> CoordTA_aux >> TargetLevelSet_aux;
			nPointTargetLevelSet++;
		}
		index_file.close();

		xCoord = new double [nPointTargetLevelSet];
		yCoord = new double [nPointTargetLevelSet];

		/*--- Store the information ---*/
		index_file.open("Target_LevelSet.dat", ios::in);
		nPointTargetLevelSet = 0;
		while (!index_file.eof()) {
			getline(index_file, text_line);
			istringstream value_line(text_line);
			value_line >> xCoord[nPointTargetLevelSet] >> yCoord[nPointTargetLevelSet];
			nPointTargetLevelSet++;
		}
		index_file.close();
	}


	/*--- Interpolate the target free surface in the original free surface ---*/
	Coord_TargetLevelSet = new double* [nVertex_LevelSet];
	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
		Coord_TargetLevelSet[iVertex] = new double [nDim];

	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)	
		for (unsigned long iVar = 0; iVar < nPointTargetLevelSet-1; iVar++) {
			if ((xCoord[iVar] <= Coord_LevelSet[iVertex][0]) && (Coord_LevelSet[iVertex][0] <= xCoord[iVar+1])) {
				Coord_TargetLevelSet[iVertex][0] = Coord_LevelSet[iVertex][0];
				Coord_TargetLevelSet[iVertex][1] = yCoord[iVar] + (Coord_LevelSet[iVertex][0]-xCoord[iVar])*(yCoord[iVar+1]-yCoord[iVar] )/(xCoord[iVar+1]-xCoord[iVar]);
			}
		}


	/*--- Get coordinates of the points and compute distances to the surface ---*/
	double FreeSurface = 0.0;
	double factor = 1.0;

	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

		coord = geometry->node[iPoint]->GetCoord();
		double volume = geometry->node[iPoint]->GetVolume();

		double LevelSetDiff_Squared = 0.0;
		double LevelSetDiff = 0.0;

		if ((coord[0] > 0) && (coord[0] < 3)){

			/*--- Compute the squared distance to the rest of points, and get the minimum ---*/
			double dist_Target = 1E20;
			for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
				dist2 = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					dist2 += (coord[iDim]-Coord_TargetLevelSet[iVertex][iDim])*(coord[iDim]-Coord_TargetLevelSet[iVertex][iDim]);
				if (dist2 < fabs(dist_Target)) { 
					if (Coord_TargetLevelSet[iVertex][nDim-1] < coord[nDim-1]) dist_Target = dist2; 
					else dist_Target = -dist2; 
				}
			}
			double NumberSign = 1.0;
			if (dist_Target != 0.0) NumberSign = dist_Target/fabs(dist_Target);
			dist_Target = sqrt(fabs(dist_Target))*NumberSign;

			LevelSetDiff = (solution_container->node[iPoint]->GetSolution()[0] - dist_Target);
			LevelSetDiff_Squared = LevelSetDiff*LevelSetDiff;
			FreeSurface += 0.5*factor*LevelSetDiff_Squared*volume;

		}

		solution_container->node[iPoint]->SetDiffLevelSet(LevelSetDiff);

	}

	if ((rank == MASTER_NODE) && (iExtIter % config->GetWrt_Sol_Freq() == 0)) {

		ofstream LevelSet_file;

		/*--- Write the Level Set distribution, the target level set---*/
		LevelSet_file.precision(15);


		/*--- Write file name with extension ---*/
		char cstr[200], buffer[50];

		strcpy (cstr, "LevelSet");
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()){
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.plt", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.plt", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.plt", int(iExtIter));
		}
		else
			sprintf (buffer, ".plt");

		strcat(cstr,buffer);

		LevelSet_file.open(cstr, ios::out);	

		LevelSet_file << "TITLE = \"SU2 Level Set computation\"" << endl;
		LevelSet_file << "VARIABLES = \"x coord\",\"y coord (original)\",\"y coord (target)\"" << endl;
		LevelSet_file << "ZONE T= \"Free Surface\"" << endl;

		for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
			LevelSet_file << scientific << Coord_LevelSet[iVertex][0] << ", " << Coord_LevelSet[iVertex][1] << ", " << Coord_TargetLevelSet[iVertex][1] << endl;
		}
		LevelSet_file.close();

	}


	/*--- Store the value of the NearField coefficient ---*/
	solution_container->SetTotal_CFreeSurface(FreeSurface);

	/*--- Deallocate vector of boundary coordinates ---*/
	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
		delete [] Coord_LevelSet[iVertex];
	delete [] Coord_LevelSet;

	delete [] xCoord;
	delete [] yCoord;


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
	M2 = config->GetParticle_Mass(1);


	nDim = geometry->GetnDim();



	/*--- PRINT OUT ELECTRIC FIELD Vector ---*/

	if (config->GetElectricSolver()) {
		/*--- PRINT OUT ELECTRIC POTENTIAL ---*/
		WriteInOutputFile(geometry, solution_container[ELEC_SOL], ConsVar_file, "Electric_Potential", scalar, 0, "GetSolution", config);

		ConsVar_file << "VECTORS Electric_Field float" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
				Vector[iDim] = -1.0* solution_container[ELEC_SOL]->node[iPoint]->GetGradient(0,iDim);
			}
			if (geometry->GetnDim() == 2) ConsVar_file << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
			if (geometry->GetnDim() == 3) ConsVar_file << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
		}
	}


#ifdef OutputSource

	/*--- PRINT OUT SOURCE TERM 1 CONTRIBUTION FIELD Vector ---*/
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src1", scalar, 0, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src2", scalar, 1, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src3", scalar, 2, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src4", scalar, 3, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src5", scalar, 4, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src6", scalar, 5, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src7", scalar, 6, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src8", scalar, 7, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src9", scalar, 8, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src10", scalar, 9, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src11", scalar, 10, "GetSource", config);
	WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src12", scalar, 11, "GetSource", config);

	if (nDim == 3) {
		WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src13", scalar, 12, "GetSource", config);
		WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src14", scalar, 13, "GetSource", config);
		WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Src15", scalar, 14, "GetSource", config);
	}

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


#endif

	/*--- PRINT OUT THE NET CHARGE IN DOMAIN ---*/
	ConsVar_file << "SCALARS Net_Charge float 1" << endl;
	ConsVar_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		ConsVar_file << ELECTRON_CHARGE*(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))- solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))*M2/M3)/M2  << endl;

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
	//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species1", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);

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
	ConsVar_file << "SCALARS Argon_Mach float 1" << endl;
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
			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0))/M2 << endl;

		/*--- PRINT OUT DENSITY* Energy OF SPECIES 2 ---*/
		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species2", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);

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
		ConsVar_file << "SCALARS Ion_Mach float 1" << endl;
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
			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0))/M3 << endl;

		/*--- PRINT OUT DENSITY * Energy OF SPECIES 3 ---*/
		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species3", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);

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
		ConsVar_file << "SCALARS Electron_Mach float 1" << endl;
		ConsVar_file << "LOOKUP_TABLE default" << endl;
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			ConsVar_file << fixed <<  sqrt(solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))/solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies) << endl;

	}




}

