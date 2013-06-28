/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
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


COutput::COutput(void) {

  /*--- Initialize point and connectivity counters to zero. ---*/
  nGlobal_Poin = 0;
  nGlobal_Elem = 0;
  nGlobal_Tria = 0;
  nGlobal_Quad = 0;
  nGlobal_Tetr = 0;
  nGlobal_Hexa = 0;
  nGlobal_Wedg = 0;
  nGlobal_Pyra = 0;
  nGlobal_Line = 0;

  /*--- Initialize CGNS write flag ---*/
  wrote_CGNS_base = false;

  /*--- Initialize Tecplot write flag ---*/
  wrote_Tecplot_base = false;

}

COutput::~COutput(void) { }

//void COutput::SetDomain_Flow(CConfig *config, CGeometry *geometry, CSolution **solution_container, 
//		unsigned long iExtIter, unsigned short val_iZone, unsigned short val_nZone) {
//	char cstr[200], buffer[50];
//	ofstream ConsVar_file;
//	ofstream ConsVar_file_mean;
//	unsigned short iSpecies, iVar, loc, iDim;
//	unsigned long iPoint, iElem;
//	double coord[3] = {0.0,0.0,0.0}, velocity[3] = {0.0,0.0,0.0}, mach_number;
//	double rel_velocity[3] = {0.0,0.0,0.0}, vorticity[3] = {0.0,0.0,0.0}, rel_mach_number = 0.0, vorticity_mag = 0.0;
//	double Volume;
//
//	bool incompressible = config->GetIncompressible();
//	bool transition = (config->GetKind_Trans_Model()==LM);
//	bool scalar = true;
//	bool rotating_frame = config->GetRotating_Frame();
//	unsigned short solver = config->GetKind_Solver();
//
//	if (solver == FLUID_STRUCTURE_EULER) {
//		if (val_iZone == 0) solver = EULER;
//		if (val_iZone == 1) solver = LINEAR_ELASTICITY;
//	}
//
//	if (solver == AEROACOUSTIC_EULER) {
//		if (val_iZone == 0) solver = EULER;
//		if (val_iZone == 1) solver = WAVE_EQUATION;
//	}
//
//	if (solver == ADJ_AEROACOUSTIC_EULER) {
//		if (val_iZone == 0) solver = ADJ_EULER;
//		if (val_iZone == 1) solver = WAVE_EQUATION;
//	}
//
//	if (solver == PLASMA_EULER) {
//		if (val_iZone == 0) solver = PLASMA_EULER;
//		if (val_iZone == 1) solver = ELECTRIC_POTENTIAL;
//	}
//
//	if (solver == PLASMA_NAVIER_STOKES) {
//		if (val_iZone == 0) solver = PLASMA_NAVIER_STOKES;
//		if (val_iZone == 1) solver = ELECTRIC_POTENTIAL;
//	}
//
//	/*--- Paraview Output ---*/
//	if (config->GetOutput_FileFormat() == PARAVIEW) {
//
//		/*--- Write file name with extension ---*/
//		strcpy (cstr, config->GetFlow_FileName().c_str());
//		if (solver == ELECTRIC_POTENTIAL) strcpy (cstr, config->GetStructure_FileName().c_str());
//
//		if (config->GetUnsteady_Simulation())
//			sprintf (buffer, "_%d.vtk", int(iExtIter));
//		else
//			sprintf (buffer, ".vtk");
//
//		strcat(cstr,buffer);
//
//		/*--- Write header and open the file ---*/
//		geometry->SetParaView(cstr);
//		ConsVar_file.precision(15);
//		ConsVar_file.open(cstr, ios::out | ios::app);
//		ConsVar_file << "POINT_DATA " << geometry->GetnPoint() << endl;
//
//		if (solver == ELECTRIC_POTENTIAL) {
//			WriteInOutputFile(geometry, solution_container[ELEC_SOL], ConsVar_file, "Electric_Potential", scalar, 0, "GetSolution", config);
//			WriteInOutputFile(geometry, solution_container[ELEC_SOL], ConsVar_file, "NegElectric_Field", !scalar, 0, "GetGradient", config);
//			return;
//		}
//
//		else if ((solver == PLASMA_EULER) || (solver == PLASMA_NAVIER_STOKES)) {
//			unsigned short nDim = geometry->GetnDim();
//			unsigned short nSpecies = solution_container[PLASMA_SOL]->GetnSpecies();
//			unsigned short nDiatomics = solution_container[PLASMA_SOL]->GetnDiatomics();
//
//			if (config->GetKind_GasModel() == ARGON || config->GetKind_GasModel() == AIR21)
//				WriteReactingOutputFile(config, geometry,solution_container, ConsVar_file);
//			else if (config->GetKind_GasModel()==AIR7 || config->GetKind_GasModel() == O2 || config->GetKind_GasModel() == N2|| config->GetKind_GasModel() == AIR5 || config-> GetKind_GasModel() == ARGON_SID) {
//				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//					if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
//					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
//					WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density", scalar,  loc+0, "GetSolution", config);
//					//					WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "", scalar,  5, "GetSolution", config);
//				}
//			}
//			return;
//		}
//
//		else if  ((solver == EULER) || (solver == NAVIER_STOKES) || (solver == RANS) ||
//				(solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//
//			if (incompressible) {
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Pressure", scalar, 0, "GetSolution", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Velocity", !scalar, 0, "GetSolution", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density", scalar,  0, "GetDensityInc", config);
//				if (config->GetKind_SlopeLimit_Flow() != NONE) 
//					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Limiter", scalar,  0, "GetLimiter", config);
//				//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Viscosity", scalar, 0, "GetLaminarViscosityInc", config);
//				//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "UndLaplacian", scalar,  0, "GetUndLaplacian", config);
//			}
//
//			else {
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density", scalar,  0, "GetSolution", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density_x_Velocity", !scalar, 0, "GetSolution", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Density_x_Energy", scalar, geometry->GetnDim()+1, "GetSolution", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Velocity", !scalar,  0, "GetVelocity", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Pressure", scalar,  0, "GetPressure", config);
//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Mach_Number", scalar,  0, "GetMach", config);
//				//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "UndLaplacian", scalar,  0, "GetUndLaplacian", config);
//				//				WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Limiter", scalar,  0, "GetLimiter", config);
//
//				if (config->GetGrid_Movement())
//					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Grid_Velocity", !scalar,  0, "GetGridVel", config);
//
//				if (config->GetRotating_Frame())
//					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Rotational_Velocity", !scalar,  0, "GetRotVel", config);
//
//				if ((solver == NAVIER_STOKES) || (solver == RANS) || 
//						(solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Temperature", scalar,  0, "GetTemperature", config);
//					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Laminar_Viscosity", scalar,  0, "GetLaminarViscosity", config);
//				}
//
//				if (solver == RANS) {
//					WriteInOutputFile(geometry, solution_container[FLOW_SOL], ConsVar_file, "Eddy_Viscosity", scalar,  0, "GetEddyViscosity", config);
//					switch (config->GetKind_Turb_Model()){
//					case SA:	WriteInOutputFile(geometry, solution_container[TURB_SOL], ConsVar_file,  "Nu_Tilde", scalar, 0, "Solution", config); break;
//					case SST:	WriteInOutputFile(geometry, solution_container[TURB_SOL], ConsVar_file,  "kine", scalar, 0, "Solution", config);
//					WriteInOutputFile(geometry, solution_container[TURB_SOL], ConsVar_file,  "omega", scalar, 1, "Solution", config); break;
//					}
//				}
//			}
//
//			if ((solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS))
//				WriteInOutputFile(geometry, solution_container[LEVELSET_SOL], ConsVar_file,   "LevelSet", scalar, 0, "GetSolution", config);
//
//			return;
//		}
//	}
//
//	/*--- Tecplot Output ---*/
//	if (config->GetOutput_FileFormat() == TECPLOT) {
//
//		/*--- Write file name with extension ---*/
//		strcpy (cstr, config->GetFlow_FileName().c_str());
//		if (solver == LINEAR_ELASTICITY) strcpy (cstr, config->GetStructure_FileName().c_str());
//		if (solver == WAVE_EQUATION) strcpy (cstr, config->GetWave_FileName().c_str());
//		if ((solver == WAVE_EQUATION) && (config->GetKind_Solver() == ADJ_AEROACOUSTIC_EULER)) strcpy (cstr, config->GetAdjWave_FileName().c_str());
//		if (solver == ELECTRIC_POTENTIAL) strcpy (cstr, config->GetStructure_FileName().c_str());
//
//		/*--- Special cases where a number needs to be appended to the file name. ---*/
//		if ((solver == EULER || solver == NAVIER_STOKES || solver == RANS) && (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
//			sprintf (buffer, "_%d", int(val_iZone));
//			strcat(cstr,buffer);
//		}
//
//		if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
//			if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.plt", int(val_iZone));
//			if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.plt", int(val_iZone));
//			if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.plt", int(val_iZone));
//			if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.plt", int(val_iZone));
//			if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.plt", int(val_iZone));
//
//		} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
//			if (int(iExtIter) < 10) sprintf (buffer, "_0000%d.plt", int(iExtIter));
//			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.plt", int(iExtIter));
//			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.plt", int(iExtIter));
//			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.plt", int(iExtIter));
//			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.plt", int(iExtIter));
//		} else {
//			sprintf (buffer, ".plt");
//		}
//
//		strcat(cstr,buffer);
//
//		/*--- Write header and open the file ---*/
//		ConsVar_file.precision(15);
//		ConsVar_file.open(cstr, ios::out);
//
//		ConsVar_file << "TITLE = \"Visualization of the volumetric grid\"" << endl;
//		ConsVar_file << "VARIABLES = \"x\", \"y\", \"z\"";
//
//		if ((solver == EULER) || (solver == NAVIER_STOKES) || (solver == RANS) || 
//				(solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//			ConsVar_file << ", \"Density\", \"Velocity x\", \"Velocity y\", \"Velocity z\", \"Pressure\", \"Mach Number\", \"1stVar Residual\"";
//			
//			if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) {
//				if (config->GetKind_SlopeLimit_Flow() != NONE)
//					ConsVar_file << ", \"1stVar Limiter\"";
//				if ((config->GetKind_Upwind_Flow() == ROE_TURKEL_1ST) || (config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND)) 
//					ConsVar_file << ", \"Preconditioner_Beta\"";
//			}
//			if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
//				ConsVar_file << ", \"Diss Switch\", \"1stVar UndLapl\"";
//			}
//			if ((geometry->GetnDim() == 3) && (config->GetMagnetic_Force() == YES))  ConsVar_file << ", \"BField_x\", \"BField_y\", \"BField_z\"";
//
//		}
//
//		if ((solver == NAVIER_STOKES) || (solver == RANS) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//			if (!incompressible) {
//				if (geometry->GetnDim() == 2) ConsVar_file << ", \"Temperature\", \"Vorticity z\", \"Distance\", \"Laminar Viscosity\"";
//				if (geometry->GetnDim() == 3) ConsVar_file << ", \"Temperature\", \"Vorticity x\", \"Vorticity y\", \"Vorticity z\", \"Distance\", \"Laminar Viscosity\"";
//			}
//			else {
//				ConsVar_file << ", \"Laminar Viscosity\"";
//            }
//		}
//
//		if ((solver == RANS) || (solver == FREE_SURFACE_RANS)) {
//			switch (config->GetKind_Turb_Model()){
//			case SA:	ConsVar_file << ", \"Eddy Viscosity\", \"Nu Tilde\""; break;
//			case SST:	ConsVar_file << ", \"Eddy Viscosity\", \"kine\", \"omega\""; break;
//			}
//		}
//        
//        if (transition)
//			ConsVar_file << ", \"Intermittency\", \"Re_th\"";
//
//		if ((solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS))
//			ConsVar_file << ", \"LevelSet\", \"LevelSet Residual\", \"LevelSet Diff\"";
//
//		if (solver == WAVE_EQUATION)
//			ConsVar_file << ", \"Wave\", \"Wave Prime\", \"Wave Residual\"";
//
//		if (solver == LINEAR_ELASTICITY) {
//			if (geometry->GetnDim() == 2) ConsVar_file << ", \"Displacement x\", \"Displacement y\", \"Pressure\", \"Stress xx\", \"Stress yy\", \"Stress xy\", \"Source x\", \"Source y\"";
//			if (geometry->GetnDim() == 3) ConsVar_file << ", \"Displacement x\", \"Displacement y\", \"Displacement z\", \"Pressure\", \"Stress xx\", \"Stress yy\", \"Stress zz\", \"Stress xy\", \"Stress xz\", \"Stress yz\", \"Source x\", \"Source y\", \"Source z\"";
//		}
//
//		if (rotating_frame)
//			ConsVar_file << ", \"Relative Vx\", \"Relative Vy\", \"Relative Vz\", \"Relative Mach\"";
//
//		if (solver == ELECTRIC_POTENTIAL) {
//			if (geometry->GetnDim() == 2) ConsVar_file << "\"Phi\", \"ElectricField_x\", \"ElectricField_y\"";
//			if (geometry->GetnDim() == 3) ConsVar_file << "\"Phi\", \"ElectricField_x\", \"ElectricField_y\", \"ElectricField_z\"";
//		}
//
//		if ((solver == PLASMA_EULER) || (solver == PLASMA_NAVIER_STOKES)) {
//
//			if (config->GetKind_GasModel() == ARGON) {
//				ConsVar_file << ", ";
//                
//				if (geometry->GetnDim() == 3)  ConsVar_file << "\"BField_x\", \"BField_y\", \"BField_z\"";
//				ConsVar_file << ", ";
//				for ( iSpecies = 0; iSpecies < solution_container[PLASMA_SOL]->GetnSpecies(); iSpecies++ ) {
//					if (geometry->GetnDim() == 2) ConsVar_file << "\"Density("<< iSpecies << ")\", \"Velocity_x("<< iSpecies << ")\", \"Velocity_y("<< iSpecies << ")\", \"MachNumber("<< iSpecies << ")\", \"Pressure("<< iSpecies << ")\", \"Temperature("<< iSpecies << ")\"";
//					if (geometry->GetnDim() == 3) {
//						ConsVar_file << "\"Density("<< iSpecies << ")\", \"Velocity_x("<< iSpecies << ")\", \"Velocity_y("<< iSpecies << ")\", \"Velocity_z("<< iSpecies << ")\", \"MachNumber("<< iSpecies << ")\", \"Pressure("<< iSpecies << ")\", \"Temperature_tr("<< iSpecies << ")\"";
//					}
//					if (iSpecies < solution_container[PLASMA_SOL]->GetnSpecies() - 1)
//						ConsVar_file << ", ";
//				}
//			}
//			else {
//				ConsVar_file << ", ";
//				for ( iSpecies = 0; iSpecies < solution_container[PLASMA_SOL]->GetnSpecies(); iSpecies++ ) {
//					//				ConsVar_file << "\"Density("<< iSpecies << ")\", \"Velocity_x("<< iSpecies << ")\", \"Velocity_y("<< iSpecies << ")\", \"DensityEnergy("<< iSpecies << ")\",";
//					if (geometry->GetnDim() == 2)  ConsVar_file << "\"Density("<< iSpecies << ")\", \"Velocity_x("<< iSpecies << ")\", \"Velocity_y("<< iSpecies << ")\", \"MachNumber("<< iSpecies << ")\", \"Pressure(" << iSpecies << ")\", \"Temperature_TR(" << iSpecies << ")\",";
//                    else ConsVar_file << "\"Density("<< iSpecies << ")\", \"Velocity_x("<< iSpecies << ")\", \"Velocity_y("<< iSpecies << ")\", \"Velocity_z(" << iSpecies << ")\", \"MachNumber("<< iSpecies << ")\", \"Pressure(" << iSpecies << ")\", \"Temperature_TR(" << iSpecies << ")\",";
//                    
//					if (iSpecies < config->GetnDiatomics()) ConsVar_file << "\"Temperature_vib("<< iSpecies << ")\",";
//					
//                    if (geometry->GetnDim() == 2) ConsVar_file << "\"Source_Density("<< iSpecies << ")\", \"Source_Velocity_x("<< iSpecies << ")\", \"Source_Velocity_y("<< iSpecies << ")\", \"Source_DensityEnergy("<< iSpecies << ")\"";
//                    else ConsVar_file << "\"Source_Density("<< iSpecies << ")\", \"Source_Velocity_x("<< iSpecies << ")\", \"Source_Velocity_y("<< iSpecies << ")\", \"Source_Velocity_z(" << iSpecies << ")\", \"Source_DensityEnergy("<< iSpecies << ")\"";
//                    
//					if (iSpecies < config->GetnDiatomics()) ConsVar_file << ", \"Source_Energy_vib("<< iSpecies << ")\"";
//					if (iSpecies < solution_container[PLASMA_SOL]->GetnSpecies() - 1)
//						ConsVar_file << ", ";
//				}
//			}
//
//			/*--- Check for writing mean solution file ---*/
//			if (config->GetWrite_Mean_Solution()) {
//				char cstr_mean[200], buffer_mean[50];
//				strcpy (cstr_mean, config->GetFlow_FileName().c_str());
//				sprintf (buffer_mean, "_mean.plt");
//				strcat(cstr_mean,buffer_mean);
//
//				/*--- Write header and open the file ---*/
//				ConsVar_file_mean.precision(15);
//				ConsVar_file_mean.open(cstr_mean, ios::out);
//
//				ConsVar_file_mean << "TITLE = \"Visualization of the volumetric grid\"" << endl;
//				ConsVar_file_mean << "VARIABLES = \"x\", \"y\", \"z\"";
//				ConsVar_file_mean << ", \"Density\", \"Velocity x\", \"Velocity y\", \"Energy\"";
//			}
//		}
//
//		ConsVar_file << endl;
//		if (config->GetWrite_Mean_Solution())
//			ConsVar_file_mean << endl;
//
//		ConsVar_file << "ZONE ";
//		if (config->GetWrite_Mean_Solution())
//			ConsVar_file_mean << "ZONE ";
//
//		if (config->GetUnsteady_Simulation() == TIME_SPECTRAL)
//			ConsVar_file << "STRANDID="<<int(val_iZone+1)<<", SOLUTIONTIME="<<int(val_iZone+1)<<", ";
//
//		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
//			ConsVar_file << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*int(iExtIter)<<", ";
//
//		ConsVar_file << "NODES="<< geometry->GetnPoint() <<" , ELEMENTS="<< geometry->GetnElem() <<", DATAPACKING=POINT";
//		if (config->GetWrite_Mean_Solution())
//			ConsVar_file_mean << "NODES="<< geometry->GetnPoint() <<" , ELEMENTS="<< geometry->GetnElem() <<", DATAPACKING=POINT";
//
//		if (geometry->GetnDim() == 2) {
//			ConsVar_file << ", ZONETYPE=FEQUADRILATERAL" << endl;
//			if (config->GetWrite_Mean_Solution())
//				ConsVar_file_mean << ", ZONETYPE=FEQUADRILATERAL" << endl;
//		}
//		if (geometry->GetnDim() == 3) {
//			ConsVar_file << ", ZONETYPE=FEBRICK"<< endl;
//			if (config->GetWrite_Mean_Solution())
//				ConsVar_file_mean << ", ZONETYPE=FEBRICK"<< endl;
//		}
//
//		/*--- Cycle through all points and write the solution. ---*/
//		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//			Volume = geometry->node[iPoint]->GetVolume();
//
//			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
//				coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
//
//			if ((solver == EULER) || (solver == NAVIER_STOKES) || (solver == RANS) ||
//					(solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//
//				/*--- Compute vectorial quantities for 2D and 3D ---*/
//				for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
//					velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, config->GetIncompressible());
//				}
//
//				if ((solver == NAVIER_STOKES) || (solver == RANS) || rotating_frame ||
//					(solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//					vorticity[0] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(0);
//					vorticity[1] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(1);
//					vorticity[2] = solution_container[FLOW_SOL]->node[iPoint]->GetVorticity(2);
//				}
//
//				/*--- Compute the Mach number ---*/
//				mach_number = sqrt(solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2())/
//						solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
//			}
//
//			/*--- Compute relative velocities for rotating frame ---*/
//			if (rotating_frame) {
//				double *rot_vel = geometry->node[iPoint]->GetRotVel();
//				double rel_velocity_mag = 0.0; vorticity_mag = 0.0;
//				for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
//					rel_velocity[iDim] = velocity[iDim] - rot_vel[iDim];
//					rel_velocity_mag  += rel_velocity[iDim]*rel_velocity[iDim];
//					vorticity_mag     += vorticity[iDim]*vorticity[iDim];
//				}
//				rel_velocity_mag = sqrt(rel_velocity_mag);
//				vorticity_mag    = sqrt(vorticity_mag);
//				rel_mach_number  = rel_velocity_mag/solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
//			}
//
//			if (solver == LINEAR_ELASTICITY) {
//				for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
//					coord[iDim] += solution_container[FEA_SOL]->node[iPoint]->GetSolution(iDim);
//			}
//
//			/*--- Write the coordinates ---*/
//			ConsVar_file << coord[0] <<" "<< coord[1] <<" "<< coord[2];
//			if (config->GetWrite_Mean_Solution())
//				ConsVar_file_mean << coord[0] << " " << coord[1] << " " << coord[2];
//
//			/*--- Write the Euler variables ---*/
//			if ((solver == EULER) || (solver == NAVIER_STOKES) || (solver == RANS) ||
//					(solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//
//				if (!incompressible) {
//					mach_number = sqrt(solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2())/
//							solution_container[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
//					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
//						velocity[iDim] = solution_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim, config->GetIncompressible());
//					ConsVar_file << " " << solution_container[FLOW_SOL]->node[iPoint]->GetSolution((unsigned short)0) <<" "<<
//							velocity[0] <<" "<< velocity[1] <<" "<< velocity[2]<<" "<<
//							solution_container[FLOW_SOL]->node[iPoint]->GetPressure(incompressible) <<" "<< mach_number ;
//				}
//				else {
//					double Density = solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref();
//					double Pressure = solution_container[FLOW_SOL]->node[iPoint]->GetSolution((unsigned short)0)*config->GetPressure_Ref();
//					double SoundSpeed = sqrt(config->GetBulk_Modulus()/Density);
//					mach_number = sqrt(solution_container[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/SoundSpeed;
//					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
//						velocity[iDim] = (solution_container[FLOW_SOL]->node[iPoint]->GetSolution(iDim+1)/solution_container[FLOW_SOL]->node[iPoint]->GetDensityInc())*config->GetVelocity_Ref();
//					ConsVar_file << " " << Density <<" "<< velocity[0] <<" "<< velocity[1] <<" "<< velocity[2]<<" "<< Pressure <<" "<< mach_number;
//				}
//
//				ConsVar_file <<" "<< solution_container[FLOW_SOL]->node[iPoint]->GetResidual((unsigned short)0);
//
//				if (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND) {
//					if (config->GetKind_SlopeLimit_Flow() != NONE)
//						ConsVar_file <<" "<< solution_container[FLOW_SOL]->node[iPoint]->GetLimiter((unsigned short)0);
//					if ((config->GetKind_Upwind_Flow() == ROE_TURKEL_1ST) || (config->GetKind_Upwind_Flow() == ROE_TURKEL_2ND))
//						ConsVar_file <<" "<< solution_container[FLOW_SOL]->node[iPoint]->GetPreconditioner_Beta();
//				}
//				
//				if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
//					ConsVar_file <<" "<< solution_container[FLOW_SOL]->node[iPoint]->GetSensor()
//					<<" "<< solution_container[FLOW_SOL]->node[iPoint]->GetUnd_Lapl()[0];
//				
//				if ( (config->GetMagnetic_Force() == YES) && (geometry->GetnDim() == 3)) {
//					ConsVar_file << " " << solution_container[FLOW_SOL]->node[iPoint]->GetMagneticField()[0];
//					ConsVar_file << " " << solution_container[FLOW_SOL]->node[iPoint]->GetMagneticField()[1];
//					ConsVar_file << " " << solution_container[FLOW_SOL]->node[iPoint]->GetMagneticField()[2];
//				}
//
//			}
//
//			/*--- Write the Navier-Stokes variables ---*/
//			if ((solver == NAVIER_STOKES) || (solver == RANS) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
//				if (!incompressible) {
//					if (geometry->GetnDim() == 2) ConsVar_file << " "<< solution_container[FLOW_SOL]->node[iPoint]->GetTemperature() <<" "<< vorticity[2]<<" "<< geometry->node[iPoint]->GetWallDistance()<<" "<<solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
//					if (geometry->GetnDim() == 3) ConsVar_file << " "<< solution_container[FLOW_SOL]->node[iPoint]->GetTemperature() <<" "<< vorticity[0] <<" "<< vorticity[1] <<" "<< vorticity[2]<<" "<< geometry->node[iPoint]->GetWallDistance()<<" "<<solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
//				}
//				else {
//					double viscosity = solution_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc()*config->GetViscosity_Ref();
//					ConsVar_file << " "<< viscosity;
//				}
//			}
//
//			/*--- Write the Turbulence variables ---*/
//			if ((solver == RANS) || (solver == FREE_SURFACE_RANS)) {
//				if (!incompressible){
//					ConsVar_file << " "<< solution_container[TURB_SOL]->node[iPoint]->GetmuT();
//					switch (config->GetKind_Turb_Model()){
//                        case SA:	ConsVar_file << " " << solution_container[TURB_SOL]->node[iPoint]->GetSolution(0); break;
//                        case SST:	ConsVar_file << " " << solution_container[TURB_SOL]->node[iPoint]->GetSolution(0) <<" "<< solution_container[TURB_SOL]->node[iPoint]->GetSolution(1); break;
//					}
//				}
//				else {
//                    ConsVar_file << " "<< solution_container[TURB_SOL]->node[iPoint]->GetmuT()*config->GetViscosity_Ref();
//                    switch (config->GetKind_Turb_Model()){
//                        case SA:	ConsVar_file << " " << solution_container[TURB_SOL]->node[iPoint]->GetSolution(0); break;
//                        case SST:	ConsVar_file << " " << solution_container[TURB_SOL]->node[iPoint]->GetSolution(0) <<" "<< solution_container[TURB_SOL]->node[iPoint]->GetSolution(1); break;
//					}
//                }
//                if (transition)
//					ConsVar_file << " "<< solution_container[TRANS_SOL]->node[iPoint]->GetSolution(0) << " "<< solution_container[TRANS_SOL]->node[iPoint]->GetSolution(1);
//			}
//        
//			
//			/*--- Write the level set solution ---*/
//			if ((solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS))
//				ConsVar_file << " "<< solution_container[LEVELSET_SOL]->node[iPoint]->GetSolution(0) << " "<< solution_container[LEVELSET_SOL]->node[iPoint]->GetResSour()[0] << " "<< solution_container[LEVELSET_SOL]->node[iPoint]->GetDiffLevelSet();
//
//			/*--- Write the wave solution ---*/
//			if (solver == WAVE_EQUATION)
//				ConsVar_file << " "<< solution_container[WAVE_SOL]->node[iPoint]->GetSolution(0)
//				<< " "<< solution_container[WAVE_SOL]->node[iPoint]->GetSolution(1)
//				<< " "<< solution_container[WAVE_SOL]->node[iPoint]->GetResVisc()[0];
//
//			if (solver == LINEAR_ELASTICITY) {
//				double Displ_x = solution_container[FEA_SOL]->node[iPoint]->GetSolution(0);
//				double Displ_y = solution_container[FEA_SOL]->node[iPoint]->GetSolution(1);
//				coord[0] += Displ_x;
//				coord[1] += Displ_y;
//
//				double Source_x = solution_container[FEA_SOL]->node[iPoint]->GetResSour()[2];
//				double Source_y = solution_container[FEA_SOL]->node[iPoint]->GetResSour()[3];
//
//				double Strain_xx = solution_container[FEA_SOL]->node[iPoint]->GetGradient(0,0);
//				double Strain_yy = solution_container[FEA_SOL]->node[iPoint]->GetGradient(1,1);
//				double Strain_xy = 0.5*(solution_container[FEA_SOL]->node[iPoint]->GetGradient(0,1) +
//						solution_container[FEA_SOL]->node[iPoint]->GetGradient(1,0));
//				double Strain_Trace = Strain_xx + Strain_yy;
//
//				double Displ_z = 0.0, Strain_zz = 0.0, Strain_xz = 0.0, Strain_yz = 0.0, Source_z = 0.0;
//				if (geometry->GetnDim() == 3) {
//					Displ_z = solution_container[FEA_SOL]->node[iPoint]->GetSolution(2);
//					coord[2] += Displ_z;
//
//					Source_x = solution_container[FEA_SOL]->node[iPoint]->GetResSour()[3];
//					Source_y = solution_container[FEA_SOL]->node[iPoint]->GetResSour()[4];
//					Source_z = solution_container[FEA_SOL]->node[iPoint]->GetResSour()[5];
//
//					Strain_zz = solution_container[FEA_SOL]->node[iPoint]->GetGradient(2,2);
//					Strain_xz = 0.5*(solution_container[FEA_SOL]->node[iPoint]->GetGradient(0,2) +
//							solution_container[FEA_SOL]->node[iPoint]->GetGradient(2,0));
//					Strain_yz = 0.5*(solution_container[FEA_SOL]->node[iPoint]->GetGradient(1,2) +
//							solution_container[FEA_SOL]->node[iPoint]->GetGradient(2,1));
//					Strain_Trace += Strain_zz;
//				}
//
//				double E = config->GetElasticyMod();
//				double Nu = config->GetPoissonRatio();
//				double Mu = E / (2.0*(1.0 + Nu));
//				double Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));		// For plane strain and 3-D
//
//				double Stress_xx = 2.0*Mu*Strain_xx + Lambda*Strain_Trace;
//				double Stress_yy = 2.0*Mu*Strain_yy + Lambda*Strain_Trace;
//				double Stress_xy = 2.0*Mu*Strain_xy;
//
//				double Pressure = 1.0; //solution_container[FEA_SOL]->node[iPoint]->GetPressure(); 
//
//				double Stress_zz = 0.0, Stress_xz = 0.0, Stress_yz = 0.0;
//				if (geometry->GetnDim() == 2) {
//					ConsVar_file << " " << Displ_x << " " << Displ_y << " " << Pressure << " " << Stress_xx << " " << Stress_yy << " " << Stress_xy <<" "<< Source_x<<" "<< Source_y;
//				}
//				else {
//					Stress_zz = 2.0*Mu*Strain_zz + Lambda*Strain_Trace;
//					Stress_xz = 2.0*Mu*Strain_xz;
//					Stress_yz = 2.0*Mu*Strain_yz;
//					ConsVar_file << " " << Displ_x << " " << Displ_y << " " << Displ_z << " " << Pressure << " " << Stress_xx << " " << Stress_yy << " " << Stress_zz << " " << Stress_xy << " " << Stress_xz << " " << Stress_yz<<" "<< Source_x<<" "<< Source_y<<" "<< Source_z;
//				}
//
//			}
//			/*--- Write the variables from electric potential---*/
//			if (solver == ELECTRIC_POTENTIAL) {
//				unsigned short nDim = geometry->GetnDim();
//				ConsVar_file << " " << solution_container[ELEC_SOL]->node[iPoint]->GetSolution(0);
//				for (iDim = 0; iDim < nDim; iDim++)
//					ConsVar_file << " " << -1.0*solution_container[ELEC_SOL]->node[iPoint]->GetGradient(0,iDim);
//			}
//
//			/*--- Write the multi-species flow variables ---*/
//			if ((solver == PLASMA_EULER ) || (solver == PLASMA_NAVIER_STOKES)) {
//				unsigned short nSpecies, nDiatomics, nDim;
//				double velocity2, SoundSpeed, MachNumber, Temperature_tr, Pressure;
//				nSpecies = solution_container[PLASMA_SOL]->GetnSpecies();
//				nDiatomics = solution_container[PLASMA_SOL]->GetnDiatomics();
//				nDim = geometry->GetnDim();
//				if ( (config->GetKind_GasModel() == ARGON) && (nDim ==3)) {
//					ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetMagneticField()[0];
//					ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetMagneticField()[1];
//					ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetMagneticField()[2];
//				}
//
//				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//					if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
//					else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
//
//					/*--- Write the Density ---*/
//					//ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
//          ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, nDim+3);
//
//					/*--- Write the Velocity and the Mach number ---*/
//					velocity2 = 0.0;
//					for (iDim = 0; iDim < nDim; iDim++) {
//						//ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+iDim)/solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
//            ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, iDim+1);
//
//						//velocity2 += pow(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+iDim)/solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0),2.0);
//            velocity2 += solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, iDim+1)*solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, iDim+1);
//					}
//					SoundSpeed = solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, nDim+5);
//					MachNumber = sqrt(velocity2)/SoundSpeed;
//					ConsVar_file << " " << MachNumber;
//
//					/*--- Write the Pressure ---*/
//					//Pressure = solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies);
//          Pressure = solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, nDim+2);
//					ConsVar_file << " " << Pressure;
//
//					/*--- Write the Translational Temperature ---*/
//					//Temperature_tr = solution_container[PLASMA_SOL]->node[iPoint]->GetTemperature_tr(iSpecies);
//          Temperature_tr = solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, 0);
//					ConsVar_file << " " << Temperature_tr;
//
//					if (iSpecies < nDiatomics) {
//						//ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetTemperature_vib(iSpecies);
//            ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetPrimVar(iSpecies, nDim+1);
//          }
//
//					if (config->GetKind_GasModel() != ARGON) {
//						/*--- Source terms ---*/
//						ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetResSour()[loc+0] / Volume;
//						for (iDim = 0; iDim < nDim; iDim++)
//							ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetResSour()[loc+1+iDim]/solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0) / Volume;
//						ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetResSour()[loc+1+nDim] / Volume;
//						if (iSpecies < nDiatomics)
//							ConsVar_file << " " << solution_container[PLASMA_SOL]->node[iPoint]->GetResSour()[loc+2+nDim] / Volume;
//					}
//				}
//
//				/*--- Write mean solution file ---*/
//				if (config->GetWrite_Mean_Solution()) {
//					double *Mass_Ave_Solution;
//					double Mixture_Density, Mass_Fraction;
//
//					Mass_Ave_Solution = new double [nDim+2];
//
//					for (iVar = 0; iVar < nDim+2; iVar++) {
//						Mass_Ave_Solution[iVar] = 0.0;
//					}
//
//					/*--- Pre-process by mass-averaging ---*/
//					Mixture_Density = 0.0;
//					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//						if (iSpecies < nDiatomics) loc = (nDim+3)*iSpecies;
//						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
//						Mixture_Density += solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
//					}
//					for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//						if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
//						else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
//						Mass_Fraction = solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0) / Mixture_Density;
//						Mass_Ave_Solution[0] += solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0);
//						for (iDim = 0; iDim < nDim; iDim++) {
//							Mass_Ave_Solution[iDim+1] += solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+iDim)*Mass_Fraction;
//						}
//						Mass_Ave_Solution[nDim+1] += solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+1+nDim)*Mass_Fraction;
//					}
//
//					ConsVar_file_mean << " " << Mass_Ave_Solution[0];
//					for (iDim = 0; iDim < nDim; iDim++)
//						ConsVar_file_mean << " " << Mass_Ave_Solution[iDim+1] / Mass_Ave_Solution[0];
//					ConsVar_file_mean << " " << Mass_Ave_Solution[nDim+1];
//
//					delete [] Mass_Ave_Solution;
//				}
//			}
//
//			/*--- Write the rotating frame variables ---*/
//			if (rotating_frame) {
//				ConsVar_file <<" "<< rel_velocity[0] <<" "<< rel_velocity[1] <<" "<< rel_velocity[2]<<" "<< rel_mach_number;
//			}
//
//			/*--- End line ---*/
//			ConsVar_file << endl;
//			if (config->GetWrite_Mean_Solution())
//				ConsVar_file_mean << endl;
//
//		}
//
//		/*--- Write the grid connectivity. ---*/
//		for(iElem = 0; iElem < geometry->GetnElem(); iElem++) {
//			if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
//				ConsVar_file <<
//						geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 << endl;
//				if (config->GetWrite_Mean_Solution()) {
//					ConsVar_file_mean <<
//							geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 << endl;
//				}
//			}
//			if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) {
//				ConsVar_file <<
//						geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
//				if (config->GetWrite_Mean_Solution()) {
//					ConsVar_file_mean <<
//							geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
//				}
//			}
//			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
//				ConsVar_file <<
//						geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
//				if (config->GetWrite_Mean_Solution()) {
//					ConsVar_file_mean <<
//							geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 << endl;
//				}
//			}
//			if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
//				ConsVar_file <<
//						geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(6)+1 <<" "<< geometry->elem[iElem]->GetNode(7)+1 << endl;
//				if (config->GetWrite_Mean_Solution()) {
//					ConsVar_file_mean <<
//							geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(6)+1 <<" "<< geometry->elem[iElem]->GetNode(7)+1 << endl;
//				}
//			}
//			if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID) {
//				ConsVar_file <<
//						geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 << endl;
//				if (config->GetWrite_Mean_Solution()) {
//					ConsVar_file_mean <<
//							geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(2)+1 <<" "<< geometry->elem[iElem]->GetNode(3)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 << endl;
//				}
//			}
//			if (geometry->elem[iElem]->GetVTK_Type() == WEDGE) {
//				ConsVar_file <<
//						geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(1)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
//						geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 << endl;
//				if (config->GetWrite_Mean_Solution()) {
//					ConsVar_file_mean <<
//							geometry->elem[iElem]->GetNode(0)+1 <<" "<< geometry->elem[iElem]->GetNode(1)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(1)+1 <<" "<< geometry->elem[iElem]->GetNode(2)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(3)+1 <<" "<< geometry->elem[iElem]->GetNode(4)+1 <<" "<<
//							geometry->elem[iElem]->GetNode(4)+1 <<" "<< geometry->elem[iElem]->GetNode(5)+1 << endl;
//				}
//			}
//		}
//		ConsVar_file.close();
//		if (config->GetWrite_Mean_Solution())
//			ConsVar_file_mean.close();
//	}
//}

void COutput::SetSurface_Flow(CConfig *config, CGeometry *geometry, CSolution *FlowSolution, unsigned long iExtIter) {
	unsigned long iPoint, iSpecies, iVertex;
	unsigned short iMarker, iDim, Species_0, Species_1, Species_2;
	double PressCoeff = 0.0, SkinFrictionCoeff, HeatTransferCoeff = 0.0,WallTemperature, WallTemperature_vib, WallNumDensity, WallPressure, WallHeatFlux, Mach, *aux_press = NULL, *aux_friction = NULL, *aux_rel_mach = NULL,
	*aux_heat_transfer = NULL, *aux_y_plus = NULL;
	double *aux_wall_temperature_Ar = NULL, *aux_wall_numdensity_Ar  = NULL, *aux_wall_temperature  = NULL;
	double *aux_wall_temperature_Ion = NULL, *aux_wall_numdensity_Ion = NULL, *aux_wall_numdensity  = NULL;
	double *aux_wall_pressure_Ar = NULL, *aux_wall_pressure_Ion = NULL, *aux_wall_pressure_Elec  = NULL;
	double *aux_wall_temperature_Elec = NULL, *aux_wall_numdensity_Elec = NULL, *aux_wall_numdensity_diff = NULL;
  double **aux_temp_tr, **aux_temp_vib, **aux_numdensity, **aux_heatflux, **aux_partpress;
	char cstr[200], buffer [50];
	ofstream SurfFlow_file;

	unsigned short solver = config->GetKind_Solver();
	bool incompressible = config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool adiabatic = config->GetAdiabaticWall();
	if (config->GetOutput_FileFormat() == PARAVIEW) {

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

		case NAVIER_STOKES: case RANS: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
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
		case PLASMA_NAVIER_STOKES:	case PLASMA_EULER:
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
						aux_wall_temperature[iPoint] = FlowSolution->node[iPoint]->GetTemperature_tr(0);
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
						aux_wall_temperature[iPoint] = FlowSolution->node[iPoint]->GetTemperature_tr(1);
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
						aux_wall_temperature[iPoint] = FlowSolution->node[iPoint]->GetTemperature_tr(2);
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

		if ((solver == EULER) || (solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) ||
            (solver == NAVIER_STOKES) || (solver == RANS) || (solver == FREE_SURFACE_RANS) ) {
			SurfFlow_file << ", \"Pressure Coefficient\"";
    }
    
		if ((solver == EULER) || (solver == FREE_SURFACE_EULER))
			SurfFlow_file << ", \"Mach Number\"";

		if ((solver == NAVIER_STOKES) || (solver == RANS)) {
			if (incompressible) SurfFlow_file << ", \"Skin Friction Coefficient\", \"y+\"";
            else SurfFlow_file << ", \"Skin Friction Coefficient\", \"Heat Transfer Coefficient\", \"y+\"";
		}
		if (rotating_frame && (solver == EULER))
			SurfFlow_file << ", \"Relative Mach\"";

		if (solver == PLASMA_NAVIER_STOKES ) {
      switch (config->GetKind_GasModel()) {
      
      case ARGON:
        SurfFlow_file <<  ", \"Skin Friction Coefficient\"";
        SurfFlow_file <<  ", \"Ar Temperature\", \"Ar NumDen\", \"Ion Temperature\", \"Ion NumDen \"";
        SurfFlow_file <<  ", \"Ele Temperature\", \"Ele NumDen\", \"NumofNegCharge\"";
        SurfFlow_file <<  ", \"Ar Pressure\", \"Ion Pressure\", \"Electron Pressure \"";
        if (!adiabatic) SurfFlow_file <<  ", \"Heat_Transfer_Coefficient\"";
        break;
        
      case N2: case O2: case AIR5: case AIR7:
        unsigned short iSpecies;
        for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
          SurfFlow_file << ", \"Temperature_tr(" << iSpecies << ")\"";
          if (iSpecies < config->GetnDiatomics())
            SurfFlow_file << ", \"Temperature_vib(" << iSpecies << ")\"";
          SurfFlow_file << ", \"Pressure(" << iSpecies << ")\"";
          SurfFlow_file << ", \"NumDensity(" << iSpecies << ")\"";
          if (!adiabatic) SurfFlow_file <<  ", \"HeatFlux(" << iSpecies << ")\"";
        }
        break;
      }
      
		}
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
		aux_press 	 	  = new double [geometry->GetnPoint()];
		aux_friction 	  = new double [geometry->GetnPoint()];
		aux_heat_transfer = new double [geometry->GetnPoint()];
		aux_y_plus = new double [geometry->GetnPoint()];

		if (rotating_frame && (solver == EULER)) aux_rel_mach = new double [geometry->GetnPoint()];

		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) aux_press[iPoint] = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					aux_press[iPoint] = FlowSolution->GetCPressure(iMarker,iVertex);
				}

		if ((solver == NAVIER_STOKES) || (solver == RANS)) {
			for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				aux_friction[iPoint] = 0.0;
				aux_heat_transfer[iPoint] = 0.0;
				aux_y_plus[iPoint] = 0.0;
			}
			for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
				if (config->GetMarker_All_Plotting(iMarker) == YES) {
					for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
						aux_friction[iPoint] = FlowSolution->GetCSkinFriction(iMarker,iVertex);
						aux_heat_transfer[iPoint] = FlowSolution->GetHeatTransferCoeff(iMarker,iVertex);
						aux_y_plus[iPoint] = FlowSolution->GetYPlus(iMarker,iVertex);
					}
				}
		}

		if (solver == PLASMA_NAVIER_STOKES) {
      unsigned short iSpecies;
      unsigned short nSpecies = FlowSolution->GetnSpecies();
      unsigned short nDim = geometry->GetnDim();
      
      switch (config->GetKind_GasModel()) {
        case ARGON:
          aux_wall_temperature_Ar 	 = new double [geometry->GetnPoint()];
          aux_wall_temperature_Ion 	 = new double [geometry->GetnPoint()];
          aux_wall_temperature_Elec 	 = new double [geometry->GetnPoint()];
          aux_wall_numdensity_Ar 	 	 = new double [geometry->GetnPoint()];
          aux_wall_numdensity_Ion 	 = new double [geometry->GetnPoint()];
          aux_wall_numdensity_Elec 	 = new double [geometry->GetnPoint()];
          aux_wall_numdensity_diff 	 = new double [geometry->GetnPoint()];
          aux_wall_pressure_Ar 		 = new double [geometry->GetnPoint()];
          aux_wall_pressure_Ion		 = new double [geometry->GetnPoint()];
          aux_wall_pressure_Elec		 = new double [geometry->GetnPoint()];
          if (!adiabatic) aux_heat_transfer 		 	 = new double [geometry->GetnPoint()];
          
          double Mass_0, Mass_1, Mass_2;
          Mass_0 = config->GetParticle_Mass(0);
          Mass_1 = config->GetParticle_Mass(1);
          Mass_2 = config->GetParticle_Mass(2);
          Species_0 = 0; Species_1 = 1; Species_2 = 2;
          for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
            aux_friction[iPoint] 			  = 0.0;
            aux_wall_temperature_Ar[iPoint]   = 0.0;
            aux_wall_numdensity_Ar[iPoint] 	  = 0.0;
            aux_wall_temperature_Ion[iPoint]  = 0.0;
            aux_wall_numdensity_Ion[iPoint]   = 0.0;
            aux_wall_temperature_Elec[iPoint] = 0.0;
            aux_wall_numdensity_Elec[iPoint]  = 0.0;
            aux_wall_numdensity_diff[iPoint]  = 0.0;
            aux_wall_pressure_Ar[iPoint] 	  = 0.0;
            aux_wall_pressure_Ion[iPoint]     = 0.0;
            aux_wall_pressure_Elec[iPoint]    = 0.0;
            if (!adiabatic) 			aux_heat_transfer[iPoint] 		  = 0.0;
          }
          for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Plotting(iMarker) == YES)
              for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                aux_friction[iPoint] 				= 	FlowSolution->GetCSkinFriction(iMarker,iVertex);
                aux_wall_temperature_Ar[iPoint] 	= 	FlowSolution->node[iPoint]->GetTemperature_tr(0);
                aux_wall_numdensity_Ar[iPoint] 		= 	FlowSolution->node[iPoint]->GetDensity(0)/Mass_0;
                aux_wall_temperature_Ion[iPoint] 	= 	FlowSolution->node[iPoint]->GetTemperature_tr(1);
                aux_wall_numdensity_Ion[iPoint] 	= 	FlowSolution->node[iPoint]->GetDensity(1)/Mass_1;
                aux_wall_temperature_Elec[iPoint] 	= 	FlowSolution->node[iPoint]->GetTemperature_tr(2);
                aux_wall_numdensity_Elec[iPoint] 	= 	FlowSolution->node[iPoint]->GetDensity(2)/Mass_2;
                aux_wall_numdensity_diff[iPoint] 	= 	FlowSolution->node[iPoint]->GetDensity(2)/Mass_2 - FlowSolution->node[iPoint]->GetDensity(1)/Mass_1;
                aux_wall_pressure_Ar[iPoint] 		=   FlowSolution->node[iPoint]->GetPressure(Species_0);
                aux_wall_pressure_Ion[iPoint] 		=   FlowSolution->node[iPoint]->GetPressure(Species_1);
                aux_wall_pressure_Elec[iPoint] 		=   FlowSolution->node[iPoint]->GetPressure(Species_2);
                if (!adiabatic) aux_heat_transfer[iPoint] 			= 	FlowSolution->GetHeatTransferCoeff(iMarker,iVertex);
              }
          break;
          
        case N2: case O2: case AIR5: case AIR7:

          /*--- Allocate arrays ---*/
          aux_temp_tr    = new double*[nSpecies];
          aux_temp_vib   = new double*[nSpecies];
          aux_numdensity = new double*[nSpecies];
          aux_heatflux   = new double*[nSpecies];
          aux_partpress  = new double*[nSpecies];
          for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
            aux_temp_tr[iSpecies]    = new double[geometry->GetnPoint()];
            aux_temp_vib[iSpecies]   = new double[geometry->GetnPoint()];
            aux_numdensity[iSpecies] = new double[geometry->GetnPoint()];
            aux_heatflux[iSpecies]   = new double[geometry->GetnPoint()];
            aux_partpress[iSpecies]  = new double[geometry->GetnPoint()];
          }
          
          /*--- Initialize arrays ---*/
          for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
              aux_temp_tr[iSpecies][iPoint]    = 0.0;
              aux_temp_vib[iSpecies][iPoint]   = 0.0;
              aux_partpress[iSpecies][iPoint]  = 0.0;
              aux_numdensity[iSpecies][iPoint] = 0.0;
              aux_heatflux[iSpecies][iPoint]   = 0.0;
            }
          }
          for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
            if (config->GetMarker_All_Plotting(iMarker) == YES) {
              for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                  aux_temp_tr[iSpecies][iPoint]    = FlowSolution->node[iPoint]->GetPrimVar(iSpecies, 0);
                  aux_temp_vib[iSpecies][iPoint]   = FlowSolution->node[iPoint]->GetPrimVar(iSpecies, nDim+1);
                  aux_partpress[iSpecies][iPoint]  = FlowSolution->node[iPoint]->GetPrimVar(iSpecies, nDim+2);
                  aux_numdensity[iSpecies][iPoint] = FlowSolution->node[iPoint]->GetPrimVar(iSpecies, nDim+3)
                                                     / (config->GetMolar_Mass(iSpecies)) * AVOGAD_CONSTANT;
                  aux_heatflux[iSpecies][iPoint]   = FlowSolution->GetHeatTransferCoeff(iMarker, iVertex);
                }
              }
            }
          }
          
          break;
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
							rel_velocity[iDim] = FlowSolution->node[iPoint]->GetVelocity(iDim, config->GetIncompressible()) - rot_vel[iDim];
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
			if ((solver == EULER) || (solver == NAVIER_STOKES) || (solver == RANS) || 
					(solver == FREE_SURFACE_EULER) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << aux_press[iPoint];
			}

			/*--- Write the Mach number ---*/
			if ((solver == EULER || (solver == FREE_SURFACE_EULER))) {
				if (incompressible) {
					double Density = FlowSolution->node[iPoint]->GetDensityInc()*config->GetDensity_Ref();
					double SoundSpeed = sqrt(config->GetBulk_Modulus()/Density);					
					Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/SoundSpeed;
				}
				else Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2())/FlowSolution->node[iPoint]->GetSoundSpeed();
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << Mach;
			}

			/*--- Write the skin friction coefficient variables ---*/
			if ((solver == NAVIER_STOKES) || (solver == RANS) || (solver == FREE_SURFACE_NAVIER_STOKES) || (solver == FREE_SURFACE_RANS)) {
				if (geometry->node[iPoint]->GetBoundary()) {
                    if (incompressible) SurfFlow_file << " " << aux_friction[iPoint] << " " << aux_y_plus[iPoint];
                    else SurfFlow_file << " " << aux_friction[iPoint] << " " << aux_heat_transfer[iPoint] << " " << aux_y_plus[iPoint];
                }
			}

			if (solver == PLASMA_NAVIER_STOKES) {
        unsigned short iSpecies;
        unsigned short nSpecies = FlowSolution->GetnSpecies();
        
        switch (config->GetKind_GasModel()) {
          case ARGON:
            if (geometry->node[iPoint]->GetBoundary()) {
              SkinFrictionCoeff 	= aux_friction[iPoint]; 				SurfFlow_file << " " << SkinFrictionCoeff;
              WallTemperature 	= aux_wall_temperature_Ar[iPoint]; 		SurfFlow_file << " " << WallTemperature;
              WallNumDensity 		= aux_wall_numdensity_Ar[iPoint]; 		SurfFlow_file << " " << WallNumDensity;
              WallTemperature 	= aux_wall_temperature_Ion[iPoint]; 	SurfFlow_file << " " << WallTemperature;
              WallNumDensity 		= aux_wall_numdensity_Ion[iPoint]; 		SurfFlow_file << " " << WallNumDensity;
              WallTemperature 	= aux_wall_temperature_Elec[iPoint]; 	SurfFlow_file << " " << WallTemperature;
              WallNumDensity 		= aux_wall_numdensity_Elec[iPoint]; 	SurfFlow_file << " " << WallNumDensity;
              WallNumDensity 		= aux_wall_numdensity_diff[iPoint]; 	SurfFlow_file << " " << WallNumDensity;
              WallPressure 		= aux_wall_pressure_Ar[iPoint]; 		SurfFlow_file << " " << WallPressure;
              WallPressure 		= aux_wall_pressure_Ion[iPoint]; 		SurfFlow_file << " " << WallPressure;
              WallPressure 		= aux_wall_pressure_Elec[iPoint]; 		SurfFlow_file << " " << WallPressure;
              if (!adiabatic) {
                HeatTransferCoeff 	= aux_heat_transfer[iPoint];
                SurfFlow_file << " " << HeatTransferCoeff;
              }
            }
            break;
            
          case N2: case O2: case AIR5: case AIR7:
            if (geometry->node[iPoint]->GetBoundary()) {
              for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                WallTemperature     = aux_temp_tr[iSpecies][iPoint];      SurfFlow_file << " " << WallTemperature;
                WallTemperature_vib = aux_temp_vib[iSpecies][iPoint];     SurfFlow_file << " " << WallTemperature_vib;
                WallPressure        = aux_partpress[iSpecies][iPoint];    SurfFlow_file << " " << WallPressure;
                WallNumDensity      = aux_numdensity[iSpecies][iPoint];   SurfFlow_file << " " << WallNumDensity;
                WallHeatFlux        = aux_heatflux[iSpecies][iPoint];     SurfFlow_file << " " << WallHeatFlux;
              }
            }
            break;
        }

			}

			/*--- Write the relative mach number ---*/
			if (rotating_frame && (solver == EULER)) {
				if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << " " << aux_rel_mach[iPoint];
			}

			/*--- End line ---*/
			if (geometry->node[iPoint]->GetBoundary()) SurfFlow_file << endl;
		}

		/*--- Write the cells using the new numbering ---*/
		unsigned long iElem;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
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
		delete [] aux_heat_transfer;
		delete [] aux_y_plus;
		if (rotating_frame && (solver == EULER)) delete[] aux_rel_mach;
		delete [] PointSurface;
		SurfFlow_file.close();

		if(solver == PLASMA_NAVIER_STOKES) {
      switch(config->GetKind_GasModel()) {
        case ARGON:
          delete[]  aux_wall_temperature_Ar;
          delete[]  aux_wall_temperature_Ion;
          delete[]  aux_wall_temperature_Elec;
          delete[]  aux_wall_numdensity_Ar;
          delete[]  aux_wall_numdensity_Ion;
          delete[]  aux_wall_numdensity_Elec;
          delete[]  aux_wall_numdensity_diff;
          delete[]  aux_wall_pressure_Ar;
          delete[]  aux_wall_pressure_Ion;
          delete[]  aux_wall_pressure_Elec;
          break;
          
        case N2: case O2: case AIR5: case AIR7:
          for (iSpecies=0; iSpecies < FlowSolution->GetnSpecies(); iSpecies++) {
            delete [] aux_temp_tr[iSpecies];
            delete [] aux_temp_vib[iSpecies];
            delete [] aux_partpress[iSpecies];
            delete [] aux_numdensity[iSpecies];
            delete [] aux_heatflux[iSpecies];
          }
          delete [] aux_temp_tr;
          delete [] aux_temp_vib;
          delete [] aux_partpress;
          delete [] aux_numdensity;
          delete [] aux_heatflux;
      }

		}

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

		SurfAdj_file << "ZONE ";
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
			SurfAdj_file << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*int(iExtIter+1)<<", ";

		SurfAdj_file << "NODES="<< nPointSurface <<" , ELEMENTS="<< nElemSurface <<", DATAPACKING=POINT";

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
	unsigned long iPoint, iVertex, Global_Index;
	unsigned short iMarker;
	double PressCoeff = 0.0, SkinFrictionCoeff, xCoord, yCoord, zCoord, Mach;
	char cstr[200], buffer [50];
	ofstream SurfFlow_file;
	unsigned short solver = config->GetKind_Solver();

	/*--- Write file name with extension if unsteady ---*/
	strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
	if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
		if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
	}
	else
		sprintf (buffer, ".csv");

	strcat (cstr, buffer);
	SurfFlow_file.precision(15);
	SurfFlow_file.open(cstr, ios::out);

	if (geometry->GetnDim() == 2) {
		switch (solver) {
		case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"Global_Index\"" << endl; break;
		case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"Global_Index\"" << endl; break;
		}
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					PressCoeff = FlowSolution->GetCPressure(iMarker,iVertex);
					Global_Index = geometry->node[iPoint]->GetGlobalIndex();
					switch (solver) {
					case EULER :
						Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2()) / FlowSolution->node[iPoint]->GetSoundSpeed();
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << Mach << "," << yCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS:
						SkinFrictionCoeff = FlowSolution->GetCSkinFriction(iMarker,iVertex);
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << SkinFrictionCoeff << "," << yCoord << ", " << Global_Index << endl;
						break;
					}
				}
	}

	if (geometry->GetnDim() == 3) {
		switch (solver) {
		case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
		case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
		}
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Plotting(iMarker) == YES)
				for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					xCoord = geometry->node[iPoint]->GetCoord(0);
					yCoord = geometry->node[iPoint]->GetCoord(1);
					zCoord = geometry->node[iPoint]->GetCoord(2);
					PressCoeff = FlowSolution->GetCPressure(iMarker,iVertex);
					Global_Index = geometry->node[iPoint]->GetGlobalIndex();
					switch (solver) {
					case EULER :
						Mach = sqrt(FlowSolution->node[iPoint]->GetVelocity2()) / FlowSolution->node[iPoint]->GetSoundSpeed();
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << Mach <<"," << yCoord << "," << zCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS:
						SkinFrictionCoeff = FlowSolution->GetCSkinFriction(iMarker,iVertex);
						SurfFlow_file << scientific << xCoord << "," << PressCoeff << "," << SkinFrictionCoeff << "," << yCoord << "," << zCoord << ", " << Global_Index << endl;
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
			MaxLocalVertex_Surface = 0, nBuffer_Scalar, *Buffer_Receive_nVertex = NULL, position, Global_Index;
	ofstream SurfFlow_file;

	/*--- Write the surface .csv file, the information of all the vertices is 
   sent to the MASTER_NODE for writing ---*/
	nLocalVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
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
	unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];

	nVertex_Surface = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Buffer_Send_Press[nVertex_Surface] = FlowSolution->GetCPressure(iMarker,iVertex);
					Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
					Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
					Buffer_Send_GlobalIndex[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
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
	unsigned long *Buffer_Receive_GlobalIndex = NULL;

	if (rank == MASTER_NODE) {
		Buffer_Receive_Coord_x = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Coord_y = new double [nProcessor*MaxLocalVertex_Surface];
		if (geometry->GetnDim() == 3) 
			Buffer_Receive_Coord_z = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Press = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_Mach = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_SkinFriction = new double [nProcessor*MaxLocalVertex_Surface];
		Buffer_Receive_GlobalIndex = new  unsigned long [nProcessor*MaxLocalVertex_Surface];
	}

	nBuffer_Scalar = MaxLocalVertex_Surface;

	/*--- Send the information to the Master node ---*/
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI::DOUBLE, 
			Buffer_Receive_Coord_x, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI::DOUBLE, 
			Buffer_Receive_Coord_y, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (geometry->GetnDim() == 3) 
		MPI::COMM_WORLD.Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI::DOUBLE,
				Buffer_Receive_Coord_z, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	MPI::COMM_WORLD.Gather(Buffer_Send_Press, nBuffer_Scalar, MPI::DOUBLE, 
			Buffer_Receive_Press, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if (config->GetKind_Solver() == EULER)
		MPI::COMM_WORLD.Gather(Buffer_Send_Mach, nBuffer_Scalar, MPI::DOUBLE, 
				Buffer_Receive_Mach, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);
	if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver() == RANS))
		MPI::COMM_WORLD.Gather(Buffer_Send_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, 
				Buffer_Receive_SkinFriction, nBuffer_Scalar, MPI::DOUBLE, MASTER_NODE);

	MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG, 
			Buffer_Receive_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG, MASTER_NODE);

	/*--- The master node writes the surface files ---*/
	if (rank == MASTER_NODE) {

		/*--- Write file name with extension if unsteady ---*/
		char buffer[50];
		string filename = config->GetSurfFlowCoeff_FileName();

		/*--- Remove the domain number from the surface csv filename ---*/
		if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());

		/*--- Write file name with extension if unsteady ---*/
		strcpy (cstr, filename.c_str());
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
		}
		else
			sprintf (buffer, ".csv");

		strcat (cstr, buffer);
		SurfFlow_file.precision(15);
		SurfFlow_file.open(cstr, ios::out);

		/*--- Write the 2D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 2) {
			switch (config->GetKind_Solver()) {
			case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"Global_Index\"" << endl; break;
			case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"Global_Index\"" << endl; break;
			}
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					PressCoeff = Buffer_Receive_Press[position];
					Global_Index = Buffer_Receive_GlobalIndex[position];
					switch (config->GetKind_Solver()) {
					case EULER :
						Mach = Buffer_Receive_Mach[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << Mach << ", " << yCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS:
						SkinFrictionCoeff = Buffer_Receive_SkinFriction[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << SkinFrictionCoeff << ", " << yCoord << ", " << Global_Index << endl;
						break;
					}
				}
		}

		/*--- Write the 3D surface flow coefficient file ---*/
		if (geometry->GetnDim() == 3) {
			switch (config->GetKind_Solver()) {
			case EULER : SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Mach_Number\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
			case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"x_coord\",\"Pressure_Coefficient\",\"Skin_Friction_Coefficient\",\"y_coord\",\"z_coord\",\"Global_Index\"" << endl; break;
			}
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
				for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
					position = iProcessor*MaxLocalVertex_Surface+iVertex;
					xCoord = Buffer_Receive_Coord_x[position];
					yCoord = Buffer_Receive_Coord_y[position];
					zCoord = Buffer_Receive_Coord_z[position];
					PressCoeff = Buffer_Receive_Press[position];
					Global_Index = Buffer_Receive_GlobalIndex[position];
					switch (config->GetKind_Solver()) {
					case EULER :
						Mach = Buffer_Receive_Mach[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << Mach << ", " << yCoord << ", " << zCoord << ", " << Global_Index << endl;
						break;
					case NAVIER_STOKES: case RANS:
						SkinFrictionCoeff = Buffer_Receive_SkinFriction[position];
						SurfFlow_file << scientific << xCoord << ", " << PressCoeff << ", " << SkinFrictionCoeff << ", " << yCoord << ", " << zCoord << ", " << Global_Index << endl;
						break;
					}
				}
		}
	}

	if (rank == MASTER_NODE) {
		delete [] Buffer_Receive_Coord_x;
		delete [] Buffer_Receive_Coord_y;
		if (geometry->GetnDim() == 3) 
			delete [] Buffer_Receive_Coord_z;
		delete [] Buffer_Receive_Press;
		delete [] Buffer_Receive_Mach;
		delete [] Buffer_Receive_SkinFriction;
		delete [] Buffer_Receive_GlobalIndex;
	}

	delete [] Buffer_Send_Coord_x;
	delete [] Buffer_Send_Coord_y;
	delete [] Buffer_Send_Coord_z;
	delete [] Buffer_Send_Press;
	delete [] Buffer_Send_Mach;
	delete [] Buffer_Send_SkinFriction;
	delete [] Buffer_Send_GlobalIndex;

	SurfFlow_file.close();

#endif
}

void COutput::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolution *AdjSolution, CSolution *FlowSolution, unsigned long iExtIter, unsigned short val_iZone) {

#ifdef NO_MPI

	unsigned long iPoint, iVertex;
	double *Solution, xCoord, yCoord, zCoord, *IntBoundary_Jump;
	unsigned short iMarker;
	char cstr[200], buffer[50];
	ofstream SurfAdj_file;

	/*--- Write file name with extension if unsteady ---*/
	strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());

	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));

	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
		if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
		if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
		if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
		if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
	}
	else
		sprintf (buffer, ".csv");

	strcat(cstr, buffer);
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
		char cstr[200], buffer[50];
		ofstream SurfAdj_file;
		string filename = config->GetSurfAdjCoeff_FileName();

		/*--- Remove the domain number from the surface csv filename ---*/
		if (nProcessor > 1) filename.erase (filename.end()-2, filename.end());

		/*--- Write file name with extension if unsteady ---*/
		strcpy (cstr, filename.c_str());

		if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
			if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.csv", int(val_iZone));
			if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.csv", int(val_iZone));
			if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.csv", int(val_iZone));
			if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.csv", int(val_iZone));
			if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.csv", int(val_iZone));

		} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
		}
		else
			sprintf (buffer, ".csv");

		strcat (cstr, buffer);
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

void COutput::MergeGeometry(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  /*--- Merge coordinates of all grid nodes (excluding ghost points). ---*/
  
  MergeCoordinates(config, geometry);
  
  /*--- Merge connectivity for each type of element (excluding halos). ---*/
  
  MergeConnectivity(config, geometry, TRIANGLE    );
  MergeConnectivity(config, geometry, RECTANGLE   );
  MergeConnectivity(config, geometry, TETRAHEDRON );
  MergeConnectivity(config, geometry, HEXAHEDRON  );
  MergeConnectivity(config, geometry, WEDGE       );
  MergeConnectivity(config, geometry, PYRAMID     );
  
  /*--- Update total number of elements after merge. ---*/
  
  nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
  nGlobal_Hexa + nGlobal_Pyra + nGlobal_Wedg;
  
}

void COutput::MergeCoordinates(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables needed on all processors ---*/
  
	unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint;
    
#ifdef NO_MPI
  
	/*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/
  
  /*--- Total number of points in the mesh (excluding halos). ---*/
  nGlobal_Poin = geometry->GetnPoint(); //geometry->GetnPointDomain(); (F.P.)
  nGlobal_Doma = geometry->GetnPointDomain();

  /*--- Allocate the coordinates data structure. ---*/
  
  Coords = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new double[nGlobal_Poin];
  }
  
  /*--- Loop over the mesh to collect the coords of the local points. ---*/

  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//    if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
      
      /*--- Retrieve the current coordinates at this node. ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][jPoint] = geometry->node[iPoint]->GetCoord(iDim);
      }
      
      /*--- Increment a counter since we may be skipping over 
       some halo nodes during this loop. ---*/
      jPoint++;
//    }  (F.P.)
  }
  
#else
  
	/*--- MPI preprocessing ---*/
  
  int iProcessor;
	int nProcessor = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  
	/*--- Local variables needed for merging the geometry with MPI. ---*/
  
	unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
	unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
	unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];
  
  /*--- Sum total number of nodes that belong to the domain ---*/
  
	Buffer_Send_nPoin[0] = geometry->GetnPointDomain();
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoin, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoin, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  if (rank == MASTER_NODE) {
		nGlobal_Doma = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
		}
	}
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  
	nLocalPoint = geometry->GetnPoint(); //geometry->GetnPointDomain();  (F.P.)
	Buffer_Send_nPoin[0] = nLocalPoint;
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint,
                            1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoin, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoin, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	nBuffer_Scalar = MaxLocalPoint;
  
	/*--- Send and Recv buffers. ---*/
  
	double *Buffer_Send_X = new double[MaxLocalPoint];
	double *Buffer_Recv_X = NULL;
  
	double *Buffer_Send_Y = new double[MaxLocalPoint];
	double *Buffer_Recv_Y = NULL;
  
  double *Buffer_Send_Z, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new double[MaxLocalPoint];
  
	unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
	unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
	/*--- Prepare the receive buffers in the master node only. ---*/
  
	if (rank == MASTER_NODE) {
    
		Buffer_Recv_X = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_Y = new double[nProcessor*MaxLocalPoint];
		if (nDim == 3) Buffer_Recv_Z = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
		/*--- Sum total number of nodes to be written and allocate arrays ---*/
		nGlobal_Poin = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
		}
		Coords = new double*[nDim];
		for (iDim = 0; iDim < nDim; iDim++) {
			Coords[iDim] = new double[nGlobal_Poin];
		}
	}
  
	/*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  /*--- Loop over this partition to collect the coords of the local points.
   Note that we are NOT including the halo nodes here. ---*/
  double *Coords_Local; jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//    if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
      
      /*--- Retrieve local coordinates at this node. ---*/
      Coords_Local = geometry->node[iPoint]->GetCoord();
      
      /*--- Load local coords into the temporary send buffer. ---*/
      Buffer_Send_X[jPoint] = Coords_Local[0];
      Buffer_Send_Y[jPoint] = Coords_Local[1];
      if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];
      
      /*--- Store the global index for this local node. ---*/
      Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Increment jPoint as the counter. We need this because iPoint
       may include halo nodes that we skip over during this loop. ---*/
      jPoint++;
//    }  (F.P.)
  }
  
  /*--- Gather the coordinate data on the master node using MPI. ---*/
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(Buffer_Send_X, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Recv_X, nBuffer_Scalar, MPI::DOUBLE,
                         MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Y, nBuffer_Scalar, MPI::DOUBLE,
                         Buffer_Recv_Y, nBuffer_Scalar, MPI::DOUBLE,
                         MASTER_NODE);
  if (nDim == 3) {
    MPI::COMM_WORLD.Gather(Buffer_Send_Z, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Z, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
  }
  MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         MASTER_NODE);
  
  /*--- The master node unpacks and sorts this variable by global index ---*/
  
  if (rank == MASTER_NODE) {
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        
        /*--- Get global index, then loop over each variable and store ---*/
        iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
        Coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
        Coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
        if (nDim == 3) Coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
        jPoint++;
      }
      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
      jPoint = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary data buffers. ---*/
  
	delete [] Buffer_Send_X;
	delete [] Buffer_Send_Y;
	if (nDim == 3) delete [] Buffer_Send_Z;
	delete [] Buffer_Send_GlobalIndex;
	if (rank == MASTER_NODE) {
		delete [] Buffer_Recv_X;
		delete [] Buffer_Recv_Y;
		if (nDim == 3)  delete [] Buffer_Recv_Z;
		delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }
  
#endif
  
}

void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables needed on all processors ---*/
  
	unsigned short NODES_PER_ELEMENT;
  
  unsigned long iPoint, iNode, jNode;
	unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  unsigned long *Conn_Elem;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nLocalElem = geometry->GetnElemTria();
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case RECTANGLE:
      nLocalElem = geometry->GetnElemQuad();
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      nLocalElem = geometry->GetnElemTetr();
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      nLocalElem = geometry->GetnElemHexa();
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case WEDGE:
      nLocalElem = geometry->GetnElemWedg();
      NODES_PER_ELEMENT = N_POINTS_WEDGE;
      break;
    case PYRAMID:
      nLocalElem = geometry->GetnElemPyra();
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(0); break;
  }
  
  /*--- Merge the connectivity in serial or parallel. ---*/
  
#ifdef NO_MPI
  
  /*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
  /*--- Allocate a temporary array for the connectivity ---*/
  Conn_Elem = new unsigned long[nLocalElem*NODES_PER_ELEMENT];
  
  /*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
  jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Check if this is a halo node. ---*/
      isHalo = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        if (!geometry->node[iPoint]->GetDomain())
          isHalo = true;
      }
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the temporary array. Do not merge any
       halo cells (periodic BC) ---*/
      if (!isHalo) {
        nElem_Total++;
        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
          Conn_Elem[jNode] = geometry->elem[iElem]->GetNode(iNode);
          
          /*--- Increment jNode as the counter. ---*/
          jNode++;
        }
      }
    }
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int iProcessor, jProcessor;
  int nProcessor = MPI::COMM_WORLD.Get_size();
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0, pElem = 0;
  unsigned long MaxLocalElem = 0;
  
  bool *Write_Elem;
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[nProcessor];
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&nLocalElem, &MaxLocalElem,
                            1, MPI::UNSIGNED_LONG, MPI::MAX);
  MPI::COMM_WORLD.Gather(&Buffer_Send_nElem, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nElem, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;
  
  int *Buffer_Send_Halo = new int[MaxLocalElem];
  int *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[nProcessor*nBuffer_Scalar];
    Buffer_Recv_Halo = new int[nProcessor*MaxLocalElem];
    Conn_Elem = new unsigned long[nProcessor*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  jNode = 0; jElem = 0;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the send buffer. ---*/
      Buffer_Send_Halo[jElem] = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        
        /*--- Store the global index values directly. ---*/
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Check if this is a halo node. If so, flag this element
         as a halo cell. We will use this later to sort and remove
         any duplicates from the connectivity list. ---*/
        if (!geometry->node[iPoint]->GetDomain())
          Buffer_Send_Halo[jElem] = true;
        
        /*--- Increment jNode as the counter. We need this because iElem
         may include other elements that we skip over during this loop. ---*/
        jNode++;
      }
      jElem++;
    }
  }
  
  /*--- Gather the element connectivity information. ---*/
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         Buffer_Recv_Elem, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                         MASTER_NODE);
  MPI::COMM_WORLD.Gather(Buffer_Send_Halo, MaxLocalElem, MPI::INT,
                         Buffer_Recv_Halo, MaxLocalElem, MPI::INT,
                         MASTER_NODE);
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    Write_Elem = new bool[nProcessor*MaxLocalElem];
    for (iElem = 0; iElem < nProcessor*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Loop for flagging duplicate elements so that they are not
     included in the final connectivity list. ---*/
    kElem = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Check if this element was originally marked as a halo. ---*/
        if (Buffer_Recv_Halo[kElem+iElem]) {
          
          /*--- Check all other elements flagged as halos on the
           remaining processors for duplicates (start at iProcessor+1). ---*/
          pElem = (iProcessor+1)*MaxLocalElem;
          for (jProcessor = iProcessor+1; jProcessor < nProcessor; jProcessor++) {
            for (jElem = 0; jElem < Buffer_Recv_nElem[jProcessor]; jElem++) {
              
              /*--- Check if this element was originally marked as a halo. ---*/
              if (Buffer_Recv_Halo[pElem+jElem]) {
                
                /*--- Check for a duplicate by comparing the index of each
                 node in the element. ---*/
                bool isDuplicate = true;
                for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
                  if (Buffer_Recv_Elem[kElem*NODES_PER_ELEMENT+iElem*NODES_PER_ELEMENT+iNode] !=
                      Buffer_Recv_Elem[pElem*NODES_PER_ELEMENT+jElem*NODES_PER_ELEMENT+iNode])
                    isDuplicate = false;
                }
                
                /*--- If we have found a duplicate element, set both the
                 original flag and "write" state booleans to false. In this
                 way, this element will not be found as we continue searching
                 and it will not be written to the connectivity list. ---*/
                if (isDuplicate) {
                  Buffer_Recv_Halo[pElem+jElem] = false;
                  Write_Elem[pElem+jElem] = false;
                }
              }
            }
            pElem = (jProcessor+1)*MaxLocalElem;
          }
        }
      }
      kElem = (iProcessor+1)*MaxLocalElem;
    }
    
    /*--- Store the unique connectivity list for this element type. ---*/
    
    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store ---*/
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode];
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
#endif
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case RECTANGLE:
        nGlobal_Quad = nElem_Total;
        if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
        break;
      case TETRAHEDRON:
        nGlobal_Tetr = nElem_Total;
        if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
        break;
      case HEXAHEDRON:
        nGlobal_Hexa = nElem_Total;
        if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
        break;
      case WEDGE:
        nGlobal_Wedg = nElem_Total;
        if (nGlobal_Wedg > 0) Conn_Wedg = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(0); break;
    }
  }
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolution **solution, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  
	unsigned short Kind_Solver  = config->GetKind_Solver();
	unsigned short iVar, jVar, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
	unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0, iVar_Eddy = 0;
	unsigned short iVar_GridVel = 0, iVar_PressMach = 0, iVar_TempLam = 0;
  
	unsigned long iPoint = 0, jPoint = 0;
  
	bool Write_Mean_Solution = false;
	bool grid_movement = ((config->GetUnsteady_Simulation() &&
                         config->GetWrt_Unsteady() &&
                         config->GetGrid_Movement()) ||
                        ((config->GetUnsteady_Simulation() == TIME_SPECTRAL) &&
                         config->GetGrid_Movement()));
    bool incompressible = config->GetIncompressible();
  
	if (Kind_Solver == AEROACOUSTIC_EULER) {
		if (val_iZone == ZONE_0) Kind_Solver = EULER;
		if (val_iZone == ZONE_1) Kind_Solver = WAVE_EQUATION;
	}
	if (Kind_Solver == PLASMA_EULER) {
		if (val_iZone == ZONE_0) Kind_Solver = PLASMA_EULER;
		if (val_iZone == ZONE_1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
	if (Kind_Solver == PLASMA_NAVIER_STOKES) {
		if (val_iZone == ZONE_0) Kind_Solver = PLASMA_NAVIER_STOKES;
		if (val_iZone == ZONE_1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
  
	/*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/
	switch (Kind_Solver) {
    case EULER : case NAVIER_STOKES:
      FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
      FirstIndex = PLASMA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      Write_Mean_Solution = config->GetWrite_Mean_Solution();
      break;
    case RANS :
      FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; ThirdIndex = NONE;
      break;
    case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES:
      FirstIndex = FLOW_SOL; SecondIndex = LEVELSET_SOL; ThirdIndex = NONE;
      break;
    case FREE_SURFACE_RANS:
      FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; ThirdIndex = LEVELSET_SOL;
      break;
    case ELECTRIC_POTENTIAL:
      FirstIndex = ELEC_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case WAVE_EQUATION:
      FirstIndex = WAVE_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case LINEAR_ELASTICITY:
      FirstIndex = FEA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES :
      FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES :
      FirstIndex = ADJPLASMA_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES:  case ADJ_FREE_SURFACE_RANS:
      FirstIndex = ADJFLOW_SOL; SecondIndex = ADJLEVELSET_SOL; ThirdIndex = NONE;
      break;
    case ADJ_RANS :
      FirstIndex = ADJFLOW_SOL; SecondIndex = ADJTURB_SOL; ThirdIndex = NONE;
      break;
    case LIN_EULER : case LIN_NAVIER_STOKES : ThirdIndex = NONE;
      FirstIndex = LINFLOW_SOL; SecondIndex = NONE;
      break;
    default: SecondIndex = NONE; ThirdIndex = NONE;
      break;
	}
	nVar_First = solution[FirstIndex]->GetnVar();
	if (SecondIndex != NONE) nVar_Second = solution[SecondIndex]->GetnVar();
    if (ThirdIndex != NONE) nVar_Third = solution[ThirdIndex]->GetnVar();
	nVar_Consv = nVar_First + nVar_Second + nVar_Third;
	nVar_Total = 2*nVar_Consv;
  
	/*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
	if (grid_movement) {
		iVar_GridVel = nVar_Total;
		if (geometry->GetnDim() == 2) nVar_Total += 2;
		else if (geometry->GetnDim() == 3) nVar_Total += 3;
	}
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) ||
      (Kind_Solver == FREE_SURFACE_RANS)) {
    /*--- Pressure and Mach ---*/
    iVar_PressMach = nVar_Total;
    nVar_Total += 2;
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
    /*--- Temperature & Laminar Viscosity ---*/
    iVar_TempLam = nVar_Total;
    nVar_Total += 2;
  }
  
  if ((Kind_Solver == RANS) || (Kind_Solver == FREE_SURFACE_RANS)) {
    /*--- Eddy Viscosity ---*/
    iVar_Eddy = nVar_Total;
    nVar_Total += 1;
  }
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
#ifdef NO_MPI
  
	/*--- In serial, the single process has access to all solution data,
	so it is simple to retrieve and store inside Solution_Data. ---*/
  nGlobal_Poin = geometry->GetnPoint(); //geometry->GetnPointDomain();  (F.P.)
	Data = new double*[nVar_Total];
	for (iVar = 0; iVar < nVar_Total; iVar++) {
		Data[iVar] = new double[nGlobal_Poin];
	}
  Volume = new double[nGlobal_Poin];
  
  /*--- In case there is grid movement ---*/
  double *Grid_Vel;
  
	/*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic/sliding halo nodes). ---*/
  jPoint = 0;
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//		if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
      
			/*--- Solution (first, and second system of equations) ---*/
			jVar = 0;
			for (iVar = 0; iVar < nVar_First; iVar++) {
				Data[jVar][jPoint] = solution[FirstIndex]->node[iPoint]->GetSolution(iVar);
                jVar++;
			}
      
			for (iVar = 0; iVar < nVar_Second; iVar++) {
				Data[jVar][jPoint] = solution[SecondIndex]->node[iPoint]->GetSolution(iVar);
                jVar++;
			}
            
            for (iVar = 0; iVar < nVar_Third; iVar++) {
				Data[jVar][jPoint] = solution[ThirdIndex]->node[iPoint]->GetSolution(iVar);
                jVar++;
			}
      
			/*--- Residual (first, and second system of equations) ---*/
			for (iVar = 0; iVar < nVar_First; iVar++) {
				Data[jVar][jPoint] = solution[FirstIndex]->node[iPoint]->GetResidual(iVar);
                jVar++;
			}
      
			for (iVar = 0; iVar < nVar_Second; iVar++) {
				Data[jVar][jPoint] = solution[SecondIndex]->node[iPoint]->GetResidual(iVar);
                jVar++;
			}
            
            for (iVar = 0; iVar < nVar_Third; iVar++) {
				Data[jVar][jPoint] = solution[ThirdIndex]->node[iPoint]->GetResidual(iVar);
                jVar++;
			}
      
            /*--- Store the volume of the current CV ---*/
            Volume[jPoint] = geometry->node[iPoint]->GetVolume();
            
			/*--- For unsteady problems with grid movement, write the mesh velocities
             also, in case we need them for the unsteady adjoint. ---*/
			if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady() && grid_movement) {
				Grid_Vel = geometry->node[iPoint]->GetGridVel();
				for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
					Data[jVar][jPoint] = Grid_Vel[iDim];
                    jVar++;
				}
			}
      
			/*--- Any extra data for output files ---*/
            switch (config->GetKind_Solver()) {
                case EULER: case FREE_SURFACE_EULER:
                    /*--- Load buffers with the pressure and mach variables. ---*/
                    if (incompressible) {
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetPressure(incompressible)*config->GetPressure_Ref(); jVar++;
                        Data[jVar][jPoint] = sqrt(solution[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solution[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
                    } else {
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
                        Data[jVar][jPoint] = sqrt(solution[FLOW_SOL]->node[iPoint]->GetVelocity2())/
                        solution[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
                    }
                    break;
                    /*--- Write pressure, mach, temperature and laminar viscosity ---*/
                case NAVIER_STOKES: case FREE_SURFACE_NAVIER_STOKES:
                    if (incompressible) {
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
                        Data[jVar][jPoint] = sqrt(solution[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solution[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
                        Data[jVar][jPoint] = 0.0; jVar++;
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
                    } else {
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
                        Data[jVar][jPoint] = sqrt(solution[FLOW_SOL]->node[iPoint]->GetVelocity2())/
                        solution[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
                    }
                    break;
                    /*--- Write  pressure, mach, temperature, laminar viscosity, eddy viscosity ---*/
                case RANS: case FREE_SURFACE_RANS:
                    if (incompressible) {
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
                        Data[jVar][jPoint] = sqrt(solution[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/sqrt(config->GetBulk_Modulus()/(solution[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref())); jVar++;
                        Data[jVar][jPoint] = 0.0; jVar++;
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(); jVar++;
                    } else {
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetPressure(incompressible); jVar++;
                        Data[jVar][jPoint] = sqrt(solution[FLOW_SOL]->node[iPoint]->GetVelocity2())/
                        solution[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
                        Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
                    }
                    Data[jVar][jPoint] = solution[FLOW_SOL]->node[iPoint]->GetEddyViscosity(); jVar++;
                    break;
//			}  (F.P.)
		}
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    
	}

#else
  
	/*--- MPI preprocessing ---*/
  
  int rank = MPI::COMM_WORLD.Get_rank();
	int nProcessor = MPI::COMM_WORLD.Get_size();
	int iProcessor;
  
	/*--- Local variables needed for merging with MPI ---*/
  unsigned short CurrentIndex;
  
	unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
	unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
	unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
    
  /*--- Each processor sends its local number of nodes to the master. ---*/
	//!
	//! TO DO: MPI I/O for writing the solution files.
	//!
  nLocalPoint = geometry->GetnPoint(); // geometry->GetnPointDomain();  (F.P.)
	Buffer_Send_nPoint[0] = nLocalPoint;
	if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Gather(&Buffer_Send_nPoint, 1, MPI::UNSIGNED_LONG,
                         Buffer_Recv_nPoint, 1, MPI::UNSIGNED_LONG, MASTER_NODE);
	nBuffer_Scalar = MaxLocalPoint;

	//!
	//! TO DO: Here is where the code can be extended to an arbitrary number
	//! of variables specified by the user (name & type), and more
	//! logic needs to be done.
	//!
  
	/*--- Send and Recv buffers. ---*/
	double *Buffer_Send_Var = new double[MaxLocalPoint];
	double *Buffer_Recv_Var = NULL;
  
	double *Buffer_Send_Res = new double[MaxLocalPoint];
	double *Buffer_Recv_Res = NULL;
  
	double *Buffer_Send_Vol = new double[MaxLocalPoint];
	double *Buffer_Recv_Vol = NULL;
  
	unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
	unsigned long *Buffer_Recv_GlobalIndex = NULL;
    
	/*--- Prepare the receive buffers in the master node only. ---*/
	if (rank == MASTER_NODE) {
    
		Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_Res = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_Vol = new double[nProcessor*MaxLocalPoint];
		Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
		/*--- Sum total number of nodes to be written and allocate arrays ---*/
		nGlobal_Poin = 0;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
		}
		Data = new double*[nVar_Total];
		for (iVar = 0; iVar < nVar_Total; iVar++) {
			Data[iVar] = new double[nGlobal_Poin];
		}
		Volume = new double[nGlobal_Poin];
	
  }
  
	/*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
	for (iVar = 0; iVar < nVar_Consv; iVar++) {
    
		/*--- Logic for which solution class to draw from. ---*/
		jVar = iVar;
		CurrentIndex = FirstIndex;
		if ((SecondIndex != NONE) && (iVar > nVar_First-1)) {
			jVar = iVar - nVar_First;
			CurrentIndex = SecondIndex;
		}
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//			if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
        
				/*--- Get this variable into the temporary send buffer. ---*/
				Buffer_Send_Var[jPoint] = solution[CurrentIndex]->node[iPoint]->GetSolution(jVar);
				Buffer_Send_Res[jPoint] = solution[CurrentIndex]->node[iPoint]->GetResidual(jVar);
        
				/*--- Only send/recv the volumes & global indices during the first loop ---*/
				if (iVar == 0) {
					Buffer_Send_Vol[jPoint]         = geometry->node[iPoint]->GetVolume();
					Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
				}
				jPoint++;
//			}  (F.P.)
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		if (iVar == 0) {
			MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
			MPI::COMM_WORLD.Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI::UNSIGNED_LONG,
                             MASTER_NODE);
		}
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
      jPoint = 0;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
					Data[iVar+nVar_Consv][iGlobal_Index] = Buffer_Recv_Res[jPoint];
					if (iVar == 0) Volume[iGlobal_Index] = Buffer_Recv_Vol[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
 
	/*--- Additional communication routine for the grid velocity. Note that
   we are reusing the same temporary buffers from above for efficiency.
   Also, in the future more routines like this could be used to write
   an arbitrary number of additional variables to the file. ---*/
	if (grid_movement) {
    
		/*--- Loop over this partition to collect the current variable ---*/
		jPoint = 0; double *Grid_Vel;
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//			if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
        
				/*--- Load buffers with the three grid velocity components. ---*/
				Grid_Vel = geometry->node[iPoint]->GetGridVel();
				Buffer_Send_Var[jPoint] = Grid_Vel[0];
				Buffer_Send_Res[jPoint] = Grid_Vel[1];
				if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
				jPoint++;
//			}  (F.P.)
		}
    
		/*--- Gather the data on the master node. ---*/
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
		if (geometry->GetnDim() == 3) {
			MPI::COMM_WORLD.Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             Buffer_Recv_Vol, nBuffer_Scalar, MPI::DOUBLE,
                             MASTER_NODE);
		}
    
		/*--- The master node unpacks and sorts this variable by global index ---*/
		if (rank == MASTER_NODE) {
			jPoint = 0; iVar = iVar_GridVel;
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
				for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
					/*--- Get global index, then loop over each variable and store ---*/
					iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
					Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
					Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
					if (geometry->GetnDim() == 3)
						Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
					jPoint++;
				}
				/*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
				jPoint = (iProcessor+1)*nBuffer_Scalar;
			}
		}
	}
  
  /*--- Communicate Pressure and Mach ---*/
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//      if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        Buffer_Send_Var[jPoint] = solution[FLOW_SOL]->node[iPoint]->GetPressure(incompressible);
        if (incompressible) {
          Buffer_Send_Res[jPoint] = 0.0;
        } else {
          Buffer_Send_Res[jPoint] = sqrt(solution[FLOW_SOL]->node[iPoint]->GetVelocity2())/
          solution[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
        }
        jPoint++;
//      }  (F.P.)
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_PressMach;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate Temperature & Laminar Viscosity ---*/
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//      if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
        
        /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
        Buffer_Send_Var[jPoint] = solution[FLOW_SOL]->node[iPoint]->GetTemperature();
        Buffer_Send_Res[jPoint] = solution[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        jPoint++;
//      }  (F.P.)
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    MPI::COMM_WORLD.Gather(Buffer_Send_Res, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Res, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_TempLam;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
          Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Communicate the Eddy Viscosity ---*/
  if (Kind_Solver == RANS) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//      if (geometry->node[iPoint]->GetDomain()) {  (F.P.)
        
        /*--- Load buffers with the pressure and mach variables. ---*/
        Buffer_Send_Var[jPoint] = solution[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
        jPoint++;
//      }  (F.P.)
    }
    
    /*--- Gather the data on the master node. ---*/
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Gather(Buffer_Send_Var, nBuffer_Scalar, MPI::DOUBLE,
                           Buffer_Recv_Var, nBuffer_Scalar, MPI::DOUBLE,
                           MASTER_NODE);
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0; iVar = iVar_Eddy;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
  }
  
	/*--- Immediately release the temporary buffers. ---*/
  
	delete [] Buffer_Send_Var;
	delete [] Buffer_Send_Res;
	delete [] Buffer_Send_Vol;
	delete [] Buffer_Send_GlobalIndex;
	if (rank == MASTER_NODE) {
		delete [] Buffer_Recv_Var;
		delete [] Buffer_Recv_Res;
		delete [] Buffer_Recv_Vol;
		delete [] Buffer_Recv_GlobalIndex;
	}
  
#endif

}

void COutput::WriteRestart(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
  unsigned short iVar;
  unsigned long iPoint;
  
  ofstream restart_file;
  string filename, AdjExt, UnstExt;
  unsigned long iExtIter = config->GetExtIter();
  char buffer[50];
  unsigned short Kind_ObjFunc = config->GetKind_ObjFunc();
  
  /*--- Retrieve filename from config ---*/
  if (config->IsAdjoint())
    filename = config->GetRestart_AdjFileName();
  else
    filename = config->GetRestart_FlowFileName();
  
  /*--- Remove restart filename extension ---*/
  filename.erase(filename.end()-4, filename.end());
  
  /*--- The adjoint problem requires a particular extension. ---*/
  if (config->IsAdjoint()) {
    switch (Kind_ObjFunc) {
      case DRAG_COEFFICIENT:      AdjExt = "_cd";   break;
      case LIFT_COEFFICIENT:      AdjExt = "_cl";   break;
      case SIDEFORCE_COEFFICIENT: AdjExt = "_csf";  break;
      case PRESSURE_COEFFICIENT:  AdjExt = "_cp";   break;
      case MOMENT_X_COEFFICIENT:  AdjExt = "_cmx";  break;
      case MOMENT_Y_COEFFICIENT:  AdjExt = "_cmy";  break;
      case MOMENT_Z_COEFFICIENT:  AdjExt = "_cmz";  break;
      case EFFICIENCY:            AdjExt = "_eff";  break;
      case EQUIVALENT_AREA:       AdjExt = "_ea";   break;
      case NEARFIELD_PRESSURE:    AdjExt = "_nfp";  break;
      case FORCE_X_COEFFICIENT:   AdjExt = "_cfx";  break;
      case FORCE_Y_COEFFICIENT:   AdjExt = "_cfy";  break;
      case FORCE_Z_COEFFICIENT:   AdjExt = "_cfz";  break;
      case THRUST_COEFFICIENT:    AdjExt = "_ct";   break;
      case TORQUE_COEFFICIENT:    AdjExt = "_cq";   break;
      case FIGURE_OF_MERIT:       AdjExt = "_merit";break;
      case FREESURFACE:           AdjExt = "_fs";   break;
      case NOISE:                 AdjExt = "_fwh";  break;
    }
    filename.append(AdjExt);
  }
  
  /*--- Unsteady problems require the physical timestep to be appended. ---*/
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    
    if (int(val_iZone) < 10) sprintf (buffer, "_0000%d", int(val_iZone));
    if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d", int(val_iZone));
    if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d", int(val_iZone));
    if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d", int(val_iZone));
    if (int(val_iZone) >= 10000) sprintf (buffer, "_%d", int(val_iZone));
    UnstExt = string(buffer);
    filename.append(UnstExt);
  } else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))
      sprintf (buffer, "_0000%d", int(iExtIter));
    if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))
      sprintf (buffer, "_000%d", int(iExtIter));
    if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))
      sprintf (buffer, "_00%d", int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000))
      sprintf (buffer, "_0%d", int(iExtIter));
    if (int(iExtIter) >= 10000)
      sprintf (buffer, "_%d", int(iExtIter));
    UnstExt = string(buffer);
    filename.append(UnstExt);
  }
  
  /*--- Lastly, add the .dat extension ---*/
  filename.append(".dat");
  
  /*--- Open the restart file and write the solution. ---*/
  restart_file.open(filename.c_str(), ios::out);
  restart_file.precision(15);
  
  
  /*--- Write the basic header based on the particular solution ----*/
  restart_file << "\"PointID\"";
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
    restart_file << ", \"Conservative_" << iVar+1<<"\"";
  }
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
    restart_file << ", \"Residual_" << iVar+1<<"\"";
  }
  restart_file << ", \"Volume\"\n";
  
  //  /*--- Add names for any extra variables (this will need to be adjusted). ---*/
  //	if (grid_movement) {
  //    if (nDim == 2) {
  //      Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\"";
  //    } else {
  //      Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\",\"Grid_Velz\"";
  //    }
  //	}
  //
  //  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
  //      (Kind_Solver == FREE_SURFACE_EULER) || (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
  //    Tecplot_File << ",\"Pressure\",\"Mach\"";
  //  }
  //
  //  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
  //      (Kind_Solver == FREE_SURFACE_NAVIER_STOKES) || (Kind_Solver == FREE_SURFACE_RANS)) {
  //    Tecplot_File << ",\"Temperature\",\"Laminar Viscosity\"";
  //  }
  //
  //  if ((Kind_Solver == RANS) || (Kind_Solver == FREE_SURFACE_RANS)) {
  //    Tecplot_File << ",\"Eddy Viscosity\"";
  //  }
  
//  for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) { (F.P.)
  for (iPoint = 0; iPoint < nGlobal_Doma; iPoint++) {
    
    /*--- Index of the point ---*/
    restart_file << iPoint << "\t";
    
    /*--- Loop over the vars/residuals and write the values to file ---*/
    for (iVar = 0; iVar < 2*nVar_Consv; iVar++) {
      restart_file << scientific << Data[iVar][iPoint] << "\t";
    }
    
    /*--- Last required column is the volume (for mesh adaptation) ---*/
    restart_file << scientific << Volume[iPoint] << "\t";
    
    restart_file << endl;
  }
  
  /*--- Binary version - eventually ---*/
  //  restart_file.close();
  //  restart_file.open("restart_file.bin", ios::binary);
  //  restart_file.precision(15);
  //  restart_file.write( (char *)&Data[0][0], sizeof(double)*(nGlobal_Poin*nVar_Consv*2) );
  //  restart_file.write( (char *)&Volume[0], sizeof(double)*(nGlobal_Poin) );
  //  restart_file.close();
  
}

void COutput::SetHistory_Header(ofstream *ConvHist_file, CConfig *config) {
	char cstr[200], buffer[50];

	bool rotating_frame = config->GetRotating_Frame();
	bool equiv_area = config->GetEquivArea();
	bool turbulent = (config->GetKind_Solver() == RANS);
	unsigned short iSpecies;

	/*--- Write file name with extension ---*/
	strcpy (cstr, config->GetConv_FileName().data());
	if (config->GetOutput_FileFormat() == PARAVIEW)  sprintf (buffer, ".csv");
	if (config->GetOutput_FileFormat() == TECPLOT)  sprintf (buffer, ".plt");
	strcat(cstr,buffer);

	ConvHist_file->open(cstr, ios::out);
	ConvHist_file->precision(15);

	char begin[]= "\"Iteration\"";

	char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
	char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
	char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
	char free_surface_coeff[]= ",\"CFreeSurface\"";
	char plasma_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\",\"CQ\"";
	char wave_coeff[]= ",\"CWave\"";
	char fea_coeff[]= ",\"CFEA\"";
	char adj_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
	char adj_plasma_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";

	char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
	char adj_flow_resid[]= ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";

	char turb_resid[1000];
	switch (config->GetKind_Turb_Model()){
	case SA:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
	case SST:	sprintf (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
	}

	char adj_turb_resid[]= ",\"Res_AdjTurb[0]\"";

	char levelset_resid[]= ",\"Res_LevelSet\"";
	char adj_levelset_resid[]= ",\"Res_AdjLevelSet\"";

	char wave_resid[]= ",\"Res_Wave[0]\",\"Res_Wave[1]\"";

	char fea_resid[]= ",\"Res_FEA\"";

	char end[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";

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

	case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		ConvHist_file[0] << begin << flow_coeff << free_surface_coeff;
		ConvHist_file[0] << flow_resid << levelset_resid << end;
		break;

	case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
		ConvHist_file[0] << begin << plasma_coeff;
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
			ConvHist_file[0] << ",\"Res_Density[" << iSpecies << "]\",\"Res_Energy[" << iSpecies << "]\"";
		}		
		ConvHist_file[0] << end;
		break;

	case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS:
		ConvHist_file[0] << begin << adj_coeff << adj_flow_resid;
		if (turbulent) ConvHist_file[0] << adj_turb_resid;
		ConvHist_file[0] << end;
		break;

	case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
		ConvHist_file[0] << begin << adj_coeff << free_surface_coeff;
		ConvHist_file[0] << adj_flow_resid << adj_levelset_resid << end;
		break;

	case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES:
		ConvHist_file[0] << begin << adj_plasma_coeff;
		for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
			ConvHist_file[0] << ",\"Res_PsiDensity[" << iSpecies << "]\",\"Res_PsiEnergy[" << iSpecies << "]\"";
		}		
		ConvHist_file[0] << end;
		break;

	case WAVE_EQUATION:
		ConvHist_file[0] << begin << wave_coeff;
		ConvHist_file[0] << wave_resid << end;
		break;

	case LINEAR_ELASTICITY:
		ConvHist_file[0] << begin << fea_coeff;
		ConvHist_file[0] << fea_resid << end;
		break;

	case AEROACOUSTIC_EULER:
		ConvHist_file[0] << begin << flow_coeff;
		ConvHist_file[0] << wave_coeff;
		ConvHist_file[0] << flow_resid;
		ConvHist_file[0] << wave_resid << end;
		break;

	case ADJ_AEROACOUSTIC_EULER:
		ConvHist_file[0] << begin << adj_coeff;
		ConvHist_file[0] << adj_flow_resid;
		ConvHist_file[0] << wave_coeff;
		ConvHist_file[0] << wave_resid << end;
		break;
	}

	if (config->GetOutput_FileFormat() == TECPLOT) {
		char zone[]= "ZONE T= \"Convergence history\"";
		ConvHist_file[0] << zone << endl;
	}


}

void COutput::SetHistory_DualTime(CGeometry ***geometry, CSolution ****solution_container, CConfig **config,
		CIntegration ***integration, unsigned long iExtIter, unsigned short val_iZone) {

#ifdef NO_MPI

	bool write_heads = ((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*100)) == 0);
	bool write_solution = ((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*1)) == 0);
	bool incompressible = config[val_iZone]->GetIncompressible();
	bool rotating_frame = config[val_iZone]->GetRotating_Frame();

	/*--- Write the screen header---*/
	if ((iExtIter == 1) || (write_heads)) {

		switch (config[val_iZone]->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			if (incompressible) cout << endl << " DT Iter" << "      Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			else if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
			else cout << endl << " DT Iter" << "      Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			if (incompressible) cout << endl << " DT Iter" << "      Res[AdjRho]" << "   Sens_Geo" << "   Sens_Mach" << endl;
			else if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[AdjRho]" << "     Res[AdjRhoE]" << " Sens_Geo" << " Sens_Mach" << endl;
			else cout << endl << " DT Iter" << "      Res[AdjRho]" << "     Res[AdjRhoE]" << "   Sens_Geo" << "   Sens_Mach" << endl;
			break;
		case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "    Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "      CLevelSet" << endl;
			break;
		case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "    Res[AdjPress]" << "     Res[AdjDist]" << "   Sens_Geo" << "   Sens_Mach" << endl;
			break;
		case RANS :
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]" << " CThrust(Total)" << " CTorque(Total)" << endl;
			else cout << endl << " DT Iter" << "      Res[Rho]" << "       Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			break;
		case PLASMA_EULER : case PLASMA_NAVIER_STOKES :
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "      Res[Rho_0]" << "       Res[E_0]" << endl;
			break;
		case WAVE_EQUATION :
			cout << endl << " Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "      Res[Wave]" << "      CWave" << endl;
			break;
		case LINEAR_ELASTICITY:
			cout << endl << " Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			if (geometry[val_iZone][MESH_0]->GetnDim() == 2) cout << endl << " DT Iter" << "      Res[Displx]" << "    Res[Disply]" << "   CFEA(Total)"<<  endl;
			if (geometry[val_iZone][MESH_0]->GetnDim() == 3) cout << endl << " DT Iter" << "      Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "   CFEA(Total)"<<  endl;
			break;
		case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "      Res[Rho]" << "   Res[Displx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
			break;
		case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES:
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "      Res[Rho]" << "   Res[Wave]" << "   CLift(Total)" << "   CDrag(Total)" << "   CWave(Total)" << endl;
			break;
		case ADJ_AEROACOUSTIC_EULER:
			cout << endl << " Min DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time()<<
			". Max DT: " << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time() <<
			". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			cout << endl << " DT Iter" << "      Res[AdjRho]" << "   Res[AdjWave]" << "   Sens_Geo" << "   Sens_Mach" << endl;
			break;
		}
	}

	/*--- Write the solution on the screen and history file ---*/
	if ((iExtIter != 0) && (write_solution)) {
		switch (config[val_iZone]->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(0));
			if (!incompressible) {
				cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(geometry[val_iZone][0]->GetnDim()+1));
			}
			if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) {
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CT();
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CQ();
			} else {
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
			}
			cout << endl;
			break;
		case ADJ_EULER : case ADJ_NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_Max(0));
			if (!incompressible) {
				cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_Max(geometry[val_iZone][0]->GetnDim()+1));
			}
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Mach();
			cout << endl;
			break;
		case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][LEVELSET_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][LEVELSET_SOL]->GetTotal_CFreeSurface();
			cout << endl;
			break;
		case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(17); cout << log10(solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_Max(0));
			cout.width(17); cout << log10(solution_container[val_iZone][MESH_0][ADJLEVELSET_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Mach();
			cout << endl;
			break;
		case RANS :
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][TURB_SOL]->GetRes_Max(0));
			if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) {
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CT();
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CQ();
			} else {
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
				cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
			}
			cout << endl;
			break;
		case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			if ((config[val_iZone]->GetKind_GasModel() == AIR7) || (config[val_iZone]->GetKind_GasModel() == O2) || (config[val_iZone]->GetKind_GasModel() == N2) || (config[val_iZone]->GetKind_GasModel() == AIR5)) {
				cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetRes_Max(0));
				cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetRes_Max(3));
			}
			cout << endl;
			break;
		case WAVE_EQUATION:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(15); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][WAVE_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[val_iZone][MESH_0][WAVE_SOL]->GetTotal_CWave();
			cout << endl;
			break;
		case LINEAR_ELASTICITY:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(17); cout << log10(solution_container[val_iZone][MESH_0][FEA_SOL]->GetRes_Max(0));
			cout.width(15); cout << log10(solution_container[val_iZone][MESH_0][FEA_SOL]->GetRes_Max(1));
			if (geometry[val_iZone][MESH_0]->GetnDim() == 3) { cout.width(15); cout << log10(solution_container[val_iZone][MESH_0][FEA_SOL]->GetRes_Max(2)); }
			cout.width(14); cout << "0.000000";
			cout << endl;
			break;
		case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][FEA_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
			cout << endl;
			break;
		case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[ZONE_1][MESH_0][WAVE_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
			cout.width(15); cout << solution_container[ZONE_1][MESH_0][WAVE_SOL]->GetTotal_CWave();
			cout << endl;
			break;
		case ADJ_AEROACOUSTIC_EULER:
			cout.precision(6);
			cout.setf(ios::fixed,ios::floatfield);
			cout.width(8); cout << iExtIter;
			cout.width(14); cout << log10(solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_Max(0));
			cout.width(14); cout << log10(solution_container[ZONE_1][MESH_0][WAVE_SOL]->GetRes_Max(0));
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Mach();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_AoA();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Press();
			cout.width(15); cout << solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Temp();
			cout.width(15); cout << solution_container[ZONE_1][MESH_0][WAVE_SOL]->GetTotal_CWave();
			cout << endl;
			break;
		}
	}
	cout.unsetf(ios::fixed);

#else

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	double *sbuf_residual_flow = NULL, *sbuf_residual_turbulent = NULL, *sbuf_residual_levelset = NULL, *sbuf_force = NULL, *sbuf_time = NULL,
			*sbuf_residual_plasma = NULL, *rbuf_residual_flow = NULL, *rbuf_residual_turbulent = NULL, *rbuf_residual_levelset = NULL, *rbuf_force = NULL, *rbuf_time = NULL,
			*sbuf_residual_adjoint = NULL, *rbuf_residual_adjoint = NULL, *sbuf_residual_linearized = NULL,
			*rbuf_residual_linearized = NULL, *rbuf_residual_plasma = NULL;
	unsigned short iVar, buf_convergence = 0, *sbuf_conv = NULL, *rbuf_conv = NULL;

	bool write_heads = ((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*100)) == 0);
	bool write_solution = ((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*1)) == 0);
	bool compressible = !config[val_iZone]->GetIncompressible();
	bool incompressible = config[val_iZone]->GetIncompressible();
	bool rotating_frame = config[val_iZone]->GetRotating_Frame();
	unsigned short nDim = geometry[val_iZone][MESH_0]->GetnDim();

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
	unsigned short nVar_Plasma = config[val_iZone]->GetnMonatomics()*(nDim+2) + config[val_iZone]->GetnDiatomics()*(nDim+3);

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
	switch (config[val_iZone]->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES: case RANS :

		for (iVar = 0; iVar < nVar_Flow; iVar++)
			sbuf_residual_flow[iVar] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(iVar);

		if (config[val_iZone]->GetKind_Solver() == RANS) {
			for (iVar = 0; iVar < nVar_Turb; iVar++)
				sbuf_residual_turbulent[iVar] = log10 (solution_container[val_iZone][MESH_0][TURB_SOL]->GetRes_Max(iVar));
		}

		sbuf_force[0]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
		sbuf_force[1]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
		sbuf_force[2]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CSideForce();
		sbuf_force[4]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CMx();
		sbuf_force[5]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CMy();
		sbuf_force[6]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CMz();
		sbuf_force[7]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CEff();
		sbuf_force[8]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CEquivArea(),
				sbuf_force[9]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CNearFieldOF();
		sbuf_force[10] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CFx();
		sbuf_force[11] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CFy();
		sbuf_force[12] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CFz();
		sbuf_force[13] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CT();
		sbuf_force[14] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CQ();
		sbuf_force[15] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CMerit();
		sbuf_conv[0] = integration[val_iZone][FLOW_SOL]->GetConvergence();
		break;

	case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:

		for (iVar = 0; iVar < nVar_Flow; iVar++)
			sbuf_residual_flow[iVar] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(iVar);

		for (iVar = 0; iVar < nVar_LevelSet; iVar++)
			sbuf_residual_levelset[iVar] = solution_container[val_iZone][MESH_0][LEVELSET_SOL]->GetRes_Max(iVar);

		sbuf_force[0]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
		sbuf_force[1]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
		sbuf_force[2]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CSideForce();
		sbuf_force[4]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CMx();
		sbuf_force[5]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CMy();
		sbuf_force[6]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CMz();
		sbuf_force[7]  = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CEff();
		sbuf_force[8]  = solution_container[val_iZone][MESH_0][LEVELSET_SOL]->GetTotal_CFreeSurface();
		sbuf_conv[0] = integration[val_iZone][FLOW_SOL]->GetConvergence();
		break;

	case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
		sbuf_force[0]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CLift();
		sbuf_force[1]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CDrag();
		sbuf_force[2]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CSideForce();
		sbuf_force[4]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CMx();
		sbuf_force[5]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CMy();
		sbuf_force[6]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CMz();
		sbuf_force[7]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CEff();
		sbuf_force[8]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CEquivArea(),
				sbuf_force[9]  = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CNearFieldOF();
		sbuf_force[10] = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CFx();
		sbuf_force[11] = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CFy();
		sbuf_force[12] = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CFz();
		sbuf_force[13] = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CT();
		sbuf_force[14] = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CQ();
		sbuf_force[15] = solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetTotal_CMerit();

		for (iVar = 0; iVar < nVar_Plasma; iVar++)
			sbuf_residual_plasma[iVar] = log10 (solution_container[val_iZone][MESH_0][PLASMA_SOL]->GetRes_Max(iVar));

		sbuf_conv[0] = integration[val_iZone][PLASMA_SOL]->GetConvergence();
		break;

	case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :

		for (iVar = 0; iVar < nVar_Adj; iVar++)
			sbuf_residual_adjoint[iVar] = log10 (solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_Max(iVar));

		sbuf_force[0] = solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo();
		sbuf_force[1] = solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Mach();
		sbuf_force[2] = solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_AoA();
		sbuf_force[3] = solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Press();
		sbuf_force[4] = solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Temp();
		sbuf_conv[0] = integration[val_iZone][ADJFLOW_SOL]->GetConvergence();
		break;

	case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:

		for (iVar = 0; iVar < nVar_Adj; iVar++)
			sbuf_residual_adjoint[iVar] = log10 (solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_Max(iVar));

		sbuf_force[0] = solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo();
		sbuf_force[1] = solution_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Mach();
		sbuf_conv[0] = integration[val_iZone][ADJFLOW_SOL]->GetConvergence();
		break;

	case LIN_EULER : case LIN_NAVIER_STOKES :

		for (iVar = 0; iVar < nVar_Lin; iVar++)
			sbuf_residual_linearized[iVar] = log10 (solution_container[val_iZone][MESH_0][LINFLOW_SOL]->GetRes_Max(iVar));

		sbuf_force[0] = solution_container[val_iZone][MESH_0][LINFLOW_SOL]->GetTotal_CDeltaLift();
		sbuf_force[1] = solution_container[val_iZone][MESH_0][LINFLOW_SOL]->GetTotal_CDeltaDrag();
		sbuf_conv[0] = integration[val_iZone][LINFLOW_SOL]->GetConvergence();
		break;
	}
	sbuf_time[0] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMin_Delta_Time();
	sbuf_time[1] = solution_container[val_iZone][MESH_0][FLOW_SOL]->GetMax_Delta_Time();

	/*--- Send/Receive information ---*/
	switch (config[val_iZone]->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES:
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		MPI::COMM_WORLD.Reduce(sbuf_residual_levelset, rbuf_residual_levelset, nVar_LevelSet, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case RANS :
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		MPI::COMM_WORLD.Reduce(sbuf_residual_turbulent, rbuf_residual_turbulent, nVar_Turb, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
		MPI::COMM_WORLD.Reduce(sbuf_residual_plasma, rbuf_residual_plasma, nVar_Plasma, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	case ADJ_EULER : case ADJ_NAVIER_STOKES :
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjoint, rbuf_residual_adjoint, nVar_Adj, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
		break;
	case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
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

	switch (config[val_iZone]->GetKind_Solver()) {
	case EULER: case NAVIER_STOKES: case RANS:
		if (buf_convergence == 1) integration[val_iZone][FLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][FLOW_SOL]->SetConvergence(false); break;
	case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		if (buf_convergence == 1) integration[val_iZone][FLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][FLOW_SOL]->SetConvergence(false);
		break;
	case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
		if (buf_convergence == 1) integration[val_iZone][ADJFLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][ADJFLOW_SOL]->SetConvergence(false); break;
	case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
		if (buf_convergence == 1) integration[val_iZone][ADJFLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][ADJFLOW_SOL]->SetConvergence(false); break;
	case LIN_EULER: case LIN_NAVIER_STOKES:
		if (buf_convergence == 1) integration[val_iZone][LINFLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][LINFLOW_SOL]->SetConvergence(false); break;
	case ELECTRIC_POTENTIAL:
		if (buf_convergence == 1) integration[val_iZone][ELEC_SOL]->SetConvergence(true);
		else integration[val_iZone][ELEC_SOL]->SetConvergence(false); break;
	case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
		if (buf_convergence == 1) integration[val_iZone][PLASMA_SOL]->SetConvergence(true);
		else integration[val_iZone][PLASMA_SOL]->SetConvergence(false); break;
	case WAVE_EQUATION:
		if (buf_convergence == 1) integration[val_iZone][WAVE_SOL]->SetConvergence(true);
		else integration[val_iZone][WAVE_SOL]->SetConvergence(false); break;
	case LINEAR_ELASTICITY:
		if (buf_convergence == 1) integration[val_iZone][FEA_SOL]->SetConvergence(true);
		else integration[val_iZone][FEA_SOL]->SetConvergence(false); break;
	}

	/*--- Write result using the master node ---*/
	if (rank == MASTER_NODE) {

		/*--- Write the screen header---*/
		if ((iExtIter == 1) || (write_heads)) {
			cout << endl << " Min DT: " << rbuf_time[0] << 
					". Max DT: " << rbuf_time[1] <<
					". Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
			switch (config[val_iZone]->GetKind_Solver()) {
			case EULER : case NAVIER_STOKES:
				if (incompressible) cout << endl << " DT Iter" << "      Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				else if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
				else cout << endl << " DT Iter" << "      Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				break;
			case RANS :
				if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) cout << endl << " DT Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]" << " CThrust(Total)" << " CTorque(Total)" << endl;
				else cout << endl << " DT Iter" << "      Res[Rho]" << "       Res[nu]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				break;
			case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
				cout << endl << " DT Iter" << "    Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "      CLevelSet" << endl;
				break;
			case ADJ_EULER : case ADJ_NAVIER_STOKES :
				if (incompressible) cout << endl << " DT Iter" << "  Res[Psi_Press]" << "     Sens_Geo" << "    Sens_Mach" << "    Sens_AoA"<< "    Sens_Press"<< "    Sens_Temp"<< endl;
				else cout << endl << " DT Iter" << "  Res[Psi_Rho]" << "    Res[Psi_E]" << "     Sens_Geo" << "    Sens_Mach" << "    Sens_AoA"<< "    Sens_Press"<< "    Sens_Temp"<< endl;
				break;
			case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES : case ADJ_FREE_SURFACE_RANS :
				if (incompressible) cout << endl << " DT Iter" << "  Res[Psi_Press]" << "     Sens_Geo" << "    Sens_Mach" << "    Sens_AoA"<< "    Sens_Press"<< "    Sens_Temp"<< endl;
				else cout << endl << " DT Iter" << "  Res[Psi_Rho]" << "    Res[Psi_E]" << "     Sens_Geo" << "    Sens_Mach" << "    Sens_AoA"<< "    Sens_Press"<< "    Sens_Temp"<< endl;
				break;
			case PLASMA_EULER : case PLASMA_NAVIER_STOKES :
				cout << endl << " DT Iter" << "      Res[Rho_0]" << "       Res[E_0]" << endl;
				break;
			}
		}
		/*--- Write the solution on the screen ---*/
		if ((iExtIter != 0) && (write_solution)) {
			switch (config[val_iZone]->GetKind_Solver()) {
			case EULER : case NAVIER_STOKES: case RANS:
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(8); cout << iExtIter;
				cout.width(14); cout << rbuf_residual_flow[0];
				if (config[val_iZone]->GetKind_Solver() == RANS) {
					cout.width(14); cout << rbuf_residual_turbulent[0];
				} else if (!incompressible) {
					cout.width(13); cout << rbuf_residual_flow[3];
				}
				if (rotating_frame && geometry[val_iZone][MESH_0]->GetnDim() == 3 ) {
					cout.width(15); cout << rbuf_force[13];
					cout.width(15); cout << rbuf_force[14];
				} else {
					cout.width(15); cout << rbuf_force[0];
					cout.width(15); cout << rbuf_force[1];
				}
				cout << endl;
				break;
			case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
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
				break;
			case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES : case ADJ_FREE_SURFACE_RANS :
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
				break;
			case PLASMA_EULER : case PLASMA_NAVIER_STOKES :
				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(8); cout << iExtIter;
				cout.width(14); cout << rbuf_residual_flow[0];
				cout.width(13); cout << rbuf_residual_flow[3];
				break;
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
		CIntegration ***integration, unsigned long iExtIter, unsigned long timeused, unsigned short val_iZone) {

#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
#else
	int rank = MASTER_NODE;
#endif

	/*--- WARNING: These buffers have hard-coded lengths. Note that you
	 may have to adjust them to be larger if adding more entries. ---*/
	char begin[1000], direct_coeff[1000], adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000],
	turb_resid[1000], adj_turb_resid[1000], plasma_resid[1000], adj_plasma_resid[1000], resid_aux[1000],
	levelset_resid[1000], adj_levelset_resid[1000], wave_coeff[1000], fea_coeff[1000], wave_resid[1000], 
	fea_resid[1000], end[1000];
	double dummy = 0.0;
	unsigned short iVar, iMarker;

	unsigned long LinSolvIter = 0;
	double timeiter = double(timeused)/double(iExtIter+1);

	unsigned short FinestMesh = config->GetFinestMesh();
	unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
	unsigned short nSpecies = config->GetnSpecies();

	bool incompressible = config->GetIncompressible();
	bool compressible = !config->GetIncompressible();
	bool rotating_frame = config->GetRotating_Frame();
	bool equiv_area = config->GetEquivArea();

	bool turbulent = (config->GetKind_Solver() == RANS);
	bool adjoint = ((config->GetKind_Solver() == ADJ_EULER) || (config->GetKind_Solver() == ADJ_NAVIER_STOKES) || 
			(config->GetKind_Solver() == ADJ_RANS) || (config->GetKind_Solver() == ADJ_FREE_SURFACE_EULER) ||
			(config->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_FREE_SURFACE_RANS) ||
			(config->GetKind_Solver() == ADJ_PLASMA_EULER) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES) ||
			(config->GetKind_Solver() == ADJ_AEROACOUSTIC_EULER));
	bool fluid_structure = ((config->GetKind_Solver() == FLUID_STRUCTURE_EULER) || (config->GetKind_Solver() == FLUID_STRUCTURE_NAVIER_STOKES) || 
			(config->GetKind_Solver() == FLUID_STRUCTURE_RANS));
	bool free_surface = ((config->GetKind_Solver() == FREE_SURFACE_EULER) || (config->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES) || (config->GetKind_Solver() == FREE_SURFACE_RANS) ||
			(config->GetKind_Solver() == ADJ_FREE_SURFACE_EULER) || (config->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES) || (config->GetKind_Solver() == ADJ_FREE_SURFACE_RANS));
	bool aeroacoustic = ((config->GetKind_Solver() == AEROACOUSTIC_EULER) || (config->GetKind_Solver() == AEROACOUSTIC_NAVIER_STOKES) ||
			(config->GetKind_Solver() == AEROACOUSTIC_RANS));
	bool wave = (config->GetKind_Solver() == WAVE_EQUATION);
	bool fea = (config->GetKind_Solver() == LINEAR_ELASTICITY);
	bool plasma = ((config->GetKind_Solver() == PLASMA_EULER) || (config->GetKind_Solver() == PLASMA_NAVIER_STOKES) ||
			(config->GetKind_Solver() == ADJ_PLASMA_EULER) || (config->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES));

	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || 
			(config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
	bool write_heads = (((iExtIter % (config->GetWrt_Con_Freq()*20)) == 0) || dual_time);

	/*--- Initialize variables to store information from all domains (direct solution) ---*/
	double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0, 
			Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
			Total_CT = 0.0, Total_CQ = 0.0, Total_CFreeSurface = 0.0, Total_CWave = 0.0, Total_CFEA = 0.0;

	/*--- Initialize variables to store information from all domains (adjoint solution) ---*/
	double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
	double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;

	double *sbuf_residual_flow = NULL, *sbuf_residual_turbulent = NULL, *sbuf_residual_levelset = NULL, *sbuf_residual_plasma = NULL;
	double *sbuf_residual_adjflow = NULL, *sbuf_residual_adjturbulent = NULL, *sbuf_residual_adjlevelset = NULL, *sbuf_residual_adjplasma = NULL;
	double *sbuf_residual_wave = NULL; double *sbuf_residual_fea = NULL;
	double *rbuf_residual_flow = NULL, *rbuf_residual_turbulent = NULL, *rbuf_residual_levelset = NULL, *rbuf_residual_plasma = NULL;
	double *rbuf_residual_adjflow = NULL, *rbuf_residual_adjturbulent = NULL, *rbuf_residual_adjlevelset = NULL, *rbuf_residual_adjplasma = NULL;
	double *rbuf_residual_wave = NULL; double *rbuf_residual_fea = NULL;
	double *sbuf_force = NULL,  *rbuf_force = NULL;
	unsigned long *sbuf_time = NULL, *rbuf_time = NULL;
	unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL; 

	/*--- Initialize number of variables ---*/
	unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0, nVar_Wave = 0, nVar_FEA = 0, nVar_Plasma = 0, 
			nVar_AdjFlow = 0, nVar_AdjPlasma = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0, nVar_Force = 0, nVar_Time = 0, nVar_Conv = 0;


	/*--- Direct problem variables ---*/
	if (compressible) nVar_Flow = nDim+2; else nVar_Flow = nDim+1;	
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
		case SA:	nVar_Turb = 1; break;
		case SST: nVar_Turb = 2; break;
		}

	if (wave) nVar_Wave = 2;
	if (fea) nVar_FEA = nDim;
	if (plasma) nVar_Plasma = config->GetnMonatomics()*(nDim+2) + config->GetnDiatomics()*(nDim+3);
	if (free_surface) nVar_LevelSet = 1;

	/*--- Adjoint problem variables ---*/
	if (compressible) nVar_AdjFlow = nDim+2; else nVar_AdjFlow = nDim+1;
	nVar_AdjTurb = 1;

	if (plasma) nVar_AdjPlasma = config->GetnMonatomics()*(nDim+2) + config->GetnDiatomics()*(nDim+3);
	if (free_surface) nVar_AdjLevelSet = 1;

	/*--- Other vectors ---*/
	nVar_Force = 16; nVar_Time = 3; nVar_Conv = 1;

	/*--- Allocate memory for send buffer ---*/
	sbuf_residual_flow = new double[nVar_Flow];
	sbuf_residual_turbulent = new double[nVar_Turb];
	sbuf_residual_plasma = new double[nVar_Plasma];
	sbuf_residual_levelset = new double[nVar_LevelSet];
	sbuf_residual_wave = new double[nVar_Wave];
	sbuf_residual_fea = new double[nVar_FEA];

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
	for (iVar = 0; iVar < nVar_FEA; iVar++) { sbuf_residual_fea[iVar] = -20; }

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
		rbuf_residual_fea = new double[nVar_FEA];

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
		for (iVar = 0; iVar < nVar_FEA; iVar++) { rbuf_residual_fea[iVar] = 0.0; }

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
	case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
	case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: case FLUID_STRUCTURE_RANS:
	case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES: case AEROACOUSTIC_RANS:
	case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
	case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
	case ADJ_AEROACOUSTIC_EULER:

		/*--- Flow solution coefficients (serial) ---*/
		Total_CLift += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
		Total_CDrag += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
		Total_CSideForce += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
		Total_CMx += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
		Total_CMy += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
		Total_CMz += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
		Total_CFx += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
		Total_CFy += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
		Total_CFz += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();
		Total_CEff += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEff();

		if (equiv_area) {
			Total_CEquivArea += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
			Total_CNearFieldOF += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
		}

		if (rotating_frame) {
			Total_CMerit += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
			Total_CT += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CT();
			Total_CQ += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CQ();
		}

		if (aeroacoustic)
			Total_CWave += solution_container[ZONE_1][FinestMesh][WAVE_SOL]->GetTotal_CWave();

		if (fluid_structure)
			Total_CFEA  += solution_container[ZONE_1][FinestMesh][FEA_SOL]->GetTotal_CFEA();

		/*--- Flow solution coefficients (parallel) ---*/
		sbuf_force[0] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
		sbuf_force[1] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
		sbuf_force[2] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
		sbuf_force[3] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
		sbuf_force[4] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
		sbuf_force[5] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
		sbuf_force[6] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
		sbuf_force[7] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
		sbuf_force[8] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();

		if (free_surface)
			sbuf_force[9] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFreeSurface();

		if (fluid_structure)
			sbuf_force[9] += solution_container[ZONE_1][FinestMesh][FEA_SOL]->GetTotal_CFEA();

		if (aeroacoustic)
			sbuf_force[9]  += solution_container[ZONE_1][FinestMesh][WAVE_SOL]->GetTotal_CWave();

		if (equiv_area) {
			sbuf_force[9] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea(),
					sbuf_force[10] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
		}

		if (rotating_frame) {
			sbuf_force[9]  += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
			sbuf_force[10] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CT();
			sbuf_force[11] += solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CQ();
		}

		/*--- Convergence criteria ---*/
		sbuf_conv[0] = integration[val_iZone][FLOW_SOL]->GetConvergence();

		/*--- Flow Residuals ---*/
		for (iVar = 0; iVar < nVar_Flow; iVar++)
			sbuf_residual_flow[iVar] = solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(iVar);

		LinSolvIter = (unsigned long) solution_container[val_iZone][FinestMesh][FLOW_SOL]->GetIterLinSolver();

		/*--- Turbulent residual ---*/
		if (turbulent) {
			for (iVar = 0; iVar < nVar_Turb; iVar++)
				sbuf_residual_turbulent[iVar] = solution_container[val_iZone][FinestMesh][TURB_SOL]->GetRes_Max(iVar);
		}

		/*--- Free Surface residual ---*/
		if (free_surface) {
			for (iVar = 0; iVar < nVar_LevelSet; iVar++)
				sbuf_residual_levelset[iVar] = solution_container[val_iZone][FinestMesh][LEVELSET_SOL]->GetRes_Max(iVar);
		}


		/*--- FEA residual ---*/
		if (fluid_structure) {
			for (iVar = 0; iVar < nVar_FEA; iVar++)
				sbuf_residual_fea[iVar] = solution_container[ZONE_1][FinestMesh][FEA_SOL]->GetRes_Max(iVar);
		}

		/*--- Aeroacoustic residual ---*/
		if (aeroacoustic) {
			for (iVar = 0; iVar < nVar_Wave; iVar++)
				sbuf_residual_wave[iVar] = solution_container[ZONE_1][FinestMesh][WAVE_SOL]->GetRes_Max(iVar);
		}

		/*--- Adjoint solver ---*/
		if (adjoint) {

			/*--- Adjoint solution coefficients (serial) ---*/
			Total_Sens_Geo  = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
			Total_Sens_Mach = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
			Total_Sens_AoA  = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA();
			Total_Sens_Press = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
			Total_Sens_Temp  = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();

			/*--- Adjoint solution coefficients (parallel) ---*/
			sbuf_force[0] = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
			sbuf_force[1] = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
			sbuf_force[2] = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA();
			sbuf_force[3] = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
			sbuf_force[4] = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();

			/*--- Convergence criteria ---*/
			sbuf_conv[0] = integration[val_iZone][ADJFLOW_SOL]->GetConvergence();

			/*--- Adjoint flow residuals ---*/
			for (iVar = 0; iVar < nVar_AdjFlow; iVar++)
				sbuf_residual_adjflow[iVar] = solution_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(iVar);

			/*--- Adjoint turbulent residuals ---*/
			if (turbulent) {
				for (iVar = 0; iVar < nVar_AdjTurb; iVar++)
					sbuf_residual_adjturbulent[iVar] = solution_container[val_iZone][FinestMesh][ADJTURB_SOL]->GetRes_Max(iVar);
			}

			/*--- Adjoint level set residuals ---*/
			if (free_surface) {
				for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++)
					sbuf_residual_adjlevelset[iVar] = solution_container[val_iZone][FinestMesh][ADJLEVELSET_SOL]->GetRes_Max(iVar);
			}
		}

		break;

	case PLASMA_EULER: case ADJ_PLASMA_EULER: case PLASMA_NAVIER_STOKES: case ADJ_PLASMA_NAVIER_STOKES:

		/*--- Plasma coefficients (serial) ---*/
		Total_CLift += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CLift();
		Total_CDrag += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CDrag();
		Total_CSideForce += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CSideForce();
		Total_CMx += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMx();
		Total_CMy += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMy();
		Total_CMz += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMz();
		Total_CFx += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFx();
		Total_CFy += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFy();
		Total_CFz += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFz();
		Total_CEff += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CEff();
    Total_CQ += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CQ();
      
		/*--- Plasma coefficients (parallel) ---*/
		sbuf_force[0] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CLift();
		sbuf_force[1] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CDrag();
		sbuf_force[2] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CSideForce();
		sbuf_force[3] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMx();
		sbuf_force[4] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMy();
		sbuf_force[5] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CMz();
		sbuf_force[6] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFx();
		sbuf_force[7] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFy();
		sbuf_force[8] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CFz();
    sbuf_force[9] += solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetTotal_CQ();
      
		/*--- Convergence criteria ---*/
		sbuf_conv[0] = integration[val_iZone][PLASMA_SOL]->GetConvergence();

		/*--- Plasma Residuals ---*/
		for (iVar = 0; iVar < nVar_Plasma; iVar++)
			sbuf_residual_plasma[iVar] = solution_container[val_iZone][FinestMesh][PLASMA_SOL]->GetRes_Max(iVar);

		if (adjoint) {

			/*--- Adjoint solution coefficients (serial) ---*/
			Total_Sens_Geo  = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Geo();
			Total_Sens_Mach = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Mach();
			Total_Sens_AoA  = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_AoA();
			Total_Sens_Press = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Press();
			Total_Sens_Temp  = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Temp();

			/*--- Adjoint solution coefficients (parallel) ---*/
			sbuf_force[0] = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Geo();
			sbuf_force[1] = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Mach();
			sbuf_force[2] = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_AoA();
			sbuf_force[3] = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Press();
			sbuf_force[4] = 0.0; //solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetTotal_Sens_Temp();

			/*--- Convergence criteria ---*/
			sbuf_conv[0] = integration[val_iZone][ADJPLASMA_SOL]->GetConvergence();

			/*--- Adjoint plasma residuals ---*/
			for (iVar = 0; iVar < nVar_AdjPlasma; iVar++)
				sbuf_residual_adjplasma[iVar] = solution_container[val_iZone][FinestMesh][ADJPLASMA_SOL]->GetRes_Max(iVar);
		}

		break;

	case WAVE_EQUATION:

		/*--- Wave coefficients (serial) ---*/
		Total_CWave += solution_container[val_iZone][FinestMesh][WAVE_SOL]->GetTotal_CWave();

		/*--- Wave coefficients (parallel) ---*/
		sbuf_force[0] += solution_container[val_iZone][FinestMesh][WAVE_SOL]->GetTotal_CWave();

		/*--- Convergence criteria ---*/
		sbuf_conv[0] = integration[val_iZone][WAVE_SOL]->GetConvergence();

		/*--- Wave Residuals ---*/
		for (iVar = 0; iVar < nVar_Wave; iVar++) {
			sbuf_residual_wave[iVar] = solution_container[val_iZone][FinestMesh][WAVE_SOL]->GetRes_Max(iVar);
		}

		break;

	case LINEAR_ELASTICITY:

		/*--- FEA coefficients (serial) ---*/
		Total_CFEA += solution_container[val_iZone][FinestMesh][FEA_SOL]->GetTotal_CFEA();

		/*--- FEA coefficients (parallel) ---*/
		sbuf_force[0] += solution_container[val_iZone][FinestMesh][FEA_SOL]->GetTotal_CFEA();

		/*--- Convergence criteria ---*/
		sbuf_conv[0] = integration[val_iZone][FEA_SOL]->GetConvergence();

		/*--- Plasma Residuals ---*/
		for (iVar = 0; iVar < nVar_FEA; iVar++) {
			sbuf_residual_fea[iVar] = solution_container[val_iZone][FinestMesh][FEA_SOL]->GetRes_Max(iVar);
		}

		break;

	}

	sbuf_time[0] = (unsigned long) timeiter;
	sbuf_time[1] = (unsigned long) timeused;
	sbuf_time[2] = LinSolvIter;

#ifndef NO_MPI

	unsigned short buf_convergence = 0;

	/*--- Send/Receive information ---*/
	switch (config->GetKind_Solver()) {
	case EULER : case NAVIER_STOKES: case RANS :
	case FREE_SURFACE_EULER : case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS :
	case FLUID_STRUCTURE_EULER : case FLUID_STRUCTURE_NAVIER_STOKES: case FLUID_STRUCTURE_RANS :
	case AEROACOUSTIC_EULER : case AEROACOUSTIC_NAVIER_STOKES: case AEROACOUSTIC_RANS :
		MPI::COMM_WORLD.Reduce(sbuf_residual_flow, rbuf_residual_flow, nVar_Flow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		if (turbulent)
			MPI::COMM_WORLD.Reduce(sbuf_residual_turbulent, rbuf_residual_turbulent, nVar_Turb, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		if (free_surface)
			MPI::COMM_WORLD.Reduce(sbuf_residual_levelset, rbuf_residual_levelset, nVar_LevelSet, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		if (fluid_structure)
			MPI::COMM_WORLD.Reduce(sbuf_residual_fea, rbuf_residual_fea, nVar_FEA, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		if (aeroacoustic)
			MPI::COMM_WORLD.Reduce(sbuf_residual_wave, rbuf_residual_wave, nVar_Wave, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
		MPI::COMM_WORLD.Reduce(sbuf_residual_plasma, rbuf_residual_plasma, nVar_Plasma, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case WAVE_EQUATION:
		MPI::COMM_WORLD.Reduce(sbuf_residual_wave, rbuf_residual_wave, nVar_Wave, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case LINEAR_ELASTICITY:
		MPI::COMM_WORLD.Reduce(sbuf_residual_fea, rbuf_residual_fea, nVar_FEA, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjflow, rbuf_residual_adjflow, nVar_AdjFlow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES : case ADJ_FREE_SURFACE_RANS :
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjflow, rbuf_residual_adjflow, nVar_AdjFlow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjlevelset, rbuf_residual_adjlevelset, nVar_AdjLevelSet, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case ADJ_AEROACOUSTIC_EULER:
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjflow, rbuf_residual_adjflow, nVar_AdjFlow, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		MPI::COMM_WORLD.Reduce(sbuf_residual_wave, rbuf_residual_wave, nVar_Wave, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES :
		MPI::COMM_WORLD.Reduce(sbuf_residual_adjplasma, rbuf_residual_adjplasma, nVar_AdjPlasma, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
		break;
	}

	MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, nVar_Force, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Reduce(sbuf_time, rbuf_time, nVar_Time, MPI::UNSIGNED_LONG, MPI::MAX, MASTER_NODE);
	MPI::COMM_WORLD.Reduce(sbuf_conv, rbuf_conv, nVar_Conv, MPI::UNSIGNED_SHORT, MPI::SUM, MASTER_NODE);
	MPI::COMM_WORLD.Barrier();

	/*-- Compute global convergence criteria --*/
	if (rank == MASTER_NODE) {
		for (iVar = 0; iVar < nVar_Flow; iVar++) rbuf_residual_flow[iVar] = sqrt(rbuf_residual_flow[iVar]);
		if (turbulent) for (iVar = 0; iVar < nVar_Turb; iVar++) rbuf_residual_turbulent[iVar] = sqrt(rbuf_residual_turbulent[iVar]);
		if (free_surface) for (iVar = 0; iVar < nVar_LevelSet; iVar++) rbuf_residual_levelset[iVar] = sqrt(rbuf_residual_levelset[iVar]);
		if (plasma) for (iVar = 0; iVar < nVar_Plasma; iVar++) rbuf_residual_plasma[iVar] = sqrt(rbuf_residual_plasma[iVar]);
		if (wave) for (iVar = 0; iVar < nVar_Wave; iVar++) rbuf_residual_wave[iVar] = sqrt(rbuf_residual_wave[iVar]);
		if (fea) for (iVar = 0; iVar < nVar_FEA; iVar++) rbuf_residual_fea[iVar] = sqrt(rbuf_residual_fea[iVar]);
		if (adjoint) {
			for (iVar = 0; iVar < nVar_AdjFlow; iVar++) rbuf_residual_adjflow[iVar] = sqrt(rbuf_residual_adjflow[iVar]);
			if (turbulent) for (iVar = 0; iVar < nVar_AdjTurb; iVar++) rbuf_residual_adjturbulent[iVar] = sqrt(rbuf_residual_adjturbulent[iVar]);
			if (free_surface) for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++) rbuf_residual_adjlevelset[iVar] = sqrt(rbuf_residual_adjlevelset[iVar]);
			if (plasma) for (iVar = 0; iVar < nVar_AdjPlasma; iVar++) rbuf_residual_adjplasma[iVar] = sqrt(rbuf_residual_adjplasma[iVar]);
		}
		if (rbuf_conv[0] == size) buf_convergence = 1;
		else buf_convergence = 0;
	}

	MPI::COMM_WORLD.Bcast(&buf_convergence, 1, MPI::UNSIGNED_SHORT, MASTER_NODE);

	switch (config->GetKind_Solver()) {
	case EULER: case NAVIER_STOKES: case RANS:
		if (buf_convergence == 1) integration[val_iZone][FLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][FLOW_SOL]->SetConvergence(false); break;
	case WAVE_EQUATION:
		if (buf_convergence == 1) integration[val_iZone][WAVE_SOL]->SetConvergence(true);
		else integration[val_iZone][WAVE_SOL]->SetConvergence(false); break;
	case LINEAR_ELASTICITY:
		if (buf_convergence == 1) integration[val_iZone][FEA_SOL]->SetConvergence(true);
		else integration[val_iZone][FEA_SOL]->SetConvergence(false); break;
	case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
		if (buf_convergence == 1) integration[val_iZone][PLASMA_SOL]->SetConvergence(true);
		else integration[val_iZone][PLASMA_SOL]->SetConvergence(false); break;
	case ELECTRIC_POTENTIAL:
		if (buf_convergence == 1) integration[val_iZone][ELEC_SOL]->SetConvergence(true);
		else integration[val_iZone][ELEC_SOL]->SetConvergence(false); break;
	case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		if (buf_convergence == 1) integration[val_iZone][FLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][FLOW_SOL]->SetConvergence(false); break;
	case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
		if (buf_convergence == 1) integration[val_iZone][FLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][FLOW_SOL]->SetConvergence(false); break;
	case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES:
		if (buf_convergence == 1) integration[val_iZone][FLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][FLOW_SOL]->SetConvergence(false); break;
	case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
		if (buf_convergence == 1) integration[val_iZone][ADJFLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][ADJFLOW_SOL]->SetConvergence(false); break;
	case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES:
		if (buf_convergence == 1) integration[val_iZone][ADJPLASMA_SOL]->SetConvergence(true);
		else integration[val_iZone][ADJPLASMA_SOL]->SetConvergence(false); break;
	case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
		if (buf_convergence == 1) integration[val_iZone][ADJFLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][ADJFLOW_SOL]->SetConvergence(false); break;
	case ADJ_AEROACOUSTIC_EULER:
		if (buf_convergence == 1) integration[val_iZone][ADJFLOW_SOL]->SetConvergence(true);
		else integration[val_iZone][ADJFLOW_SOL]->SetConvergence(false); break;
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

		if (free_surface) Total_CFreeSurface = rbuf_force[9];
		if (fluid_structure) Total_CFEA = rbuf_force[9];
		if (aeroacoustic) Total_CWave = rbuf_force[9];
		if (wave) Total_CWave = rbuf_force[9];
		if (fea) Total_CFEA = rbuf_force[9];

		/*--- Note that the nearfield based functionals should be redefined using the Cd weight ---*/
		if (equiv_area) {
			Total_CEquivArea  = config->GetWeightCd()*Total_CDrag + (1.0-config->GetWeightCd())*(rbuf_force[9]/double(size));
			Total_CNearFieldOF  = config->GetWeightCd()*Total_CDrag + (1.0-config->GetWeightCd())*rbuf_force[10];
		}

		if (rotating_frame) {
			Total_CMerit = rbuf_force[9];
			Total_CT = rbuf_force[10];
			Total_CQ = rbuf_force[11];
		}
    
    if (plasma) {
      Total_CQ = rbuf_force[9];
    }

		if (adjoint) {
			Total_Sens_Geo  = rbuf_force[0];
			Total_Sens_Mach = rbuf_force[1];
			Total_Sens_AoA  = rbuf_force[2];
			Total_Sens_Press = rbuf_force[3];
			Total_Sens_Temp  = rbuf_force[4];
		}

		timeused = rbuf_time[1];
		LinSolvIter = rbuf_time[2];

	}

#else

	/*--- Update the residual without MPI stuff ---*/
	for (iVar = 0; iVar < nVar_Flow; iVar++) rbuf_residual_flow[iVar] = sbuf_residual_flow[iVar];
	if (turbulent) for (iVar = 0; iVar < nVar_Turb; iVar++) rbuf_residual_turbulent[iVar] = sbuf_residual_turbulent[iVar];
	if (free_surface) for (iVar = 0; iVar < nVar_LevelSet; iVar++) rbuf_residual_levelset[iVar] = sbuf_residual_levelset[iVar];
	if (wave) for (iVar = 0; iVar < nVar_Wave; iVar++) rbuf_residual_wave[iVar] = sbuf_residual_wave[iVar];
	if (fea) for (iVar = 0; iVar < nVar_FEA; iVar++) rbuf_residual_fea[iVar] = sbuf_residual_fea[iVar];
	if (fluid_structure) for (iVar = 0; iVar < nVar_FEA; iVar++) rbuf_residual_fea[iVar] = sbuf_residual_fea[iVar];
	if (aeroacoustic) for (iVar = 0; iVar < nVar_Wave; iVar++) rbuf_residual_wave[iVar] = sbuf_residual_wave[iVar];
	if (plasma) for (iVar = 0; iVar < nVar_Plasma; iVar++) rbuf_residual_plasma[iVar] = sbuf_residual_plasma[iVar];

	if (adjoint) {
		for (iVar = 0; iVar < nVar_AdjFlow; iVar++) rbuf_residual_adjflow[iVar] = sbuf_residual_adjflow[iVar];
		if (turbulent) for (iVar = 0; iVar < nVar_AdjTurb; iVar++) rbuf_residual_adjturbulent[iVar] = sbuf_residual_adjturbulent[iVar];
		if (free_surface) for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++) rbuf_residual_adjlevelset[iVar] = sbuf_residual_adjlevelset[iVar];
		if (plasma) for (iVar = 0; iVar < nVar_AdjPlasma; iVar++) rbuf_residual_adjplasma[iVar] = sbuf_residual_adjplasma[iVar];
	}

	/*--- Note that the nearfield based functionals should be redefined using the Cd weight ---*/
	if (equiv_area) {
		Total_CEquivArea  = config->GetWeightCd()*Total_CDrag + (1.0-config->GetWeightCd())*Total_CEquivArea;
		Total_CNearFieldOF  = config->GetWeightCd()*Total_CDrag + (1.0-config->GetWeightCd())*Total_CNearFieldOF;
	}	

#endif	

	/*--- Output using all the processors ---*/
#ifndef NO_MPI
	MPI::COMM_WORLD.Barrier();
#endif	
	if (write_heads) {
		switch (config->GetKind_Solver()) {
		case EULER : case NAVIER_STOKES: case RANS:
			for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
				switch (config->GetMarker_All_Boundary(iMarker)) {
				case NACELLE_EXHAUST:
					cout << "Exhaust surface: "<< config->GetMarker_All_Tag(iMarker) << ". Mass flow rate: "
					<< solution_container[ZONE_0][MESH_0][FLOW_SOL]->GetMassFlow_Rate(iMarker) << "."<<endl;
					break;
				case NACELLE_INFLOW:
					cout << "Inflow surface: "<< config->GetMarker_All_Tag(iMarker) << ". Mass flow rate: "
					<< solution_container[ZONE_0][MESH_0][FLOW_SOL]->GetMassFlow_Rate(iMarker)
					<< ". Fan face Mach: " << solution_container[ZONE_0][MESH_0][FLOW_SOL]->GetFanFace_Mach(iMarker) <<"."<<endl;
					break;
				}
			}
			break;
		}
	}
#ifndef NO_MPI
	MPI::COMM_WORLD.Barrier();
#endif	

	if (rank == MASTER_NODE) {

		/*--- Write the begining of the history file ---*/
		sprintf (begin, "%12d", int(iExtIter));

		/*--- Write the end of the history file ---*/
		sprintf (end, ", %12.10f, %12.10f\n", double(LinSolvIter), double(timeused)/(CLOCKS_PER_SEC*60.0));
		/*--- Write the solution and residual of the history file ---*/
		switch (config->GetKind_Solver()) {

		case EULER : case NAVIER_STOKES: case RANS:
		case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: case FLUID_STRUCTURE_RANS:
		case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES: case AEROACOUSTIC_RANS:
		case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
		case ADJ_FREE_SURFACE_EULER: case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
		case ADJ_AEROACOUSTIC_EULER:
			sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
					Total_CFz, Total_CEff);
			if (equiv_area)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
						Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldOF);
			if (rotating_frame)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
						Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CMerit, Total_CT, Total_CQ);
			if (free_surface)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
						Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFreeSurface);
			if (fluid_structure)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
						Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFEA);
			if (aeroacoustic)
				sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
						Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CWave);

			if (nDim == 2) {
				if (incompressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), dummy, dummy );
				else sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), log10 (rbuf_residual_flow[3]), dummy );
			}
			else {
				if (incompressible) sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), log10 (rbuf_residual_flow[3]), dummy );
				else sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_flow[0]), log10 (rbuf_residual_flow[1]), log10 (rbuf_residual_flow[2]), log10 (rbuf_residual_flow[3]), log10 (rbuf_residual_flow[4]) );
			}

			if (turbulent){
				switch(nVar_Turb) {
				case 1: sprintf (turb_resid, ", %12.10f", log10 (rbuf_residual_turbulent[0])); break;
				case 2: sprintf (turb_resid, ", %12.10f\t%12.10f", log10(rbuf_residual_turbulent[0]), log10(rbuf_residual_turbulent[1])); break;
				}
			}

			if (free_surface) sprintf (levelset_resid, ", %12.10f", log10 (rbuf_residual_levelset[0]));

			if (fluid_structure) {
				if (nDim == 2) sprintf (levelset_resid, ", %12.10f, %12.10f, 0.0", log10 (rbuf_residual_fea[0]), log10 (rbuf_residual_fea[1]));
				else sprintf (levelset_resid, ", %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_fea[0]), log10 (rbuf_residual_fea[1]), log10 (rbuf_residual_fea[2]));
			}

			if (aeroacoustic) sprintf (levelset_resid, ", %12.10f, %12.10f", log10 (rbuf_residual_wave[0]), log10 (rbuf_residual_wave[1]));

			if (adjoint) {
				sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
				if (nDim == 2) sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (rbuf_residual_adjflow[0]),log10 (rbuf_residual_adjflow[1]),log10 (rbuf_residual_adjflow[2]),log10 (rbuf_residual_adjflow[3]) );
				else sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_adjflow[0]),log10 (rbuf_residual_adjflow[1]),log10 (rbuf_residual_adjflow[2]),log10 (rbuf_residual_adjflow[3]), log10 (rbuf_residual_adjflow[4]) );
				if (turbulent) sprintf (adj_turb_resid, ", %12.10f", log10 (rbuf_residual_adjturbulent[0]));
				if (free_surface) sprintf (adj_levelset_resid, ", %12.10f", log10 (rbuf_residual_adjlevelset[0]));
			}

			break;

		case PLASMA_EULER : case ADJ_PLASMA_EULER : case PLASMA_NAVIER_STOKES: case ADJ_PLASMA_NAVIER_STOKES:
			unsigned short iSpecies, loc;

			/*--- Direct problem coefficients ---*/
			sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
					Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
					Total_CFz, Total_CEff, Total_CQ);

			/*--- Direct problem residual ---*/
			for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
				if ( iSpecies < config->GetnDiatomics() ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*config->GetnDiatomics() + (nDim+2)*(iSpecies-config->GetnDiatomics());
				sprintf (resid_aux, ", %12.10f, %12.10f", log10 (rbuf_residual_plasma[loc+0]), log10 (rbuf_residual_plasma[loc+nDim+1]));
				if (iSpecies == 0) strcpy(plasma_resid, resid_aux);
				else strcat(plasma_resid, resid_aux);
			}

			/*--- Adjoint problem coefficients ---*/
			if (adjoint) {
				sprintf (adjoint_coeff, ", 0.0, 0.0, 0.0, 0.0");
				for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
					if ( iSpecies < config->GetnDiatomics() ) loc = (nDim+3)*iSpecies;
					else loc = (nDim+3)*config->GetnDiatomics() + (nDim+2)*(iSpecies-config->GetnDiatomics());
					sprintf (resid_aux, ", %12.10f, %12.10f", log10 (rbuf_residual_adjplasma[loc+0]),log10 (rbuf_residual_adjplasma[loc+nDim+1]));
					if (iSpecies == 0) strcpy(adj_plasma_resid, resid_aux);
					else strcat(adj_plasma_resid, resid_aux);
				}
			}

			break;

		case WAVE_EQUATION:

			sprintf (direct_coeff, ", %12.10f", Total_CWave);
			sprintf (wave_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_wave[0]), log10 (rbuf_residual_wave[1]), dummy, dummy, dummy );

			break;

		case LINEAR_ELASTICITY:

			sprintf (direct_coeff, ", %12.10f", Total_CFEA);
			sprintf (fea_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (rbuf_residual_fea[0]), dummy, dummy, dummy, dummy );

			break;

		}

		/*--- Write headers and some important other information ---*/
		if (write_heads) {

			switch (config->GetKind_Solver()) {
			case EULER : case NAVIER_STOKES: case FREE_SURFACE_EULER : case FREE_SURFACE_NAVIER_STOKES:
			case FLUID_STRUCTURE_EULER : case FLUID_STRUCTURE_NAVIER_STOKES: case AEROACOUSTIC_EULER : case AEROACOUSTIC_NAVIER_STOKES:
				cout << endl << " Iter" << "    Time(s)" << "     Res[Rho]";
				if (!fluid_structure && !aeroacoustic) {
					if (incompressible) cout << "   CLift(Total)" << "   CDrag(Total)" << endl;
					else if (rotating_frame && nDim == 3) cout << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
					else if (equiv_area) cout << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
					else cout << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				}
				else if (fluid_structure) cout << "   Res[Displx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
				else if (aeroacoustic) cout << "   Res[Wave]" << "   CLift(Total)" << "   CDrag(Total)" << "   CWave(Total)" << endl;
				break;

			case RANS : case FREE_SURFACE_RANS: case FLUID_STRUCTURE_RANS: case AEROACOUSTIC_RANS:
				switch (config->GetKind_Turb_Model()){
				case SA:	cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "       Res[nu]"; break;
				case SST:	cout << endl << " Iter" << "   Time(s)" << "      Res[Rho]" << "     Res[kine]" << "     Res[omega]"; break;
				}
				if (rotating_frame && nDim == 3 )
					cout << "   CThrust(Total)" << "   CTorque(Total)" << endl;
				else
					cout << "   CLift(Total)"   << "   CDrag(Total)"   << endl;
				break;

				case PLASMA_EULER : case PLASMA_NAVIER_STOKES :
					if (config->GetKind_GasModel() == ARGON)
						cout << endl << " Iter" << "   Time(s)" << "      Res[r1]" << "       Res[r2]" << "   Res(r3)"<<  endl;
					if (config->GetKind_GasModel() == AIR21)
						cout << endl << " Iter" << "   Time(s)" << "      Res[r1]" << "       Res[r2]" << "   Res(r3)"<<  endl;
					if ((config->GetKind_GasModel() == ARGON_SID) || (config->GetKind_GasModel() == AIR7) || (config->GetKind_GasModel() == AIR5) || (config->GetKind_GasModel() == N2) || (config->GetKind_GasModel() == O2))
						cout << endl << " Iter" << "   Time(s)" << "      Res[Rho0]" << "      Res[E0]" << "	   CQ(Total)" << "	  CDrag(Total)" << endl;
					break;

				case WAVE_EQUATION :
					cout << endl << " Iter" << "   Time(s)" << "      Res[Wave]" << "   CWave(Total)"<<  endl;
					break;

				case LINEAR_ELASTICITY :
					if (nDim == 2) cout << endl << " Iter" << "   Time(s)" << "    Res[Displx]" << "    Res[Disply]" << "   CFEA(Total)"<<  endl;
					if (nDim == 3) cout << endl << " Iter" << "   Time(s)" << "    Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "   CFEA(Total)"<<  endl;
					break;

				case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES :
					if ((config->GetKind_GasModel() == ARGON_SID) || (config->GetKind_GasModel() == AIR7) || (config->GetKind_GasModel() == AIR5) || (config->GetKind_GasModel() == N2) || (config->GetKind_GasModel() == O2))
						cout << endl << " Iter" << "   Time(s)" << "        Res[Psi_Rho0]" << "       Res[Psi_E0]" << "	      Sens_Geo" << endl;
					break;

				case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES : case ADJ_AEROACOUSTIC_EULER:
					if (incompressible) cout << endl << " Iter" << "   Time(s)" << "   Res[Psi_Rho]" << "     Sens_Geo" << "    Sens_Mach" << endl;
					else cout << endl << " Iter" << "   Time(s)" << "   Res[Psi_Rho]" << "     Res[Psi_E]" << "     Sens_Geo" << "    Sens_Mach" << endl;
					break;

				case ADJ_RANS : case ADJ_FREE_SURFACE_RANS :
					cout << endl << "  Iter" << "    Time(s)" << "     Res[Psi_Rho]" << "      Res[Psi_nu]" << "       Sens_Geo" << endl;
					break;
			}

		}

		/*--- Write the solution on the screen and history file ---*/
		if ((iExtIter != 0)
				|| ((config->IsAdjoint()) && (config->GetKind_Adjoint() == DISCRETE))) {
			switch (config->GetKind_Solver()) {

			case EULER : case NAVIER_STOKES:
			case FREE_SURFACE_EULER: case FREE_SURFACE_NAVIER_STOKES:
			case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
			case AEROACOUSTIC_EULER: case AEROACOUSTIC_NAVIER_STOKES:
				ConvHist_file[0] << begin << direct_coeff << flow_resid;
				if (free_surface) ConvHist_file[0] << levelset_resid;
				if (fluid_structure) ConvHist_file[0] << fea_resid;
				if (aeroacoustic) ConvHist_file[0] << levelset_resid;
				ConvHist_file[0] << end;

				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(11); cout << double(timeiter)/CLOCKS_PER_SEC;
				cout.width(13); cout << log10(rbuf_residual_flow[0]);
				if (!fluid_structure && !aeroacoustic && !equiv_area) {
					if (!incompressible ) {
						if (nDim == 2 ) { cout.width(14); cout << log10(rbuf_residual_flow[3]); }
						else { cout.width(14); cout << log10(rbuf_residual_flow[4]); }
					}
				}
				else if (fluid_structure) { cout.width(14); cout << log10(rbuf_residual_fea[0]); }
				else if (aeroacoustic) { cout.width(14); cout << log10(rbuf_residual_wave[0]); }

				if (rotating_frame && nDim == 3 ) { cout.width(15); cout << Total_CT; cout.width(15); cout << Total_CQ; }
				else if (aeroacoustic) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; cout.width(15); cout << Total_CWave; }
				else if (equiv_area) { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; cout.width(15);
				cout.precision(4);
				cout.setf(ios::scientific,ios::floatfield);
				cout << Total_CNearFieldOF; }
				else { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; }
				cout << endl;
				break;

			case RANS : case FREE_SURFACE_RANS:
				ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid << end;

				cout.precision(6);
				cout.setf(ios::fixed,ios::floatfield);
				cout.width(5); cout << iExtIter;
				cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
				cout.width(14); cout << log10(rbuf_residual_flow[0]);
				switch(nVar_Turb) {
				case 1: cout.width(14); cout << log10(rbuf_residual_turbulent[0]); break;
				case 2: cout.width(14); cout << log10(rbuf_residual_turbulent[0]);
				cout.width(14); cout << log10(rbuf_residual_turbulent[1]); break;
				}
				if (rotating_frame  && nDim == 3 ) { cout.width(15); cout << Total_CT; cout.width(15); cout << Total_CQ; }
				else { cout.width(15); cout << Total_CLift; cout.width(15); cout << Total_CDrag; }
				cout << endl;
				break;

				case PLASMA_EULER : case PLASMA_NAVIER_STOKES:
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
					if ((config->GetKind_GasModel() == ARGON_SID) || config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == O2 || config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {
						cout.width(14); cout << log10(rbuf_residual_plasma[0]);
						cout.width(14); cout << log10(rbuf_residual_plasma[nDim+1]);
						cout.width(14); cout << Total_CQ;
						cout.width(14); cout << Total_CDrag;
					}
					cout << endl;
					break;

				case WAVE_EQUATION:
					ConvHist_file[0] << begin << wave_coeff << wave_resid << end;

					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
					cout.width(14); cout << log10(rbuf_residual_wave[0]);
					cout.width(14); cout << Total_CWave;
					cout << endl;
					break;

				case LINEAR_ELASTICITY:
					ConvHist_file[0] << begin << fea_coeff << fea_resid << end;

					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
					cout.width(15); cout << log10(rbuf_residual_fea[0]);
					cout.width(15); cout << log10(rbuf_residual_fea[1]);
					if (nDim == 3) { cout.width(15); cout << log10(rbuf_residual_fea[2]); }
					cout.width(14); cout << Total_CFEA;
					cout << endl;
					break;

				case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES : case ADJ_AEROACOUSTIC_EULER:
					ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
					cout.width(15); cout << log10(rbuf_residual_adjflow[0]);
					if (!incompressible) {
						if (geometry[val_iZone][FinestMesh]->GetnDim() == 2 ) { cout.width(15); cout << log10(rbuf_residual_adjflow[3]); }
						else { cout.width(15); cout << log10(rbuf_residual_adjflow[4]); }
					}
					cout.precision(4);
					cout.setf(ios::scientific,ios::floatfield);
					cout.width(14); cout << Total_Sens_Geo;
					cout.width(14); cout << Total_Sens_Mach;
					cout << endl;
					break;

				case ADJ_RANS : case ADJ_FREE_SURFACE_RANS :
					ConvHist_file[0] << begin << adjoint_coeff << direct_coeff << adj_flow_resid << adj_turb_resid << end;
					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(6); cout << iExtIter;
					cout.width(11); cout << double(timeiter)/CLOCKS_PER_SEC;
					cout.width(17); cout << log10(rbuf_residual_adjflow[0]);
					cout.width(17); cout << log10(rbuf_residual_adjturbulent[0]);
					cout.width(16); cout << Total_Sens_Geo;
					cout << endl;
					break;

				case ADJ_PLASMA_EULER : case ADJ_PLASMA_NAVIER_STOKES:
					ConvHist_file[0] << begin << adjoint_coeff << adj_plasma_resid << end;

					cout.precision(6);
					cout.setf(ios::fixed,ios::floatfield);
					cout.width(5); cout << iExtIter;
					cout.width(10); cout << double(timeiter)/CLOCKS_PER_SEC;
					if ((config->GetKind_GasModel() == ARGON_SID) || config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == O2 || config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {
						cout.width(19); cout << log10(rbuf_residual_adjplasma[0]);
						cout.width(19); cout << log10(rbuf_residual_adjplasma[nDim+1]);
						cout.width(19); cout << Total_Sens_Geo;
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
	delete [] sbuf_residual_fea;

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
		delete [] rbuf_residual_fea;

		delete [] rbuf_residual_adjflow;
		delete [] rbuf_residual_adjlevelset;
		delete [] rbuf_residual_adjplasma;
		delete [] rbuf_residual_adjturbulent;

		delete [] rbuf_force;
		delete [] rbuf_time;
		delete [] rbuf_conv;

	}
}

void COutput::SetResult_Files(CSolution ****solution_container, CGeometry ***geometry, CConfig **config,
		unsigned long iExtIter, unsigned short val_nZone) {

  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	unsigned short iZone;

	for (iZone = 0; iZone < val_nZone; iZone++) {

		/*--- Flags identifying the types of files to be written. ---*/
		bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
		bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
		bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();

		switch (config[iZone]->GetKind_Solver()) {

		case EULER : case NAVIER_STOKES : case RANS :
		case FREE_SURFACE_EULER : case FREE_SURFACE_NAVIER_STOKES: case FREE_SURFACE_RANS:
		case FLUID_STRUCTURE_EULER : case FLUID_STRUCTURE_NAVIER_STOKES : case FLUID_STRUCTURE_RANS:
			if (Wrt_Srf && iZone == ZONE_0) SetSurface_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			if (Wrt_Csv && iZone == ZONE_0) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			break;

		case PLASMA_EULER : case PLASMA_NAVIER_STOKES :
			if (Wrt_Srf) SetSurface_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[0][MESH_0][PLASMA_SOL], iExtIter);
			break;

		case WAVE_EQUATION:
			if (Wrt_Srf) SetSurface_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][WAVE_SOL], iExtIter);
			break;

		case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS : case ADJ_FREE_SURFACE_EULER : case ADJ_FREE_SURFACE_NAVIER_STOKES: case ADJ_FREE_SURFACE_RANS:
			if (Wrt_Srf && iZone == ZONE_0) SetSurface_Adjoint(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][ADJFLOW_SOL], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][ADJFLOW_SOL], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
			break;

		case LIN_EULER : case LIN_NAVIER_STOKES :
			if (Wrt_Srf) SetSurface_Linearized(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][LINFLOW_SOL], config[iZone]->GetSurfLinCoeff_FileName(), iExtIter);
			if (Wrt_Csv) SetSurfaceCSV_Linearized(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][LINFLOW_SOL], config[iZone]->GetSurfLinCoeff_FileName(), iExtIter);
			break;

		case AEROACOUSTIC_EULER : case AEROACOUSTIC_NAVIER_STOKES : case AEROACOUSTIC_RANS:
			if (Wrt_Srf) SetSurface_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			break;
		case ADJ_AEROACOUSTIC_EULER:
			if (iZone == ZONE_0) {
				if (Wrt_Srf) SetSurface_Adjoint(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][ADJFLOW_SOL], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
				if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][ADJFLOW_SOL], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter, iZone);
			} else if (iZone == ZONE_1) {
				if (Wrt_Srf) SetSurface_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
				if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0][FLOW_SOL], iExtIter);
			}
			break;
		}
	
    /*--- Construction/Testing of new I/O routines. The conditional guards will
     come off soon and these will be the default options. ---*/
    if (Wrt_Vol) {
      /*--- Merge the node coordinates and connectivity. ---*/
      if (config[iZone]->GetWrt_Sol_Tec_ASCII()  || config[iZone]->GetWrt_Sol_Tec_Binary() || config[iZone]->GetWrt_Sol_CGNS())
        MergeGeometry(config[iZone], geometry[iZone][MESH_0], iZone);
    }
    
    if (config[iZone]->GetWrt_Sol_Tec_ASCII()  || config[iZone]->GetWrt_Sol_Tec_Binary() || config[iZone]->GetWrt_Sol_CGNS() || config[iZone]->GetWrt_Restart())
    MergeSolution(config[iZone], geometry[iZone][MESH_0], solution_container[iZone][MESH_0], iZone);
    
    /*--- Write restart, CGNS, or Tecplot files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- Write a native restart file ---*/
      if (config[iZone]->GetWrt_Restart())
        WriteRestart(config[iZone], geometry[iZone][MESH_0], iZone);
      
      if (Wrt_Vol) {
        /*--- Write a Tecplot ASCII file ---*/
        if (config[iZone]->GetWrt_Sol_Tec_ASCII())
          WriteTecplotASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone);
        
        /*--- Write a Tecplot binary file ---*/
        if (config[iZone]->GetWrt_Sol_Tec_Binary()) WriteTecplotBinary(config[iZone], geometry[iZone][MESH_0], iZone);
        
        /*--- Write a CGNS file ---*/
        if (config[iZone]->GetWrt_Sol_CGNS()) WriteCGNS(config[iZone], geometry[iZone][MESH_0], iZone);
      }
      
      /*--- Clean up memory after writing all output ---*/
      CleanUp(config[iZone], geometry[iZone][MESH_0]);
      
    }
    
  }
}

void COutput::SetFlowRate(CSolution *solution_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {

	ifstream index_file;
	int integration_node[100];
	int nPointFlowRate = 0;

	/*--- Read list of nodes ---*/
	index_file.open("flowrate_nodes.dat", ios::in);
	while (index_file.good()) {
		index_file >> integration_node[nPointFlowRate];
		nPointFlowRate++;
	}
	nPointFlowRate--;
	index_file.close();


	/*--- Perform trapezoid integration ---*/
	double y1,y2,q1,q2;
	double integral = 0;
	for (int j=0;j<nPointFlowRate-1; j++) {
		y1 = geometry->node[integration_node[j]]->GetCoord(1);
		y2 = geometry->node[integration_node[j+1]]->GetCoord(1);
		q1 = solution_container->node[integration_node[j]]->GetVelocity(0, config->GetIncompressible());
		q2 = solution_container->node[integration_node[j+1]]->GetVelocity(0, config->GetIncompressible());

		integral = integral + 0.5*(q1+q2)*(y2-y1);
	}

	/*--- Store integral in solution_container ---*/
	solution_container->SetTotal_CEquivArea(integral);  // integral shows up as EquivArea in history
}

void COutput::SetEquivalentArea(CSolution *solution_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {

	ofstream EquivArea_file, FuncGrad_file;
	unsigned short iMarker = 0, iDim;
	short *AzimuthalAngle = NULL;
	double Gamma, auxXCoord, auxYCoord, auxZCoord, InverseDesign, DeltaX, Coord_i, Coord_j, jp1Coord, *Coord = NULL, MeanFuntion, 
			*Face_Normal = NULL, auxArea, auxPress, Mach, Beta, R_Plane, Pressure_Inf, Density_Inf,
			RefAreaCoeff, ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL,
			*Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, *NearFieldWeight = NULL,
			*Weight = NULL, jFunction, jp1Function;
	unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint, 
			*IdPoint = NULL, *IdDomain = NULL, auxDomain;
	unsigned short iPhiAngle;
	ofstream NearFieldEA_file; ifstream TargetEA_file;

	bool incompressible = config->GetIncompressible();
	double XCoordBegin_OF = config->GetEA_IntLimit(0);
	double XCoordEnd_OF = config->GetEA_IntLimit(1);

	unsigned short nDim = geometry->GetnDim();
	double AoA = -(config->GetAoA()*PI_NUMBER/180.0);

	int rank = MESH_0;

	Mach  = config->GetMach_FreeStreamND();
	Gamma = config->GetGamma();
	Beta = sqrt(Mach*Mach-1.0);
	R_Plane = fabs(config->GetEA_IntLimit(2));
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

	/*--- Compute the total number of points on the near-field ---*/
	nVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();
				/*--- Using Face_Normal(z), and Coord(z) we identify only a surface, 
				 note that there are 2 NEARFIELD_BOUNDARY surfaces ---*/
				if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) nVertex_NearField ++;
			}

	/*--- Create an array with all the coordinates, points, pressures, face area,
	 equivalent area, and nearfield weight ---*/
	Xcoord = new double[nVertex_NearField]; 
	Ycoord = new double[nVertex_NearField];
	Zcoord = new double[nVertex_NearField];
	AzimuthalAngle = new short[nVertex_NearField];
	IdPoint = new unsigned long[nVertex_NearField];
	IdDomain = new unsigned long[nVertex_NearField];
	Pressure = new double[nVertex_NearField];
	FaceArea = new double[nVertex_NearField];
	EquivArea = new double[nVertex_NearField];
	TargetArea = new double[nVertex_NearField];
	NearFieldWeight = new double[nVertex_NearField];
	Weight = new double[nVertex_NearField];

	/*--- Copy the boundary information to an array ---*/
	nVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {

					IdPoint[nVertex_NearField] = iPoint;
					Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
					Ycoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(1);

					if (nDim ==2) {
						AzimuthalAngle[nVertex_NearField] = 0;
					}

					if (nDim == 3) {
						Zcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(2);

						/*--- Rotate the nearfield cylinder (AoA) only 3D ---*/
						double YcoordRot = Ycoord[nVertex_NearField];
						double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA); 

						/* Compute the Azimuthal angle (resolution of degress in the Azimuthal angle)---*/
						double AngleDouble; short AngleInt;
						AngleDouble = atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER;
						AngleInt = (short) floor(AngleDouble + 0.5);
						if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
						else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
					}

					if (AzimuthalAngle[nVertex_NearField] <= 60) {						
						Pressure[nVertex_NearField] = solution_container->node[iPoint]->GetPressure(incompressible);
						FaceArea[nVertex_NearField] = fabs(Face_Normal[nDim-1]);
						nVertex_NearField ++;
					}

				}
			}

#else

	int nProcessor = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
	int iProcessor;

	unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
	unsigned long *Buffer_Send_nVertex = new unsigned long [1];

	/*--- Compute the total number of points of the near-field ghost nodes ---*/
	nLocalVertex_NearField = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Coord = geometry->node[iPoint]->GetCoord();

				if (geometry->node[iPoint]->GetDomain())
					if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0))
						nLocalVertex_NearField ++;
			}

	Buffer_Send_nVertex[0] = nLocalVertex_NearField;

	/*--- Send Near-Field vertex information --*/
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI::UNSIGNED_LONG, MPI::MAX);
	MPI::COMM_WORLD.Allgather(Buffer_Send_nVertex, 1, MPI::UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI::UNSIGNED_LONG);

	double *Buffer_Send_Xcoord = new double[MaxLocalVertex_NearField];
	double *Buffer_Send_Ycoord = new double[MaxLocalVertex_NearField];
	double *Buffer_Send_Zcoord = new double[MaxLocalVertex_NearField];
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
		Buffer_Send_Ycoord[iVertex] = 0.0; Buffer_Send_Zcoord[iVertex] = 0.0;
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
					if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
						Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
						Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
						Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
						Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
						Buffer_Send_Pressure[nLocalVertex_NearField] = solution_container->node[iPoint]->GetPressure(incompressible);
						Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim-1]);
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
		AzimuthalAngle = new short[nVertex_NearField];
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

				if (nDim == 2) {
					AzimuthalAngle[nVertex_NearField] = 0;
				}

				if (nDim == 3) {
					Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField+iVertex];

					/*--- Rotate the nearfield cylinder  ---*/
					double YcoordRot = Ycoord[nVertex_NearField];
					double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA); 

					/*--- Compute the Azimuthal angle ---*/
					double AngleDouble; short AngleInt;
					AngleDouble = atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER;
					AngleInt = (short) floor(AngleDouble + 0.5);
					if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
					else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
				}

				if (AzimuthalAngle[nVertex_NearField] <= 60) {
					IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
					Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
					FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
					IdDomain[nVertex_NearField] = iProcessor;
					nVertex_NearField++;
				}

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

		vector<short> PhiAngleList;
		vector<short>::iterator IterPhiAngleList;

		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			PhiAngleList.push_back(AzimuthalAngle[iVertex]);

		sort( PhiAngleList.begin(), PhiAngleList.end());
		IterPhiAngleList = unique( PhiAngleList.begin(), PhiAngleList.end());
		PhiAngleList.resize( IterPhiAngleList - PhiAngleList.begin() );

		/*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
		vector<vector<double> > Xcoord_PhiAngle; Xcoord_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > Ycoord_PhiAngle; Ycoord_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > Zcoord_PhiAngle; Zcoord_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<unsigned long> > IdPoint_PhiAngle; IdPoint_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<unsigned long> > IdDomain_PhiAngle; IdDomain_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > Pressure_PhiAngle; Pressure_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > FaceArea_PhiAngle; FaceArea_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > EquivArea_PhiAngle; EquivArea_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > TargetArea_PhiAngle; TargetArea_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > NearFieldWeight_PhiAngle; NearFieldWeight_PhiAngle.resize(PhiAngleList.size()); 
		vector<vector<double> > Weight_PhiAngle; Weight_PhiAngle.resize(PhiAngleList.size());

		/*--- Distribute the values among the different PhiAngles ---*/
		for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				if (AzimuthalAngle[iVertex] == PhiAngleList[iPhiAngle]) {
					Xcoord_PhiAngle[iPhiAngle].push_back(Xcoord[iVertex]);
					Ycoord_PhiAngle[iPhiAngle].push_back(Ycoord[iVertex]);
					Zcoord_PhiAngle[iPhiAngle].push_back(Zcoord[iVertex]);
					IdPoint_PhiAngle[iPhiAngle].push_back(IdPoint[iVertex]);
					IdDomain_PhiAngle[iPhiAngle].push_back(IdDomain[iVertex]);
					Pressure_PhiAngle[iPhiAngle].push_back(Pressure[iVertex]);
					FaceArea_PhiAngle[iPhiAngle].push_back(FaceArea[iVertex]);
					EquivArea_PhiAngle[iPhiAngle].push_back(EquivArea[iVertex]);
					TargetArea_PhiAngle[iPhiAngle].push_back(TargetArea[iVertex]);
					NearFieldWeight_PhiAngle[iPhiAngle].push_back(NearFieldWeight[iVertex]);
					Weight_PhiAngle[iPhiAngle].push_back(Weight[iVertex]);
				}

		/*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++)
				for (jVertex = 0; jVertex < Xcoord_PhiAngle[iPhiAngle].size() - 1 - iVertex; jVertex++)
					if (Xcoord_PhiAngle[iPhiAngle][jVertex] > Xcoord_PhiAngle[iPhiAngle][jVertex+1]) {
						auxXCoord = Xcoord_PhiAngle[iPhiAngle][jVertex]; Xcoord_PhiAngle[iPhiAngle][jVertex] = Xcoord_PhiAngle[iPhiAngle][jVertex+1]; Xcoord_PhiAngle[iPhiAngle][jVertex+1] = auxXCoord;
						auxYCoord = Ycoord_PhiAngle[iPhiAngle][jVertex]; Ycoord_PhiAngle[iPhiAngle][jVertex] = Ycoord_PhiAngle[iPhiAngle][jVertex+1]; Ycoord_PhiAngle[iPhiAngle][jVertex+1] = auxYCoord;
						auxZCoord = Zcoord_PhiAngle[iPhiAngle][jVertex]; Zcoord_PhiAngle[iPhiAngle][jVertex] = Zcoord_PhiAngle[iPhiAngle][jVertex+1]; Zcoord_PhiAngle[iPhiAngle][jVertex+1] = auxZCoord;
						auxPress = Pressure_PhiAngle[iPhiAngle][jVertex]; Pressure_PhiAngle[iPhiAngle][jVertex] = Pressure_PhiAngle[iPhiAngle][jVertex+1]; Pressure_PhiAngle[iPhiAngle][jVertex+1] = auxPress;
						auxArea = FaceArea_PhiAngle[iPhiAngle][jVertex]; FaceArea_PhiAngle[iPhiAngle][jVertex] = FaceArea_PhiAngle[iPhiAngle][jVertex+1]; FaceArea_PhiAngle[iPhiAngle][jVertex+1] = auxArea;
						auxPoint = IdPoint_PhiAngle[iPhiAngle][jVertex]; IdPoint_PhiAngle[iPhiAngle][jVertex] = IdPoint_PhiAngle[iPhiAngle][jVertex+1]; IdPoint_PhiAngle[iPhiAngle][jVertex+1] = auxPoint;
						auxDomain = IdDomain_PhiAngle[iPhiAngle][jVertex]; IdDomain_PhiAngle[iPhiAngle][jVertex] = IdDomain_PhiAngle[iPhiAngle][jVertex+1]; IdDomain_PhiAngle[iPhiAngle][jVertex+1] = auxDomain;
					}


		/*--- Check that all the azimuth lists have the same size ---*/
		unsigned short nVertex = Xcoord_PhiAngle[0].size();
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			unsigned short nVertex_aux = Xcoord_PhiAngle[iPhiAngle].size();
			if (nVertex_aux != nVertex) cout <<"Be careful!!! one azimuth list is shorter than the other"<< endl;
			nVertex = min(nVertex, nVertex_aux);
		}

		/*--- Compute equivalent area distribution at each azimuth angle ---*/
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			EquivArea_PhiAngle[iPhiAngle][0] = 0.0;
			for (iVertex = 1; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
				EquivArea_PhiAngle[iPhiAngle][iVertex] = 0.0;

				Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][iVertex]*sin(AoA);

				for (jVertex = 0; jVertex < iVertex-1; jVertex++) {

					Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex]*sin(AoA);
					jp1Coord = Xcoord_PhiAngle[iPhiAngle][jVertex+1]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex+1]*sin(AoA);

					jFunction = factor*(Pressure_PhiAngle[iPhiAngle][jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
					jp1Function = factor*(Pressure_PhiAngle[iPhiAngle][jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);

					DeltaX = (jp1Coord-Coord_j);
					MeanFuntion = 0.5*(jp1Function + jFunction);
					EquivArea_PhiAngle[iPhiAngle][iVertex] += DeltaX * MeanFuntion;
				}
			}
		}

		/*--- Create a file with the equivalent area distribution at each azimuthal angle ---*/
		NearFieldEA_file.precision(15);
		NearFieldEA_file.open("NearFieldEA.plt", ios::out); 
		NearFieldEA_file << "TITLE = \"Nearfield Equivalent Area at each azimuthal angle \"" << endl;
		NearFieldEA_file << "VARIABLES = \"Coord (local to the near-field cylinder)\"";

		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			NearFieldEA_file << ", \"Equivalent Area, Phi= " << PhiAngleList[iPhiAngle] << " deg.\"";
		}

		NearFieldEA_file << endl;
		for (iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {
			double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
			double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);
			NearFieldEA_file << scientific << XcoordRot - XcoordRot_init;
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
				NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
			}
			NearFieldEA_file << endl;
		}
		NearFieldEA_file.close();

		/*--- Read target equivalent area from the configuration file, 
		 this first implementation requires a complete table (same as the original 
		 EA table). so... no interpolation. ---*/

		vector<vector<double> > TargetArea_PhiAngle_Trans;
		TargetEA_file.open("TargetEA.dat", ios::in);

		if (TargetEA_file.fail()) {
			if (iExtIter == 0) { cout << "There is no Target Equivalent Area file (TargetEA.dat)!!"<< endl;
			cout << "Using default parameters (Target Equiv Area = 0.0)" << endl;
			}
			/*--- Set the table to 0 ---*/
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				for (iVertex = 0; iVertex < TargetArea_PhiAngle[iPhiAngle].size(); iVertex++)
					TargetArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
		}
		else {
		
    		/*--- skip header lines ---*/
            string line;		
    		getline(TargetEA_file, line);
    		getline(TargetEA_file, line);
    		
			while (TargetEA_file) {
			
				string line;
				getline(TargetEA_file, line);
				istringstream is(line);
				vector<double> row;
				unsigned short iter = 0;
				
				while (is.good()) {
					string token;
					getline(is,token,',');

					istringstream js(token);
					
					double data;
                    js >> data;

					/*--- The first element in the table is the coordinate ---*/
					if (iter != 0) row.push_back(data); 
					iter++;
				}
				TargetArea_PhiAngle_Trans.push_back(row);
			}

			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++)
					TargetArea_PhiAngle[iPhiAngle][iVertex] = TargetArea_PhiAngle_Trans[iVertex][iPhiAngle];
		}

		/*--- Divide by the number of Phi angles in the nearfield ---*/
		double PhiFactor = 1.0/double(PhiAngleList.size());

		/*--- Evaluate the objective function ---*/
		InverseDesign = 0;
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
				Weight_PhiAngle[iPhiAngle][iVertex] = 1.0; 
				Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];

				double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex]-TargetArea_PhiAngle[iPhiAngle][iVertex];
				if ((Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;

				InverseDesign += PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*Difference*Difference;

			}

		/*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
				Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
				NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
				for (jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
					Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
					Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;

					double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex]-TargetArea_PhiAngle[iPhiAngle][jVertex];
					if ((Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;

					NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*2.0*Difference*factor*sqrt(Coord_j-Coord_i);
				}
			}		

		/*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/
		EquivArea_file.precision(15);
		EquivArea_file.open("nearfield_flow.plt", ios::out);
		EquivArea_file << "TITLE = \"SU2 Equivalent area computation at each azimuthal angle \"" << endl;
		EquivArea_file << "VARIABLES = \"Coord (local to the near-field cylinder)\",\"Equivalent area\",\"Target equivalent area\",\"NearField weight\",\"Pressure coefficient\"" << endl;

		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
			EquivArea_file << fixed << "ZONE T= \"Azimuthal angle " << PhiAngleList[iPhiAngle] << " deg.\"" << endl;
			for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {

				double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
				double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);

				EquivArea_file << scientific << XcoordRot - XcoordRot_init << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
				                                                                                                    << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << NearFieldWeight_PhiAngle[iPhiAngle][iVertex] << ", " <<
				                                                                                                    (Pressure_PhiAngle[iPhiAngle][iVertex]-Pressure_Inf)/Pressure_Inf << endl;
			}
		}

		EquivArea_file.close();

		/*--- Write Weight file for adjoint computation ---*/
		FuncGrad_file.precision(15);
		FuncGrad_file.open("WeightNF.dat", ios::out);

		FuncGrad_file << scientific << "-1.0";
		for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
			FuncGrad_file << scientific << "\t" << PhiAngleList[iPhiAngle];
		FuncGrad_file << endl;

		for (iVertex = 0; iVertex < NearFieldWeight_PhiAngle[0].size(); iVertex++) {
			double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
			FuncGrad_file << scientific << XcoordRot;
			for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
				FuncGrad_file << scientific << "\t" << NearFieldWeight_PhiAngle[iPhiAngle][iVertex];
			FuncGrad_file << endl;
		}
		FuncGrad_file.close();		

		/*--- Delete structures ---*/
		delete [] Xcoord; delete [] Ycoord; delete [] Zcoord; 
		delete [] AzimuthalAngle; delete [] IdPoint; delete [] IdDomain;
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

void COutput::SetFreeSurface(CSolution *solution_container, CGeometry *geometry,
		CConfig *config, unsigned long iExtIter) {
	double *coord = NULL, dist2, *iCoord = NULL, *jCoord = NULL, *U_i = NULL, *U_j = NULL, 
			**Coord_TargetLevelSet = NULL, **Coord_LevelSet = NULL, CoordTA_aux, TargetLevelSet_aux,
			*xCoord = NULL, *yCoord = NULL, auxCoordx, auxCoordy, FreeSurface, factor,
            volume, LevelSetDiff_Squared, LevelSetDiff, dist_Target, NumberSign;
	unsigned short iDim;
	unsigned long iPoint, jPoint, iVertex, jVertex, nVertex_LevelSet, iEdge, nPointTargetLevelSet, iVar;
	ifstream index_file;
    ofstream LevelSet_file;
	string text_line;
	int rank = MASTER_NODE;    
    char cstr[200], buffer[50];

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

	index_file.open("LevelSetTarget.dat", ios::in);
	if (index_file.fail()) {
		if ((iExtIter == 0) && (rank == MASTER_NODE)) { 
			cout << "There is no Target Level Set file (LevelSetTarget.dat)!!"<< endl;
			cout << "Using default parameters (Target Level Set = 0.0)" << endl;
		}
		nPointTargetLevelSet = 2;

		xCoord = new double [nPointTargetLevelSet];
		yCoord = new double [nPointTargetLevelSet];

		xCoord[0] = -10E6; yCoord[0] = config->GetFreeSurface_Zero();
		xCoord[1] = 10E6; yCoord[1] = config->GetFreeSurface_Zero();
	}
    
	else {

		/*--- Dimensionalization loop ---*/
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
		index_file.open("LevelSetTarget.dat", ios::in);
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

	for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
		for (iVar = 0; iVar < nPointTargetLevelSet-1; iVar++) {
			if ((xCoord[iVar] <= Coord_LevelSet[iVertex][0]) && (Coord_LevelSet[iVertex][0] <= xCoord[iVar+1])) {
				Coord_TargetLevelSet[iVertex][0] = Coord_LevelSet[iVertex][0];
				Coord_TargetLevelSet[iVertex][1] = yCoord[iVar] + (Coord_LevelSet[iVertex][0]-xCoord[iVar])*(yCoord[iVar+1]-yCoord[iVar] )/(xCoord[iVar+1]-xCoord[iVar]);
			}
		}
    }

	/*--- Get coordinates of the points and compute distances to the surface ---*/
	FreeSurface = 0.0, factor = 1.0;
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		coord = geometry->node[iPoint]->GetCoord();
		volume = geometry->node[iPoint]->GetVolume();
		LevelSetDiff_Squared = 0.0; LevelSetDiff = 0.0;
        
        /*--- Only in a limited region of the domain ---*/
		if ((coord[0] > 0.0) && (coord[0] < 5.0)){

			/*--- Compute the squared distance to the rest of points, and get the minimum ---*/
			dist_Target = 1E20;
			for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
				dist2 = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					dist2 += (coord[iDim]-Coord_TargetLevelSet[iVertex][iDim])*(coord[iDim]-Coord_TargetLevelSet[iVertex][iDim]);
				if (dist2 < fabs(dist_Target)) { 
					if (Coord_TargetLevelSet[iVertex][nDim-1] < coord[nDim-1]) dist_Target = dist2; 
					else dist_Target = -dist2; 
				}
			}
			NumberSign = 1.0;
			if (dist_Target != 0.0) NumberSign = dist_Target/fabs(dist_Target);
			dist_Target = sqrt(fabs(dist_Target))*NumberSign;
			LevelSetDiff = (solution_container->node[iPoint]->GetSolution()[0] - dist_Target);
			LevelSetDiff_Squared = LevelSetDiff*LevelSetDiff;
			FreeSurface += 0.5*factor*LevelSetDiff_Squared*volume;
		}
		solution_container->node[iPoint]->SetDiffLevelSet(LevelSetDiff);
	}

	if ((rank == MASTER_NODE) && (iExtIter % config->GetWrt_Sol_Freq() == 0)) {

		/*--- Write the Level Set distribution, the target level set---*/
		LevelSet_file.precision(15);

		/*--- Write file name with extension ---*/
		strcpy (cstr, "LevelSet");
		if (config->GetUnsteady_Simulation()){
			if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.dat", int(iExtIter));
			if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.dat", int(iExtIter));
			if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.dat", int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.dat", int(iExtIter));
			if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.dat", int(iExtIter));
		}
		else {
			sprintf (buffer, ".dat");
        }

		strcat(cstr,buffer);

		LevelSet_file.open(cstr, ios::out);	
		LevelSet_file << "TITLE = \"SU2 Free surface simulation\"" << endl;
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

//
///*!
// * \method WriteReactingOutputFile
// * \brief Output of the Plasma solver in Paraview Format
// * \author A. Lonkar
// */
//
//void COutput::WriteReactingOutputFile(CConfig *config, CGeometry *geometry,CSolution **solution_container, ofstream & ConsVar_file) {
//
//	unsigned short iPoint, iDim, nDim;
//	double M1, M2, M3;
//	double Vector[3];
//	ConsVar_file.precision(15);
//
//	M1 = config->GetMolar_Mass(0)/AVOGAD_CONSTANT;
//	M2 = config->GetMolar_Mass(1)/AVOGAD_CONSTANT;
//	M3 = config->GetMolar_Mass(2)/AVOGAD_CONSTANT;
//
//
//	nDim = geometry->GetnDim();
//
//
//
//	/*--- PRINT OUT ELECTRIC FIELD Vector ---*/
//
//#ifdef OutputSource
//	bool scalar = true;
//
//	if (config->GetElectricSolver()) {
//		/*--- PRINT OUT ELECTRIC POTENTIAL ---*/
//		WriteInOutputFile(geometry, solution_container[ELEC_SOL], ConsVar_file, "Electric_Potential", scalar, 0, "GetSolution", config);
//
//		ConsVar_file << "VECTORS Electric_Field float" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
//				Vector[iDim] = -1.0* solution_container[ELEC_SOL]->node[iPoint]->GetGradient(0,iDim);
//			}
//			if (geometry->GetnDim() == 2) ConsVar_file << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
//			if (geometry->GetnDim() == 3) ConsVar_file << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
//		}
//	}
//
//	/*--- PRINT OUT THE Total Density IN DOMAIN ---*/
//	ConsVar_file << "SCALARS Total_Density float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(0*(nDim+2)) + solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))+solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))  << endl;
//
//	/*--- PRINT OUT Total Pressure  ---*/
//	ConsVar_file << "SCALARS Total_Pressure float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(0) +solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(1)+solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(2) << endl;
//
//	/*--- PRINT OUT THE NEGATIVE CHARGE IN DOMAIN ---*/
//	ConsVar_file << "SCALARS Negative_Charge float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << ELECTRON_CHARGE*solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))/M3  << endl;
//
//	/*--- PRINT OUT THE POSITIVE CHARGE IN DOMAIN ---*/
//	ConsVar_file << "SCALARS Positive_Charge float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << ELECTRON_CHARGE*solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))/M2  << endl;
//
//
//
//	/*--- PRINT OUT THE NET CHARGE IN DOMAIN ---*/
//	ConsVar_file << "SCALARS Net_Charge float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << ELECTRON_CHARGE*(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))- solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))*M2/M3)/M2  << endl;
//
//
//	/*--- PRINT OUT THE degree of ionization  IN DOMAIN ---*/
//	ConsVar_file << "SCALARS alpha float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))/M3)/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))/M3 + solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(0*(nDim+2))/M1) << endl;
//
//
//#endif
//
//
//	/*--- PRINT OUT THE NET CHARGE IN DOMAIN ---*/
//	/*		ConsVar_file << "SCALARS Limiter float 1" << endl;
//	 ConsVar_file << "LOOKUP_TABLE default" << endl;
//	 for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//	 ConsVar_file << solution_container[PLASMA_SOL]->node[iPoint]->GetLimiter(0) << endl;
//	 */
//
//	/*--- PRINT OUT THE NET CHARGE IN DOMAIN ---*/
//	ConsVar_file << "SCALARS Net_Charge float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << ELECTRON_CHARGE*(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(1*(nDim+2))- solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(2*(nDim+2))*M2/M3)/M2  << endl;
//
//	unsigned short loc, iSpecies;
//
//	/* ***************************************** ---*/
//	/* SPECIES 1 DATA OUTPUT ---*/
//	/* ***************************************** ---*/
//
//	iSpecies = 0;
//	loc = (nDim+2)*iSpecies;
//
//	/*--- PRINT OUT THE DENSITY OF SPECIES 1 ---*/
//	//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Argon_Density", scalar, loc + 0, "GetSolution");
//
//	ConsVar_file << "SCALARS Argon_Density float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0)) << endl;
//
//	/*--- PRINT OUT DENSITY* Energy OF SPECIES 1 ---*/
//	//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species1", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);
//
//	/*--- PRINT OUT Velocity vector OF SPECIES 1 ---*/
//	ConsVar_file << "VECTORS Argon_Velocity float" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//		for (iDim = 0; iDim < nDim; iDim++)
//			Vector[iDim] = solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity(iDim, iSpecies);
//		if (nDim== 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
//		if (nDim == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
//	}
//
//	/*--- PRINT OUT Pressure OF SPECIES 1 ---*/
//	ConsVar_file << "SCALARS Argon_Pressure float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies) << endl;
//
//	/*--- PRINT OUT Mach number OF SPECIES 1 ---*/
//	ConsVar_file << "SCALARS Argon_Mach float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << fixed <<  sqrt(solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))/solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies) << endl;
//
//
//	/*--- PRINT OUT Temperature OF SPECIES 1 ---*/
//	ConsVar_file << "SCALARS Argon_Temperature float 1" << endl;
//	ConsVar_file << "LOOKUP_TABLE default" << endl;
//	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//		ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies))/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(0)*8314.0/(M1*AVOGAD_CONSTANT)) << endl;
//
//	int nSpecies = config->GetnSpecies();
//	if (nSpecies > 1) {
//		/* ***************************************** ---*/
//		/* SPECIES 2 DATA OUTPUT ---*/
//		/* ***************************************** ---*/
//
//		iSpecies = 1;
//		loc = (nDim+2)*iSpecies;
//
//		/*--- PRINT OUT THE DENSITY OF SPECIES 2 ---*/
//		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Ion_Density", scalar,loc + 0, "GetSolution");
//		ConsVar_file << "SCALARS Ion_Density float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0)) << endl;
//
//		/*--- PRINT OUT DENSITY* Energy OF SPECIES 2 ---*/
//		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species2", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);
//
//		/*--- PRINT OUT Velocity vector OF SPECIES 2 ---*/
//		ConsVar_file << "VECTORS Ion_Velocity float" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
//				Vector[iDim] = solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity(iDim, iSpecies);
//			if (nDim == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
//			if (nDim == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
//		}
//
//		/*--- PRINT OUT Pressure OF SPECIES 2 ---*/
//		ConsVar_file << "SCALARS Ion_Pressure float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies) << endl;
//
//		/*--- PRINT OUT Mach number OF SPECIES 2 ---*/
//		ConsVar_file << "SCALARS Ion_Mach float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed <<  sqrt(solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))/solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies) << endl;
//
//		/*--- PRINT OUT Temperature OF SPECIES 1 ---*/
//		ConsVar_file << "SCALARS Ion_Temperature float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies))/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc + 0)*8314.0/(M2*AVOGAD_CONSTANT)) << endl;
//
//	}
//
//	if (nSpecies > 2) {
//
//		/* ***************************************** ---*/
//		/* SPECIES 3 DATA OUTPUT ---*/
//		/* ***************************************** ---*/
//		iSpecies = 2;
//		loc = (nDim+2)*iSpecies;
//
//		/*--- PRINT OUT THE DENSITY OF SPECIES 3 ---*/
//		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Electron_Density", scalar, loc + 0, "GetSolution");
//		ConsVar_file << "SCALARS Electron_Density float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc+0)) << endl;
//
//		/*--- PRINT OUT DENSITY * Energy OF SPECIES 3 ---*/
//		//WriteInOutputFile(geometry, solution_container[PLASMA_SOL], ConsVar_file, "Density_x_Energy_Species3", scalar, loc + geometry->GetnDim()+1, "GetSolution", config);
//
//		/*--- PRINT OUT Velocity vector OF SPECIES 3 ---*/
//		ConsVar_file << "VECTORS Electron_Velocity float" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
//				Vector[iDim] = solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity(iDim, iSpecies);
//			if (geometry->GetnDim() == 2) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << "0.0" << endl;
//			if (geometry->GetnDim() == 3) ConsVar_file << fixed << Vector[0] << "\t" << Vector[1] << "\t" << Vector[2] << endl;
//		}
//
//		/*--- PRINT OUT Pressure OF SPECIES 3 ---*/
//		ConsVar_file << "SCALARS Electron_Pressure float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed << solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies) << endl;
//
//		/*--- PRINT OUT Temperature OF SPECIES 3 ---*/
//		ConsVar_file << "SCALARS Electron_Temperature float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed << (solution_container[PLASMA_SOL]->node[iPoint]->GetPressure(iSpecies))/(solution_container[PLASMA_SOL]->node[iPoint]->GetSolution(loc + 0)*8314.0/(M3*AVOGAD_CONSTANT)) << endl;
//
//
//		/*--- PRINT OUT Mach number OF SPECIES 3 ---*/
//		ConsVar_file << "SCALARS Electron_Mach float 1" << endl;
//		ConsVar_file << "LOOKUP_TABLE default" << endl;
//		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
//			ConsVar_file << fixed <<  sqrt(solution_container[PLASMA_SOL]->node[iPoint]->GetVelocity2(iSpecies))/solution_container[PLASMA_SOL]->node[iPoint]->GetSoundSpeed(iSpecies) << endl;
//
//	}
//
//
//
//
//}

