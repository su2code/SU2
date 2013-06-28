/*!
 * \file SU2_CFD.cpp
 * \brief Main file of Computational Fluid Dynamics Code (SU2_CFD).
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

#include "../include/SU2_CFD.hpp"

int main(int argc, char *argv[]) {
	
	/*--- Variable definitions ---*/
	bool StopCalc = false;
	unsigned long StartTime, StopTime, TimeUsed = 0, ExtIter = 0, IntIter = 0;
	double Physical_dt, Physical_t;
	unsigned short iMesh;
	ofstream ConvHist_file;
	char file_name[200];
	int rank = MASTER_NODE;
	
#ifndef NO_MPI
	/*--- MPI initialization, and buffer setting ---*/
	void *buffer, *old_buffer;
	int size, bufsize;
	bufsize = MAX_MPI_BUFFER;
	buffer = new char[bufsize];
	MPI::Init(argc,argv);
	MPI::Attach_buffer(buffer, bufsize);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
#endif
	
	/*--- Pointer to different structures, that will be used 
	 in the entire code ---*/
	CConfig *config;
	CIntegration **integration_container;
	COutput *output;
	CGeometry **geometry;
	CSolution ***solution_container;
	CNumerics ****solver_container;
	
	/*--- Definition of the configuration class, and open the config file ---*/
	if (argc == 2) config = new CConfig(argv[1], SU2_CFD);
	else {
		strcpy (file_name, "default.cfg");
		config = new CConfig(file_name, SU2_CFD);
	}
	
	/*--- Definition of the output class ---*/
	output = new COutput(config);
	
#ifndef NO_MPI
	if (rank == MASTER_NODE) {
		/*--- Open the convergence history file ---*/
		output->SetHistory_file(&ConvHist_file, config);
	}
	else {
		/*--- Change the name of the input-output files for the 
		 parallel computation ---*/
		config->SetFileNameDomain(rank);
	}
#else
	/*--- Open the convergence history file ---*/
	output->SetHistory_file(&ConvHist_file, config);
#endif
	
	/*--- Definition of the geometry class and open the mesh file ---*/
	geometry = new CGeometry *[config->GetMGLevels()+1];
	geometry[MESH_0] = new CPhysicalGeometry(config, config->GetMesh_FileName(), 
																					 config->GetMesh_FileFormat());
	
	/*--- Perform the non-dimensionalization for the flow equations ---*/
	config->SetNondimensionalization(geometry[MESH_0]->GetnDim(), rank);
	
#ifndef NO_MPI
  if (rank == MASTER_NODE) {
		/*--- Set the master node to the NO_SOLVER type. ---*/
		config->SetKind_Solver(NO_SOLVER);
	}
  
	/*--- Synchronization point after reading the grid ---*/
	MPI::COMM_WORLD.Barrier();
#endif
	
	/*--- Definition of the geometry (edge based structure) ---*/
	Geometrical_Definition(geometry, config);
	
	/*--- Definition of the solution class (solution_container[#GRIDS][#EQ_SYSTEMS]) ---*/
	solution_container = new CSolution** [config->GetMGLevels()+1];
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
		solution_container[iMesh] = new CSolution* [MAX_SOLS];
	
	/*--- Definition of the numerical method class (solver_container[#GRIDS][#EQ_SYSTEMS][#EQ_TERMS]) ---*/
	solver_container = new CNumerics*** [config->GetMGLevels()+1];
	
	/*--- Definition of the integration class (integration_container[#EQ_SYSTEMS]) ---*/
	integration_container = new CIntegration* [MAX_SOLS];
	
	/*--- Solver and integration definition ---*/
	Solver_Definition(solver_container, solution_container, integration_container, geometry, config);
	
	/*--- Computation of wall distance ---*/
	if ((config->GetMPI_Kind_Solver() == RANS) || (config->GetMPI_Kind_Solver() == ADJ_RANS))
		geometry[MESH_0]->SetWall_Distance(config);
	
	/*--- Computation of positive area in the z-plane ---*/
	geometry[MESH_0]->SetPositive_ZArea(config);
	
	/*--- Store the wall distance as the Eikonal solution ---*/
	if ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS))
		integration_container[EIKONAL_SOL]->SetEikonal_Solver(geometry, solution_container, config);	
	
	/*--- Set the near-field boundary conditions (be careful there is some 
	 MPI stuff in these subroutines) ---*/
	for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
#ifndef NO_MPI
		if ((rank != MASTER_NODE) || (iMesh == MESH_0))
#endif
		{ 
			geometry[iMesh]->MachNearField(config);
			geometry[iMesh]->MachInterface(config); 
		}
	
	
	/*--- External loop of the solver ---*/
	if (rank == MASTER_NODE) 
		cout << endl <<"------------------------------ Start solver -----------------------------" << endl;
	
	
	while (ExtIter < config->GetnExtIter()) {
		
		StartTime = clock();
		config->UpdateCFL(ExtIter);
		
		switch (config->GetKind_Solver()) {
			case NO_SOLVER:
				/*--- No solution (parallel mode) ---*/
				break;
				
			case POTENTIAL_FLOW:
				/*--- Potential flow equation ---*/
				config->SetGlobalParam(POTENTIAL_FLOW, RUNTIME_POT_SYS);
				integration_container[FLOW_SOL]->SetPotential_Solver(geometry, solution_container, solver_container,
																														 config, RUNTIME_POT_SYS, MESH_0);
				config->SetnExtIter(0);
				break;
				
			case EULER:
				/*--- Set the initial condition at the first iteration ---*/
				integration_container[FLOW_SOL]->SetInitialCondition(geometry, solution_container, solver_container, 
																														 config, ExtIter);
				
				/*--- Euler equations ---*/
				config->SetGlobalParam(EULER, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																														 config, RUNTIME_FLOW_SYS, ExtIter);
				
				/*--- Dual time stepping strategy ---*/
				if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
					for(IntIter = 1; IntIter < config->GetnIntIter(); IntIter++) {
						integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																																 config, RUNTIME_FLOW_SYS, IntIter);
						output->SetHistoryDT_Serial(geometry, config, solution_container, IntIter);
						if (integration_container[FLOW_SOL]->GetConvergence()) break;
					}
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[FLOW_SOL]->SetDualTime_Solver(geometry[iMesh], solution_container[iMesh][FLOW_SOL]);
						integration_container[FLOW_SOL]->SetConvergence(false);
					}
					Physical_dt = config->GetUnst_TimeStep();
					Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time())
						integration_container[FLOW_SOL]->SetConvergence(true);
				}
				
				break;
				
			case NAVIER_STOKES:
				/*--- Navier-Stokes equations ---*/
				config->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																														 config, RUNTIME_FLOW_SYS, ExtIter);
				
				/*--- Dual time stepping strategy ---*/
				if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
					for(IntIter = 1; IntIter < config->GetnIntIter(); IntIter++) {
						integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																																 config, RUNTIME_FLOW_SYS, IntIter);
						output->SetHistoryDT_Serial(geometry, config, solution_container, IntIter);
						if (integration_container[FLOW_SOL]->GetConvergence()) break;
					}
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[FLOW_SOL]->SetDualTime_Solver(geometry[iMesh], solution_container[iMesh][FLOW_SOL]);
						integration_container[FLOW_SOL]->SetConvergence(false);
					}
					Physical_dt = config->GetUnst_TimeStep();
					Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time())
						integration_container[FLOW_SOL]->SetConvergence(true);
				}
				
				break;
				
			case RANS:
				/*--- Reynolds-averaged Navier-Stokes equations ---*/
				config->SetGlobalParam(RANS, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																														 config, RUNTIME_FLOW_SYS, ExtIter);
				
				/*--- Turbulent model solution ---*/
				config->SetGlobalParam(RANS, RUNTIME_TURB_SYS);
				integration_container[TURB_SOL]->SetSingleGrid_Solver(geometry, solution_container, solver_container,
																															config, RUNTIME_TURB_SYS, ExtIter);
				
				/*--- The turbulent solution is only computed on the finest mesh, so it is necessary
				 to project the solution into higher mesh levels ---*/
				for (iMesh = 0; iMesh < config->GetMGLevels(); iMesh++) {
					integration_container[TURB_SOL]->SetRestricted_Solution(RUNTIME_TURB_SYS, solution_container[iMesh],
																																	solution_container[iMesh+1], geometry[iMesh],
																																	geometry[iMesh+1], config, iMesh, false);
					integration_container[TURB_SOL]->SetRestricted_Gradient(RUNTIME_TURB_SYS, solution_container[iMesh],
																																	solution_container[iMesh+1], geometry[iMesh],
																																	geometry[iMesh+1], config);
				}
				break;
				
			case ADJ_EULER:
				/*--- Continuous adjoint Euler equations ---*/
				if (ExtIter == 0) {
					cout << "Iteration over the direct problem to store all flow information." << endl;
					config->SetGlobalParam(ADJ_EULER, RUNTIME_FLOW_SYS);
					integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, 
																															 config, RUNTIME_FLOW_SYS, ExtIter);
					
					/*--- Compute gradiens of the flow variables, this is necessary
					 for sensitivity computation, note that in the direct problem 
					 we are not computing the gradients ---*/	
					if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
						solution_container[MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry[MESH_0], config);
					if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
							(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES))
						solution_container[MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry[MESH_0], config);
					
					/*--- Set force projection vector for boundary conditions ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						solution_container[iMesh][FLOW_SOL]->SetTotal_CDrag(solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag());
						solution_container[iMesh][FLOW_SOL]->SetTotal_CLift(solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift());
						solution_container[iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry[iMesh], solution_container[iMesh], config);
						if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) ||
								(config->GetKind_ObjFunc() == NEARFIELD_PRESSURE))
							solution_container[iMesh][ADJFLOW_SOL]->SetIntBoundary_Jump(geometry[iMesh], solution_container[iMesh], config);
					}
					
				}
				config->SetGlobalParam(ADJ_EULER, RUNTIME_ADJFLOW_SYS);
				integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_ADJFLOW_SYS, ExtIter);
				break;
				
			case ADJ_NAVIER_STOKES:
				/*--- Continuous adjoint Navier-Stokes equations (using frozen viscosity) ---*/
				if (ExtIter == 0) {
					cout << "Iteration over the direct problem to store all flow information." << endl;
					config->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_FLOW_SYS);
					integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																															 RUNTIME_FLOW_SYS, ExtIter);
					
					/*--- Set force projection vector for adjoint boundary conditions ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						solution_container[iMesh][FLOW_SOL]->SetTotal_CDrag(solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag());
						solution_container[iMesh][FLOW_SOL]->SetTotal_CLift(solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift());
						solution_container[iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry[iMesh], solution_container[iMesh], config);
					}
				}
				
				config->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS);
				integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_ADJFLOW_SYS, ExtIter);
				break;
				
			case ADJ_RANS:
				/*--- Continuous adjoint RANS equations (using adjoint turbulence model) ---*/
				if (ExtIter == 0) {
					cout << "Computation of distance field" << endl;
					integration_container[EIKONAL_SOL]->SetEikonal_Solver(geometry, solution_container, config);
					
					cout << "Iteration over the direct problem to store all flow information." << endl;
					config->SetGlobalParam(ADJ_RANS, RUNTIME_FLOW_SYS);
					integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																															 RUNTIME_FLOW_SYS, ExtIter);
					
					/*--- Compute gradiens of flow variables (only on the finest grid) ---*/
					if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
						solution_container[MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry[MESH_0], config);
						solution_container[MESH_0][TURB_SOL]->SetSolution_Gradient_GG(geometry[MESH_0]);
					}
					if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
							(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)) {
						solution_container[MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry[MESH_0], config);
						solution_container[MESH_0][TURB_SOL]->SetSolution_Gradient_LS(geometry[MESH_0], config);
					}
					
					/*--- Set force projection vector for boundary conditions ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						solution_container[iMesh][FLOW_SOL]->SetTotal_CDrag(solution_container[MESH_0][FLOW_SOL]->GetTotal_CDrag());
						solution_container[iMesh][FLOW_SOL]->SetTotal_CLift(solution_container[MESH_0][FLOW_SOL]->GetTotal_CLift());
						solution_container[iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry[iMesh], solution_container[iMesh], config);
					}
				}
				
				/*--- Iteration over the flow-adjoint equations ---*/
				config->SetGlobalParam(ADJ_RANS, RUNTIME_ADJFLOW_SYS);
				integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_ADJFLOW_SYS, ExtIter);
				
				/*--- Iteration over the turbulent-adjoint equations ---*/
				config->SetGlobalParam(ADJ_RANS, RUNTIME_ADJTURB_SYS);
				integration_container[ADJTURB_SOL]->SetSingleGrid_Solver(geometry, solution_container, solver_container, config,
																																 RUNTIME_ADJTURB_SYS, ExtIter);
				
				/*--- The turbulent solution is only computed in the finest mesh,
				 projection of the solution to higher mesh levels ---*/
				for (iMesh = 0; iMesh < config->GetMGLevels(); iMesh++) {
					integration_container[ADJTURB_SOL]->SetRestricted_Solution(RUNTIME_ADJTURB_SYS, solution_container[iMesh],
																																		 solution_container[iMesh+1], geometry[iMesh],
																																		 geometry[iMesh+1], config, iMesh, false);
					
					integration_container[ADJTURB_SOL]->SetRestricted_Gradient(RUNTIME_ADJTURB_SYS, solution_container[iMesh],
																																		 solution_container[iMesh+1], geometry[iMesh],
																																		 geometry[iMesh+1], config);
				}
				break;
				
			case LIN_EULER:
				/*--- Linearized Euler equations ---*/
				if (ExtIter == 0) {
					cout << "Iteration over the flow variables to store all the information." << endl;
					config->SetGlobalParam(LIN_EULER, RUNTIME_FLOW_SYS);
					integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																															 RUNTIME_FLOW_SYS, ExtIter);
					
					/*--- Compute gradiens of flow variables (only on the finest grid) ---*/
					if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
						solution_container[MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry[MESH_0], config);
					if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
							(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES))
						solution_container[MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry[MESH_0], config);
				}
				
				config->SetGlobalParam(LIN_EULER, RUNTIME_LINFLOW_SYS);
				integration_container[LINFLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_LINFLOW_SYS, ExtIter);
				break;
				
			case ELECTRIC_POTENTIAL:
				/*--- Potential flow equation ---*/
				config->SetGlobalParam(ELECTRIC_POTENTIAL, RUNTIME_ELEC_SYS);
				integration_container[ELEC_SOL]->SetPotential_Solver(geometry, solution_container, solver_container, config,
																														 RUNTIME_ELEC_SYS, MESH_0);
				config->SetnExtIter(0);
				break;
				
			case ADJ_ELECTRIC_POTENTIAL:
				/*--- Continuous adjoint electric potential equation ---*/
				config->SetGlobalParam(ADJ_ELECTRIC_POTENTIAL, RUNTIME_ADJELEC_SYS);
				integration_container[ADJELEC_SOL]->SetPotential_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_ADJELEC_SYS, MESH_0);
				config->SetnExtIter(0);
				break;
				
			case LIN_ELECTRIC_POTENTIAL:
				/*--- Linearized electric potential equation ---*/
				config->SetGlobalParam(LIN_ELECTRIC_POTENTIAL, RUNTIME_LINELEC_SYS);
				integration_container[LINELEC_SOL]->SetPotential_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_LINELEC_SYS, MESH_0);
				config->SetnExtIter(0);
				break;
				
			case ONE_SHOT_ADJ_EULER:
				/*--- One Shot continuous adjoint Euler equations ---*/
				config->SetGlobalParam(ONE_SHOT_ADJ_EULER, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																														 RUNTIME_FLOW_SYS, ExtIter);
				
				config->SetGlobalParam(ONE_SHOT_ADJ_EULER, RUNTIME_FLOW_ADJFLOW_SYS);
				integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_FLOW_ADJFLOW_SYS, ExtIter);
				break;
				
			case ONE_SHOT_ADJ_NAVIER_STOKES:
				/*--- One Shot continuous adjoint Navier-Stokes equations ---*/
				config->SetGlobalParam(ONE_SHOT_ADJ_NAVIER_STOKES, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																														 RUNTIME_FLOW_SYS, ExtIter);
				
				config->SetGlobalParam(ONE_SHOT_ADJ_NAVIER_STOKES, RUNTIME_FLOW_ADJFLOW_SYS);
				integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, config,
																																RUNTIME_FLOW_ADJFLOW_SYS, ExtIter);
				break;
				
			case MULTI_SPECIES_NAVIER_STOKES:
				/*--- Plasma equations ---*/
				config->SetGlobalParam(MULTI_SPECIES_NAVIER_STOKES, RUNTIME_PLASMA_SYS);
				integration_container[PLASMA_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																															 config, RUNTIME_PLASMA_SYS, ExtIter);
				
				/*--- Dual time stepping strategy ---*/
				if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
					for(IntIter = 1; IntIter < config->GetnIntIter(); IntIter++) {
						integration_container[PLASMA_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																																	 config, RUNTIME_PLASMA_SYS, IntIter);
						output->SetHistoryDT_Serial(geometry, config, solution_container, IntIter);
						if (integration_container[PLASMA_SOL]->GetConvergence()) break;
					}
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[PLASMA_SOL]->SetDualTime_Solver(geometry[iMesh], solution_container[iMesh][FLOW_SOL]);
						integration_container[PLASMA_SOL]->SetConvergence(false);
					}
				}

			/*--- Electric potential equation ---*/
			config->SetGlobalParam(MULTI_SPECIES_NAVIER_STOKES, RUNTIME_ELEC_SYS);
			integration_container[ELEC_SOL]->SetPotential_Solver(geometry, solution_container, solver_container, config,
					RUNTIME_ELEC_SYS, MESH_0);
			break;
	
			case TWO_PHASE_FLOW:
				/*--- Update density and viscosity, using level set ---*/
				integration_container[FLOW_SOL]->SetTwoPhase_Solver(geometry, solution_container, solver_container, 
																														config, ExtIter);
				
				/*--- Navier-Stokes equations ---*/	
				config->SetGlobalParam(TWO_PHASE_FLOW, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container, 
																														 config, RUNTIME_FLOW_SYS, ExtIter);
		    
				/*--- Dual time stepping strategy ---*/
				if (config->GetUnsteady_Simulation() == DUAL_TIME_STEPPING) {
					for(IntIter = 1; IntIter < config->GetnIntIter(); IntIter++) {
						config->SetGlobalParam(TWO_PHASE_FLOW, RUNTIME_FLOW_SYS);
						integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																																 config, RUNTIME_FLOW_SYS, IntIter);
						output->SetHistoryDT_Serial(geometry, config, solution_container, IntIter);
						if (integration_container[FLOW_SOL]->GetConvergence()) break;
					}
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[FLOW_SOL]->SetDualTime_Solver(geometry[iMesh], solution_container[iMesh][FLOW_SOL]);
						integration_container[FLOW_SOL]->SetConvergence(false);
					}
					Physical_dt = config->GetUnst_TimeStep();
					Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time())
						integration_container[FLOW_SOL]->SetConvergence(true);
				}
				
				/*--- Level-Set model solution ---*/
/*				if (config->GetFreeSurface_Movement()) {
					config->SetGlobalParam(TWO_PHASE_FLOW, RUNTIME_LEVELSET_SYS);
					integration_container[LEVELSET_SOL]->SetSingleGrid_Solver(geometry, solution_container, solver_container,
																																		config, RUNTIME_LEVELSET_SYS, ExtIter);
				}*/
				
				break;
				
				
			case COMBUSTION:
				/*--- Euler equations ---*/
				config->SetGlobalParam(COMBUSTION, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry, solution_container, solver_container,
																														 config, RUNTIME_FLOW_SYS, ExtIter);
				
				/*--- Combustion equations ---*/
				config->SetGlobalParam(COMBUSTION, RUNTIME_COMBUSTION_SYS);
				integration_container[COMBUSTION_SOL]->SetSingleGrid_Solver(geometry, solution_container, solver_container,
																																		config, RUNTIME_COMBUSTION_SYS, ExtIter);
				
				break;
				
		}
		
#ifndef NO_MPI
		MPI::COMM_WORLD.Barrier();
#endif
		StopTime = clock(); TimeUsed += (StopTime - StartTime);
		
		/*--- Convergence history for serial and parallel computation ---*/	
		if ( ExtIter % config->GetWrt_Con_Freq() == 0 ) {
			/*--- Evaluate and plot the equivalent area ---*/			
			if ((config->GetMPI_Kind_Solver() == EULER) && (config->GetEquivArea() == YES))
				output->SetEquivalentArea(solution_container[MESH_0][FLOW_SOL], geometry[MESH_0], config, ExtIter);
			
#ifdef NO_MPI
			output->SetHistory_Serial(&ConvHist_file, geometry, config, solution_container, 
																ExtIter, TimeUsed);
#else
			output->SetHistory_Parallel(&ConvHist_file, geometry, solution_container, config, 
																	integration_container, ExtIter, rank, size, TimeUsed);
#endif
		}
		/*--- Review convergence criteria ---*/
		switch (config->GetKind_Solver()) {
			case NO_SOLVER:
				StopCalc = integration_container[NO_SOLVER]->GetConvergence(); break;
			case POTENTIAL_FLOW: case EULER: case NAVIER_STOKES: case RANS: case TWO_PHASE_FLOW:
				StopCalc = integration_container[FLOW_SOL]->GetConvergence(); break;
			case MULTI_SPECIES_NAVIER_STOKES:
				StopCalc = integration_container[PLASMA_SOL]->GetConvergence(); break;
			case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
				StopCalc = integration_container[ADJFLOW_SOL]->GetConvergence(); break;
			case LIN_EULER: case LIN_NAVIER_STOKES:
				StopCalc = integration_container[LINFLOW_SOL]->GetConvergence(); break;
			case ELECTRIC_POTENTIAL:
				StopCalc = integration_container[ELEC_SOL]->GetConvergence(); break;
			case ADJ_ELECTRIC_POTENTIAL:
				StopCalc = integration_container[ADJELEC_SOL]->GetConvergence(); break;
			case ONE_SHOT_ADJ_EULER: case ONE_SHOT_ADJ_NAVIER_STOKES:
				StopCalc = ( integration_container[FLOW_SOL]->GetConvergence() &&
										integration_container[ADJFLOW_SOL]->GetConvergence()); break;
		}
		
		
		/*--- Solution output (be careful there is some MPI stuff in these subroutines) ---*/
		if ((((ExtIter+1 == config->GetnExtIter()) || 
					(ExtIter % config->GetWrt_Sol_Freq() == 0) || StopCalc ) && (ExtIter != 0)
				 ) || (config->GetnExtIter() == 1)) {
			output->SetResult_Files(solution_container, geometry[MESH_0], config, ExtIter);
		}
		
		/*--- Stop criteria	---*/	
		if (StopCalc) break;
		
		ExtIter++;
	}
	ConvHist_file.close();
	
	/*--- Memory deallocation ---*/
	//Solver_Deallocation(solver_container, solution_container, integration_container, output, geometry, config);
	//Geometrical_Deallocation(geometry, config);

#ifndef NO_MPI
	/*--- Finalize MPI parallelization ---*/	
	old_buffer = buffer;
	MPI::Detach_buffer(old_buffer);
	//	delete [] buffer;
	MPI::Finalize();
#endif
	
	return 0;
}
