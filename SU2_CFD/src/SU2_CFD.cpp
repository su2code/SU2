/*!
 * \file SU2_CFD.cpp
 * \brief Main file of Computational Fluid Dynamics Code (SU2_CFD).
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

#include "../include/SU2_CFD.hpp"
using namespace std;

int main(int argc, char *argv[]) {
	/*--- Variable definitions ---*/
	bool StopCalc = false;
	unsigned long StartTime, StopTime, TimeUsed = 0, ExtIter = 0, IntIter = 0;
	double Physical_dt, Physical_t;
	unsigned short iMesh, iDomain, nDomain;
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

	/*--- Pointer to different structures that will be used throughout the entire code ---*/
	CConfig *config = NULL;
	CIntegration **integration_container = NULL;
	COutput *output = NULL;
	CGeometry ***geometry = NULL;
	CSolution ****solution_container = NULL;
	CNumerics ****solver_container = NULL;
  CSurfaceMovement *surface_mov = NULL;
  CVolumetricMovement *grid_movement = NULL;
	CFreeFormChunk** chunk = NULL;

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
		output->SetHistory_Header(&ConvHist_file, config);
	}
	/*--- Change the name of the input-output files for the 
	 parallel computation ---*/
	config->SetFileNameDomain(rank+1);
#else
	/*--- Open the convergence history file ---*/
	output->SetHistory_Header(&ConvHist_file, config);
#endif

	/*--- Check if there are multiple domains in the mesh file ---*/
	nDomain = GetnDomain(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());

	/*--- Definition of the geometry class and open the mesh file ---*/
	geometry = new CGeometry **[nDomain];
	for (iDomain = 0; iDomain < nDomain; iDomain++) {
		geometry[iDomain] = new CGeometry *[config->GetMGLevels()+1];
		geometry[iDomain][MESH_0] = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat(), iDomain+1, nDomain);
	}

	/*--- Perform the non-dimensionalization for the flow equations ---*/
	for (iDomain = 0; iDomain < nDomain; iDomain++)
		config->SetNondimensionalization(geometry[iDomain][MESH_0]->GetnDim(), rank, iDomain);

#ifndef NO_MPI
	/*--- Synchronization point after reading the grid ---*/
	MPI::COMM_WORLD.Barrier();
#endif

	/*--- Definition of the geometry (edge based structure) for all domains ---*/
	Geometrical_Definition(geometry, config, nDomain);

	/*--- Definition of the solution class (solution_container[#DOMAINS][#GRIDS][#EQ_SYSTEMS]) ---*/
	solution_container = new CSolution***[nDomain];

	for (iDomain = 0; iDomain < nDomain; iDomain++) {
		solution_container[iDomain] = new CSolution** [config->GetMGLevels()+1];
		for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
			solution_container[iDomain][iMesh] = new CSolution* [MAX_SOLS];
		Solution_Definition(solution_container[iDomain], geometry[iDomain], config);
	}

	/*--- Definition of the integration class (integration_container[#EQ_SYSTEMS]) ---*/
	/*--- Note that this assumes the same integration on all domains ---*/
	integration_container = new CIntegration*[MAX_SOLS];
	Integration_Definition(integration_container, geometry[DOMAIN_0], config);

	/*--- Definition of the numerical method class (solver_container[#GRIDS][#EQ_SYSTEMS][#EQ_TERMS]) ---*/
	/*--- Note that this assumes the same numerics on all domains ---*/
	solver_container = new CNumerics***[config->GetMGLevels()+1];
	Solver_Definition(solver_container, solution_container[DOMAIN_0], geometry[DOMAIN_0], config);

	/*--- Computation of wall distance (note that there is MPI in this subroutine) ---*/
	for (iDomain = 0; iDomain < nDomain; iDomain++)
		if ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS))
			geometry[iDomain][MESH_0]->SetWall_Distance(config);

	/*--- Computation of positive area in the z-plane (note that there is MPI in this subroutine) ---*/
	for (iDomain = 0; iDomain < nDomain; iDomain++)
		geometry[iDomain][MESH_0]->SetPositive_ZArea(config);

	/*--- Set the near-field boundary conditions (note that there is MPI in this subroutine) ---*/
	for (iDomain = 0; iDomain < nDomain; iDomain++)
		for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
			geometry[iDomain][iMesh]->MachNearField(config);
			geometry[iDomain][iMesh]->MachInterface(config); 
		}
	
  /*--- If there is unsteady mesh movement, instatiate the geometry classes ---*/
  if (config->GetGrid_Movement()) { 
		
		if (rank == MASTER_NODE) 
			cout << endl <<"---------------------- Grid movement preprocessing ----------------------" << endl;
		
    /*--- Definition and initialization of the surface deformation class ---*/
		if (rank == MASTER_NODE) cout << "Set surface movement." << endl;
    surface_mov = new CSurfaceMovement();
    surface_mov->CopyBoundary(geometry[DOMAIN_0][MESH_0], config);
		
		/*--- Definition of the FFD deformation class (3D case) ---*/
		unsigned short nChunk = MAX_NUMBER_CHUNCK;
		chunk = new CFreeFormChunk*[nChunk];
		
		/*--- Read the FFD information fron the grid file (3D problems)---*/
		if (geometry[DOMAIN_0][MESH_0]->GetnDim() == 3)
			surface_mov->ReadFFDInfo(config, geometry[DOMAIN_0][MESH_0], chunk, config->GetMesh_FileName());
		
    /*--- Definition of the Class for grid movement ---*/
		if (rank == MASTER_NODE) cout << "Set volumetric grid movement." << endl;
    grid_movement = new CVolumetricMovement(geometry[DOMAIN_0][MESH_0]);		
		
  }
  
	/*--- External loop of the solver ---*/
	if (rank == MASTER_NODE) 
		cout << endl <<"------------------------------ Start solver -----------------------------" << endl;

	while (ExtIter < config->GetnExtIter()) {

		StartTime = clock();
		config->UpdateCFL(ExtIter);

		switch (config->GetKind_Solver()) {
				
			case EULER: case NAVIER_STOKES: case RANS:
        
        for (iDomain = 0; iDomain < nDomain; iDomain++) {
					
					/*--- Set the value of the internal iteration ---*/
					IntIter = ExtIter;
					if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || 
							(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
					
          /*--- Set the initial condition at the first iteration ---*/
          integration_container[FLOW_SOL]->SetInitialCondition(geometry[iDomain], solution_container[iDomain], solver_container, 
                                                               config, ExtIter);
					
          /*--- Solve the Euler's, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/
					if (config->GetKind_Solver() == EULER) config->SetGlobalParam(EULER, RUNTIME_FLOW_SYS);
					if (config->GetKind_Solver() == NAVIER_STOKES) config->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS);
					if (config->GetKind_Solver() == RANS) config->SetGlobalParam(RANS, RUNTIME_FLOW_SYS);
					
          integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
                                                               config, RUNTIME_FLOW_SYS, IntIter);
					
          /*--- Dual time stepping strategy ---*/
          if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
            for(IntIter = 1; IntIter < config->GetUnst_nIntIter(); IntIter++) {
              integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
                                                                   config, RUNTIME_FLOW_SYS, IntIter);
							output->SetHistory_DualTime(geometry[iDomain], solution_container[iDomain], config, integration_container, IntIter);
              if (integration_container[FLOW_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout<<endl; break;}
            }
						
						/*--- Update dual time solver ---*/
            for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
              integration_container[FLOW_SOL]->SetDualTime_Solver(geometry[iDomain][iMesh], solution_container[iDomain][iMesh][FLOW_SOL]);
              integration_container[FLOW_SOL]->SetConvergence(false);
            }
            /*--- Perform mesh motion, if necessary ---*/
            if (config->GetGrid_Movement()) {
							
							if (rank == MASTER_NODE) 
								cout << endl << "Performing the deformation of the volumetric grid." << endl;
							
							/*--- Surface grid deformation ---*/
							if (geometry[iDomain][MESH_0]->GetnDim() == 2)
								surface_mov->SetBoundary_Flutter2D(geometry[iDomain][MESH_0], config, ExtIter);
							else
								surface_mov->SetBoundary_Flutter3D(geometry[iDomain][MESH_0], config, chunk, ExtIter);
							
							/*--- Volumetric grid deformation ---*/
              grid_movement->SpringMethod(geometry[DOMAIN_0][MESH_0], config);
            }
            Physical_dt = config->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
            if (Physical_t >=  config->GetUnst_Time()) integration_container[FLOW_SOL]->SetConvergence(true);
          }
					
					/*--- Solve the turbulence model ---*/
					if (config->GetKind_Solver() == RANS) {
						/*--- Turbulent model solution ---*/
						config->SetGlobalParam(RANS, RUNTIME_TURB_SYS);
						integration_container[TURB_SOL]->SetSingleGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
																																	config, RUNTIME_TURB_SYS, ExtIter);
					}
				}
				
				break;
				
		case ADJ_EULER: case ADJ_NAVIER_STOKES:

			for (iDomain = 0; iDomain < nDomain; iDomain++) {
				/*--- Continuous adjoint Euler equations ---*/

				if (ExtIter == 0) {

					if (rank == MASTER_NODE) cout << "Iteration over the direct problem to store all flow information." << endl;

					if (config->GetKind_Solver() == ADJ_EULER) config->SetGlobalParam(ADJ_EULER, RUNTIME_FLOW_SYS);
					if (config->GetKind_Solver() == ADJ_NAVIER_STOKES) config->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_FLOW_SYS);

					integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
							config, RUNTIME_FLOW_SYS, ExtIter);

					/*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
						 note that in the direct problem we are not computing the gradients ---*/	
					if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
						solution_container[iDomain][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry[iDomain][MESH_0], config);
					if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
							(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES))
						solution_container[iDomain][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry[iDomain][MESH_0], config);

					/*--- Set force projection vector for boundary conditions ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						solution_container[iDomain][iMesh][FLOW_SOL]->SetTotal_CDrag(solution_container[iDomain][MESH_0][FLOW_SOL]->GetTotal_CDrag());
						solution_container[iDomain][iMesh][FLOW_SOL]->SetTotal_CLift(solution_container[iDomain][MESH_0][FLOW_SOL]->GetTotal_CLift());
						solution_container[iDomain][iMesh][FLOW_SOL]->SetTotal_CT(solution_container[iDomain][MESH_0][FLOW_SOL]->GetTotal_CT());
						solution_container[iDomain][iMesh][FLOW_SOL]->SetTotal_CQ(solution_container[iDomain][MESH_0][FLOW_SOL]->GetTotal_CQ());
						solution_container[iDomain][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry[iDomain][iMesh], solution_container[iDomain][iMesh], config);
						if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) || (config->GetKind_ObjFunc() == NEARFIELD_PRESSURE))
							solution_container[iDomain][iMesh][ADJFLOW_SOL]->SetIntBoundary_Jump(geometry[iDomain][iMesh], solution_container[iDomain][iMesh], config);
					}
				}

				/*--- Set the value of the internal iteration ---*/
				IntIter = ExtIter;
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
						(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
					/*--- Read the flow solution from a restart file, and interpolate in the coarse grid ---*/


					/*--- Set the iner iteration to zero ---*/
					IntIter = 0;
				}

				if (config->GetKind_Solver() == ADJ_EULER) config->SetGlobalParam(ADJ_EULER, RUNTIME_ADJFLOW_SYS);
				if (config->GetKind_Solver() == ADJ_NAVIER_STOKES) config->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS);

				integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config,
						RUNTIME_ADJFLOW_SYS, IntIter);

				/*--- Dual time stepping strategy ---*/
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
						(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

					for(IntIter = 1; IntIter < config->GetUnst_nIntIter(); IntIter++) {
						integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config,
								RUNTIME_ADJFLOW_SYS, IntIter);

						output->SetHistory_DualTime(geometry[iDomain], solution_container[iDomain], config, integration_container, IntIter);

						if (integration_container[ADJFLOW_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout<<endl; break;}
					}

					/*--- Update dual time solver ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[ADJFLOW_SOL]->SetDualTime_Solver(geometry[iDomain][iMesh], solution_container[iDomain][iMesh][ADJFLOW_SOL]);
						integration_container[ADJFLOW_SOL]->SetConvergence(false);
					}

					Physical_dt = config->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time()) integration_container[ADJFLOW_SOL]->SetConvergence(true);
				}

			}
			break;

		case NS_PLASMA:

			for (iDomain = 0; iDomain < nDomain; iDomain++) {

				/*--- Set the value of the internal iteration ---*/
				IntIter = ExtIter;
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
						(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

				/*--- Plasma equations ---*/
				config->SetGlobalParam(NS_PLASMA, RUNTIME_PLASMA_SYS);
				integration_container[PLASMA_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config, RUNTIME_PLASMA_SYS, IntIter);

				/*--- Dual time stepping strategy ---*/
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
					for(IntIter = 1; IntIter < config->GetUnst_nIntIter(); IntIter++) {
						integration_container[PLASMA_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
								config, RUNTIME_PLASMA_SYS, IntIter);
						output->SetHistory_DualTime(geometry[iDomain], solution_container[iDomain], config, integration_container, IntIter);
						if (integration_container[PLASMA_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout<<endl; break;}
					}

					/*--- Update dual time solver ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[PLASMA_SOL]->SetDualTime_Solver(geometry[iDomain][iMesh], solution_container[iDomain][iMesh][PLASMA_SOL]);
						integration_container[PLASMA_SOL]->SetConvergence(false);
					}

					Physical_dt = config->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time()) integration_container[PLASMA_SOL]->SetConvergence(true);
				}

				/*--- Electric potential equation ---*/
				if (config->GetKind_GasModel() == ARGON || config->GetKind_GasModel() == AIR21) {

					if (config->GetElectricSolver()) {
						config->SetGlobalParam(NS_PLASMA, RUNTIME_ELEC_SYS);
						integration_container[ELEC_SOL]->SetPotential_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config, RUNTIME_ELEC_SYS, MESH_0);
					}
				}
				if (config->GetKind_GasModel() == AIR7) {
					config->SetGlobalParam(NS_PLASMA, RUNTIME_ELEC_SYS);
				}
			}
			break;
		case ADJ_NS_PLASMA:

			for (iDomain = 0; iDomain < nDomain; iDomain++) {

				/*--- Plasma equations ---*/
				config->SetGlobalParam(ADJ_NS_PLASMA, RUNTIME_PLASMA_SYS);
				integration_container[PLASMA_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config, RUNTIME_PLASMA_SYS, ExtIter);

				/*--- Adjoint Plasma equations ---*/
				config->SetGlobalParam(ADJ_NS_PLASMA, RUNTIME_ADJPLASMA_SYS);
				integration_container[ADJPLASMA_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config, RUNTIME_ADJPLASMA_SYS, ExtIter);

			}
			break;

		case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES:

			for (iDomain = 0; iDomain < nDomain; iDomain++) {

				/*--- Set the value of the internal iteration ---*/
				IntIter = ExtIter;
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
						(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

				/*--- Update density and viscosity, using level set, including dual time strategy (first iteration) ---*/
				solution_container[iDomain][MESH_0][LEVELSET_SOL]->SetLevelSet_Distance(geometry[iDomain][MESH_0], config);
				output->SetFreeSurface(solution_container[iDomain][MESH_0][LEVELSET_SOL], geometry[iDomain][MESH_0], config, ExtIter);

				integration_container[FLOW_SOL]->SetFreeSurface_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
						config, ExtIter);

				/*--- Navier-Stokes equations ---*/
				if (config->GetKind_Solver() == FREE_SURF_EULER) config->SetGlobalParam(FREE_SURF_EULER, RUNTIME_FLOW_SYS);
				if (config->GetKind_Solver() == FREE_SURF_NAVIER_STOKES) config->SetGlobalParam(FREE_SURF_NAVIER_STOKES, RUNTIME_FLOW_SYS);
				integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
						config, RUNTIME_FLOW_SYS, IntIter);

				/*--- Level-Set model solution ---*/
				if (config->GetKind_Solver() == FREE_SURF_EULER) config->SetGlobalParam(FREE_SURF_EULER, RUNTIME_LEVELSET_SYS);
				if (config->GetKind_Solver() == FREE_SURF_NAVIER_STOKES) config->SetGlobalParam(FREE_SURF_NAVIER_STOKES, RUNTIME_LEVELSET_SYS);
				integration_container[LEVELSET_SOL]->SetSingleGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
						config, RUNTIME_LEVELSET_SYS, IntIter);

				/*--- Dual time stepping strategy for the flow equations ---*/
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
					for(IntIter = 1; IntIter < config->GetUnst_nIntIter(); IntIter++) {

						/*--- Update density and viscosity, using level set, including dual time strategy (first iteration) ---*/
						integration_container[FLOW_SOL]->SetFreeSurface_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
								config, IntIter);

						/*--- Navier-Stokes equations ---*/
						if (config->GetKind_Solver() == FREE_SURF_EULER) config->SetGlobalParam(FREE_SURF_EULER, RUNTIME_FLOW_SYS);
						if (config->GetKind_Solver() == FREE_SURF_NAVIER_STOKES) config->SetGlobalParam(FREE_SURF_NAVIER_STOKES, RUNTIME_FLOW_SYS);
						integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
								config, RUNTIME_FLOW_SYS, IntIter);

						/*--- Level-Set model solution ---*/
						if (config->GetKind_Solver() == FREE_SURF_EULER) config->SetGlobalParam(FREE_SURF_EULER, RUNTIME_LEVELSET_SYS);
						if (config->GetKind_Solver() == FREE_SURF_NAVIER_STOKES) config->SetGlobalParam(FREE_SURF_NAVIER_STOKES, RUNTIME_LEVELSET_SYS);
						integration_container[LEVELSET_SOL]->SetSingleGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
								config, RUNTIME_LEVELSET_SYS, IntIter);

						output->SetHistory_DualTime(geometry[iDomain], solution_container[iDomain], config, integration_container, IntIter);

						if (integration_container[FLOW_SOL]->GetConvergence()) break;
					}

					/*--- Set convergence the global convergence criteria to false, and dual time solution ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[FLOW_SOL]->SetDualTime_Solver(geometry[iDomain][iMesh], solution_container[iDomain][iMesh][FLOW_SOL]);
						integration_container[FLOW_SOL]->SetConvergence(false);
					}

					integration_container[LEVELSET_SOL]->SetDualTime_Solver(geometry[iDomain][MESH_0], solution_container[iDomain][MESH_0][LEVELSET_SOL]);
					integration_container[LEVELSET_SOL]->SetConvergence(false);

					/*--- Set the value of the global convergence criteria ---*/
					Physical_dt = config->GetDelta_UnstTime();
					Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time())
						integration_container[FLOW_SOL]->SetConvergence(true);
				}

			}

			break;


		case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES:

			for (iDomain = 0; iDomain < nDomain; iDomain++) {
				/*--- Continuous adjoint free surface equations ---*/

				if (ExtIter == 0) {

					if (rank == MASTER_NODE) cout << "Iteration over the direct problem to store all flow information." << endl;

					/*--- Update density and viscosity, using level set, including dual time strategy (first iteration) ---*/
					solution_container[iDomain][MESH_0][LEVELSET_SOL]->SetLevelSet_Distance(geometry[iDomain][MESH_0], config);
					integration_container[FLOW_SOL]->SetFreeSurface_Solver(geometry[iDomain], solution_container[iDomain], solver_container, 
							config, ExtIter);

					if (config->GetKind_Solver() == ADJ_FREE_SURF_EULER) config->SetGlobalParam(ADJ_FREE_SURF_EULER, RUNTIME_FLOW_SYS);
					if (config->GetKind_Solver() == ADJ_FREE_SURF_NAVIER_STOKES) config->SetGlobalParam(ADJ_FREE_SURF_NAVIER_STOKES, RUNTIME_FLOW_SYS);

					integration_container[FLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
							config, RUNTIME_FLOW_SYS, ExtIter);

					config->SetGlobalParam(ADJ_FREE_SURF_EULER, RUNTIME_LEVELSET_SYS);
					integration_container[LEVELSET_SOL]->SetSingleGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
							config, RUNTIME_LEVELSET_SYS, ExtIter);

					/*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
						 note that in the direct problem we are not computing the gradients ---*/	
					if (config->GetKind_Gradient_Method() == GREEN_GAUSS)
						solution_container[iDomain][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry[iDomain][MESH_0], config);
					if ((config->GetKind_Gradient_Method() == LEAST_SQUARES) ||
							(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES))
						solution_container[iDomain][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry[iDomain][MESH_0], config);

					/*--- Set force projection vector for boundary conditions ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						solution_container[iDomain][iMesh][FLOW_SOL]->SetTotal_CDrag(solution_container[iDomain][MESH_0][FLOW_SOL]->GetTotal_CDrag());
						solution_container[iDomain][iMesh][FLOW_SOL]->SetTotal_CLift(solution_container[iDomain][MESH_0][FLOW_SOL]->GetTotal_CLift());
						solution_container[iDomain][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry[iDomain][iMesh], solution_container[iDomain][iMesh], config);
					}

					solution_container[DOMAIN_0][MESH_0][LEVELSET_SOL]->SetLevelSet_Distance(geometry[DOMAIN_0][MESH_0], config);
					output->SetFreeSurface(solution_container[DOMAIN_0][MESH_0][LEVELSET_SOL], geometry[DOMAIN_0][MESH_0], config, ExtIter);

				}

				/*--- Set the value of the internal iteration ---*/
				IntIter = ExtIter;
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
						(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

				/*--- Iteration over the flow-adjoint equations ---*/
				config->SetGlobalParam(ADJ_FREE_SURF_EULER, RUNTIME_ADJFLOW_SYS);
				integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config,
						RUNTIME_ADJFLOW_SYS, IntIter);

				/*--- Iteration over the level-set-adjoint equations ---*/
				config->SetGlobalParam(ADJ_FREE_SURF_EULER, RUNTIME_ADJLEVELSET_SYS);
				integration_container[ADJLEVELSET_SOL]->SetSingleGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config,
						RUNTIME_ADJLEVELSET_SYS, IntIter);

				/*--- Dual time stepping strategy for the flow equations ---*/
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
					for(IntIter = 1; IntIter < config->GetUnst_nIntIter(); IntIter++) {

						/*--- Navier-Stokes equations ---*/
						config->SetGlobalParam(ADJ_FREE_SURF_EULER, RUNTIME_ADJFLOW_SYS);
						integration_container[ADJFLOW_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
								config, RUNTIME_ADJFLOW_SYS, IntIter);

						/*--- Level-Set model solution ---*/
						config->SetGlobalParam(ADJ_FREE_SURF_EULER, RUNTIME_ADJLEVELSET_SYS);
						integration_container[ADJLEVELSET_SOL]->SetSingleGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
								config, RUNTIME_ADJLEVELSET_SYS, IntIter);

						output->SetHistory_DualTime(geometry[iDomain], solution_container[iDomain], config, integration_container, IntIter);

						if (integration_container[ADJFLOW_SOL]->GetConvergence()) break;
					}

					/*--- Set convergence the global convergence criteria to false, and dual time solution ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[ADJFLOW_SOL]->SetDualTime_Solver(geometry[iDomain][iMesh], solution_container[iDomain][iMesh][FLOW_SOL]);
						integration_container[ADJFLOW_SOL]->SetConvergence(false);
					}

					integration_container[ADJLEVELSET_SOL]->SetDualTime_Solver(geometry[iDomain][MESH_0], solution_container[iDomain][MESH_0][LEVELSET_SOL]);
					integration_container[ADJLEVELSET_SOL]->SetConvergence(false);

					/*--- Set the value of the global convergence criteria ---*/
					Physical_dt = config->GetDelta_UnstTime();
					Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time())
						integration_container[ADJFLOW_SOL]->SetConvergence(true);
				}
			}
			break;


		case WAVE:

			for (iDomain = 0; iDomain < nDomain; iDomain++) {

				/*--- Set the value of the internal iteration ---*/
				IntIter = ExtIter;
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
						(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

				/*--- Wave equations ---*/
				config->SetGlobalParam(WAVE, RUNTIME_WAVE_SYS);
				integration_container[WAVE_SOL]->SetSingleGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container, config, RUNTIME_WAVE_SYS, IntIter);

				/*--- Dual time stepping strategy ---*/
				if ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
					for(IntIter = 1; IntIter < config->GetUnst_nIntIter(); IntIter++) {
						integration_container[WAVE_SOL]->SetMultiGrid_Solver(geometry[iDomain], solution_container[iDomain], solver_container,
								config, RUNTIME_WAVE_SYS, IntIter);
						output->SetHistory_DualTime(geometry[iDomain], solution_container[iDomain], config, integration_container, IntIter);
						if (integration_container[WAVE_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout<<endl; break;}
					}

					/*--- Update dual time solver ---*/
					for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
						integration_container[WAVE_SOL]->SetDualTime_Solver(geometry[iDomain][iMesh], solution_container[iDomain][iMesh][WAVE_SOL]);
						integration_container[WAVE_SOL]->SetConvergence(false);
					}

					Physical_dt = config->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
					if (Physical_t >=  config->GetUnst_Time()) integration_container[WAVE_SOL]->SetConvergence(true);
				}

			}
			break;

		}

#ifndef NO_MPI
		MPI::COMM_WORLD.Barrier();
#endif
		StopTime = clock(); TimeUsed += (StopTime - StartTime);

		/*--- Convergence history for serial and parallel computation ---*/	
		if (ExtIter % config->GetWrt_Con_Freq() == 0) {

			/*--- Evaluate and plot the equivalent area ---*/
			if ((config->GetKind_Solver() == EULER) && (config->GetEquivArea() == YES))
				output->SetEquivalentArea(solution_container[DOMAIN_0][MESH_0][FLOW_SOL], geometry[DOMAIN_0][MESH_0], config, ExtIter);

			/*--- Print the history file --*/
			output->SetHistory_MainIter(&ConvHist_file, geometry, solution_container, config, integration_container, ExtIter, TimeUsed, nDomain);

		}

		/*--- Review convergence criteria ---*/
		switch (config->GetKind_Solver()) {
		case EULER: case NAVIER_STOKES: case RANS:
			StopCalc = integration_container[FLOW_SOL]->GetConvergence(); break;
		case FREE_SURF_EULER: case FREE_SURF_NAVIER_STOKES: case FREE_SURF_RANS:
			StopCalc = integration_container[FLOW_SOL]->GetConvergence(); break;
		case NS_PLASMA:
			StopCalc = integration_container[PLASMA_SOL]->GetConvergence(); break;
		case ELECTRIC_POTENTIAL:
			StopCalc = integration_container[ELEC_SOL]->GetConvergence(); break;
		case WAVE:
			StopCalc = integration_container[WAVE_SOL]->GetConvergence(); break;
		case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
			StopCalc = integration_container[ADJFLOW_SOL]->GetConvergence(); break;
		case ADJ_FREE_SURF_EULER: case ADJ_FREE_SURF_NAVIER_STOKES: case ADJ_FREE_SURF_RANS:
			StopCalc = integration_container[ADJFLOW_SOL]->GetConvergence(); break;
		case LIN_EULER: case LIN_NAVIER_STOKES:
			StopCalc = integration_container[LINFLOW_SOL]->GetConvergence(); break;
		}

		/*--- Solution output (be careful there is some MPI stuff in these subroutines) ---*/
		if (
				(((ExtIter+1 == config->GetnExtIter()) || (ExtIter % config->GetWrt_Sol_Freq() == 0) || StopCalc ) && (ExtIter != 0)) 
				|| (config->GetnExtIter() == 1) 
				|| (((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) && (ExtIter == 0))
		) {
			output->SetResult_Files(solution_container, geometry, config, ExtIter, nDomain);
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
