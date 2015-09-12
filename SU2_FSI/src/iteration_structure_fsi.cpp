/*!
 * \file iteration_structure_fsi.cpp
 * \brief Main subroutines used by SU2_FSI for iteration
 * \author R. Sanchez
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/iteration_structure_fsi.hpp"

void FSI_BGS_Iteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                             CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                             CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
							 CTransfer*** transfer_container, unsigned long iFluidIt, unsigned long nFluidIt) {

	su2double Physical_dt, Physical_t;
	unsigned short iMesh;
	unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
	unsigned long IntIter_Struct = 0; config_container[ZONE_1]->SetIntIter(IntIter_Struct);
	unsigned long iFSIIter = 0;
	unsigned long nFSIIter = config_container[ZONE_0]->GetnIterFSI();
	unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

	bool fem_solver = (config_container[ZONE_1]->GetKind_Solver() == FEM_ELASTICITY);

	int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/*-----------------------------------------------------------------*/
	/*---------------- Predict structural displacements ---------------*/
	/*-----------------------------------------------------------------*/

	FSI_Disp_Predictor(output, integration_container, geometry_container,
                          solver_container, numerics_container, config_container,
                          surface_movement, grid_movement, FFDBox);

	while (iFSIIter<nFSIIter){

		/*-----------------------------------------------------------------*/
		/*------------------------ Update mesh ----------------------------*/
		/*-----------------------------------------------------------------*/

        FSI_Disp_Transfer(output, integration_container, geometry_container,
                          solver_container, numerics_container, config_container,
                          surface_movement, grid_movement, FFDBox, transfer_container);

		/*-----------------------------------------------------------------*/
		/*-------------------- Fluid subiteration -------------------------*/
		/*-----------------------------------------------------------------*/

        Flow_Subiteration(output, integration_container, geometry_container,
                          solver_container, numerics_container, config_container,
                          surface_movement, grid_movement, FFDBox);

		/*--- Write the convergence history for the fluid (only screen output) ---*/

		output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

		/*-----------------------------------------------------------------*/
		/*------------------- Set FEA loads from fluid --------------------*/
		/*-----------------------------------------------------------------*/

        FSI_Load_Transfer(output, integration_container, geometry_container,
	                 	 solver_container, numerics_container, config_container,
	                 	 surface_movement, grid_movement, FFDBox, transfer_container, ExtIter);

		/*-----------------------------------------------------------------*/
		/*------------------ Structural subiteration ----------------------*/
		/*-----------------------------------------------------------------*/

        if (fem_solver){
    	    FEM_Subiteration(output, integration_container, geometry_container,
    	                 	 solver_container, numerics_container, config_container,
    	                 	 surface_movement, grid_movement, FFDBox);
        }
        else{
    	    FEA_Subiteration(output, integration_container, geometry_container,
    	                 	 solver_container, numerics_container, config_container,
    	                 	 surface_movement, grid_movement, FFDBox);
        }

		/*--- Write the convergence history for the structure (only screen output) ---*/

		output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_1);

		/*-----------------------------------------------------------------*/
		/*----------------- Displacements relaxation ----------------------*/
		/*-----------------------------------------------------------------*/

	    FSI_Disp_Relaxation(output, geometry_container, solver_container, config_container, iFSIIter);

		/*-----------------------------------------------------------------*/
		/*-------------------- Check convergence --------------------------*/
		/*-----------------------------------------------------------------*/

		integration_container[ZONE_1][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[ZONE_1][MESH_0], config_container[ZONE_1],
																solver_container[ZONE_1][MESH_0][FEA_SOL], iFSIIter);

		if (integration_container[ZONE_1][FEA_SOL]->GetConvergence_FSI()) break;

		/*-----------------------------------------------------------------*/
		/*--------------------- Update iFSIIter ---------------------------*/
		/*-----------------------------------------------------------------*/

		iFSIIter++;

	}

	/*-----------------------------------------------------------------*/
  	/*-------------------- Update fluid solver ------------------------*/
	/*-----------------------------------------------------------------*/

    Flow_Update(output, integration_container, geometry_container,
                solver_container, numerics_container, config_container,
                surface_movement, grid_movement, FFDBox, ExtIter);

	/*-----------------------------------------------------------------*/
  	/*----------------- Update structural solver ----------------------*/
	/*-----------------------------------------------------------------*/

    if (fem_solver){
    	FEM_Update(output, integration_container, geometry_container,
                	solver_container, numerics_container, config_container, ExtIter);
    }
    else{
    	FEA_Update(output, integration_container, geometry_container,
                	solver_container, config_container, ExtIter);
    }


	/*-----------------------------------------------------------------*/
	/*--------------- Update convergence parameter --------------------*/
	/*-----------------------------------------------------------------*/
	integration_container[ZONE_1][FEA_SOL]->SetConvergence_FSI(false);

}

void Flow_Subiteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                       CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                       CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {

	unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
	unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

	unsigned short iMesh, iZone;

	/*--- Only one zone allowed for the fluid as for now ---*/
	unsigned short nFluidZone = 1;

	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	/*--- Set the initial condition ---*/
	for (iZone = 0; iZone < nFluidZone; iZone++)
		solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

	/*--- Apply a Wind Gust ---*/
	/*--- Initial set up for unsteady problems with gusts - Not enabled yet. ---*/

//	for (iZone = 0; iZone < nFluidZone; iZone++) {
//		if (config_container[ZONE_0]->GetWind_Gust()){
//			  SetWind_GustField(config_container[iZone],geometry_container[iZone],solver_container[iZone]);
//		}
//	}

	for (iZone = 0; iZone < nFluidZone; iZone++) {

		/*--- Set the value of the internal iteration ---*/

		IntIter = ExtIter;

		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
		   (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

		/*--- Update global parameters ---*/

		if (config_container[iZone]->GetKind_Solver() == EULER){
			config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
		}
		if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES){
			config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
		}
		if (config_container[iZone]->GetKind_Solver() == RANS){
			config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
		 }

		/*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/

		integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
																	config_container, RUNTIME_FLOW_SYS, IntIter, iZone);

		if (config_container[iZone]->GetKind_Solver() == RANS) {

			/*--- Solve the turbulence model ---*/

			config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
			integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																		   config_container, RUNTIME_TURB_SYS, IntIter, iZone);

		/*--- Solve transition model ---*/

			if (config_container[iZone]->GetKind_Trans_Model() == LM) {
				config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
				integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																			  config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
			}

		}

	}

	/*--- Dual time stepping strategy ---*/

	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
		   (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

			for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {

			/*--- Write the convergence history (only screen output) ---*/

			output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

			/*--- Set the value of the internal iteration ---*/

			config_container[ZONE_0]->SetIntIter(IntIter);

			for (iZone = 0; iZone < nFluidZone; iZone++) {

				/*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/

				if (config_container[iZone]->GetKind_Solver() == EULER)
					config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
				if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES)
					config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
				if (config_container[iZone]->GetKind_Solver() == RANS)
					config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);

				  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/

				integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
																				config_container, RUNTIME_FLOW_SYS, IntIter, iZone);

				/*--- Pseudo-timestepping the turbulence model ---*/

				if (config_container[iZone]->GetKind_Solver() == RANS) {

					/*--- Solve the turbulence model ---*/

					config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
					integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																				  config_container, RUNTIME_TURB_SYS, IntIter, iZone);

					/*--- Solve transition model ---*/

						if (config_container[iZone]->GetKind_Trans_Model() == LM) {
							config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
							integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																					  config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
						}
				}
			}
			if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) break;
		}
	}
}


void Flow_Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                       CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                       CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
                       unsigned long ExtIter) {

	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif


	su2double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;

	/*--- Only one zone allowed for the fluid as for now ---*/
	unsigned short nFluidZone = 1;

	for (iZone = 0; iZone < nFluidZone; iZone++) {

		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
		   (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

			/*--- Update dual time solver on all mesh levels ---*/

			for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
				integration_container[iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][FLOW_SOL], config_container[iZone], iMesh);
				integration_container[iZone][FLOW_SOL]->SetConvergence(false);
			}

			/*--- Update dual time solver for the turbulence model ---*/

			if (config_container[iZone]->GetKind_Solver() == RANS) {
				integration_container[iZone][TURB_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][TURB_SOL], config_container[iZone], MESH_0);
				integration_container[iZone][TURB_SOL]->SetConvergence(false);
			}

			/*--- Update dual time solver for the transition model ---*/

			if (config_container[iZone]->GetKind_Trans_Model() == LM) {
				integration_container[iZone][TRANS_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][TRANS_SOL], config_container[iZone], MESH_0);
				integration_container[iZone][TRANS_SOL]->SetConvergence(false);
			}

			/*--- Verify convergence criteria (based on total time) ---*/

			Physical_dt = config_container[iZone]->GetDelta_UnstTime();
			Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime())
				integration_container[iZone][FLOW_SOL]->SetConvergence(true);

		}

	}


}


void FEA_Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                  CSolver ****solver_container, CConfig **config_container,
                  unsigned long ExtIter) {

	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif


	/*--- Only one zone allowed for the structure as for now ---*/
	unsigned short nFluidZone = 1, nStrucZone=1, nTotalZone;
	unsigned short iZone;

	nTotalZone = nFluidZone + nStrucZone;

	for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

		/*----------------- Update structural solver ----------------------*/

		integration_container[iZone][FEA_SOL]->SetStructural_Solver(geometry_container[iZone][MESH_0],
												solver_container[iZone][MESH_0][FEA_SOL], config_container[iZone], MESH_0);

	}


}

void FEM_Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                  CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                  unsigned long ExtIter) {

	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	/*--- Only one zone allowed for the structure as for now ---*/
	su2double Physical_dt, Physical_t;
	unsigned short nFluidZone = 1, nStrucZone=1, nTotalZone;
	unsigned short iZone;
	unsigned int ZONE_STRUC = nFluidZone;

	nTotalZone = nFluidZone + nStrucZone;

	bool dynamic = (config_container[ZONE_STRUC]->GetDynamic_Analysis() == DYNAMIC);					// Dynamic problems
	bool nonlinear = (config_container[ZONE_STRUC]->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Geometrically non-linear problems


	/*----------------- Compute averaged nodal stress ------------------------*/

	for (iZone = nFluidZone; iZone < nTotalZone; iZone++)
		solver_container[iZone][MESH_0][FEA_SOL]->Compute_NodalStress(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], numerics_container[iZone][MESH_0][FEA_SOL][VISC_TERM], config_container[iZone]);

	/*----------------- Update structural solver ----------------------*/

	if (dynamic){
		for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {
			integration_container[iZone][FEA_SOL]->SetFEM_StructuralSolver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][FEA_SOL], config_container[iZone], MESH_0);
			integration_container[iZone][FEA_SOL]->SetConvergence(false);
		}

	    /*--- Verify convergence criteria (based on total time) ---*/
		Physical_dt = config_container[ZONE_STRUC]->GetDelta_DynTime();
		Physical_t  = (ExtIter+1)*Physical_dt;
		if (Physical_t >=  config_container[ZONE_STRUC]->GetTotal_DynTime())
			integration_container[ZONE_STRUC][FEA_SOL]->SetConvergence(true);
	}


}


void FEA_Subiteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                  CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                  CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	su2double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	unsigned long IntIter_Struct = 0; //config_container[ZONE_0]->SetIntIter(IntIter);
  	unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

#ifdef HAVE_MPI
	int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/*--- Only one zone allowed for the fluid as for now ---*/
	unsigned short nFluidZone = 1, nStrucZone=1, nTotalZone;

	nTotalZone = nFluidZone + nStrucZone;

	for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

		/*--- Set the initial condition at the first iteration ---*/

	//	solver_container[ZONE_1][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[ZONE_1], solver_container[ZONE_1], config_container[ZONE_1], ExtIter);

		/*--- Set the value of the internal iteration ---*/

		IntIter_Struct = ExtIter;

		/*--- FEA equations ---*/

		config_container[iZone]->SetGlobalParam(LINEAR_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

		/*--- Run the iteration ---*/

		integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
																	 config_container, RUNTIME_FEA_SYS, IntIter_Struct, iZone);

	}

}

void FEM_Subiteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                  CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                  CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {


	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	/*--- Only one zone allowed for the fluid as for now ---*/
	unsigned short nFluidZone = 1, nStrucZone=1, nTotalZone;

	nTotalZone = nFluidZone + nStrucZone;

	unsigned int ZONE_STRUC = nFluidZone;

	su2double Physical_dt, Physical_t;
	su2double loadIncrement;
	unsigned short iMesh, iZone;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	unsigned long IntIter = 0; config_container[ZONE_STRUC]->SetIntIter(IntIter);
  	unsigned long ExtIter = config_container[ZONE_STRUC]->GetExtIter();

  	unsigned long iIncrement;
  	unsigned long nIncrements = config_container[ZONE_STRUC]->GetNumberIncrements();

	bool dynamic = (config_container[ZONE_STRUC]->GetDynamic_Analysis() == DYNAMIC);					// Dynamic problems
	bool nonlinear = (config_container[ZONE_STRUC]->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Geometrically non-linear problems
	bool linear = (config_container[ZONE_STRUC]->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Geometrically non-linear problems

	bool initial_calc = config_container[ZONE_STRUC]->GetExtIter() == 0;				// Checks if it is the first calculation.
	bool first_iter = config_container[ZONE_STRUC]->GetIntIter() == 0;				// Checks if it is the first iteration
	bool restart = config_container[ZONE_STRUC]->GetRestart();												// Restart analysis
	bool initial_calc_restart = (config_container[ZONE_STRUC]->GetExtIter() == config_container[ZONE_STRUC]->GetDyn_RestartIter()); // Initial calculation for restart


	su2double CurrentTime = config_container[ZONE_STRUC]->GetCurrent_DynTime();
	su2double Static_Time = config_container[ZONE_STRUC]->GetStatic_Time();

	bool statTime = (CurrentTime <= Static_Time);

	bool incremental_load = config_container[ZONE_STRUC]->GetIncrementalLoad();							// If an incremental load is applied

	/*--- Set the convergence monitor to false, to prevent the solver to stop in intermediate FSI subiterations ---*/
	integration_container[ZONE_STRUC][FEA_SOL]->SetConvergence(false);

	if (linear){

		for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

			/*--- Set the value of the internal iteration ---*/

			IntIter = ExtIter;

			/*--- FEA equations ---*/

			config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

			/*--- Run the iteration ---*/

			integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
																			config_container, RUNTIME_FEA_SYS, IntIter, iZone);

		}

	}
	/*--- If the structure is held static and the solver is nonlinear, we don't need to solve for static time, but we need to compute Mass Matrix and Integration constants ---*/
	else if ((nonlinear) && (!statTime)){

		/*--- THIS IS THE DIRECT APPROACH (NO INCREMENTAL LOAD APPLIED) ---*/

		if (!incremental_load){

			for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

				/*--- Set the value of the internal iteration ---*/

				IntIter = 0;

				/*--- FEA equations ---*/

				config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

				/*--- Run the iteration ---*/

				integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
		                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);



			}

			/*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

			for (IntIter = 1; IntIter < config_container[ZONE_STRUC]->GetDyn_nIntIter(); IntIter++){

				for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

					/*--- Write the convergence history (only screen output) ---*/

					output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUC);

					config_container[iZone]->SetIntIter(IntIter);

					integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
																					config_container, RUNTIME_FEA_SYS, IntIter, iZone);

				}

				if (integration_container[ZONE_STRUC][FEA_SOL]->GetConvergence()) break;

			}

		}
		/*--- The incremental load is only used in nonlinear cases ---*/
		else if (incremental_load){

			/*--- Set the initial condition: store the current solution as Solution_Old ---*/

			for (iZone = nFluidZone; iZone < nTotalZone; iZone++)
				solver_container[iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

			for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

					/*--- The load increment is 1.0: all the load applied ---*/
					loadIncrement = 1.0;
					solver_container[iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);

					/*--- Set the value of the internal iteration ---*/

					IntIter = 0;

					/*--- FEA equations ---*/

					config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

					/*--- Run the first iteration ---*/

					integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
			                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);


					/*--- Write the convergence history (only screen output) ---*/

					output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUC);

					/*--- Run the second iteration ---*/

					IntIter = 1;

					config_container[iZone]->SetIntIter(IntIter);

					integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
			                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);

			}

			bool meetCriteria;
			su2double Residual_UTOL, Residual_RTOL, Residual_ETOL;
			su2double Criteria_UTOL, Criteria_RTOL, Criteria_ETOL;

			Criteria_UTOL = config_container[ZONE_STRUC]->GetIncLoad_Criteria(0);
			Criteria_RTOL = config_container[ZONE_STRUC]->GetIncLoad_Criteria(1);
			Criteria_ETOL = config_container[ZONE_STRUC]->GetIncLoad_Criteria(2);

			Residual_UTOL = log10(solver_container[ZONE_STRUC][MESH_0][FEA_SOL]->GetRes_FEM(0));
			Residual_RTOL = log10(solver_container[ZONE_STRUC][MESH_0][FEA_SOL]->GetRes_FEM(1));
			Residual_ETOL = log10(solver_container[ZONE_STRUC][MESH_0][FEA_SOL]->GetRes_FEM(2));

			meetCriteria = ( ( Residual_UTOL <  Criteria_UTOL ) &&
					 	 	 ( Residual_RTOL <  Criteria_RTOL ) &&
							 ( Residual_ETOL <  Criteria_ETOL ) );

			/*--- If the criteria is met and the load is not "too big", do the regular calculation ---*/
			if (meetCriteria){

				for (IntIter = 2; IntIter < config_container[ZONE_STRUC]->GetDyn_nIntIter(); IntIter++){

					for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

					/*--- Write the convergence history (only screen output) ---*/

					output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUC);

					config_container[iZone]->SetIntIter(IntIter);

					integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
																					config_container, RUNTIME_FEA_SYS, IntIter, iZone);

					}

					if (integration_container[ZONE_STRUC][FEA_SOL]->GetConvergence()) break;

				}

			}

			/*--- If the criteria is not met, a whole set of subiterations for the different loads must be done ---*/

			else {

				/*--- Here we have to restart the solution to the original one of the iteration ---*/
				/*--- Retrieve the Solution_Old as the current solution before subiterating ---*/

				for (iZone = nFluidZone; iZone < nTotalZone; iZone++)
					solver_container[iZone][MESH_0][FEA_SOL]->ResetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

				/*--- For the number of increments ---*/
				for (iIncrement = 0; iIncrement < nIncrements; iIncrement++){

					loadIncrement = (iIncrement + 1.0) * (1.0 / nIncrements);

					/*--- Set the load increment and the initial condition, and output the parameters of UTOL, RTOL, ETOL for the previous iteration ---*/

					for (iZone = nFluidZone; iZone < nTotalZone; iZone++){

						/*--- Set the convergence monitor to false, to force se solver to converge every subiteration ---*/
						integration_container[iZone][FEA_SOL]->SetConvergence(false);

						output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUC);

						/*--- FEA equations ---*/

						config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);


						solver_container[iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);
					}

					cout << endl;
					cout << "-- Incremental load: increment " << iIncrement + 1 << " ------------------------------------------" << endl;

					for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

						/*--- Set the value of the internal iteration ---*/
						IntIter = 0;
						config_container[iZone]->SetIntIter(IntIter);

						/*--- FEA equations ---*/

						config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

						/*--- Run the iteration ---*/

						integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
				                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);



					}

					/*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

					for (IntIter = 1; IntIter < config_container[ZONE_STRUC]->GetDyn_nIntIter(); IntIter++){

						for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

							/*--- Write the convergence history (only screen output) ---*/

							output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUC);

							config_container[iZone]->SetIntIter(IntIter);

							integration_container[iZone][FEA_SOL]->Structural_Iteration_FEM(geometry_container, solver_container, numerics_container,
						                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);

						}

						if (integration_container[ZONE_STRUC][FEA_SOL]->GetConvergence()) break;

					}

				}

			}

		}


	}
	else if ( (nonlinear && statTime) &&
			  ((first_iter && initial_calc) || (restart && initial_calc_restart))
			){

		/*--- We need to do the preprocessing to compute the Mass Matrix and integration constants ---*/
		for (iZone = nFluidZone; iZone < nTotalZone; iZone++)
			  solver_container[iZone][MESH_0][FEA_SOL]->Preprocessing(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0],
					  config_container[iZone], numerics_container[iZone][MESH_0][FEA_SOL], MESH_0, 0, RUNTIME_FEA_SYS, false);
	}

}


void FSI_Disp_Transfer(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
						 CTransfer*** transfer_container){


	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	bool MatchingMesh = config_container[ZONE_0]->GetMatchingMesh();

	/*--- Displacement transfer --  This will have to be modified for non-matching meshes ---*/

	if (MatchingMesh){
		solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetFlow_Displacement(geometry_container[ZONE_0], grid_movement[ZONE_0],
																   	   config_container[ZONE_0], config_container[ZONE_1],
																   	   geometry_container[ZONE_1], solver_container[ZONE_1]);
	}
	else{

		/*--- Transfer the information scattered (this is, each processor only receives the information it needs ---*/
		transfer_container[ZONE_1][ZONE_0]->Broadcast_InterfaceData_Matching(solver_container[ZONE_1][MESH_0][FEA_SOL],solver_container[ZONE_0][MESH_0][FLOW_SOL],
																		   geometry_container[ZONE_1][MESH_0],geometry_container[ZONE_0][MESH_0],
																		   config_container[ZONE_1], config_container[ZONE_0]);

		/*--- Set the volume deformation for the fluid zone ---*/
		grid_movement[ZONE_0]->SetVolume_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], true);

	}


}


void FSI_Load_Transfer(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
						 CTransfer*** transfer_container, unsigned long ExtIter){

	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	bool MatchingMesh = config_container[ZONE_0]->GetMatchingMesh();

	/*--- Load transfer --  This will have to be modified for non-matching meshes ---*/

	unsigned short SolContainer_Position_fea = config_container[ZONE_1]->GetContainerPosition(RUNTIME_FEA_SYS);

	/*--- FEA equations -- Necessary as the SetFEA_Load routine is as of now contained in the structural solver ---*/

	config_container[ZONE_1]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

	if (MatchingMesh){

		solver_container[ZONE_1][MESH_0][FEA_SOL]->SetFEA_Load(solver_container[ZONE_0], geometry_container[ZONE_1], geometry_container[ZONE_0],
															   config_container[ZONE_1], config_container[ZONE_0], numerics_container[ZONE_1][MESH_0][SolContainer_Position_fea][VISC_TERM]);

	}
	else{
//		solver_container[ZONE_1][MESH_0][FEA_SOL]->SetFEA_Load_Int(solver_container[ZONE_0], geometry_container[ZONE_1], geometry_container[ZONE_0],
//															   config_container[ZONE_1], config_container[ZONE_0], numerics_container[ZONE_1][MESH_0][SolContainer_Position_fea][VISC_TERM]);

		transfer_container[ZONE_0][ZONE_1]->Broadcast_InterfaceData_Matching(solver_container[ZONE_0][MESH_0][FLOW_SOL],solver_container[ZONE_1][MESH_0][FEA_SOL],
																		   geometry_container[ZONE_0][MESH_0],geometry_container[ZONE_1][MESH_0],
																		   config_container[ZONE_0], config_container[ZONE_1]);

	}


}


void FSI_Disp_Relaxation(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
							CConfig **config_container, unsigned long iFSIIter) {

	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	/*-------------------- Aitken's relaxation ------------------------*/

	/*--- Only one zone allowed for the fluid as for now ---*/
	unsigned short nFluidZone = 1, nStrucZone=1, nTotalZone;
	unsigned short iZone;

	nTotalZone = nFluidZone + nStrucZone;

	for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

		/*------------------- Compute the coefficient ---------------------*/

		solver_container[iZone][MESH_0][FEA_SOL]->ComputeAitken_Coefficient(geometry_container[iZone], config_container[iZone],
																		    solver_container[iZone], iFSIIter);

		/*----------------- Set the relaxation parameter ------------------*/

		solver_container[iZone][MESH_0][FEA_SOL]->SetAitken_Relaxation(geometry_container[iZone], config_container[iZone],
																	   solver_container[iZone]);

	}

	for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {
		solver_container[iZone][MESH_0][FEA_SOL]->Set_MPI_Solution_Pred_Old(geometry_container[iZone][MESH_0], config_container[iZone]);
	}

}

void FSI_Load_Relaxation(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
	     	 	 	 	 	CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
	     	 	 	 	 	CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox){


}


void FSI_Disp_Predictor(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox){

	#ifdef HAVE_MPI
		int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	/*--- Only one zone allowed for the fluid as for now ---*/
	unsigned short nFluidZone = 1, nStrucZone=1, nTotalZone;
	unsigned short iZone;

	nTotalZone = nFluidZone + nStrucZone;

	for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {

		solver_container[iZone][MESH_0][FEA_SOL]->PredictStruct_Displacement(geometry_container[iZone], config_container[iZone],
																		  solver_container[iZone]);

	}

	/*--- For parallel simulations we need to communicate the predicted solution before updating the fluid mesh ---*/

	for (iZone = nFluidZone; iZone < nTotalZone; iZone++) {
		solver_container[iZone][MESH_0][FEA_SOL]->Set_MPI_Solution_Pred(geometry_container[iZone][MESH_0], config_container[iZone]);
	}

}

void FSI_Load_Predictor(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
					     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
						 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
						 unsigned long ExtIter){


}
