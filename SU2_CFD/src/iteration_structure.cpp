/*!
 * \file iteration_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
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

#include "../include/iteration_structure.hpp"

void MeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {

	double Physical_dt, Physical_t;
	unsigned short iMesh; // Index for multi-grid level
    unsigned short iZone; // Index for zone of the mesh
	bool time_spectral = (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL);
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	if (time_spectral) nZone = config_container[ZONE_0]->GetnTimeInstances();
	bool relative_motion = config_container[ZONE_0]->GetRelative_Motion();
    unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
    unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#endif
    
  /*--- Initial set up for unsteady problems with dynamic meshes. ---*/
	for (iZone = 0; iZone < nZone; iZone++) {
		/*--- Dynamic mesh update ---*/
		if ((config_container[iZone]->GetGrid_Movement()) && (!time_spectral) && (config_container[ZONE_0]->GetKind_GridMovement(ZONE_0) != AEROELASTIC))
			SetGrid_Movement(geometry_container[iZone], surface_movement[iZone],
					grid_movement[iZone], FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, ExtIter);
	}

	/*--- If any relative motion between zones was found, perform a search
   and interpolation for any sliding interfaces before the next timestep. ---*/
	if (relative_motion) {
		SetSliding_Interfaces(geometry_container, solver_container, config_container, nZone);
    }

	for (iZone = 0; iZone < nZone; iZone++) {

		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

		/*--- Set the initial condition ---*/
		solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

		/*--- Update global parameters ---*/
		if (config_container[iZone]->GetKind_Solver() == EULER) config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == RANS) config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);

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
				integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
			}
		}

		/*--- Compute & store time-spectral source terms across all zones ---*/
		if (time_spectral)
			SetTimeSpectral(geometry_container, solver_container, config_container, nZone, (iZone+1)%nZone);

	}

	/*--- Dual time stepping strategy ---*/
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
        (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
        
		for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {
            
            /*--- Write the convergence history (only screen output) ---*/
            output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, ZONE_0);
            
            /*--- Set the value of the internal iteration ---*/
            config_container[ZONE_0]->SetIntIter(IntIter);
            
			/*--- All zones must be advanced and coupled with each pseudo timestep. ---*/
			for (iZone = 0; iZone < nZone; iZone++) {
                
				/*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/
				if (config_container[iZone]->GetKind_Solver() == EULER) config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
				if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
				if (config_container[iZone]->GetKind_Solver() == RANS) config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
				integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                            config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
                
				/*--- Pseudo-timestepping the turbulence model ---*/
				if (config_container[iZone]->GetKind_Solver() == RANS) {
					/*--- Turbulent model solution ---*/
					config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
					integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                                 config_container, RUNTIME_TURB_SYS, IntIter, iZone);
					if (config_container[iZone]->GetKind_Trans_Model() == LM) {
						/*--- Transition model solution ---*/
						config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
						integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                                      config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
					}
				}
                
				/*--- Call if AEROELASTIC motion was specified ---*/
				if ((config_container[ZONE_0]->GetGrid_Movement()) && (config_container[ZONE_0]->GetKind_GridMovement(ZONE_0) == AEROELASTIC)) {
					SetGrid_Movement(geometry_container[iZone], surface_movement[iZone],
                                     grid_movement[iZone], FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, IntIter);
                    /* If unsteady step converged, write out the plunge and pitch for that step */
                    int rank = MASTER_NODE;
                    #ifndef NO_MPI
                    rank = MPI::COMM_WORLD.Get_rank();
                    #endif
                    if (rank == MASTER_NODE && integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) {
                        std::fstream output_file;
                        output_file.open("plunging_pitching2.txt", std::fstream::in | std::fstream::out | std::fstream::app);
                        
                        output_file << std::setprecision(15) << config_container[ZONE_0]->GetAeroelastic_plunge() << "    " << config_container[ZONE_0]->GetAeroelastic_pitch() << "\n";
                        output_file.close();
                        
                        std::cout.precision(15);
                        std::cout << "plunge = " << config_container[ZONE_0]->GetAeroelastic_plunge() << std::endl;
                        std::cout << "pitch = " << config_container[ZONE_0]->GetAeroelastic_pitch() << std::endl;
                    }

				}
                
			}
            
            if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) break;
      
		}

		for (iZone = 0; iZone < nZone; iZone++) {
			/*--- Update dual time solver on all mesh levels ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
				integration_container[iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][FLOW_SOL], config_container[iZone]);
				integration_container[iZone][FLOW_SOL]->SetConvergence(false);
			}

			/*--- Update dual time solver for the turbulence model ---*/
			if (config_container[iZone]->GetKind_Solver() == RANS) {
				integration_container[iZone][TURB_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][TURB_SOL], config_container[iZone]);
				integration_container[iZone][TURB_SOL]->SetConvergence(false);
			}

			Physical_dt = config_container[iZone]->GetDelta_UnstTime();
			Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime())
				integration_container[iZone][FLOW_SOL]->SetConvergence(true);
		}
	}

}

void AdjMeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
  
	bool time_spectral = (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL);
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	if (time_spectral) nZone = config_container[ZONE_0]->GetnTimeInstances();
	bool relative_motion = config_container[ZONE_0]->GetRelative_Motion();
	int rank = MASTER_NODE;
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Initial set up for unsteady problems with dynamic meshes. ---*/
	for (iZone = 0; iZone < nZone; iZone++) {
		/*--- Dynamic mesh update ---*/
		if (config_container[iZone]->GetGrid_Movement() && !time_spectral)
			SetGrid_Movement(geometry_container[iZone], surface_movement[iZone],
					grid_movement[iZone], FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, ExtIter);
	}
	/*--- If any relative motion between zones was found, perform a search
   and interpolation for any sliding interfaces before the next timestep. ---*/
	if (relative_motion) {
		SetSliding_Interfaces(geometry_container, solver_container, config_container, nZone);
  }

	for (iZone = 0; iZone < nZone; iZone++) {
		if ((ExtIter == 0) || (config_container[iZone]->GetUnsteady_Simulation() && !time_spectral)) {
			if (rank == MASTER_NODE && iZone == ZONE_0)
				cout << "Single iteration of the direct solver to store flow data." << endl;
			if (config_container[iZone]->GetUnsteady_Simulation())
				solver_container[iZone][MESH_0][FLOW_SOL]->GetRestart(geometry_container[iZone][MESH_0], config_container[iZone], iZone);
		}
	}

	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Continuous adjoint Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations ---*/
		if ((ExtIter == 0) || config_container[iZone]->GetUnsteady_Simulation()) {
      
			if (config_container[iZone]->GetKind_Solver() == ADJ_EULER)         config_container[iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_FLOW_SYS, ExtIter);
			if (config_container[iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) config_container[iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
			if (config_container[iZone]->GetKind_Solver() == ADJ_RANS)          config_container[iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_FLOW_SYS, ExtIter);
      
      /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
			integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                  config_container, RUNTIME_FLOW_SYS, 0, iZone);
      
			if (config_container[iZone]->GetKind_Solver() == ADJ_RANS) {
        
        /*--- Solve the turbulence model ---*/
				config_container[iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_TURB_SYS, ExtIter);
				integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                     config_container, RUNTIME_TURB_SYS, IntIter, iZone);
        
        /*--- Solve transition model ---*/
        if (config_container[iZone]->GetKind_Trans_Model() == LM) {
          config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
          integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                        config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
        }
			}
      
			/*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
			 note that in the direct Euler problem we are not computing the gradients of the primitive variables ---*/
			if (config_container[iZone]->GetKind_Gradient_Method() == GREEN_GAUSS)
				solver_container[iZone][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry_container[iZone][MESH_0], config_container[iZone]);
			if (config_container[iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
				solver_container[iZone][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry_container[iZone][MESH_0], config_container[iZone]);
      
			/*--- Set contribution from cost function for boundary conditions ---*/
      for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
        
        /*--- Set the value of the non-dimensional coefficients in the coarse levels, using the fine level solution ---*/
        solver_container[iZone][iMesh][FLOW_SOL]->SetTotal_CDrag(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag());
        solver_container[iZone][iMesh][FLOW_SOL]->SetTotal_CLift(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CLift());
        solver_container[iZone][iMesh][FLOW_SOL]->SetTotal_CT(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CT());
        solver_container[iZone][iMesh][FLOW_SOL]->SetTotal_CQ(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CQ());
        
        /*--- Compute the adjoint boundary condition on Euler walls ---*/
        solver_container[iZone][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry_container[iZone][iMesh], solver_container[iZone][iMesh], config_container[iZone]);
        
        /*--- Set the internal boundary condition on nearfield surfaces ---*/
        if ((config_container[iZone]->GetKind_ObjFunc() == EQUIVALENT_AREA) || (config_container[iZone]->GetKind_ObjFunc() == NEARFIELD_PRESSURE))
          solver_container[iZone][iMesh][ADJFLOW_SOL]->SetIntBoundary_Jump(geometry_container[iZone][iMesh], solver_container[iZone][iMesh], config_container[iZone]);
      }
		}

		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
			IntIter = 0;
		}

		if (config_container[iZone]->GetKind_Solver() == ADJ_EULER) config_container[iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_ADJFLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) config_container[iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == ADJ_RANS) config_container[iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJFLOW_SYS, ExtIter);

		/*--- Iteration of the flow adjoint problem ---*/
		integration_container[iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
				config_container, RUNTIME_ADJFLOW_SYS, IntIter, iZone);

		/*--- Iteration of the turbulence model adjoint ---*/
		if ((config_container[iZone]->GetKind_Solver() == ADJ_RANS) && (((config_container[iZone]->GetKind_Adjoint() != HYBRID) && (!config_container[iZone]->GetFrozen_Visc())) || (config_container[iZone]->GetKind_Adjoint() == HYBRID))) {
			/*--- Turbulent model solution ---*/
			config_container[iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJTURB_SYS, ExtIter);
			integration_container[iZone][ADJTURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_ADJTURB_SYS, IntIter, iZone);

		}
	}

	/*--- Dual time stepping strategy ---*/
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
		for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      if (iZone == ZONE_0) output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, iZone);
      
      /*--- Set the value of the internal iteration ---*/
      config_container[ZONE_0]->SetIntIter(IntIter);

			/*--- All zones must be advanced and coupled with each pseudo timestep. ---*/
			for (iZone = 0; iZone < nZone; iZone++) {
				integration_container[iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
						config_container, RUNTIME_ADJFLOW_SYS, IntIter, iZone);

				if (integration_container[iZone][ADJFLOW_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout<<endl; break;}
			}
		}

		for (iZone = 0; iZone < nZone; iZone++) {
			/*--- Update dual time solver ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
				integration_container[iZone][ADJFLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][ADJFLOW_SOL], config_container[iZone]);
				integration_container[iZone][ADJFLOW_SOL]->SetConvergence(false);
			}

			Physical_dt = config_container[iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime()) integration_container[iZone][ADJFLOW_SOL]->SetConvergence(true);
		}
	}
}

void PlasmaIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {

	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

	/*--- Set the value of the internal iteration ---*/
	IntIter = ExtIter;
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

	/*--- Plasma solver with electric potential ---*/
	if (nZone > 1)
		solver_container[ZONE_1][MESH_0][ELEC_SOL]->Copy_Zone_Solution(solver_container[ZONE_1], geometry_container[ZONE_1], config_container[ZONE_1],
				solver_container[ZONE_0], geometry_container[ZONE_0], config_container[ZONE_0]);
  
	/*--- Plasma solver ---*/
	if (config_container[ZONE_0]->GetKind_Solver() == PLASMA_EULER) config_container[ZONE_0]->SetGlobalParam(PLASMA_EULER, RUNTIME_PLASMA_SYS, ExtIter);
	if (config_container[ZONE_0]->GetKind_Solver() == PLASMA_NAVIER_STOKES) config_container[ZONE_0]->SetGlobalParam(PLASMA_NAVIER_STOKES, RUNTIME_PLASMA_SYS, ExtIter);
	integration_container[ZONE_0][PLASMA_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_PLASMA_SYS, IntIter, ZONE_0);

	/*--- Electric potential solver ---*/
	if (nZone > 1) {

		if (config_container[ZONE_1]->GetKind_GasModel() == ARGON || config_container[ZONE_1]->GetKind_GasModel() == AIR21) {

			if (config_container[ZONE_1]->GetElectricSolver()) {
				solver_container[ZONE_0][MESH_0][PLASMA_SOL]->Copy_Zone_Solution(solver_container[ZONE_0], geometry_container[ZONE_0], config_container[ZONE_0],
						solver_container[ZONE_1], geometry_container[ZONE_1], config_container[ZONE_1]);
				config_container[ZONE_1]->SetGlobalParam(PLASMA_NAVIER_STOKES, RUNTIME_ELEC_SYS, ExtIter);
				integration_container[ZONE_1][ELEC_SOL]->SetPotential_Solver(geometry_container, solver_container, numerics_container,
						config_container, RUNTIME_ELEC_SYS, MESH_0, ZONE_1);

			}

		}
	}
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

    for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Set the value of the internal iteration ---*/
      config_container[ZONE_0]->SetIntIter(IntIter);

      /*--- Plasma solver ---*/
      if (config_container[ZONE_0]->GetKind_Solver() == PLASMA_EULER) config_container[ZONE_0]->SetGlobalParam(PLASMA_EULER, RUNTIME_PLASMA_SYS, ExtIter);
      if (config_container[ZONE_0]->GetKind_Solver() == PLASMA_NAVIER_STOKES) config_container[ZONE_0]->SetGlobalParam(PLASMA_NAVIER_STOKES, RUNTIME_PLASMA_SYS, ExtIter);
      integration_container[ZONE_0][PLASMA_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,config_container, RUNTIME_PLASMA_SYS, IntIter, ZONE_0);      

      if (integration_container[ZONE_0][PLASMA_SOL]->GetConvergence()) break;
    }
    
    /*--- Update dual time solver on all mesh levels ---*/
    integration_container[ZONE_0][PLASMA_SOL]->SetDualTime_Solver(geometry_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_0][PLASMA_SOL], config_container[ZONE_0]);
    integration_container[ZONE_0][PLASMA_SOL]->SetConvergence(false);
    
/*    Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
    Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
      integration_container[ZONE_0][FLOW_SOL]->SetConvergence(true);*/
    
  }
}

void AdjPlasmaIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
	int rank = MASTER_NODE;
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Continuous adjoint Euler equations ---*/
	if (ExtIter == 0 || config_container[ZONE_0]->GetUnsteady_Simulation()) {

		if (rank == MASTER_NODE) cout << "Iteration over the direct problem to store all flow information." << endl;

		/*--- Plasma equations ---*/
		if (config_container[ZONE_0]->GetKind_Solver() == ADJ_PLASMA_EULER) config_container[ZONE_0]->SetGlobalParam(ADJ_PLASMA_EULER, RUNTIME_PLASMA_SYS, ExtIter);
		if (config_container[ZONE_0]->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES) config_container[ZONE_0]->SetGlobalParam(ADJ_PLASMA_NAVIER_STOKES, RUNTIME_PLASMA_SYS, ExtIter);
		integration_container[ZONE_0][PLASMA_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                   config_container, RUNTIME_PLASMA_SYS, ExtIter, ZONE_0);
    solver_container[ZONE_0][MESH_0][ADJPLASMA_SOL]->SetForceProj_Vector(geometry_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_0], config_container[ZONE_0]);
    
	}

	/*--- Adjoint Plasma equations ---*/
	if (config_container[ZONE_0]->GetKind_Solver() == ADJ_PLASMA_EULER) config_container[ZONE_0]->SetGlobalParam(ADJ_PLASMA_EULER, RUNTIME_ADJPLASMA_SYS, ExtIter);
	if (config_container[ZONE_0]->GetKind_Solver() == ADJ_PLASMA_NAVIER_STOKES) config_container[ZONE_0]->SetGlobalParam(ADJ_PLASMA_NAVIER_STOKES, RUNTIME_ADJPLASMA_SYS, ExtIter);
	integration_container[ZONE_0][ADJPLASMA_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container, 
			config_container, RUNTIME_ADJPLASMA_SYS, ExtIter, ZONE_0);
}

void FreeSurfaceIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                          CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                          CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
  
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
        
        /*--- Set the initial condition ---*/
		solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
		
		/*--- Navier-Stokes equations ---*/
		if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(FREE_SURFACE_EULER, RUNTIME_FLOW_SYS, IntIter);
		if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(FREE_SURFACE_NAVIER_STOKES, RUNTIME_FLOW_SYS, IntIter);
		if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(FREE_SURFACE_RANS, RUNTIME_FLOW_SYS, IntIter);
		integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
		
        /*--- Turbulence model ---*/
		if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_RANS) {
			config_container[iZone]->SetGlobalParam(FREE_SURFACE_RANS, RUNTIME_TURB_SYS, IntIter);
			integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_TURB_SYS, IntIter, iZone);
		}
    
		/*--- Level-Set model solution ---*/
		if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(FREE_SURFACE_EULER, RUNTIME_LEVELSET_SYS, IntIter);
		if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(FREE_SURFACE_NAVIER_STOKES, RUNTIME_LEVELSET_SYS, IntIter);
		if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(FREE_SURFACE_RANS, RUNTIME_LEVELSET_SYS, IntIter);
		integration_container[iZone][LEVELSET_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_LEVELSET_SYS, IntIter, iZone);
		
		/*--- Dual time stepping strategy for the flow equations ---*/
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
      /*--- Loop over the internal itarations ---*/
			for(IntIter = 1; IntIter < config_container[iZone]->GetUnst_nIntIter(); IntIter++) {
        
        /*--- Write the convergence history (only screen output) ---*/
				output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, iZone);
        
        /*--- Set the value of the internal iteration ---*/
        config_container[iZone]->SetIntIter(IntIter);
        
				/*--- Euler, Navier-Stokes, and RANS equations solution ---*/
				if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(FREE_SURFACE_EULER, RUNTIME_FLOW_SYS, IntIter);
				if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(FREE_SURFACE_NAVIER_STOKES, RUNTIME_FLOW_SYS, IntIter);
        if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(FREE_SURFACE_RANS, RUNTIME_FLOW_SYS, IntIter);
				integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
        
        /*--- Turbulence model solution ---*/
        if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_RANS) {
          config_container[iZone]->SetGlobalParam(FREE_SURFACE_RANS, RUNTIME_TURB_SYS, IntIter);
          integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_TURB_SYS, IntIter, iZone);
        }
        
				/*--- Level-Set model solution ---*/
        if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(FREE_SURFACE_EULER, RUNTIME_LEVELSET_SYS, IntIter);
        if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(FREE_SURFACE_NAVIER_STOKES, RUNTIME_LEVELSET_SYS, IntIter);
        if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(FREE_SURFACE_RANS, RUNTIME_LEVELSET_SYS, IntIter);
        integration_container[iZone][LEVELSET_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_LEVELSET_SYS, IntIter, iZone);
        
				if (integration_container[iZone][FLOW_SOL]->GetConvergence()) break;
			}
      
			/*--- Set convergence the global convergence criteria to false, and dual time solution ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
				integration_container[iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][FLOW_SOL], config_container[iZone]);
				integration_container[iZone][FLOW_SOL]->SetConvergence(false);
			}
            
			/*--- Update dual time solver for the turbulence model ---*/
      if (config_container[iZone]->GetKind_Solver() == FREE_SURFACE_RANS) {
				integration_container[iZone][TURB_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][TURB_SOL], config_container[iZone]);
				integration_container[iZone][TURB_SOL]->SetConvergence(false);
			}
      
			integration_container[iZone][LEVELSET_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][LEVELSET_SOL], config_container[iZone]);
			integration_container[iZone][LEVELSET_SOL]->SetConvergence(false);
      
			/*--- Set the value of the global convergence criteria ---*/
			Physical_dt = config_container[iZone]->GetDelta_UnstTime();
			Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime())
				integration_container[iZone][FLOW_SOL]->SetConvergence(true);
		}
    
	}
  
}

void AdjFreeSurfaceIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {

	double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
	int rank = MASTER_NODE;
  
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	for (iZone = 0; iZone < nZone; iZone++) {

		/*--- Continuous adjoint free surface equations ---*/

		if (ExtIter == 0) {

			if (rank == MASTER_NODE) cout << "Iteration over the direct problem to store all flow information." << endl;

      /*--- Compute level set function using the distance to the free surface ---*/
      solver_container[iZone][MESH_0][LEVELSET_SOL]->SetLevelSet_Distance(geometry_container[iZone][MESH_0], config_container[iZone], true, true);
      
      /*--- Set the initial condition ---*/
      solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

      /*--- Navier-Stokes equations ---*/
      if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_EULER, RUNTIME_FLOW_SYS, IntIter);
      if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_NAVIER_STOKES, RUNTIME_FLOW_SYS, IntIter);
      if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_RANS, RUNTIME_FLOW_SYS, IntIter);
      integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_FLOW_SYS, IntIter, iZone);

      /*--- Turbulence model ---*/
      if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS) {
        config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_RANS, RUNTIME_TURB_SYS, IntIter);
        integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_TURB_SYS, IntIter, iZone);
      }
      
      /*--- Level-Set model solution ---*/
      if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_EULER, RUNTIME_LEVELSET_SYS, IntIter);
      if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_NAVIER_STOKES, RUNTIME_LEVELSET_SYS, IntIter);
      if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_RANS, RUNTIME_LEVELSET_SYS, IntIter);
      integration_container[iZone][LEVELSET_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_LEVELSET_SYS, IntIter, iZone);

			/*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
			 note that in the direct problem we are not computing the gradients ---*/
			if (config_container[iZone]->GetKind_Gradient_Method() == GREEN_GAUSS)
				solver_container[iZone][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry_container[iZone][MESH_0], config_container[iZone]);
			if (config_container[iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
				solver_container[iZone][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry_container[iZone][MESH_0], config_container[iZone]);

			/*--- Set contribution from cost function for boundary conditions ---*/
      for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
        solver_container[iZone][iMesh][FLOW_SOL]->SetTotal_CDrag(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag());
        solver_container[iZone][iMesh][FLOW_SOL]->SetTotal_CLift(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CLift());
        solver_container[iZone][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry_container[iZone][iMesh], solver_container[iZone][iMesh], config_container[iZone]);
      }
      
		}

		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
    
    /*--- Set the initial condition ---*/
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
    
    /*--- Euler, Navier-Stokes, and RANS adjoint equations solution ---*/
    if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_EULER, RUNTIME_ADJFLOW_SYS, IntIter);
    if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS, IntIter);
    if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_RANS, RUNTIME_ADJFLOW_SYS, IntIter);
		integration_container[iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_ADJFLOW_SYS, IntIter, iZone);

		/*--- Iteration over the level-set-adjoint equations ---*/
    if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_EULER, RUNTIME_ADJLEVELSET_SYS, IntIter);
    if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_NAVIER_STOKES, RUNTIME_ADJLEVELSET_SYS, IntIter);
    if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_RANS, RUNTIME_ADJLEVELSET_SYS, IntIter);
		integration_container[iZone][ADJLEVELSET_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container, config_container, RUNTIME_ADJLEVELSET_SYS, IntIter, iZone);

		/*--- Dual time stepping strategy for the flow equations ---*/
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
			for(IntIter = 1; IntIter < config_container[iZone]->GetUnst_nIntIter(); IntIter++) {

        /*--- Write the convergence history (only screen output) ---*/
        output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, iZone);
        
        /*--- Set the value of the internal iteration ---*/
        config_container[iZone]->SetIntIter(IntIter);
        
				/*--- Euler, Navier-Stokes, and RANS adjoint equations solution ---*/
        if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_EULER, RUNTIME_ADJFLOW_SYS, IntIter);
        if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS, IntIter);
        if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_RANS, RUNTIME_ADJFLOW_SYS, IntIter);
				integration_container[iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
						config_container, RUNTIME_ADJFLOW_SYS, IntIter, iZone);

				/*--- Level-Set model solution ---*/
        if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_EULER)          config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_EULER, RUNTIME_ADJLEVELSET_SYS, IntIter);
        if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_NAVIER_STOKES)  config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_NAVIER_STOKES, RUNTIME_ADJLEVELSET_SYS, IntIter);
        if (config_container[iZone]->GetKind_Solver() == ADJ_FREE_SURFACE_RANS)           config_container[iZone]->SetGlobalParam(ADJ_FREE_SURFACE_RANS, RUNTIME_ADJLEVELSET_SYS, IntIter);
				integration_container[iZone][ADJLEVELSET_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
						config_container, RUNTIME_ADJLEVELSET_SYS, IntIter, iZone);


				if (integration_container[iZone][ADJFLOW_SOL]->GetConvergence()) break;
			}

			/*--- Set convergence the global convergence criteria to false, and dual time solution ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
				integration_container[iZone][ADJFLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][ADJFLOW_SOL], config_container[iZone]);
				integration_container[iZone][ADJFLOW_SOL]->SetConvergence(false);
			}

			integration_container[iZone][ADJLEVELSET_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][ADJLEVELSET_SOL], config_container[iZone]);
			integration_container[iZone][ADJLEVELSET_SOL]->SetConvergence(false);

			/*--- Set the value of the global convergence criteria ---*/
			Physical_dt = config_container[iZone]->GetDelta_UnstTime();
			Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime())
				integration_container[iZone][ADJFLOW_SOL]->SetConvergence(true);
		}

	}	

}

void WaveIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {

	double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
	int rank = MASTER_NODE;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	for (iZone = 0; iZone < nZone; iZone++) {

		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

		/*--- Wave equations ---*/
		config_container[iZone]->SetGlobalParam(WAVE_EQUATION, RUNTIME_WAVE_SYS, ExtIter);
		integration_container[iZone][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
				config_container, RUNTIME_WAVE_SYS, IntIter, iZone);

		/*--- Dual time stepping strategy ---*/
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
			for(IntIter = 1; IntIter < config_container[iZone]->GetUnst_nIntIter(); IntIter++) {
        output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, iZone);
        config_container[iZone]->SetIntIter(IntIter);
				integration_container[iZone][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
						config_container, RUNTIME_WAVE_SYS, IntIter, iZone);
				if (integration_container[iZone][WAVE_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout<<endl; break;}
			}

			/*--- Update dual time solver ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
				integration_container[iZone][WAVE_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][WAVE_SOL], config_container[iZone]);
				integration_container[iZone][WAVE_SOL]->SetConvergence(false);
			}

			Physical_dt = config_container[iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime()) integration_container[iZone][WAVE_SOL]->SetConvergence(true);
		}

	}

}

void FEAIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
	int rank = MASTER_NODE;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	for (iZone = 0; iZone < nZone; iZone++) {

		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

		/*--- Set the initial condition at the first iteration ---*/
		solver_container[iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

		/*--- FEA equations ---*/
		config_container[iZone]->SetGlobalParam(LINEAR_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);
		integration_container[iZone][FEA_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
				config_container, RUNTIME_FEA_SYS, IntIter, iZone);

		/*--- Dual time stepping strategy ---*/
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

			for(IntIter = 1; IntIter < config_container[iZone]->GetUnst_nIntIter(); IntIter++) {
        output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, iZone);
        config_container[iZone]->SetIntIter(IntIter);
				integration_container[iZone][FEA_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
						config_container, RUNTIME_FEA_SYS, IntIter, iZone);
				if (integration_container[iZone][FEA_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout << endl; break;}
			}

			/*--- Update dual time solver ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
				integration_container[iZone][FEA_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][FEA_SOL], config_container[iZone]);
				integration_container[iZone][FEA_SOL]->SetConvergence(false);
			}

			Physical_dt = config_container[iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime()) integration_container[iZone][FEA_SOL]->SetConvergence(true);
		}

	}	

}

void FluidStructureIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	double Physical_dt, Physical_t;
	unsigned short iMesh;
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

	/*--- Set the value of the internal iteration ---*/
	IntIter = ExtIter;
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

	/*--- Euler equations ---*/
	config_container[ZONE_0]->SetGlobalParam(FLUID_STRUCTURE_EULER, RUNTIME_FLOW_SYS, ExtIter);
	integration_container[ZONE_0][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_FLOW_SYS, IntIter, ZONE_0);

	/*--- Update loads for the FEA model ---*/
	solver_container[ZONE_1][MESH_0][FEA_SOL]->SetFEA_Load(solver_container[ZONE_0], geometry_container[ZONE_1], geometry_container[ZONE_0], config_container[ZONE_1], config_container[ZONE_0]);

	/*--- FEA model solution ---*/
	config_container[ZONE_1]->SetGlobalParam(FLUID_STRUCTURE_EULER, RUNTIME_FEA_SYS, ExtIter);
	integration_container[ZONE_1][FEA_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_FEA_SYS, IntIter, ZONE_1);

	/*--- Update the FEA geometry (ZONE 1), and the flow geometry (ZONE 0) --*/
	solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetFlow_Displacement(geometry_container[ZONE_0], grid_movement[ZONE_0],
			config_container[ZONE_0], config_container[ZONE_1],
			geometry_container[ZONE_1], solver_container[ZONE_1]);

	/*--- Dual time stepping strategy for the coupled system ---*/
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

		for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {

      /*--- Write the convergence history (only screen output) ---*/
			output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, ZONE_0);
      
      /*--- Set the value of the internal iteration ---*/
      config_container[ZONE_0]->SetIntIter(IntIter);
      
			/*--- Euler equations ---*/
			config_container[ZONE_0]->SetGlobalParam(FLUID_STRUCTURE_EULER, RUNTIME_FLOW_SYS, ExtIter);
			integration_container[ZONE_0][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_FLOW_SYS, IntIter, ZONE_0);

			/*--- Update loads for the FEA model ---*/
			solver_container[ZONE_1][MESH_0][FEA_SOL]->SetFEA_Load(solver_container[ZONE_0], geometry_container[ZONE_1], geometry_container[ZONE_0], config_container[ZONE_1], config_container[ZONE_0]);

			/*--- FEA model solution ---*/
			config_container[ZONE_1]->SetGlobalParam(FLUID_STRUCTURE_EULER, RUNTIME_FEA_SYS, ExtIter);
			integration_container[ZONE_1][FEA_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_FEA_SYS, IntIter, ZONE_1);

			/*--- Update the FEA geometry (ZONE 1), and the flow geometry (ZONE 0) --*/
			solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetFlow_Displacement(geometry_container[ZONE_0], grid_movement[ZONE_0],
					config_container[ZONE_0], config_container[ZONE_1],
					geometry_container[ZONE_1], solver_container[ZONE_1]);

			if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) break;
		}

		/*--- Set convergence the global convergence criteria to false, and dual time solution ---*/
		for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
			integration_container[ZONE_0][FLOW_SOL]->SetDualTime_Solver(geometry_container[ZONE_0][iMesh], solver_container[ZONE_0][iMesh][FLOW_SOL], config_container[ZONE_0]);
			integration_container[ZONE_0][FLOW_SOL]->SetConvergence(false);
		}

		integration_container[ZONE_1][FEA_SOL]->SetDualTime_Solver(geometry_container[ZONE_1][MESH_0], solver_container[ZONE_1][MESH_0][FEA_SOL], config_container[ZONE_1]);
		integration_container[ZONE_1][FEA_SOL]->SetConvergence(false);

		/*--- Set the value of the global convergence criteria ---*/
		Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
		Physical_t  = (ExtIter+1)*Physical_dt;
		if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
			integration_container[ZONE_0][FLOW_SOL]->SetConvergence(true);
	}

}

void AeroacousticIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	double Physical_dt, Physical_t;
	unsigned short iMesh;
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
	/*--- Set the value of the internal iteration ---*/
	IntIter = ExtIter;
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

	/*--- Euler equations ---*/
	config_container[ZONE_0]->SetGlobalParam(AEROACOUSTIC_EULER, RUNTIME_FLOW_SYS, ExtIter);
	integration_container[ZONE_0][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_FLOW_SYS, IntIter, ZONE_0);

	/*--- Update the noise source terms for the wave model ---*/
	solver_container[ZONE_1][MESH_0][WAVE_SOL]->SetNoise_Source(solver_container[ZONE_0], geometry_container[ZONE_1], config_container[ZONE_1]);

	/*--- Wave model solution ---*/
	config_container[ZONE_1]->SetGlobalParam(AEROACOUSTIC_EULER, RUNTIME_WAVE_SYS, ExtIter);
	integration_container[ZONE_1][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_WAVE_SYS, IntIter, ZONE_1);

	/*--- Dual time stepping strategy for the coupled system ---*/
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

		for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {

      /*--- Write the convergence history (only screen output) ---*/
			output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, ZONE_0);
      
      /*--- Set the value of the internal iteration ---*/
      config_container[ZONE_0]->SetIntIter(IntIter);
      
			/*--- Euler equations ---*/
			config_container[ZONE_0]->SetGlobalParam(AEROACOUSTIC_EULER, RUNTIME_FLOW_SYS, ExtIter);
			integration_container[ZONE_0][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_FLOW_SYS, IntIter, ZONE_0);

			/*--- Update source terms for the wave model ---*/
			solver_container[ZONE_1][MESH_0][WAVE_SOL]->SetNoise_Source(solver_container[ZONE_0], geometry_container[ZONE_1], config_container[ZONE_1]);

			/*--- Wave model solution ---*/
			config_container[ZONE_1]->SetGlobalParam(AEROACOUSTIC_EULER, RUNTIME_WAVE_SYS, ExtIter);
			integration_container[ZONE_1][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_WAVE_SYS, IntIter, ZONE_1);

			if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) break;
		}

		/*--- Set convergence the global convergence criteria to false, and dual time solution ---*/
		for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
			integration_container[ZONE_0][FLOW_SOL]->SetDualTime_Solver(geometry_container[ZONE_0][iMesh], solver_container[ZONE_0][iMesh][FLOW_SOL], config_container[ZONE_0]);
			integration_container[ZONE_0][FLOW_SOL]->SetConvergence(false);
		}

		integration_container[ZONE_1][WAVE_SOL]->SetDualTime_Solver(geometry_container[ZONE_1][MESH_0], solver_container[ZONE_1][MESH_0][WAVE_SOL], config_container[ZONE_1]);
		integration_container[ZONE_1][WAVE_SOL]->SetConvergence(false);

		/*--- Perform mesh motion for flow problem only, if necessary ---*/
		if (config_container[ZONE_0]->GetGrid_Movement())
			SetGrid_Movement(geometry_container[ZONE_0], surface_movement[ZONE_0],
					grid_movement[ZONE_0], FFDBox[ZONE_0], solver_container[ZONE_0], config_container[ZONE_0], ZONE_0, ExtIter);

		/*--- Set the value of the global convergence criteria ---*/
		Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
		Physical_t  = (ExtIter+1)*Physical_dt;
		if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
			integration_container[ZONE_0][FLOW_SOL]->SetConvergence(true);
	}	

}

void AdjAeroacousticIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container, 
		CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
		CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	double Physical_dt, Physical_t;
	unsigned short iMesh;
	int rank = MASTER_NODE;
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Continuous adjoint Euler equations + continuous adjoint wave equation ---*/

	if (rank == MASTER_NODE) cout << "Iteration over the direct problem to store all flow information." << endl;

	/*--- Load direct solutions for the flow and wave problems from file in reverse time ---*/
	solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRestart(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ZONE_0);
	solver_container[ZONE_1][MESH_0][WAVE_SOL]->GetRestart(geometry_container[ZONE_1][MESH_0], config_container[ZONE_1], ZONE_1);

	if (config_container[ZONE_0]->GetKind_Solver() == ADJ_AEROACOUSTIC_EULER)
		config_container[ZONE_0]->SetGlobalParam(ADJ_AEROACOUSTIC_EULER, RUNTIME_FLOW_SYS, ExtIter);

	/*--- Run one iteration of the direct flow problem to store needed flow variables ---*/
	integration_container[ZONE_0][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_FLOW_SYS, 0, ZONE_0);

	/*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
	 note that in the direct problem we are not computing the gradients ---*/	
	if (config_container[ZONE_0]->GetKind_Gradient_Method() == GREEN_GAUSS)
		solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_GG(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);
	if (config_container[ZONE_0]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
		solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimVar_Gradient_LS(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);

	/*--- Set contribution from cost function for boundary conditions ---*/
  for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
    solver_container[ZONE_0][iMesh][FLOW_SOL]->SetTotal_CDrag(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag());
    solver_container[ZONE_0][iMesh][FLOW_SOL]->SetTotal_CLift(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift());
    solver_container[ZONE_0][iMesh][FLOW_SOL]->SetTotal_CT(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CT());
    solver_container[ZONE_0][iMesh][FLOW_SOL]->SetTotal_CQ(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CQ());
    solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry_container[ZONE_0][iMesh], solver_container[ZONE_0][iMesh], config_container[ZONE_0]);
    if ((config_container[ZONE_0]->GetKind_ObjFunc() == EQUIVALENT_AREA) || (config_container[ZONE_0]->GetKind_ObjFunc() == NEARFIELD_PRESSURE))
      solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetIntBoundary_Jump(geometry_container[ZONE_0][iMesh], solver_container[ZONE_0][iMesh], config_container[ZONE_0]);
  }

	/*--- Set the value of the internal iteration ---*/
	IntIter = ExtIter;
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
		IntIter = 0;
	}

	if (config_container[ZONE_0]->GetKind_Solver() == ADJ_AEROACOUSTIC_EULER)
		config_container[ZONE_0]->SetGlobalParam(ADJ_AEROACOUSTIC_EULER, RUNTIME_ADJFLOW_SYS, ExtIter);

	/*--- Adjoint Wave Solver (Note that we use the same solver for the direct and adjoint problems) ---*/
	config_container[ZONE_1]->SetGlobalParam(ADJ_AEROACOUSTIC_EULER, RUNTIME_WAVE_SYS, ExtIter);
	integration_container[ZONE_1][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_WAVE_SYS, IntIter, ZONE_1);

	/*--- Update aeroacoustic adjoint coupling terms ---*/
	solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetAeroacoustic_Coupling(solver_container[ZONE_1], solver_container[ZONE_0], numerics_container[ZONE_0][MESH_0][ADJFLOW_SOL][CONV_TERM], geometry_container[ZONE_0], config_container[ZONE_0]);

	/*--- Adjoint Flow Solver ---*/
	integration_container[ZONE_0][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
			config_container, RUNTIME_ADJFLOW_SYS, IntIter, ZONE_0);

	/*--- Dual time stepping strategy ---*/
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

		for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {

      /*--- Write the convergence history (only screen output) ---*/
			output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, ZONE_0);
      
      /*--- Set the value of the internal iteration ---*/
      config_container[ZONE_0]->SetIntIter(IntIter);
      
			/*--- Adjoint Wave Solver ---*/
			config_container[ZONE_1]->SetGlobalParam(ADJ_AEROACOUSTIC_EULER, RUNTIME_WAVE_SYS, ExtIter);
			integration_container[ZONE_1][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_WAVE_SYS, IntIter, ZONE_1);

			/*--- Update aeroacoustic adjoint coupling terms ---*/
			solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetAeroacoustic_Coupling(solver_container[ZONE_1], solver_container[ZONE_0], numerics_container[ZONE_0][MESH_0][ADJFLOW_SOL][CONV_TERM], geometry_container[ZONE_0], config_container[ZONE_0]);

			/*--- Adjoint Flow Solver ---*/
			integration_container[ZONE_0][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_ADJFLOW_SYS, IntIter, ZONE_0);

			if (integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence()) {if (rank == MASTER_NODE) cout<<endl; break;}
		}

		/*--- Update dual time solver ---*/
		for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetMGLevels(); iMesh++) {
			integration_container[ZONE_0][ADJFLOW_SOL]->SetDualTime_Solver(geometry_container[ZONE_0][iMesh], solver_container[ZONE_0][iMesh][ADJFLOW_SOL], config_container[ZONE_0]);
			integration_container[ZONE_0][ADJFLOW_SOL]->SetConvergence(false);
		}

		integration_container[ZONE_1][WAVE_SOL]->SetDualTime_Solver(geometry_container[ZONE_1][MESH_0], solver_container[ZONE_1][MESH_0][WAVE_SOL], config_container[ZONE_1]);
		integration_container[ZONE_1][WAVE_SOL]->SetConvergence(false);

		/*--- Perform mesh motion, if necessary ---*/
		if (config_container[ZONE_0]->GetGrid_Movement())
			SetGrid_Movement(geometry_container[ZONE_0], surface_movement[ZONE_0],
					grid_movement[ZONE_0], FFDBox[ZONE_0], solver_container[ZONE_0],config_container[ZONE_0], ZONE_0, ExtIter);

		Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
		if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime()) integration_container[ZONE_0][ADJFLOW_SOL]->SetConvergence(true);
	}

}

void SetGrid_Movement(CGeometry **geometry_container, CSurfaceMovement *surface_movement,
                      CVolumetricMovement *grid_movement, CFreeFormDefBox **FFDBox,
                      CSolver ***solver_container, CConfig *config_container, unsigned short iZone, unsigned long ExtIter)   {
  
	unsigned short Kind_Grid_Movement = config_container->GetKind_GridMovement(iZone);
  bool adjoint = config_container->GetAdjoint();
	bool time_spectral = (config_container->GetUnsteady_Simulation() == TIME_SPECTRAL);
  
	/*--- For a time-spectral case, set "iteration number" to the zone number,
   so that the meshes are positioned correctly for each instance. ---*/
	if (time_spectral) {
		ExtIter = iZone;
		Kind_Grid_Movement = config_container->GetKind_GridMovement(ZONE_0);
	}
  
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	if ((rank == MASTER_NODE) && (Kind_Grid_Movement != AEROELASTIC)) {
		if (geometry_container[MESH_0]->GetnZone() > 1)
			cout << endl << "Performing the dynamic mesh update for Zone " << iZone << "." << endl;
		else
			cout << endl << "Performing the dynamic mesh update." << endl;
	}
  
	/*--- Perform mesh movement depending on specified type ---*/
	switch (Kind_Grid_Movement) {
      
    case RIGID_MOTION:
      
      /*--- Move each point in the volume mesh using the specified type of rigid
       mesh motion. These routines also compute the analytic grid velocities. ---*/
      if (rank == MASTER_NODE)
        cout << "Updating vertex locations through rigid mesh motion." << endl;
      grid_movement->SetRigidTranslation(geometry_container[MESH_0], config_container, iZone, ExtIter);
      grid_movement->SetRigidPlunging(geometry_container[MESH_0], config_container, iZone, ExtIter);
      grid_movement->SetRigidPitching(geometry_container[MESH_0], config_container, iZone, ExtIter);
      grid_movement->SetRigidRotation(geometry_container[MESH_0], config_container, iZone, ExtIter);
      
      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/
      for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
        /*--- Create the control volume structures ---*/
        geometry_container[iMGlevel]->SetControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
        geometry_container[iMGlevel]->SetBoundControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
        geometry_container[iMGlevel]->SetCoord(geometry_container[iMGlevel-1]);
        geometry_container[iMGlevel]->SetRestricted_GridVelocity(geometry_container[iMGlevel-1], config_container, ExtIter);
      }
      
      /*--- Adjoint mesh velocities: compute an approximation of the
       mesh velocities on the coarser levels starting from the values
       read in from the restart file for the finest mesh. ---*/
      if (adjoint) {
        for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
          /*---Set the grid velocity on the coarser levels ---*/
          geometry_container[iMGlevel]->SetRestricted_GridVelocity(geometry_container[iMGlevel-1], config_container, ExtIter);
        }
      }
      
      break;
      
    case EXTERNAL: case EXTERNAL_ROTATION:
      
      /*--- Apply rigid rotation first, if necessary ---*/
      if (Kind_Grid_Movement == EXTERNAL_ROTATION) {
        /*--- Rotate each point in the volume mesh (no deformation) ---*/
        if (rank == MASTER_NODE)
          cout << "Updating vertex locations by rigid rotation." << endl;
        grid_movement->SetRigidRotation(geometry_container[MESH_0], config_container, iZone, ExtIter);
      }
      
      /*--- Surface grid deformation read from external files ---*/
      if (rank == MASTER_NODE)
        cout << "Updating surface locations from the mesh motion file." << endl;
      surface_movement->SetExternal_Deformation(geometry_container[MESH_0], config_container, iZone, ExtIter);
      
      /*--- Volume grid deformation ---*/
      if (rank == MASTER_NODE)
        cout << "Deforming the volume grid using the spring analogy." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0], config_container, true);
      
      /*--- Update the multigrid structure after moving the finest grid ---*/
      for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
        /*--- Create the control volume structures ---*/
        geometry_container[iMGlevel]->SetControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
        geometry_container[iMGlevel]->SetBoundControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
        geometry_container[iMGlevel]->SetCoord(geometry_container[iMGlevel-1]);
      }
      
      if (!adjoint) {
        /*--- Update grid node velocities on all mesh levels ---*/
        if (rank == MASTER_NODE)
          cout << "Computing mesh velocities using a finite diff. approx." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
        for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
          /*---Set the grid velocity on the coarser levels ---*/
          geometry_container[iMGlevel]->SetGridVelocity(config_container,ExtIter);
        }
      } else {
        /*--- Adjoint mesh velocities: compute an approximation of the
         mesh velocities on the coarser levels starting from the values
         read in from the restart file for the finest mesh. ---*/
        for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
          /*---Set the grid velocity on the coarser levels ---*/
          geometry_container[iMGlevel]->SetRestricted_GridVelocity(geometry_container[iMGlevel-1], config_container, ExtIter);
        }
      }
      
      break;
      
    case FLUTTER:
      
      /*--- Surface grid deformation ---*/
      if (rank == MASTER_NODE)
        cout << "Updating flutter surfaces locations." << endl;
      if (geometry_container[MESH_0]->GetnDim() == 2)
        surface_movement->SetBoundary_Flutter2D(geometry_container[MESH_0], config_container, ExtIter);
      else
        surface_movement->SetBoundary_Flutter3D(geometry_container[MESH_0], config_container, FFDBox, ExtIter);
      
      /*--- Volume grid deformation ---*/
      if (rank == MASTER_NODE)
        cout << "Deforming the volume grid using the spring analogy." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0], config_container, true);
      
      /*--- Update the multigrid structure after moving the finest grid ---*/
      for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
        /*--- Create the control volume structures ---*/
        geometry_container[iMGlevel]->SetControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
        geometry_container[iMGlevel]->SetBoundControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
        geometry_container[iMGlevel]->SetCoord(geometry_container[iMGlevel-1]);
      }
      
      if (!adjoint) {
        /*--- Update grid node velocities on all mesh levels ---*/
        if (rank == MASTER_NODE)
          cout << "Computing mesh velocities using a finite diff. approx." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
        for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
          /*---Set the grid velocity on the coarser levels ---*/
          geometry_container[iMGlevel]->SetGridVelocity(config_container,ExtIter);
        }
      } else {
        /*--- Adjoint mesh velocities: compute an approximation of the
         mesh velocities on the coarser levels starting from the values
         read in from the restart file for the finest mesh. ---*/
        for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
          /*---Set the grid velocity on the coarser levels ---*/
          geometry_container[iMGlevel]->SetRestricted_GridVelocity(geometry_container[iMGlevel-1], config_container, ExtIter);
        }
      }
      
      break;
      
    case AEROELASTIC:
    
            /*--- Variables used for Aeroelastic case ---*/
            double Cl, Cm, Cl_proc, Cm_proc, pitch, plunge;
            double structural_solution[4]; //contains solution of typical section wing model.
            unsigned short Aeroelastic_Type, Aeroelastic_Grid_Vel;
            Aeroelastic_Type = config_container->GetType_Aeroelastic();
            Aeroelastic_Grid_Vel = config_container->GetAeroelastic_GridVelocity();
            
            /*--- Reset value of plunge and pitch for the new unsteady step ---*/
            if (ExtIter == 1) {
                pitch = 0.0;
                plunge = 0.0;
                config_container->SetAeroelastic_pitch(pitch);
                config_container->SetAeroelastic_plunge(plunge);
            }
            
            /*--- Forces per processor ---*/
            Cl_proc = solver_container[MESH_0][FLOW_SOL]->GetTotal_CLift();
            Cm_proc = -1.0*solver_container[MESH_0][FLOW_SOL]->GetTotal_CMz();
            
#ifndef NO_MPI
            /*--- Add the forces over all the processors ---*/
            MPI::COMM_WORLD.Allreduce(&Cl_proc, &Cl, 1, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Allreduce(&Cm_proc, &Cm, 1, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Barrier();
#else
            /*--- Set the forces to the forces on the sole processor ---*/
            Cl = Cl_proc;
            Cm = Cm_proc;
#endif
            
            /*--- Use the if statement to move the grid only at selected dual time step iterations. ---*/
            if (ExtIter % 3 ==0) {
                
                /*--- Solve typical section wing model, solution returned in structural solution ---*/
                grid_movement->SolveTypicalSectionWingModel(geometry_container[MESH_0], Cl, Cm, config_container, iZone, ExtIter, structural_solution);
                
                /*--- Check type of Aeroelastic mesh motion ---*/
                if (Aeroelastic_Type == RIGID)
                    grid_movement->AeroelasticRigid(geometry_container[MESH_0], config_container, iZone, structural_solution);
                if (Aeroelastic_Type == DEFORM) {
                    Aeroelastic_Grid_Vel = FD;  // For deforming mesh can't do analytic mesh velocities.
                    grid_movement->AeroelasticDeform(geometry_container[MESH_0], config_container, iZone, structural_solution);
                }
                
                /*--- Update the multigrid structure after moving the finest grid ---*/
                for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
                    /*--- Create the control volume structures ---*/
                    geometry_container[iMGlevel]->SetControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
                    geometry_container[iMGlevel]->SetBoundControlVolume(config_container,geometry_container[iMGlevel-1], UPDATE);
                    geometry_container[iMGlevel]->SetCoord(geometry_container[iMGlevel-1]);
                    if (Aeroelastic_Grid_Vel == ANALYTIC) {
                        /*---Set the grid velocity on the coarser levels ---*/
                        geometry_container[iMGlevel]->SetRestricted_GridVelocity(geometry_container[iMGlevel-1], config_container, ExtIter);
                    }
                }
                
                if (Aeroelastic_Grid_Vel == FD) {
                    /* Grid Velocity by finite Difference */
                    geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
                    for (unsigned short iMGlevel = 1; iMGlevel <= config_container->GetMGLevels(); iMGlevel++) {
                        /*---Set the grid velocity on the coarser levels ---*/
                        geometry_container[iMGlevel]->SetGridVelocity(config_container,ExtIter);
                    }
                }
            }
      break;
      
    case NONE: default:
      
      /*--- There is no mesh motion specified for this zone. ---*/
      if (rank == MASTER_NODE)
        cout << "No mesh motion specified." << endl;
      
      break;
  }
  
}

void SetTimeSpectral(CGeometry ***geometry_container, CSolver ****solver_container,
		CConfig **config_container, unsigned short nZone, unsigned short iZone) {

	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif

	/*--- Local variables and initialization ---*/
	unsigned short iVar, kZone, jZone, iMGlevel;
	unsigned short nVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
	unsigned long iPoint;
	bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool adjoint = (config_container[ZONE_0]->GetAdjoint());
	if (adjoint) {
		implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	}

	/*--- Retrieve values from the config file ---*/
	double *U      	= new double[nVar];
	double *U_old  	= new double[nVar];
	double *Psi	   	= new double[nVar];
	double *Psi_old	= new double[nVar];
	double *Source = new double[nVar];
	double deltaU, deltaPsi;
	//	double Omega[3], Center[3], Ampl[3];
	//	double Lref   = config_container[ZONE_0]->GetLength_Ref();
	//	double current_time;
	//	double alpha_dot[3];
	//  double DEG2RAD = PI_NUMBER/180.0;
	//	double GridVel[3];

//	/*--- rotational velocity ---*/
//	Omega[0]  = config_container[ZONE_0]->GetPitching_Omega_X(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	Omega[1]  = config_container[ZONE_0]->GetPitching_Omega_Y(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	Omega[2]  = config_container[ZONE_0]->GetPitching_Omega_Z(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	double Omega_mag = sqrt(pow(Omega[0],2)+pow(Omega[1],2)+pow(Omega[2],2));
//
//	/*--- amplitude of motion ---*/
//	Ampl[0]   = config_container[ZONE_0]->GetPitching_Ampl_X(ZONE_0)*DEG2RAD;
//	Ampl[1]   = config_container[ZONE_0]->GetPitching_Ampl_Y(ZONE_0)*DEG2RAD;
//	Ampl[2]   = config_container[ZONE_0]->GetPitching_Ampl_Z(ZONE_0)*DEG2RAD;

	/*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
	double period = config_container[ZONE_0]->GetTimeSpectral_Period();
	//	double deltaT = period/(double)(config_container[ZONE_0]->GetnTimeInstances());

//	/*--- For now, hard code the origin for the pitching airfoil. ---*/
//	Center[0] = config_container[ZONE_0]->GetMotion_Origin_X(ZONE_0);
//	Center[1] = config_container[ZONE_0]->GetMotion_Origin_Y(ZONE_0);
//	Center[2] = config_container[ZONE_0]->GetMotion_Origin_Z(ZONE_0);

	//	/*--- Compute the Mesh Velocities ---*/
	//	if (config_container[ZONE_0]->GetExtIter() == 0) {
	//
	////
	////		/*--- Loop over all grid levels ---*/
	////		for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {
	////
	////			/*--- Loop over each zone ---*/
	////			for (iZone = 0; iZone < nZone; iZone++) {
	////
	////				/*--- Time at current instance, i.e. at the current zone ---*/
	////				current_time = (static_cast<double>(iZone))*deltaT;
	////
	////				/*--- Angular velocity at the given time instance ---*/
	////				alpha_dot[0] = -Omega[0]*Ampl[0]*cos(Omega[0]*current_time);
	////				alpha_dot[1] = -Omega[1]*Ampl[1]*cos(Omega[1]*current_time);
	////				alpha_dot[2] = -Omega[2]*Ampl[2]*cos(Omega[2]*current_time);
	////
	////				/*--- Loop over each node in the volume mesh ---*/
	////				for (iPoint = 0; iPoint < geometry_container[iZone][iMGlevel]->GetnPoint(); iPoint++) {
	////
	////					/*--- Coordinates of the current point ---*/
	////					Coord = geometry_container[iZone][iMGlevel]->node[iPoint]->GetCoord();
	////
	////					/*--- Calculate non-dim. position from rotation center ---*/
	////					double r[3] = {0,0,0};
	////					for (iDim = 0; iDim < nDim; iDim++)
	////						r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
	////					if (nDim == 2) r[nDim] = 0.0;
	////
	////					/*--- Cross Product of angular velocity and distance from center ---*/
	////					GridVel[0] = alpha_dot[1]*r[2] - alpha_dot[2]*r[1];
	////					GridVel[1] = alpha_dot[2]*r[0] - alpha_dot[0]*r[2];
	////					GridVel[2] = alpha_dot[0]*r[1] - alpha_dot[1]*r[0];
	////
	////					// PRINT STATEMENTS FOR DEBUGGING!
	////					if (rank == MASTER_NODE && iMGlevel == 0 && iZone == 0) {
	////						if (iPoint == 440) {
	////							cout << "Analytical Mesh Velocity = (" <<  GridVel[0] << ", " << GridVel[1] << ", " << GridVel[2] << ")" << endl;
	////						}
	////					}
	////					/*--- Set Grid Velocity for the point in the given zone ---*/
	////					for(iDim = 0; iDim < nDim; iDim++) {
	////
	////						/*--- Store grid velocity for this point ---*/
	////						geometry_container[iZone][iMGlevel]->node[iPoint]->SetGridVel(iDim, GridVel[iDim]);
	////
	////					}
	////				}
	////			}
	////		}
	//
	//		/*--- by fourier interpolation ---*/
	//
	//	}




	/*--- allocate dynamic memory for D ---*/
	double **D = new double*[nZone];
	for (kZone = 0; kZone < nZone; kZone++) {
		D[kZone] = new double[nZone];
	}

	/*--- Build the time-spectral operator matrix ---*/
	for (kZone = 0; kZone < nZone; kZone++) {
		for (jZone = 0; jZone < nZone; jZone++) {

			/*--- For an even number of time instances ---*/
			if (nZone%2 == 0) {
				if (kZone == jZone) {
					D[kZone][jZone] = 0.0;
				} else {
					D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/tan(PI_NUMBER*(kZone-jZone)/nZone));
				}
			} else {
				/*--- For an odd number of time instances ---*/
				if (kZone == jZone) {
					D[kZone][jZone] = 0.0;
				} else {
					D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/sin(PI_NUMBER*(kZone-jZone)/nZone));
				}
			}
		}
	}

	/*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
	/*--- Loop over all grid levels ---*/
	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {

		/*--- Loop over each node in the volume mesh ---*/
		for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {

			for (iVar = 0; iVar < nVar; iVar++) {
				Source[iVar] = 0.0;
			}

			/*--- Step across the columns ---*/
			for (jZone = 0; jZone < nZone; jZone++) {

				/*--- Retrieve solution at this node in current zone ---*/
				for (iVar = 0; iVar < nVar; iVar++) {

					if (!adjoint) {

						U[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution(iVar);
						Source[iVar] += U[iVar]*D[iZone][jZone];

						if (implicit) {
							U_old[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution_Old(iVar);
							deltaU = U[iVar] - U_old[iVar];
							Source[iVar] += deltaU*D[iZone][jZone];
						}

					} else {

						Psi[iVar] = solver_container[jZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->GetSolution(iVar);
						Source[iVar] += Psi[iVar]*D[jZone][iZone];

						if (implicit) {
							Psi_old[iVar] = solver_container[jZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->GetSolution_Old(iVar);
							deltaPsi = Psi[iVar] - Psi_old[iVar];
							Source[iVar] += deltaPsi*D[jZone][iZone];
						}
					}
				}


				/*--- Store sources for current row ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					if (!adjoint) {
						solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source[iVar]);
					} else {
						solver_container[iZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source[iVar]);
					}
				}

			}
		}
	}



	//	/*--- Loop over all grid levels ---*/
	//	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {
	//
	//		/*--- Loop over each node in the volume mesh ---*/
	//		for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {
	//
	//			for (iZone = 0; iZone < nZone; iZone++) {
	//				for (iVar = 0; iVar < nVar; iVar++) Source[iVar] = 0.0;
	//				for (jZone = 0; jZone < nZone; jZone++) {
	//
	//					/*--- Retrieve solution at this node in current zone ---*/
	//					for (iVar = 0; iVar < nVar; iVar++) {
	//						U[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution(iVar);
	//						Source[iVar] += U[iVar]*D[iZone][jZone];
	//					}
	//				}
	//				/*--- Store sources for current iZone ---*/
	//				for (iVar = 0; iVar < nVar; iVar++)
	//					solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source[iVar]);
	//			}
	//		}
	//	}

	/*--- Source term for a turbulence model ---*/
	if (config_container[ZONE_0]->GetKind_Solver() == RANS) {

		/*--- Extra variables needed if we have a turbulence model. ---*/
		unsigned short nVar_Turb = solver_container[ZONE_0][MESH_0][TURB_SOL]->GetnVar();
		double *U_Turb      = new double[nVar_Turb];
		double *Source_Turb = new double[nVar_Turb];

		/*--- Loop over only the finest mesh level (turbulence is always solved
     on the original grid only). ---*/
		for (iPoint = 0; iPoint < geometry_container[ZONE_0][MESH_0]->GetnPoint(); iPoint++) {
			//for (iZone = 0; iZone < nZone; iZone++) {
			for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
			for (jZone = 0; jZone < nZone; jZone++) {

				/*--- Retrieve solution at this node in current zone ---*/
				for (iVar = 0; iVar < nVar_Turb; iVar++) {
					U_Turb[iVar] = solver_container[jZone][MESH_0][TURB_SOL]->node[iPoint]->GetSolution(iVar);
					Source_Turb[iVar] += U_Turb[iVar]*D[iZone][jZone];
				}
			}
			/*--- Store sources for current iZone ---*/
			for (iVar = 0; iVar < nVar_Turb; iVar++)
				solver_container[iZone][MESH_0][TURB_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source_Turb[iVar]);
			//}
		}

		delete [] U_Turb;
		delete [] Source_Turb;
  }

  /*--- delete dynamic memory for D ---*/
  for (kZone = 0; kZone < nZone; kZone++) {
	  delete [] D[kZone];
  }
  delete [] D;
  delete [] U;
  delete [] U_old;
  delete [] Psi;
  delete [] Psi_old;
  delete [] Source;

  /*--- Write file with force coefficients ---*/
	ofstream TS_Flow_file;
  ofstream mean_TS_Flow_file;

  /*--- MPI Send/Recv buffers ---*/
  double *sbuf_force = NULL,  *rbuf_force = NULL;

  /*--- Other variables ---*/
  unsigned short nVar_Force = 8;
  unsigned long current_iter = config_container[ZONE_0]->GetExtIter();

  /*--- Allocate memory for send buffer ---*/
  sbuf_force = new double[nVar_Force];

  double *averages = new double[nVar_Force];
  for (iVar = 0; iVar < nVar_Force; iVar++)
	  averages[iVar] = 0;

  /*--- Allocate memory for receive buffer ---*/
  if (rank == MASTER_NODE) {
	  rbuf_force = new double[nVar_Force];

	  TS_Flow_file.precision(15);
	  TS_Flow_file.open("TS_force_coefficients.csv", ios::out);
	  TS_Flow_file <<  "\"time_instance\",\"lift_coeff\",\"drag_coeff\",\"moment_coeff_x\",\"moment_coeff_y\",\"moment_coeff_z\"" << endl;

	  mean_TS_Flow_file.precision(15);
	  if (current_iter == 0 && iZone == 1) {
		  mean_TS_Flow_file.open("history_TS_forces.plt", ios::trunc);
		  mean_TS_Flow_file << "TITLE = \"SU2 TIME-SPECTRAL SIMULATION\"" << endl;
		  mean_TS_Flow_file <<  "VARIABLES = \"Iteration\",\"CLift\",\"CDrag\",\"CMx\",\"CMy\",\"CMz\",\"CT\",\"CQ\",\"CMerit\"" << endl;
		  mean_TS_Flow_file << "ZONE T= \"Average Convergence History\"" << endl;
	  }
	  else
		  mean_TS_Flow_file.open("history_TS_forces.plt", ios::out | ios::app);
  }

  for (kZone = 0; kZone < nZone; kZone++) {

	  /*--- Flow solution coefficients (parallel) ---*/
	  sbuf_force[0] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
	  sbuf_force[1] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
	  sbuf_force[2] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMx();
	  sbuf_force[3] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMy();
	  sbuf_force[4] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMz();
	  sbuf_force[5] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CT();
	  sbuf_force[6] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CQ();
	  sbuf_force[7] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMerit();

#ifndef NO_MPI

	  /*--- for a given zone, sum the coefficients across the processors ---*/
	  MPI::COMM_WORLD.Reduce(sbuf_force, rbuf_force, nVar_Force, MPI::DOUBLE, MPI::SUM, MASTER_NODE);
	  MPI::COMM_WORLD.Barrier();
#else
	  for (iVar = 0; iVar < nVar_Force; iVar++) {
		  rbuf_force[iVar] = sbuf_force[iVar];
	  }
#endif

	  if (rank == MASTER_NODE) {
		  TS_Flow_file << kZone << ", ";
		  for (iVar = 0; iVar < nVar_Force; iVar++)
			  TS_Flow_file << rbuf_force[iVar] << ", ";
		  TS_Flow_file << endl;

		  /*--- Increment the total contributions from each zone, dividing by nZone as you go ---*/
		  for (iVar = 0; iVar < nVar_Force; iVar++) {
			  averages[iVar] += (1.0/double(nZone))*rbuf_force[iVar];
		  }
	  }
  }


  if (rank == MASTER_NODE && iZone == ZONE_0) {

	  mean_TS_Flow_file << current_iter << ", ";
	  for (iVar = 0; iVar < nVar_Force; iVar++) {
		  mean_TS_Flow_file << averages[iVar];
		  if (iVar < nVar_Force-1)
			  mean_TS_Flow_file << ", ";
	  }
	  mean_TS_Flow_file << endl;
  }

  if (rank == MASTER_NODE) {
	  TS_Flow_file.close();
	  mean_TS_Flow_file.close();
	  delete [] rbuf_force;
  }

  delete [] sbuf_force;
  delete [] averages;

}

void SetSliding_Interfaces(CGeometry ***geometry_container, CSolver ****solver_container,
		CConfig **config_container, unsigned short nZone) {

#ifndef NO_MPI
	cout << "!!! Error: Sliding mesh interfaces not yet supported in parallel. !!!" << endl;
	cout << "Press any key to exit..." << endl;
	cin.get();
	MPI::COMM_WORLD.Abort(1);
	MPI::Finalize();
#endif

	unsigned short nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();
	unsigned short iMarker, iZone, donorZone;
	unsigned long iVertex, iPoint;
	double *Coord_i, N_0, N_1, N_2;
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	double a[4], b[4], c[4], Area;
	double eps = 1e-10;

	cout << endl << "Searching and interpolating across sliding interfaces." << endl;

	/*--- Check all zones for a sliding interface. ---*/
	for (iZone = 0; iZone < nZone; iZone++) {

		/*--- Check all markers for any SEND/RECV boundaries that might be sliding. ---*/
		for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
			if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

				short SendRecv = config_container[iZone]->GetMarker_All_SendRecv(iMarker);
				/*--- We are only interested in receive information  ---*/
				if (SendRecv < 0) {

					/*--- Loop over all points on this boundary. ---*/
					for (iVertex = 0; iVertex < geometry_container[iZone][MESH_0]->GetnVertex(iMarker); iVertex++) {

						/*--- If the donor is not in a different zone, do nothing (periodic). ---*/
						donorZone = geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->GetMatching_Zone();

						if (donorZone != iZone) {
							/*--- Get the index and coordinates of the point for which we need to interpolate a solution. ---*/
							iPoint  = geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->GetNode();
							Coord_i = geometry_container[iZone][MESH_0]->node[iPoint]->GetCoord();

							/*--- Loop over the donor marker in the donor zone and search
               the neighboring elements to find the element that contains the
               point to be interpolated. This simply MUST be improved before
               going to any large cases: for now we are doing a brute force
               search, and there are much more efficient alternatives. ---*/

							for (iElem = 0; iElem < geometry_container[donorZone][MESH_0]->GetnElem(); iElem++) {

								/*--- Compute the basis functions for triangular elements. ---*/
								if (nDim == 2) {

									Point_0 = geometry_container[donorZone][MESH_0]->elem[iElem]->GetNode(0);
									Point_1 = geometry_container[donorZone][MESH_0]->elem[iElem]->GetNode(1);
									Point_2 = geometry_container[donorZone][MESH_0]->elem[iElem]->GetNode(2);

									Coord_0 = geometry_container[donorZone][MESH_0]->node[Point_0]->GetCoord();
									Coord_1 = geometry_container[donorZone][MESH_0]->node[Point_1]->GetCoord();
									Coord_2 = geometry_container[donorZone][MESH_0]->node[Point_2]->GetCoord();

									for (unsigned short iDim = 0; iDim < nDim; iDim++) {
										a[iDim] = Coord_0[iDim]-Coord_2[iDim];
										b[iDim] = Coord_1[iDim]-Coord_2[iDim];
									}

									/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/
											Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

									a[0] = 0.5 * (Coord_1[0]*Coord_2[1]-Coord_2[0]*Coord_1[1]) / Area;
									a[1] = 0.5 * (Coord_2[0]*Coord_0[1]-Coord_0[0]*Coord_2[1]) / Area;
									a[2] = 0.5 * (Coord_0[0]*Coord_1[1]-Coord_1[0]*Coord_0[1]) / Area;

									b[0] = 0.5 * (Coord_1[1]-Coord_2[1]) / Area;
									b[1] = 0.5 * (Coord_2[1]-Coord_0[1]) / Area;
									b[2] = 0.5 * (Coord_0[1]-Coord_1[1]) / Area;

									c[0] = 0.5 * (Coord_2[0]-Coord_1[0]) / Area;
									c[1] = 0.5 * (Coord_0[0]-Coord_2[0]) / Area;
									c[2] = 0.5 * (Coord_1[0]-Coord_0[0]) / Area;

									/*--- Basis functions ---*/
									N_0 = a[0] + b[0]*Coord_i[0] + c[0]*Coord_i[1];
									N_1 = a[1] + b[1]*Coord_i[0] + c[1]*Coord_i[1];
									N_2 = a[2] + b[2]*Coord_i[0] + c[2]*Coord_i[1];

									/*--- If 0 <= N <= 1 for all three basis functions, then
                   iPoint is found within this element. Store the element index
                   and the value for the basis functions to be used during the
                   send receive along the sliding interfaces. ---*/
									if ((N_0 >= -eps && N_0 <= 1.0+eps) &&
											(N_1 >= -eps && N_1 <= 1.0+eps) &&
											(N_2 >= -eps && N_2 <= 1.0+eps)) {

										/*--- Store this element as the donor for the SEND/RECV ---*/
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetDonorElem(iElem);
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetBasisFunction(0, N_0);
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetBasisFunction(1, N_1);
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetBasisFunction(2, N_2);

									}

								} else {
									cout << "!!! Error: Sliding mesh interfaces not yet supported in 3-D. !!!" << endl;
									cout << "Press any key to exit..." << endl;
									cin.get();
									exit(1);
								}
							}
						}
					}
				}
			}
		}
	}
}

void SetTimeSpectral_Velocities(CGeometry ***geometry_container,
		CConfig **config_container, unsigned short nZone) {

	unsigned short iZone, jDegree, iDim, iMGlevel;
	unsigned short nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();

	double angular_interval = 2.0*PI_NUMBER/(double)(nZone);
	double *Coord;
//	double Omega[3], Ampl[3];
	unsigned long iPoint;
	// double DEG2RAD = PI_NUMBER/180.0;
	//	int rank = MASTER_NODE;

//	/*--- rotational velocity ---*/
//	Omega[0]  = config_container[ZONE_0]->GetPitching_Omega_X(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	Omega[1]  = config_container[ZONE_0]->GetPitching_Omega_Y(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	Omega[2]  = config_container[ZONE_0]->GetPitching_Omega_Z(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	double Omega_mag = sqrt(pow(Omega[0],2)+pow(Omega[1],2)+pow(Omega[2],2));
//
//	/*--- amplitude of motion ---*/
//	Ampl[0]   = config_container[ZONE_0]->GetPitching_Ampl_X(ZONE_0)*DEG2RAD;
//	Ampl[1]   = config_container[ZONE_0]->GetPitching_Ampl_Y(ZONE_0)*DEG2RAD;
//	Ampl[2]   = config_container[ZONE_0]->GetPitching_Ampl_Z(ZONE_0)*DEG2RAD;

	/*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
	double period = config_container[ZONE_0]->GetTimeSpectral_Period();
	double deltaT = period/(double)(config_container[ZONE_0]->GetnTimeInstances());

	/*--- allocate dynamic memory for angular positions (these are the abscissas) ---*/
	double *angular_positions = new double [nZone];
	for (iZone = 0; iZone < nZone; iZone++) {
		angular_positions[iZone] = iZone*angular_interval;
	}

	/*--- find the highest-degree trigonometric polynomial allowed by the Nyquist criterion---*/
	double high_degree = (nZone-1)/2.0;
	int highest_degree = (int)(high_degree);

	/*--- allocate dynamic memory for a given point's coordinates ---*/
	double **coords = new double *[nZone];
	for (iZone = 0; iZone < nZone; iZone++) {
		coords[iZone] = new double [nDim];
	}

	/*--- allocate dynamic memory for vectors of Fourier coefficients ---*/
	double *a_coeffs = new double [highest_degree+1];
	double *b_coeffs = new double [highest_degree+1];

	/*--- allocate dynamic memory for the interpolated positions and velocities ---*/
	double *fitted_coords = new double [nZone];
	double *fitted_velocities = new double [nZone];

	/*--- Loop over all grid levels ---*/
	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {

		/*--- Loop over each node in the volume mesh ---*/
		for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {

			/*--- Populate the 2D coords array with the
			coordinates of a given mesh point across
			the time instances (i.e. the Zones) ---*/
			/*--- Loop over each zone ---*/
			for (iZone = 0; iZone < nZone; iZone++) {
				/*--- get the coordinates of the given point ---*/
				Coord = geometry_container[iZone][iMGlevel]->node[iPoint]->GetCoord();
				for (iDim = 0; iDim < nDim; iDim++) {
					/*--- add them to the appropriate place in the 2D coords array ---*/
					coords[iZone][iDim] = Coord[iDim];
				}
			}

			/*--- Loop over each Dimension ---*/
			for (iDim = 0; iDim < nDim; iDim++) {

				/*--- compute the Fourier coefficients ---*/
				for (jDegree = 0; jDegree < highest_degree+1; jDegree++) {
					a_coeffs[jDegree] = 0;
					b_coeffs[jDegree] = 0;
					for (iZone = 0; iZone < nZone; iZone++) {
						a_coeffs[jDegree] = a_coeffs[jDegree] + (2.0/(double)nZone)*cos(jDegree*angular_positions[iZone])*coords[iZone][iDim];
						b_coeffs[jDegree] = b_coeffs[jDegree] + (2.0/(double)nZone)*sin(jDegree*angular_positions[iZone])*coords[iZone][iDim];
					}
				}

				/*--- find the interpolation of the coordinates and its derivative (the velocities) ---*/
				for (iZone = 0; iZone < nZone; iZone++) {
					fitted_coords[iZone] = a_coeffs[0]/2.0;
					fitted_velocities[iZone] = 0.0;
					for (jDegree = 1; jDegree < highest_degree+1; jDegree++) {
						fitted_coords[iZone] = fitted_coords[iZone] + a_coeffs[jDegree]*cos(jDegree*angular_positions[iZone]) + b_coeffs[jDegree]*sin(jDegree*angular_positions[iZone]);
						fitted_velocities[iZone] = fitted_velocities[iZone] + (angular_interval/deltaT)*jDegree*(b_coeffs[jDegree]*cos(jDegree*angular_positions[iZone]) - a_coeffs[jDegree]*sin(jDegree*angular_positions[iZone]));
					}
				}

				/*--- Store grid velocity for this point, at this given dimension, across the Zones ---*/
				for (iZone = 0; iZone < nZone; iZone++) {
					geometry_container[iZone][iMGlevel]->node[iPoint]->SetGridVel(iDim, fitted_velocities[iZone]);
				}



			}
		}
	}

	/*--- delete dynamic memory for the abscissas, coefficients, et cetera ---*/
	delete [] angular_positions;
	delete [] a_coeffs;
	delete [] b_coeffs;
	delete [] fitted_coords;
	delete [] fitted_velocities;
	for (iZone = 0; iZone < nDim; iZone++) {
		delete [] coords[iZone];
	}
	delete [] coords;

}
