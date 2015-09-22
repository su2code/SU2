/*!
 * \file iteration_structure.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
 * \version 4.0.0 "Cardinal"
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

#include "../include/iteration_structure.hpp"

void MeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                       CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                       CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
	su2double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
  
	bool spectral_method = (config_container[ZONE_0]->GetUnsteady_Simulation() == SPECTRAL_METHOD);
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	if (spectral_method) {
    nZone = config_container[ZONE_0]->GetnTimeInstances();
  }
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
#ifdef HAVE_MPI
  int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Set the initial condition ---*/
  
  for (iZone = 0; iZone < nZone; iZone++)
    solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
  
  /*--- Initial set up for unsteady problems with dynamic meshes. ---*/
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
    /*--- Dynamic mesh update ---*/
    
		if ((config_container[iZone]->GetGrid_Movement()) && (!spectral_method)) {
			SetGrid_Movement(geometry_container[iZone], surface_movement[iZone], grid_movement[iZone], FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, IntIter, ExtIter);
    }
    
    /*--- Apply a Wind Gust ---*/
    
    if (config_container[ZONE_0]->GetWind_Gust()) {
      SetWind_GustField(config_container[iZone], geometry_container[iZone], solver_container[iZone]);
    }
	}
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Set the value of the internal iteration ---*/
    
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
        (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
    
		/*--- Update global parameters ---*/
    
    if ((config_container[iZone]->GetKind_Solver() == EULER) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_EULER)) {
      config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
    }
    if ((config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)) {
      config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
    }
    if ((config_container[iZone]->GetKind_Solver() == RANS) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
      config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
    }
    
		/*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
    
		integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
    
    if ((config_container[iZone]->GetKind_Solver() == RANS) ||
        (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
      
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
    
		/*--- Compute & store time-spectral source terms across all zones ---*/
    
		if (spectral_method)
			SetSpectralMethod(geometry_container, solver_container, config_container, nZone, (iZone+1)%nZone);
    
	}
  
	/*--- Dual time stepping strategy ---*/
  
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
		for (IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);
      
      /*--- Set the value of the internal iteration ---*/
      
      config_container[ZONE_0]->SetIntIter(IntIter);
      
			/*--- All zones must be advanced and coupled with each pseudo timestep. ---*/
      
			for (iZone = 0; iZone < nZone; iZone++) {
        
				/*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/
        
        if ((config_container[iZone]->GetKind_Solver() == EULER) ||
            (config_container[iZone]->GetKind_Solver() == DISC_ADJ_EULER)) {
          config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
        }
        if ((config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) ||
            (config_container[iZone]->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)) {
          config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
        }
        if ((config_container[iZone]->GetKind_Solver() == RANS) ||
            (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
          config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
        }
        
        /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
        
				integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                    config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
        
				/*--- Pseudo-timestepping the turbulence model ---*/
        
        if ((config_container[iZone]->GetKind_Solver() == RANS) ||
            (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
          
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
        
				/*--- Call Dynamic mesh update if AEROELASTIC motion was specified ---*/
        if ((config_container[ZONE_0]->GetGrid_Movement()) && (config_container[ZONE_0]->GetAeroelastic_Simulation())) {
					SetGrid_Movement(geometry_container[iZone], surface_movement[iZone], grid_movement[iZone], FFDBox[iZone],
                           solver_container[iZone], config_container[iZone], iZone, IntIter, ExtIter);
          /*--- Apply a Wind Gust ---*/
          if (config_container[ZONE_0]->GetWind_Gust()) {
            if (IntIter % config_container[iZone]->GetAeroelasticIter() ==0)
              SetWind_GustField(config_container[iZone], geometry_container[iZone], solver_container[iZone]);
          }
        }
      }
      
      if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) break;
      
    }
    
		for (iZone = 0; iZone < nZone; iZone++) {
      
			/*--- Update dual time solver on all mesh levels ---*/
      
			for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
				integration_container[iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][FLOW_SOL], config_container[iZone], iMesh);
				integration_container[iZone][FLOW_SOL]->SetConvergence(false);
			}
      
			/*--- Update dual time solver for the turbulence model ---*/
      
      if ((config_container[iZone]->GetKind_Solver() == RANS) ||
          (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
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

void AdjMeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                          CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                          CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox) {
  
	su2double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
  
	bool spectral_method = (config_container[ZONE_0]->GetUnsteady_Simulation() == SPECTRAL_METHOD);
  bool grid_movement = config_container[ZONE_0]->GetGrid_Movement();
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	if (spectral_method) nZone = config_container[ZONE_0]->GetnTimeInstances();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- For the unsteady adjoint, load a new direct solution from a restart file. ---*/
  
	for (iZone = 0; iZone < nZone; iZone++) {
		if (((grid_movement && ExtIter == 0) || config_container[ZONE_0]->GetUnsteady_Simulation()) && !spectral_method) {
      int Direct_Iter = SU2_TYPE::Int(config_container[iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 1;
      if (rank == MASTER_NODE && iZone == ZONE_0 && config_container[iZone]->GetUnsteady_Simulation())
        cout << endl << " Loading flow solution from direct iteration " << Direct_Iter << "." << endl;
      solver_container[iZone][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[iZone], solver_container[iZone], config_container[iZone], Direct_Iter);
    }
	}
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Continuous adjoint Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations ---*/
    
		if ((ExtIter == 0) || config_container[iZone]->GetUnsteady_Simulation()) {
      
			if (config_container[iZone]->GetKind_Solver() == ADJ_EULER)         config_container[iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_FLOW_SYS, ExtIter);
			if (config_container[iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) config_container[iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
			if (config_container[iZone]->GetKind_Solver() == ADJ_RANS)          config_container[iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_FLOW_SYS, ExtIter);
      
      /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
      
      if (rank == MASTER_NODE && iZone == ZONE_0)
				cout << "Begin direct solver to store flow data (single iteration)." << endl;
      
      if (rank == MASTER_NODE && iZone == ZONE_0)
        cout << "Compute residuals to check the convergence of the direct problem." << endl;

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
      
      /*--- Output the residual (visualization purpouses to identify if
       the direct solution is converged)---*/
      if (rank == MASTER_NODE && iZone == ZONE_0)
        cout << "log10[Maximum residual]: " << log10(solver_container[iZone][MESH_0][FLOW_SOL]->GetRes_Max(0))
        <<", located at point "<< solver_container[iZone][MESH_0][FLOW_SOL]->GetPoint_Max(0) << "." << endl;
      
			/*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
			 note that in the direct Euler problem we are not computing the gradients of the primitive variables ---*/
      
			if (config_container[iZone]->GetKind_Gradient_Method() == GREEN_GAUSS)
				solver_container[iZone][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_GG(geometry_container[iZone][MESH_0], config_container[iZone]);
			if (config_container[iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
				solver_container[iZone][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_LS(geometry_container[iZone][MESH_0], config_container[iZone]);
      
			/*--- Set contribution from cost function for boundary conditions ---*/
      
      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
        
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
      
      if (rank == MASTER_NODE && iZone == ZONE_0)
        cout << "End direct solver, begin adjoint problem." << endl;

		}
    

		/*--- Set the value of the internal iteration ---*/
    
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
			IntIter = 0;
		}
    
		if (config_container[iZone]->GetKind_Solver() == ADJ_EULER)         config_container[iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_ADJFLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) config_container[iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == ADJ_RANS)          config_container[iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJFLOW_SYS, ExtIter);
    
		/*--- Iteration of the flow adjoint problem ---*/
    
		integration_container[iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                   config_container, RUNTIME_ADJFLOW_SYS, IntIter, iZone);
    
		/*--- Iteration of the turbulence model adjoint ---*/
    
		if ((config_container[iZone]->GetKind_Solver() == ADJ_RANS) && (!config_container[iZone]->GetFrozen_Visc())) {
      
			/*--- Turbulent model solution ---*/
      
			config_container[iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJTURB_SYS, ExtIter);
			integration_container[iZone][ADJTURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                      config_container, RUNTIME_ADJTURB_SYS, IntIter, iZone);
      
		}
    
	}
  
	/*--- Dual time stepping strategy ---*/
  
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
			(config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
		for (IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);
      
      /*--- Set the value of the internal iteration ---*/
      
      config_container[ZONE_0]->SetIntIter(IntIter);
      
			/*--- All zones must be advanced and coupled with each pseudo timestep ---*/
      
			for (iZone = 0; iZone < nZone; iZone++)
				integration_container[iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_ADJFLOW_SYS, IntIter, iZone);
      
      /*--- Check to see if the convergence criteria has been met ---*/
      
      if (integration_container[ZONE_0][ADJFLOW_SOL]->GetConvergence()) break;
		}
    
		for (iZone = 0; iZone < nZone; iZone++) {
      
			/*--- Update dual time solver ---*/
      
			for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
				integration_container[iZone][ADJFLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][ADJFLOW_SOL], config_container[iZone], iMesh);
				integration_container[iZone][ADJFLOW_SOL]->SetConvergence(false);
			}
      
			Physical_dt = config_container[iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime()) integration_container[iZone][ADJFLOW_SOL]->SetConvergence(true);
      
		}
    
	}
  
}


void TNE2Iteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                   CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                   CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
  unsigned short iZone; // Index for zone of the mesh
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();  
  
#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
    
		/*--- Set the initial condition ---*/
		solver_container[iZone][MESH_0][TNE2_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
    
		/*--- Update global parameters ---*/
		if (config_container[iZone]->GetKind_Solver() == TNE2_EULER)
      config_container[iZone]->SetGlobalParam(TNE2_EULER, RUNTIME_TNE2_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES)
      config_container[iZone]->SetGlobalParam(TNE2_NAVIER_STOKES, RUNTIME_TNE2_SYS, ExtIter);
    
		/*--- Solve the inviscid or viscous two-temperature flow equations (one iteration) ---*/
		integration_container[iZone][TNE2_SOL]->MultiGrid_Iteration(geometry_container, solver_container,
                                                                numerics_container, config_container,
                                                                RUNTIME_TNE2_SYS, IntIter, iZone);
	}
  
}

void DiscAdjMeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                          CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                          CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox) {

  unsigned short iZone, iMesh, ExtIter = config_container[ZONE_0]->GetExtIter();
  unsigned short nZone = config_container[ZONE_0]->GetnZone();

  bool turbulent = false, flow = false;

  switch(config_container[ZONE_0]->GetKind_Solver()){
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: flow = true; break;
    case DISC_ADJ_RANS: flow = true; turbulent = true; break;
  }

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (ExtIter == 0){

    /*--- Start the recording of all operations ---*/

    AD::StartRecording();

    /*--- Register all necessary variables on the tape ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Register the node coordinates as input of the iteration ---*/

      geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);

      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

        /*--- Register the conservative variables as input of the iteration ---*/
        if (flow){
          solver_container[iZone][iMesh][ADJFLOW_SOL]->RegisterInput(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
        }
        if (turbulent){
          solver_container[iZone][iMesh][ADJTURB_SOL]->RegisterInput(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
        }
      }
    }

    /*--- Update the mesh structure to get the influence of node coordinates ---*/

    for (iZone = 0; iZone < nZone; iZone++){
      geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
    }

    if (rank == MASTER_NODE){
      cout << "Begin direct solver to store computational graph (single iteration)." << endl;
    }


    /* --- Preprocessing of flow solver and postprocessing of turbulent solver.
     * We need this to get the dependency of the flow variables on the eddy viscosity. ---*/
    if (turbulent){
      for (iZone = 0; iZone < nZone; iZone++){
        for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
          solver_container[iZone][iMesh][FLOW_SOL]->Preprocessing(geometry_container[iZone][iMesh], solver_container[iZone][iMesh], config_container[iZone], iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
          solver_container[iZone][iMesh][TURB_SOL]->Postprocessing(geometry_container[iZone][iMesh],solver_container[iZone][iMesh], config_container[iZone], iMesh);
        }
      }
    }

    /*--- One iteration of the flow solver ---*/

    if (flow){
      MeanFlowIteration(output, integration_container, geometry_container,
                      solver_container, numerics_container, config_container,
                      surface_movement, volume_grid_movement, FFDBox);
    }

    if (rank == MASTER_NODE){
      if(flow){
        cout << "log10[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))
             <<", Drag: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()
            <<", Lift: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift() << "." << endl;
      }
      if (turbulent){
        cout << "log10[RMS k]: " << log10(solver_container[ZONE_0][MESH_0][TURB_SOL]->GetRes_RMS(0)) << std::endl;
      }
    }

    for (iZone = 0; iZone < nZone; iZone++){
      /*--- Register objective function as output of the iteration ---*/
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
      }
      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

        /*--- Register conservative variables as output of the iteration ---*/
        if (flow){
          solver_container[iZone][iMesh][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
        }
        if (turbulent){
          solver_container[iZone][iMesh][ADJTURB_SOL]->RegisterOutput(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
        }
      }
    }

    /*--- Stop the recording ---*/

    AD::StopRecording();

  }
  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Initialize the adjoint of the objective function (typically with 1.0) ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone]);
    }
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

      /*--- Initialize the adjoints the conservative variables ---*/
      if (flow){
        solver_container[iZone][iMesh][ADJFLOW_SOL]->SetAdjointOutput(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
      }
      if (turbulent){
        solver_container[iZone][iMesh][ADJTURB_SOL]->SetAdjointOutput(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
      }
    }
  }

  /*--- Run the adjoint computation ---*/

  AD::ComputeAdjoint();

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the sensitivities by extracting the adjoints of the node coordinates ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]);
    }
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

      /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/
      if (flow){
        solver_container[iZone][iMesh][ADJFLOW_SOL]->SetAdjointInput(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
      }
      if (turbulent){
        solver_container[iZone][iMesh][ADJTURB_SOL]->SetAdjointInput(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
      }
    }
  }

  /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/

  AD::ClearAdjoints();

  /*--- Set the convergence criteria ---*/

  if (config_container[ZONE_0]->GetConvCriteria() == RESIDUAL){
    if (flow){
      integration_container[ZONE_0][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[ZONE_0][MESH_0],config_container[ZONE_0],
                                                                         ExtIter,log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0)), MESH_0);
    }
  }
  else if (config_container[ZONE_0]->GetConvCriteria() == CAUCHY){
    if (flow){
      integration_container[ZONE_0][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[ZONE_0][MESH_0],config_container[ZONE_0],
                                                                         ExtIter,solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo(), MESH_0);
    }
  }
}

void AdjTNE2Iteration(COutput *output, CIntegration ***integration_container,
                      CGeometry ***geometry_container,
                      CSolver ****solver_container,
                      CNumerics *****numerics_container,
                      CConfig **config_container,
                      CSurfaceMovement **surface_movement,
                      CVolumetricMovement **grid_movement,
                      CFreeFormDefBox*** FFDBox) {
  
	unsigned short iMesh, iZone, nZone;
  unsigned long IntIter, ExtIter;
  int rank;

  /*--- Initialize parameters ---*/
  nZone   = geometry_container[ZONE_0][MESH_0]->GetnZone();
  IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  ExtIter = config_container[ZONE_0]->GetExtIter();
  rank    = MASTER_NODE;

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Continuous adjoint two-temperature equations  ---*/
		if (ExtIter == 0) {
      
      /*--- Set the value of the internal iteration ---*/
      IntIter = ExtIter;
      
      /*--- Set the initial condition ---*/
      solver_container[iZone][MESH_0][TNE2_SOL]->SetInitialCondition(geometry_container[iZone],
                                                                     solver_container[iZone],
                                                                     config_container[iZone], ExtIter);
      
      /*--- Update global parameters ---*/
			if (config_container[iZone]->GetKind_Solver() == ADJ_TNE2_EULER)
        config_container[iZone]->SetGlobalParam(TNE2_EULER,
                                                RUNTIME_TNE2_SYS, ExtIter);
			else if (config_container[iZone]->GetKind_Solver() == ADJ_TNE2_NAVIER_STOKES)
        config_container[iZone]->SetGlobalParam(ADJ_TNE2_NAVIER_STOKES,
                                                RUNTIME_TNE2_SYS, ExtIter);
      
      /*--- Perform one iteration of the gov. eqns. to store data ---*/
      if (rank == MASTER_NODE && iZone == ZONE_0)
				cout << "Single iteration of the direct solver to store flow data...";
			integration_container[iZone][TNE2_SOL]->MultiGrid_Iteration(geometry_container,
                                                                  solver_container,
                                                                  numerics_container,
                                                                  config_container,
                                                                  RUNTIME_TNE2_SYS,
                                                                  IntIter, iZone);
      if (rank == MASTER_NODE && iZone == ZONE_0)
        cout << " Done." << endl;
      
			/*--- Compute gradients of the flow variables, this is necessary for
       sensitivity computation, note that in the direct Euler problem we
       are not computing the gradients of the primitive variables ---*/
			if (config_container[iZone]->GetKind_Gradient_Method() == GREEN_GAUSS)
				solver_container[iZone][MESH_0][TNE2_SOL]->SetPrimitive_Gradient_GG(geometry_container[iZone][MESH_0],
                                                                          config_container[iZone]);
			if (config_container[iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
				solver_container[iZone][MESH_0][TNE2_SOL]->SetPrimitive_Gradient_LS(geometry_container[iZone][MESH_0],
                                                                          config_container[iZone]);
      
			/*--- Set contribution from cost function for boundary conditions ---*/
      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
        
        /*--- Set the value of the non-dimensional coefficients in the coarse
         levels, using the fine level solution ---*/
        solver_container[iZone][iMesh][TNE2_SOL]->SetTotal_CDrag(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CDrag());
        solver_container[iZone][iMesh][TNE2_SOL]->SetTotal_CLift(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CLift());
        solver_container[iZone][iMesh][TNE2_SOL]->SetTotal_CT(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CT());
        solver_container[iZone][iMesh][TNE2_SOL]->SetTotal_CQ(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CQ());
        
        /*--- Compute the adjoint boundary condition on Euler walls ---*/
        solver_container[iZone][iMesh][ADJTNE2_SOL]->SetForceProj_Vector(geometry_container[iZone][iMesh],
                                                                         solver_container[iZone][iMesh],
                                                                         config_container[iZone]);
      }
		}
    
		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
    
		if (config_container[iZone]->GetKind_Solver() == ADJ_TNE2_EULER)
      config_container[iZone]->SetGlobalParam(ADJ_TNE2_EULER,
                                              RUNTIME_ADJTNE2_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == ADJ_TNE2_NAVIER_STOKES)
      config_container[iZone]->SetGlobalParam(ADJ_TNE2_NAVIER_STOKES,
                                              RUNTIME_ADJTNE2_SYS, ExtIter);
    
		/*--- Iteration of the flow adjoint problem ---*/
		integration_container[iZone][ADJTNE2_SOL]->MultiGrid_Iteration(geometry_container,
                                                                   solver_container,
                                                                   numerics_container,
                                                                   config_container,
                                                                   RUNTIME_ADJTNE2_SYS,
                                                                   IntIter, iZone);
	}
}

void WaveIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                   CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                   CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
	su2double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
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
      
			for (IntIter = 1; IntIter < config_container[iZone]->GetUnst_nIntIter(); IntIter++) {
        output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, iZone);
        config_container[iZone]->SetIntIter(IntIter);
				integration_container[iZone][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                     config_container, RUNTIME_WAVE_SYS, IntIter, iZone);
				if (integration_container[iZone][WAVE_SOL]->GetConvergence()) break;
			}
      
			/*--- Update dual time solver ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
				integration_container[iZone][WAVE_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][WAVE_SOL], config_container[iZone], iMesh);
				integration_container[iZone][WAVE_SOL]->SetConvergence(false);
			}
      
			Physical_dt = config_container[iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime()) integration_container[iZone][WAVE_SOL]->SetConvergence(true);
		}
    
	}
  
}

void HeatIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                   CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                   CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
	su2double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
    
		/*--- Wave equations ---*/
		config_container[iZone]->SetGlobalParam(HEAT_EQUATION, RUNTIME_HEAT_SYS, ExtIter);
		integration_container[iZone][HEAT_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                 config_container, RUNTIME_HEAT_SYS, IntIter, iZone);
    
		/*--- Dual time stepping strategy ---*/
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
			for (IntIter = 1; IntIter < config_container[iZone]->GetUnst_nIntIter(); IntIter++) {
        output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, iZone);
        config_container[iZone]->SetIntIter(IntIter);
				integration_container[iZone][HEAT_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                     config_container, RUNTIME_HEAT_SYS, IntIter, iZone);
				if (integration_container[iZone][HEAT_SOL]->GetConvergence()) break;
			}
      
			/*--- Update dual time solver ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
				integration_container[iZone][HEAT_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][HEAT_SOL], config_container[iZone], iMesh);
				integration_container[iZone][HEAT_SOL]->SetConvergence(false);
			}
      
			Physical_dt = config_container[iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime()) integration_container[iZone][HEAT_SOL]->SetConvergence(true);
		}
    
	}
  
}

void PoissonIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                      CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                      CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
  
	unsigned short iZone;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
    
		/*--- Wave equations ---*/
		config_container[iZone]->SetGlobalParam(POISSON_EQUATION, RUNTIME_POISSON_SYS, ExtIter);
		integration_container[iZone][POISSON_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                    config_container, RUNTIME_POISSON_SYS, IntIter, iZone);
    
	}
  
}

void FEAIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                  CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                  CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {
	su2double Physical_dt, Physical_t;
	unsigned short iMesh, iZone;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  	unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
    if (config_container[iZone]->GetGrid_Movement())
      SetGrid_Movement(geometry_container[iZone], surface_movement[iZone],
                       grid_movement[iZone], FFDBox[iZone], solver_container[iZone], config_container[iZone], iZone, IntIter, ExtIter);
    
		/*--- Set the value of the internal iteration ---*/
    
		IntIter = ExtIter;
    
		/*--- Set the initial condition at the first iteration ---*/
    
		solver_container[iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
    
		/*--- FEA equations ---*/
    
		config_container[iZone]->SetGlobalParam(LINEAR_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

		/*--- Run the iteration ---*/

		integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
                                                                config_container, RUNTIME_FEA_SYS, IntIter, iZone);

		/*----------------- Update structural solver ----------------------*/

		bool dynamic = (config_container[iZone]->GetDynamic_Analysis() == DYNAMIC);

		if (dynamic){
			integration_container[iZone][FEA_SOL]->SetStructural_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][FEA_SOL], config_container[iZone], MESH_0);
		}

	}
  
}

void FluidStructureIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                             CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                             CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
                             unsigned long iFluidIt, unsigned long nFluidIt) {

	su2double Physical_dt, Physical_t;
	unsigned short iMesh;
	unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
	unsigned long IntIter_Struct = 0; config_container[ZONE_1]->SetIntIter(IntIter_Struct);
	unsigned long iFSIIter = 0;
	unsigned long nFSIIter = config_container[ZONE_0]->GetnIterFSI();
	unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

	unsigned short SolContainer_Position_fea = config_container[ZONE_1]->GetContainerPosition(RUNTIME_FEA_SYS);

	/*------------- Structural predictor for displacements ------------*/

	/*--- Predict structural displacement --*/
	solver_container[ZONE_1][MESH_0][FEA_SOL]->PredictStruct_Displacement(geometry_container[ZONE_1], config_container[ZONE_1],
																		  solver_container[ZONE_1]);


	while (iFSIIter<nFSIIter){

		/*------------------------ Mesh movement --------------------------*/

		/*--- Update the the flow geometry (ZONE 0) --*/

		solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetFlow_Displacement(geometry_container[ZONE_0], grid_movement[ZONE_0],
																	   config_container[ZONE_0], config_container[ZONE_1],
																	   geometry_container[ZONE_1], solver_container[ZONE_1]);

		/*---------------------- Fluid iteration --------------------------*/

		/*--- Set the initial condition ---*/

		solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[ZONE_0], solver_container[ZONE_0], config_container[ZONE_0], ExtIter);

		/*--- Apply a Wind Gust ---*/

		if (config_container[ZONE_0]->GetWind_Gust()){
			  SetWind_GustField(config_container[ZONE_0],geometry_container[ZONE_0],solver_container[ZONE_0]);
		}

		/*--- Set the value of the internal iteration ---*/

		IntIter = ExtIter;

		if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
		   (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;

		/*--- Update global parameters ---*/

		if (config_container[ZONE_0]->GetKind_Solver() == EULER){
			config_container[ZONE_0]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
		}
		if (config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES){
			config_container[ZONE_0]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
		}
		if (config_container[ZONE_0]->GetKind_Solver() == RANS){
			config_container[ZONE_0]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
		 }

		/*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/

		integration_container[ZONE_0][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
																	config_container, RUNTIME_FLOW_SYS, IntIter, ZONE_0);

		if (config_container[ZONE_0]->GetKind_Solver() == RANS) {

			/*--- Solve the turbulence model ---*/

			config_container[ZONE_0]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
			integration_container[ZONE_0][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																		   config_container, RUNTIME_TURB_SYS, IntIter, ZONE_0);

		/*--- Solve transition model ---*/

			if (config_container[ZONE_0]->GetKind_Trans_Model() == LM) {
				config_container[ZONE_0]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
				integration_container[ZONE_0][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																			  config_container, RUNTIME_TRANS_SYS, IntIter, ZONE_0);
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

			/*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/

			if (config_container[ZONE_0]->GetKind_Solver() == EULER)
				config_container[ZONE_0]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
			if (config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES)
				config_container[ZONE_0]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
			if (config_container[ZONE_0]->GetKind_Solver() == RANS)
				config_container[ZONE_0]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);

			  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/

			integration_container[ZONE_0][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
																			config_container, RUNTIME_FLOW_SYS, IntIter, ZONE_0);

			/*--- Pseudo-timestepping the turbulence model ---*/

			if (config_container[ZONE_0]->GetKind_Solver() == RANS) {

				/*--- Solve the turbulence model ---*/

				config_container[ZONE_0]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
				integration_container[ZONE_0][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																			  config_container, RUNTIME_TURB_SYS, IntIter, ZONE_0);

				/*--- Solve transition model ---*/

					if (config_container[ZONE_0]->GetKind_Trans_Model() == LM) {
						config_container[ZONE_0]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
						integration_container[ZONE_0][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
																				  config_container, RUNTIME_TRANS_SYS, IntIter, ZONE_0);
					}
			}

		    if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) break;

		}

		}

		/*-------------------- Structural iteration -----------------------*/

		/*--- Set the initial condition at the first iteration ---*/

		// solver_container[ZONE_1][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[ZONE_1], solver_container[ZONE_1], config_container[ZONE_1], ExtIter);

		/*--- Set the value of the internal iteration ---*/

		IntIter_Struct = ExtIter;

		/*--- FEA equations ---*/

		config_container[ZONE_1]->SetGlobalParam(LINEAR_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

		/*--- Update loads for the FEA model ---*/

		solver_container[ZONE_1][MESH_0][FEA_SOL]->SetFEA_Load(solver_container[ZONE_0], geometry_container[ZONE_1], geometry_container[ZONE_0],
															   config_container[ZONE_1], config_container[ZONE_0], numerics_container[ZONE_1][MESH_0][SolContainer_Position_fea][VISC_TERM]);

		/*--- Run the iteration ---*/

		integration_container[ZONE_1][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
																	 config_container, RUNTIME_FEA_SYS, IntIter_Struct, ZONE_1);

		/*-------------------- Aitken's relaxation ------------------------*/

		solver_container[ZONE_1][MESH_0][FEA_SOL]->ComputeAitken_Coefficient(geometry_container[ZONE_1], config_container[ZONE_1],
																		solver_container[ZONE_1], iFSIIter);


		solver_container[ZONE_1][MESH_0][FEA_SOL]->SetAitken_Relaxation(geometry_container[ZONE_1], config_container[ZONE_1],
																		solver_container[ZONE_1]);

		/*-------------------- Check convergence --------------------------*/

		integration_container[ZONE_1][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[ZONE_1][MESH_0], config_container[ZONE_1],
																solver_container[ZONE_1][MESH_0][FEA_SOL], iFSIIter);

		if (integration_container[ZONE_1][FEA_SOL]->GetConvergence_FSI()) break;

		/*--------------------- Update iFSIIter ---------------------------*/

		iFSIIter++;

	}

	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
	   (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

	  	/*-------------------- Update fluid solver ------------------------*/

		/*--- Update dual time solver on all mesh levels ---*/

		for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
			integration_container[ZONE_0][FLOW_SOL]->SetDualTime_Solver(geometry_container[ZONE_0][iMesh], solver_container[ZONE_0][iMesh][FLOW_SOL], config_container[ZONE_0], iMesh);
			integration_container[ZONE_0][FLOW_SOL]->SetConvergence(false);
		}

		/*--- Update dual time solver for the turbulence model ---*/

		if (config_container[ZONE_0]->GetKind_Solver() == RANS) {
			integration_container[ZONE_0][TURB_SOL]->SetDualTime_Solver(geometry_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_0][TURB_SOL], config_container[ZONE_0], MESH_0);
			integration_container[ZONE_0][TURB_SOL]->SetConvergence(false);
		}

  /*--- Update dual time solver for the transition model ---*/

		if (config_container[ZONE_0]->GetKind_Trans_Model() == LM) {
			integration_container[ZONE_0][TRANS_SOL]->SetDualTime_Solver(geometry_container[ZONE_0][MESH_0], solver_container[ZONE_0][MESH_0][TRANS_SOL], config_container[ZONE_0], MESH_0);
			integration_container[ZONE_0][TRANS_SOL]->SetConvergence(false);
		}

  /*--- Verify convergence criteria (based on total time) ---*/

		Physical_dt = config_container[ZONE_0]->GetDelta_UnstTime();
		Physical_t  = (ExtIter+1)*Physical_dt;
		if (Physical_t >=  config_container[ZONE_0]->GetTotal_UnstTime())
			integration_container[ZONE_0][FLOW_SOL]->SetConvergence(true);

	}

	/*----------------- Update structural solver ----------------------*/

	integration_container[ZONE_1][FEA_SOL]->SetStructural_Solver(geometry_container[ZONE_1][MESH_0], solver_container[ZONE_1][MESH_0][FEA_SOL], config_container[ZONE_1], MESH_0);

	/*-----------------------------------------------------------------*/
	/*--------------- Update convergence parameter --------------------*/
	/*-----------------------------------------------------------------*/

	integration_container[ZONE_1][FEA_SOL]->SetConvergence_FSI(false);

  
}

void SetWind_GustField(CConfig *config_container, CGeometry **geometry_container, CSolver ***solver_container) {
  // The gust is imposed on the flow field via the grid velocities. This method called the Field Velocity Method is described in the
  // NASA TM–2012-217771 - Development, Verification and Use of Gust Modeling in the NASA Computational Fluid Dynamics Code FUN3D
  // the desired gust is prescribed as the negative of the grid velocity.
  
  // If a source term is included to account for the gust field, the method is described by Jones et al. as the Split Velocity Method in
  // Simulation of Airfoil Gust Responses Using Prescribed Velocities.
  // In this routine the gust derivatives needed for the source term are calculated when applicable.
  // If the gust derivatives are zero the source term is also zero.
  // The source term itself is implemented in the class CSourceWindGust
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl << "Running simulation with a Wind Gust." << endl;
  unsigned short iDim, nDim = geometry_container[MESH_0]->GetnDim(); //We assume nDim = 2
  if (nDim != 2) {
    if (rank == MASTER_NODE) {
      cout << endl << "WARNING - Wind Gust capability is only verified for 2 dimensional simulations." << endl;
    }
  }
  
  /*--- Gust Parameters from config ---*/
  unsigned short Gust_Type = config_container->GetGust_Type();
  su2double xbegin = config_container->GetGust_Begin_Loc();    // Location at which the gust begins.
  su2double L = config_container->GetGust_WaveLength();        // Gust size
  su2double tbegin = config_container->GetGust_Begin_Time();   // Physical time at which the gust begins.
  su2double gust_amp = config_container->GetGust_Ampl();       // Gust amplitude
  su2double n = config_container->GetGust_Periods();           // Number of gust periods
  unsigned short GustDir = config_container->GetGust_Dir(); // Gust direction

  /*--- Variables needed to compute the gust ---*/
  unsigned short Kind_Grid_Movement = config_container->GetKind_GridMovement(ZONE_0);
  unsigned long iPoint;
  unsigned short iMGlevel, nMGlevel = config_container->GetnMGLevels();

  su2double x, y, x_gust, dgust_dx, dgust_dy, dgust_dt;
  su2double *Gust, *GridVel;
  su2double NewGridVel[2] = {0.0,0.0};
  su2double GustDer[3] = {0.0,0.0,0.0};

  su2double Physical_dt = config_container->GetDelta_UnstTime();
  unsigned long ExtIter = config_container->GetExtIter();
  su2double Physical_t = ExtIter*Physical_dt;
  
  su2double Uinf = solver_container[MESH_0][FLOW_SOL]->GetVelocity_Inf(0); // Assumption gust moves at infinity velocity
  
  Gust = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Gust[iDim] = 0.0;
  }

  // Vortex variables
  unsigned long nVortex = 0;
  vector<su2double> x0, y0, vort_strenth, r_core; //vortex is positive in clockwise direction.
  if (Gust_Type == VORTEX) {
    InitializeVortexDistribution(nVortex, x0, y0, vort_strenth, r_core);
  }
  
  /*--- Check to make sure gust lenght is not zero or negative (vortex gust doesn't use this). ---*/
  if (L <= 0.0 && Gust_Type != VORTEX) {
    if (rank == MASTER_NODE) cout << "ERROR: The gust length needs to be positive" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- Loop over all multigrid levels ---*/
  
  for (iMGlevel = 0; iMGlevel <= nMGlevel; iMGlevel++) {
    
    /*--- Loop over each node in the volume mesh ---*/
    
    for (iPoint = 0; iPoint < geometry_container[iMGlevel]->GetnPoint(); iPoint++) {
      
      /*--- Reset the Grid Velocity to zero if there is no grid movement ---*/
      if (Kind_Grid_Movement == GUST) {
        for (iDim = 0; iDim < nDim; iDim++)
          geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, 0.0);
      }
      
      /*--- initialize the gust and derivatives to zero everywhere ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {Gust[iDim]=0.0;}
      dgust_dx = 0.0; dgust_dy = 0.0; dgust_dt = 0.0;
      
      /*--- Begin applying the gust ---*/
      
      if (Physical_t >= tbegin) {
        
        x = geometry_container[iMGlevel]->node[iPoint]->GetCoord()[0]; // x-location of the node.
        y = geometry_container[iMGlevel]->node[iPoint]->GetCoord()[1]; // y-location of the node.
        
        // Gust coordinate
        x_gust = (x - xbegin - Uinf*(Physical_t-tbegin))/L;
        
        /*--- Calculate the specified gust ---*/
        switch (Gust_Type) {

          case TOP_HAT:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp;
              // Still need to put the gust derivatives. Think about this.
            }
            break;

          case SINE:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp*(sin(2*PI_NUMBER*x_gust));

              // Gust derivatives
              //dgust_dx = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))/L;
              //dgust_dy = 0;
              //dgust_dt = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))*(-Uinf)/L;
            }
            break;

          case ONE_M_COSINE:
             // Check if we are in the region where the gust is active
             if (x_gust > 0 && x_gust < n) {
               Gust[GustDir] = gust_amp*(1-cos(2*PI_NUMBER*x_gust));

               // Gust derivatives
               //dgust_dx = gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))/L;
               //dgust_dy = 0;
               //dgust_dt = gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))*(-Uinf)/L;
             }
             break;

          case EOG:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = -0.37*gust_amp*sin(3*PI_NUMBER*x_gust)*(1-cos(2*PI_NUMBER*x_gust));
            }
            break;
            
          case VORTEX:

             /*--- Use vortex distribution ---*/
             // Algebraic vortex equation.
             for (unsigned long i=0; i<nVortex; i++) {
               su2double r2 = pow(x-(x0[i]+Uinf*(Physical_t-tbegin)), 2) + pow(y-y0[i], 2);
               su2double r = sqrt(r2);
               su2double v_theta = vort_strenth[i]/(2*PI_NUMBER) * r/(r2+pow(r_core[i],2));
               Gust[0] = Gust[0] + v_theta*(y-y0[i])/r;
               Gust[1] = Gust[1] - v_theta*(x-(x0[i]+Uinf*(Physical_t-tbegin)))/r;
             }
             break;

           case NONE: default:

             /*--- There is no wind gust specified. ---*/
             if (rank == MASTER_NODE) {
               cout << "No wind gust specified." << endl;
             }
             break;

        }
      }
      
      /*--- Set the Wind Gust, Wind Gust Derivatives and the Grid Velocities ---*/
      
      GustDer[0] = dgust_dx;
      GustDer[1] = dgust_dy;
      GustDer[2] = dgust_dt;
      
      solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGust(Gust);
      solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGustDer(GustDer);
      
      GridVel = geometry_container[iMGlevel]->node[iPoint]->GetGridVel();
      
      /*--- Store new grid velocity ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {
        NewGridVel[iDim] = GridVel[iDim] - Gust[iDim];
        geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, NewGridVel[iDim]);
      }
      
    }
  }
  
  delete [] Gust;
  
}

void InitializeVortexDistribution(unsigned long &nVortex, vector<su2double>& x0, vector<su2double>& y0, vector<su2double>& vort_strength, vector<su2double>& r_core) {
  /*--- Read in Vortex Distribution ---*/
  std::string line;
  std::ifstream file;
  su2double x_temp, y_temp, vort_strength_temp, r_core_temp;
  file.open("vortex_distribution.txt");
  /*--- In case there is no vortex file ---*/
  if (file.fail()) {
    cout << "There is no vortex data file!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(EXIT_FAILURE);
  }
  
  // Ignore line containing the header
  getline(file, line);
  // Read in the information of the vortices (xloc, yloc, lambda(strength), eta(size, gradient))
  while (file.good())
  {
    getline(file, line);
    std::stringstream ss(line);
    if (line.size() != 0) { //ignore blank lines if they exist.
      ss >> x_temp;
      ss >> y_temp;
      ss >> vort_strength_temp;
      ss >> r_core_temp;
      x0.push_back(x_temp);
      y0.push_back(y_temp);
      vort_strength.push_back(vort_strength_temp);
      r_core.push_back(r_core_temp);
    }
  }
  file.close();
  // number of vortices
  nVortex = x0.size();

}

void SetGrid_Movement(CGeometry **geometry_container, CSurfaceMovement *surface_movement,
                      CVolumetricMovement *grid_movement, CFreeFormDefBox **FFDBox,
                      CSolver ***solver_container, CConfig *config_container, unsigned short iZone, unsigned long IntIter, unsigned long ExtIter)   {
  
  unsigned short iDim, iMGlevel, nMGlevels = config_container->GetnMGLevels();
	unsigned short Kind_Grid_Movement = config_container->GetKind_GridMovement(iZone);
  unsigned long iPoint;
  bool adjoint = config_container->GetAdjoint();
	bool spectral_method = (config_container->GetUnsteady_Simulation() == SPECTRAL_METHOD);
  
	/*--- For a time-spectral case, set "iteration number" to the zone number,
   so that the meshes are positioned correctly for each instance. ---*/
	if (spectral_method) {
		ExtIter = iZone;
		Kind_Grid_Movement = config_container->GetKind_GridMovement(ZONE_0);
	}
  
	int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
	/*--- Perform mesh movement depending on specified type ---*/
	switch (Kind_Grid_Movement) {
      
    case MOVING_WALL:
      
      /*--- Fixed wall velocities: set the grid velocities only one time
       before the first iteration flow solver. ---*/
      
      if (ExtIter == 0) {
        
        if (rank == MASTER_NODE)
          cout << endl << " Setting the moving wall velocities." << endl;
        
        surface_movement->Moving_Walls(geometry_container[MESH_0],
                                       config_container, iZone, ExtIter);
        
        /*--- Update the grid velocities on the coarser multigrid levels after
         setting the moving wall velocities for the finest mesh. ---*/
        
        grid_movement->UpdateMultiGrid(geometry_container, config_container);
        
      }
      
      break;
      
      
    case ROTATING_FRAME:
      
      /*--- Steadily rotating frame: set the grid velocities just once
       before the first iteration flow solver. ---*/
      
      if (ExtIter == 0) {
        
        if (rank == MASTER_NODE) {
          cout << endl << " Setting rotating frame grid velocities";
          cout << " for zone " << iZone << "." << endl;
        }
        
        /*--- Set the grid velocities on all multigrid levels for a steadily
         rotating reference frame. ---*/
        
        for (iMGlevel = 0; iMGlevel <= nMGlevels; iMGlevel++)
          geometry_container[iMGlevel]->SetRotationalVelocity(config_container, iZone);
        
      }
      
      break;
            
    case STEADY_TRANSLATION:
      
      /*--- Set the translational velocity and hold the grid fixed during
       the calculation (similar to rotating frame, but there is no extra
       source term for translation). ---*/
      
      if (ExtIter == 0) {
        
        if (rank == MASTER_NODE)
          cout << endl << " Setting translational grid velocities." << endl;
        
          /*--- Set the translational velocity on all grid levels. ---*/
          
          for (iMGlevel = 0; iMGlevel <= nMGlevels; iMGlevel++)
              geometry_container[iMGlevel]->SetTranslationalVelocity(config_container);
        
      }
      
      break;
      
    case RIGID_MOTION:
      
      if (rank == MASTER_NODE) {
        cout << endl << " Performing rigid mesh transformation." << endl;
      }
      
      /*--- Move each node in the volume mesh using the specified type
       of rigid mesh motion. These routines also compute analytic grid
       velocities for the fine mesh. ---*/
      
      grid_movement->Rigid_Translation(geometry_container[MESH_0],
                                       config_container, iZone, ExtIter);
      grid_movement->Rigid_Plunging(geometry_container[MESH_0],
                                    config_container, iZone, ExtIter);
      grid_movement->Rigid_Pitching(geometry_container[MESH_0],
                                    config_container, iZone, ExtIter);
      grid_movement->Rigid_Rotation(geometry_container[MESH_0],
                                    config_container, iZone, ExtIter);
      
      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/
      
      grid_movement->UpdateMultiGrid(geometry_container, config_container);
      
      break;
      
    case DEFORMING:
      
      if (rank == MASTER_NODE)
        cout << endl << " Updating surface positions." << endl;
      
      /*--- Translating ---*/

      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Translating(geometry_container[MESH_0],
                                         config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Plunging ---*/
      
      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Plunging(geometry_container[MESH_0],
                                            config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Pitching ---*/
      
      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Pitching(geometry_container[MESH_0],
                                            config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Rotating ---*/
      
      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Rotating(geometry_container[MESH_0],
                                            config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/
      
      if (!adjoint) {
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
      }
      
      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/
      
      grid_movement->UpdateMultiGrid(geometry_container, config_container);
      
      break;
      
    case EXTERNAL: case EXTERNAL_ROTATION:
      
      /*--- Apply rigid rotation to entire grid first, if necessary ---*/
      
      if (Kind_Grid_Movement == EXTERNAL_ROTATION) {
        if (rank == MASTER_NODE)
          cout << " Updating node locations by rigid rotation." << endl;
        grid_movement->Rigid_Rotation(geometry_container[MESH_0],
                                      config_container, iZone, ExtIter);
      }
      
      /*--- Load new surface node locations from external files ---*/
      
      if (rank == MASTER_NODE)
        cout << " Updating surface locations from file." << endl;
      surface_movement->SetExternal_Deformation(geometry_container[MESH_0],
                                                config_container, iZone, ExtIter);
      
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/
      
      if (!adjoint) {
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
      }
      
      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/
      
      grid_movement->UpdateMultiGrid(geometry_container, config_container);
      
      break;
      
    case AEROELASTIC: case AEROELASTIC_RIGID_MOTION:
      
      /*--- Apply rigid mesh transformation to entire grid first, if necessary ---*/
      if (IntIter == 0) {
        if (Kind_Grid_Movement == AEROELASTIC_RIGID_MOTION) {
          
          if (rank == MASTER_NODE) {
            cout << endl << " Performing rigid mesh transformation." << endl;
          }
          
          /*--- Move each node in the volume mesh using the specified type
           of rigid mesh motion. These routines also compute analytic grid
           velocities for the fine mesh. ---*/
          
          grid_movement->Rigid_Translation(geometry_container[MESH_0],
                                           config_container, iZone, ExtIter);
          grid_movement->Rigid_Plunging(geometry_container[MESH_0],
                                        config_container, iZone, ExtIter);
          grid_movement->Rigid_Pitching(geometry_container[MESH_0],
                                        config_container, iZone, ExtIter);
          grid_movement->Rigid_Rotation(geometry_container[MESH_0],
                                        config_container, iZone, ExtIter);
          
          /*--- Update the multigrid structure after moving the finest grid,
           including computing the grid velocities on the coarser levels. ---*/
          
          grid_movement->UpdateMultiGrid(geometry_container, config_container);
        }
        
      }
      
      /*--- Use the if statement to move the grid only at selected dual time step iterations. ---*/
      else if (IntIter % config_container->GetAeroelasticIter() ==0) {

        if (rank == MASTER_NODE)
          cout << endl << " Solving aeroelastic equations and updating surface positions." << endl;
        
        /*--- Solve the aeroelastic equations for the new node locations of the moving markers(surfaces) ---*/
        
        solver_container[MESH_0][FLOW_SOL]->Aeroelastic(surface_movement, geometry_container[MESH_0], config_container, ExtIter);
        
        /*--- Deform the volume grid around the new boundary locations ---*/
        
        if (rank == MASTER_NODE)
          cout << " Deforming the volume grid due to the aeroelastic movement." << endl;
        grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                             config_container, true);
        
        /*--- Update the grid velocities on the fine mesh using finite
         differencing based on node coordinates at previous times. ---*/
        
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
        
        /*--- Update the multigrid structure after moving the finest grid,
         including computing the grid velocities on the coarser levels. ---*/
        
        grid_movement->UpdateMultiGrid(geometry_container, config_container);
      }
      
      break;
      
    case ELASTICITY:
      
      if (ExtIter != 0) {
        
        if (rank == MASTER_NODE)
          cout << " Deforming the grid using the Linear Elasticity solution." << endl;
        
        /*--- Update the coordinates of the grid using the linear elasticity solution. ---*/
        for (iPoint = 0; iPoint < geometry_container[MESH_0]->GetnPoint(); iPoint++) {
          
          su2double *U_time_nM1 = solver_container[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_time_n1();
          su2double *U_time_n   = solver_container[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_time_n();
          
          for (iDim = 0; iDim < geometry_container[MESH_0]->GetnDim(); iDim++)
            geometry_container[MESH_0]->node[iPoint]->AddCoord(iDim, U_time_n[iDim] - U_time_nM1[iDim]);
          
        }
        
      }
      
      break;
      
    case NO_MOVEMENT: case GUST: default:
      
      /*--- There is no mesh motion specified for this zone. ---*/
      if (rank == MASTER_NODE)
        cout << "No mesh motion specified." << endl;
      
      break;
  }
  
}

void SetSpectralMethod(CGeometry ***geometry_container, CSolver ****solver_container,
                     CConfig **config_container, unsigned short nZone, unsigned short iZone) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Local variables and initialization ---*/
  //unsigned short iVar, kZone, jZone, iMGlevel;
  unsigned short iVar, jZone, kZone, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetAdjoint());
  if (adjoint) {
    implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  }
  
  /*--- Retrieve values from the config file ---*/
  su2double *U = new su2double[nVar];
  su2double *U_old = new su2double[nVar];
  su2double *Psi = new su2double[nVar];
  su2double *Psi_old = new su2double[nVar];
  su2double *Source = new su2double[nVar];
  su2double deltaU, deltaPsi;
  
  /*--- Compute period of oscillation ---*/
  su2double period = config_container[ZONE_0]->GetSpectralMethod_Period();
  
  /*--- allocate dynamic memory for D ---*/
  su2double **D = new su2double*[nZone];
  for (kZone = 0; kZone < nZone; kZone++) {
    D[kZone] = new su2double[nZone];
  }
  
  /*--- Build the time-spectral operator matrix ---*/
  //ComputeSpectralMethod_Operator(D, period, nZone);
    
  /* Build operator matrix for Harmonic Balance method */
  /* frequency values - hardcoded for now */
  //su2double *Omega_HB = new su2double[nZone];
  su2double *Omega_HB = config_container[ZONE_0]->GetOmega_HB();
    
  ComputeHarmonicBalance_Operator(D, Omega_HB, period, nZone);
  //delete [] Omega_HB;
  
  /*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {
    
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
            
          } 
          
          else {
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
            solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetSpectralMethod_Source(iVar, Source[iVar]);
          } 
          else {
            solver_container[iZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->SetSpectralMethod_Source(iVar, Source[iVar]);
          }
        }
        
      }
    }
  }
  
  
  /*--- Source term for a turbulence model ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == RANS) {
    
    /*--- Extra variables needed if we have a turbulence model. ---*/
    unsigned short nVar_Turb = solver_container[ZONE_0][MESH_0][TURB_SOL]->GetnVar();
    su2double *U_Turb = new su2double[nVar_Turb];
    su2double *Source_Turb = new su2double[nVar_Turb];
    
    /*--- Loop over only the finest mesh level (turbulence is always solved
          on the original grid only). ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][MESH_0]->GetnPoint(); iPoint++) {
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
          solver_container[iZone][MESH_0][TURB_SOL]->node[iPoint]->SetSpectralMethod_Source(iVar, Source_Turb[iVar]);
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
  su2double *sbuf_force = NULL,  *rbuf_force = NULL;
  
  /*--- Other variables ---*/
  unsigned short nVar_Force = 8;
  unsigned long current_iter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Allocate memory for send buffer ---*/
  sbuf_force = new su2double[nVar_Force];
  
  su2double *averages = new su2double[nVar_Force];
  for (iVar = 0; iVar < nVar_Force; iVar++)
    averages[iVar] = 0;
  
  /*--- Allocate memory for receive buffer ---*/
  if (rank == MASTER_NODE) {
    rbuf_force = new su2double[nVar_Force];
    
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
 
  if (rank == MASTER_NODE) { 
    
    /*--- Run through the zones, collecting the forces coefficients
          N.B. Summing across processors within a given zone is being done
          elsewhere. ---*/
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
      
      for (iVar = 0; iVar < nVar_Force; iVar++) {
        rbuf_force[iVar] = sbuf_force[iVar];
      }
      
      TS_Flow_file << kZone << ", ";
      for (iVar = 0; iVar < nVar_Force; iVar++)
        TS_Flow_file << rbuf_force[iVar] << ", ";
      TS_Flow_file << endl;
      
      /*--- Increment the total contributions from each zone, dividing by nZone as you go ---*/
      for (iVar = 0; iVar < nVar_Force; iVar++) {
        averages[iVar] += (1.0/su2double(nZone))*rbuf_force[iVar];
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

void ComputeSpectralMethod_Operator(su2double **D, su2double period, unsigned short nZone) {
  
  unsigned short kZone, jZone;

  /*--- Build the time-spectral operator matrix ---*/
  for (kZone = 0; kZone < nZone; kZone++) {
    for (jZone = 0; jZone < nZone; jZone++) {
      
      if (nZone%2 == 0) {
        
        /*--- For an even number of time instances ---*/
        if (kZone == jZone) {
          D[kZone][jZone] = 0.0;
        } 
        else {
          D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/tan(PI_NUMBER*(kZone-jZone)/nZone));
        }
      }
      else {
        
        /*--- For an odd number of time instances ---*/
        if (kZone == jZone) {
          D[kZone][jZone] = 0.0;
        } 
        else {
          D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/sin(PI_NUMBER*(kZone-jZone)/nZone));
        }
      }

    }
  }

}


/* Matrix - matrix product */
void MatrixMatrixProduct(unsigned short nRows_prod, unsigned short nCols_prod, su2double *matrix_a, su2double *matrix_b, su2double *product)  {
    
    unsigned short iVar, jVar, kVar;
    
    for (iVar = 0; iVar < nRows_prod; iVar++) {
        for (jVar = 0; jVar < nCols_prod; jVar++) {
            product[iVar*nCols_prod+jVar] = 0.0;
            for (kVar = 0; kVar < nRows_prod; kVar++) {
                product[iVar*nCols_prod+jVar] += matrix_a[iVar*nRows_prod+kVar]*matrix_b[kVar*nCols_prod+jVar];
            }
        }
    }
    
}

/* Matrix inverse using Gauss-Jordan elimination */
void InverseBlock(unsigned short nVar_mat, su2double *block, su2double *invBlock) {
    
    unsigned short i, j, k;
    unsigned short temp;
    su2double temporary, r;
    su2double **augmentedmatrix = new su2double*[nVar_mat];
    
    for (i = 0; i < nVar_mat; i++) {
        augmentedmatrix[i] = new su2double[2*nVar_mat];
    }
    
    for(i=0; i<nVar_mat; i++)
        for(j=0; j<nVar_mat; j++)
            augmentedmatrix[i][j] = block[i*nVar_mat+j] ;
    
    /* augmenting with identity matrix of similar dimensions */
    
    for(i=0;i<nVar_mat; i++)
        for(j=nVar_mat; j<2*nVar_mat; j++)
            if(i==j%nVar_mat)
                augmentedmatrix[i][j]=1;
            else
                augmentedmatrix[i][j]=0;
    
    /* using gauss-jordan elimination */
    
    for(j=0; j<nVar_mat; j++)
    {
        temp=j;
        
        /* swapping row if a[j][j]=0 with first non-zero row below jth row */
        
        if (augmentedmatrix[j][j]==0)
        {
            /* Finding first non-zero row */
            for(i=j+1; i<nVar_mat; i++)
                if(augmentedmatrix[i][j]!=0)
                {
                    temp=i;
                    break;
                }
            
            /* Swapping current row with temp row */
            for(k=0; k<2*nVar_mat; k++)
            {
                temporary=augmentedmatrix[j][k] ;
                augmentedmatrix[j][k]=augmentedmatrix[temp][k] ;
                augmentedmatrix[temp][k]=temporary ;
            }
        }
        
        
        
        /* performing row operations to form required identity matrix out of the input matrix */
        
        for(i=0; i<nVar_mat; i++)
            if(i!=j)
            {
                r=augmentedmatrix[i][j];
                for(k=0; k<2*nVar_mat; k++)
                {
                    if (augmentedmatrix[j][j]!=0) {
                        augmentedmatrix[i][k]-=(augmentedmatrix[j][k]/augmentedmatrix[j][j])*r ;
                        
                    }
                }
            }
            else
            {
                r=augmentedmatrix[i][j];
                if (r!=0) {
                    for(k=0; k<2*nVar_mat; k++)
                        augmentedmatrix[i][k]/=r ;
                }
                
                
            }
        
    }
    
    for(i=0; i<nVar_mat; i++)
    {
        for(j=nVar_mat; j<2*nVar_mat; j++)
            invBlock[i*nVar_mat+j-nVar_mat] = augmentedmatrix[i][j];
    }
    
    /*--- delete dynamic memory for augmented matrix ---*/
    for (k = 0; k < nVar_mat; k++) {
        delete [] augmentedmatrix[k];
    }
    delete [] augmentedmatrix;
    
}

/* ----- Computing harmonic balance operator --------- */
void ComputeHarmonicBalance_Operator(su2double **D, su2double *Omega_HB, su2double period, unsigned short nZone) {
    
    unsigned short iVar, jVar;
    
    su2double *Einv_Re  = new su2double[nZone*nZone];
    su2double *Einv_Im  = new su2double[nZone*nZone];
    su2double *E_Re     = new su2double[nZone*nZone];
    su2double *E_Im     = new su2double[nZone*nZone];
    su2double *D_diag   = new su2double[nZone*nZone];
    
    for (iVar = 0; iVar < nZone; iVar++) {
        for (jVar = 0; jVar < nZone; jVar++) {
            E_Re[iVar*nZone+jVar] = cos(Omega_HB[jVar]*(iVar*period/nZone));
            E_Im[iVar*nZone+jVar] = sin(Omega_HB[jVar]*(iVar*period/nZone));
            D_diag[iVar*nZone+jVar] = 0.;
        }
        D_diag[iVar*nZone+iVar] = Omega_HB[iVar];
    }
    
    
    /*---- Finding inverse of matrix E = E_Re + 1i*E_Im -----*/
    su2double *L = new su2double[(2*nZone)*(2*nZone)];
    su2double *L_inv = new su2double[(2*nZone)*(2*nZone)];
    su2double *M = new su2double[(2*nZone)*nZone];
    
    for (iVar = 0; iVar < nZone; iVar++) {
        for (jVar = 0; jVar < nZone; jVar++) {
            L[iVar*(2*nZone)+jVar] = E_Re[iVar*nZone+jVar];
            L[iVar*(2*nZone)+(nZone+jVar)] = -E_Im[iVar*nZone+jVar];
            L[(iVar+nZone)*(2*nZone)+jVar] = E_Im[iVar*nZone+jVar];
            L[(iVar+nZone)*(2*nZone)+(nZone+jVar)] = E_Re[iVar*nZone+jVar];
            M[iVar*nZone+jVar] = 0.;
            M[(iVar+nZone)*nZone+jVar] = 0.;
        }
        M[iVar*nZone+iVar] = 1.;
    }
    
    /* Inverse of E ( 2nZone*nZone )- inv(L)*M */
    InverseBlock(2*nZone, L, L_inv);
    
    su2double *Linv_M = new su2double[(2*nZone)*nZone];
    MatrixMatrixProduct(2*nZone, nZone, L_inv, M, Linv_M);
    
    
    for (iVar = 0; iVar < nZone; iVar++) {
        for (jVar = 0; jVar < nZone; jVar++) {
            Einv_Re[iVar*nZone+jVar] = Linv_M[iVar*nZone+jVar];
            Einv_Im[iVar*nZone+jVar] = Linv_M[(iVar+nZone)*nZone+jVar];
        }
    }
    delete [] L;
    delete [] L_inv;
    delete [] M;
    delete [] Linv_M;
    
    
    /* Spectral operator - inv(E)*D_diag*E= -E_Im*D_diag*Einv_Re - E_Re*D_diag*Einv_Im */
    su2double *prod1 = new su2double[nZone*nZone];
    su2double *prod2 = new su2double[nZone*nZone];
    
    /* E_Im*D_diag*Einv_Re stored in prod1 */
    MatrixMatrixProduct(nZone, nZone, E_Im, D_diag, prod2);
    MatrixMatrixProduct(nZone, nZone, prod2, Einv_Re, prod1);
    
    /* E_Re*D_diag*Einv_Im stored in prod2 */
    MatrixMatrixProduct(nZone, nZone, E_Re, D_diag, E_Im);
    MatrixMatrixProduct(nZone, nZone, E_Im, Einv_Im, prod2);
    
    /* Harmonic Balance spectral operator */
    for (iVar = 0; iVar < nZone; iVar++) {
        for (jVar = 0; jVar < nZone; jVar++) {
            D[iVar][jVar] = -prod1[iVar*nZone+jVar] - prod2[iVar*nZone+jVar];
        }
    }
    
    delete [] prod1;
    delete [] prod2;
    delete [] E_Re;
    delete [] E_Im;
    delete [] Einv_Re;
    delete [] Einv_Im;
    delete [] D_diag;
    
}

void SetSpectralMethod_Velocities(CGeometry ***geometry_container,
                                CConfig **config_container, unsigned short nZone) {
  
	unsigned short iZone, jDegree, iDim, iMGlevel;
	unsigned short nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();
  
	su2double angular_interval = 2.0*PI_NUMBER/(su2double)(nZone);
	su2double *Coord;
	unsigned long iPoint;
  
  
	/*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
	su2double period = config_container[ZONE_0]->GetSpectralMethod_Period();
	su2double deltaT = period/(su2double)(config_container[ZONE_0]->GetnTimeInstances());
  
	/*--- allocate dynamic memory for angular positions (these are the abscissas) ---*/
	su2double *angular_positions = new su2double [nZone];
	for (iZone = 0; iZone < nZone; iZone++) {
		angular_positions[iZone] = iZone*angular_interval;
	}
  
	/*--- find the highest-degree trigonometric polynomial allowed by the Nyquist criterion---*/
	su2double high_degree = (nZone-1)/2.0;
  int highest_degree = SU2_TYPE::Int(high_degree);
  
	/*--- allocate dynamic memory for a given point's coordinates ---*/
	su2double **coords = new su2double *[nZone];
	for (iZone = 0; iZone < nZone; iZone++) {
		coords[iZone] = new su2double [nDim];
	}
  
	/*--- allocate dynamic memory for vectors of Fourier coefficients ---*/
	su2double *a_coeffs = new su2double [highest_degree+1];
	su2double *b_coeffs = new su2double [highest_degree+1];
  
	/*--- allocate dynamic memory for the interpolated positions and velocities ---*/
	su2double *fitted_coords = new su2double [nZone];
	su2double *fitted_velocities = new su2double [nZone];
  
	/*--- Loop over all grid levels ---*/
	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {
    
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
						a_coeffs[jDegree] = a_coeffs[jDegree] + (2.0/(su2double)nZone)*cos(jDegree*angular_positions[iZone])*coords[iZone][iDim];
						b_coeffs[jDegree] = b_coeffs[jDegree] + (2.0/(su2double)nZone)*sin(jDegree*angular_positions[iZone])*coords[iZone][iDim];
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
	for (iZone = 0; iZone < nZone; iZone++) {
		delete [] coords[iZone];
	}
	delete [] coords;
  
}
