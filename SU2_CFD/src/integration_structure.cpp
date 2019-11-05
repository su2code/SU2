/*!
 * \file integration_structure.cpp
 * \brief This subroutine includes the space and time integration structure
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../include/integration_structure.hpp"

CIntegration::CIntegration(CConfig *config) {
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Convergence = false;
  Convergence_FSI = false;
  Convergence_FullMG = false;
  Cauchy_Serie = new su2double [config->GetCauchy_Elems()+1];
  InitResidual = 0.0;
}

CIntegration::~CIntegration(void) {
  delete [] Cauchy_Serie;
}

void CIntegration::Space_Integration(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics **numerics,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep,
                                     unsigned short RunTime_EqSystem) {
  unsigned short iMarker, KindBC;
  
  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));

  /*--- Compute inviscid residuals ---*/
  
  switch (config->GetKind_ConvNumScheme()) {
    case SPACE_CENTERED:
      solver_container[MainSolver]->Centered_Residual(geometry, solver_container, numerics[CONV_TERM], config, iMesh, iRKStep);
      break;
    case SPACE_UPWIND:
      solver_container[MainSolver]->Upwind_Residual(geometry, solver_container, numerics[CONV_TERM], config, iMesh);
      break;
    case FINITE_ELEMENT:
      solver_container[MainSolver]->Convective_Residual(geometry, solver_container, numerics[CONV_TERM], config, iMesh, iRKStep);
      break;
  }
  
  /*--- Compute viscous residuals ---*/
  
  solver_container[MainSolver]->Viscous_Residual(geometry, solver_container, numerics[VISC_TERM], config, iMesh, iRKStep);
  
  /*--- Compute source term residuals ---*/

  solver_container[MainSolver]->Source_Residual(geometry, solver_container, numerics[SOURCE_FIRST_TERM], numerics[SOURCE_SECOND_TERM], config, iMesh);
  
  /*--- Add viscous and convective residuals, and compute the Dual Time Source term ---*/
  
  if (dual_time)
    solver_container[MainSolver]->SetResidual_DualTime(geometry, solver_container, config, iRKStep, iMesh, RunTime_EqSystem);
  
  /*--- Boundary conditions that depend on other boundaries (they require MPI sincronization)---*/

  solver_container[MainSolver]->BC_Fluid_Interface(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config);

  /*--- Compute Fourier Transformations for markers where NRBC_BOUNDARY is applied---*/

  if (config->GetBoolGiles() && config->GetSpatialFourier()){
    solver_container[MainSolver]->PreprocessBC_Giles(geometry, config, numerics[CONV_BOUND_TERM], INFLOW);

    solver_container[MainSolver]->PreprocessBC_Giles(geometry, config, numerics[CONV_BOUND_TERM], OUTFLOW);
  }

  /*--- Weak boundary conditions ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    KindBC = config->GetMarker_All_KindBC(iMarker);
    switch (KindBC) {
      case EULER_WALL:
        solver_container[MainSolver]->BC_Euler_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ACTDISK_INLET:
        solver_container[MainSolver]->BC_ActDisk_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ENGINE_INFLOW:
        solver_container[MainSolver]->BC_Engine_Inflow(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case INLET_FLOW:
        solver_container[MainSolver]->BC_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ACTDISK_OUTLET:
        solver_container[MainSolver]->BC_ActDisk_Outlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ENGINE_EXHAUST:
        solver_container[MainSolver]->BC_Engine_Exhaust(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case SUPERSONIC_INLET:
        solver_container[MainSolver]->BC_Supersonic_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case OUTLET_FLOW:
        solver_container[MainSolver]->BC_Outlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case SUPERSONIC_OUTLET:
        solver_container[MainSolver]->BC_Supersonic_Outlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case GILES_BOUNDARY:
        solver_container[MainSolver]->BC_Giles(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
      	break;
      case RIEMANN_BOUNDARY:
      	if (config->GetBoolTurbomachinery()){
      		solver_container[MainSolver]->BC_TurboRiemann(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
      	}
      	else{
          solver_container[MainSolver]->BC_Riemann(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
      	}
      	break;
      case FAR_FIELD:
        solver_container[MainSolver]->BC_Far_Field(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case SYMMETRY_PLANE:
        solver_container[MainSolver]->BC_Sym_Plane(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ELECTRODE_BOUNDARY:
        solver_container[MainSolver]->BC_Electrode(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case DIELEC_BOUNDARY:
        solver_container[MainSolver]->BC_Dielec(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case NEUMANN:
        solver_container[MainSolver]->BC_Neumann(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
    }
  }
  
  /*--- Strong boundary conditions (Navier-Stokes and Dirichlet type BCs) ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
        solver_container[MainSolver]->BC_Isothermal_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case HEAT_FLUX:
        solver_container[MainSolver]->BC_HeatFlux_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case DIRICHLET:
        solver_container[MainSolver]->BC_Dirichlet(geometry, solver_container, config, iMarker);
        break;
      case CUSTOM_BOUNDARY:
        solver_container[MainSolver]->BC_Custom(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case CHT_WALL_INTERFACE: 
        if ((MainSolver == HEAT_SOL) || (MainSolver == FLOW_SOL && ((config->GetKind_Regime() == COMPRESSIBLE) || config->GetEnergy_Equation()))) {
          solver_container[MainSolver]->BC_ConjugateHeat_Interface(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        }
        else {
          solver_container[MainSolver]->BC_HeatFlux_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        }
        break;
    }
  
  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated corectly during the communications. ---*/
  
  if (config->GetnMarker_Periodic() > 0) {
    solver_container[MainSolver]->BC_Periodic(geometry, solver_container, numerics[CONV_BOUND_TERM], config);
  }
  
}

void CIntegration::Space_Integration_FEM(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics **numerics,
                                     CConfig *config,
                                     unsigned short RunTime_EqSystem) {

    unsigned short iMarker;

    bool initial_calc = (config->GetTimeIter() == 0);                  // Checks if it is the first calculation.
    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool first_iter = (config->GetInnerIter() == 0);                  // Checks if it is the first iteration
    unsigned short IterativeScheme = config->GetKind_SpaceIteScheme_FEA();       // Iterative schemes: NEWTON_RAPHSON, MODIFIED_NEWTON_RAPHSON
    unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

    bool restart = config->GetRestart();                                  // Restart solution
    bool initial_calc_restart = (SU2_TYPE::Int(config->GetTimeIter()) ==  SU2_TYPE::Int(config->GetRestart_Iter()));  // Restart iteration

    /*--- Compute Mass Matrix ---*/
    /*--- The mass matrix is computed only once, at the beginning of the calculation, no matter whether the ---*/
    /*--- problem is linear or nonlinear. This is done in the preprocessing step. ---*/

    /*--- If the analysis is linear, only a the constitutive term of the stiffness matrix has to be computed ---*/
    /*--- This is done only once, at the beginning of the calculation. From then on, K is constant ---*/
    if ((linear_analysis && (initial_calc && first_iter)) ||
      (linear_analysis && restart && initial_calc_restart)) {
      solver_container[MainSolver]->Compute_StiffMatrix(geometry, numerics, config);
    }
    else if (!linear_analysis) {
      /*--- If the analysis is nonlinear, also the stress terms need to be computed ---*/
      /*--- If the method is full Newton-Raphson, the stiffness matrix and the nodal term are updated every time ---*/
      /*--- They are calculated together to avoid looping twice over the elements ---*/
      if (IterativeScheme == NEWTON_RAPHSON) {
        /*--- The Jacobian is reinitialized every time in Preprocessing (before calling Space_Integration_FEM) */
        solver_container[MainSolver]->Compute_StiffMatrix_NodalStressRes(geometry, numerics, config);
      }

      /*--- If the method is modified Newton-Raphson, the stiffness matrix is only computed once at the beginning of the time-step ---*/
      /*--- Nevertheless, the Nodal Stress Term has to be computed for each iteration ---*/
      else if (IterativeScheme == MODIFIED_NEWTON_RAPHSON) {

        if (first_iter) {
          solver_container[MainSolver]->Compute_StiffMatrix_NodalStressRes(geometry, numerics, config);
        }

        else {
          solver_container[MainSolver]->Compute_NodalStressRes(geometry, numerics, config);
        }

      }

    }

    /*--- Apply the NATURAL BOUNDARY CONDITIONS (loads). ---*/
    /*--- If there are FSI loads, they have to be previously applied at other level involving both zones. ---*/

    /*--- Some external loads may be considered constant over the time step ---*/
    if (first_iter) {
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        switch (config->GetMarker_All_KindBC(iMarker)) {
          case LOAD_DIR_BOUNDARY:
        solver_container[MainSolver]->BC_Dir_Load(geometry, numerics[FEA_TERM], config, iMarker);
        break;
          case LOAD_SINE_BOUNDARY:
        solver_container[MainSolver]->BC_Sine_Load(geometry, numerics[FEA_TERM], config, iMarker);
        break;
        }
      }
    }

    /*--- Others are not, because they depend on the geometry ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case LOAD_BOUNDARY:
          solver_container[MainSolver]->BC_Normal_Load(geometry, numerics[FEA_TERM], config, iMarker);
          break;
        case DAMPER_BOUNDARY:
          solver_container[MainSolver]->BC_Damper(geometry, numerics[FEA_TERM], config, iMarker);
        break;
      }
    }

}

void CIntegration::Adjoint_Setup(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                                 unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) {
  
  unsigned short iMGLevel;
  
  if ( ( (RunTime_EqSystem == RUNTIME_ADJFLOW_SYS) && (Iteration == 0) ) ) {
    for (iMGLevel = 0; iMGLevel <= config[iZone]->GetnMGLevels(); iMGLevel++) {
      
      /*--- Set the time step in all the MG levels ---*/
      
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTime_Step(geometry[iZone][INST_0][iMGLevel], solver_container[iZone][INST_0][iMGLevel], config[iZone], iMGLevel, Iteration);
      
      /*--- Set the force coefficients ---*/
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CD(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CD());
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CL(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CL());
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CT(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CT());
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CQ(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CQ());
      
      /*--- Restrict solution and gradients to the coarse levels ---*/
      
      if (iMGLevel != config[iZone]->GetnMGLevels()) {
        SetRestricted_Solution(RUNTIME_FLOW_SYS, solver_container[iZone][INST_0][iMGLevel][FLOW_SOL], solver_container[iZone][INST_0][iMGLevel+1][FLOW_SOL],
                               geometry[iZone][INST_0][iMGLevel], geometry[iZone][INST_0][iMGLevel+1], config[iZone]);
//        ToDo: The flow solvers do not use the conservative variable gradients
//        SetRestricted_Gradient(RUNTIME_FLOW_SYS, solver_container[iZone][INST_0][iMGLevel][FLOW_SOL], solver_container[iZone][INST_0][iMGLevel+1][FLOW_SOL],
//                               geometry[iZone][INST_0][iMGLevel], geometry[iZone][INST_0][iMGLevel+1], config[iZone]);
      }
      
    }
  }
  
}

void CIntegration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                    unsigned short RunTime_EqSystem) {
  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
  unsigned short KindSolver = config->GetKind_Solver();
  
  /*--- Perform the time integration ---*/

  /*--- Fluid time integration schemes ---*/

  if (KindSolver != FEM_ELASTICITY) {

    switch (config->GetKind_TimeIntScheme()) {
      case (RUNGE_KUTTA_EXPLICIT):
        solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iRKStep);
        break;
      case (CLASSICAL_RK4_EXPLICIT):
        solver_container[MainSolver]->ClassicalRK4_Iteration(geometry, solver_container, config, iRKStep);
        break;
      case (EULER_EXPLICIT):
        solver_container[MainSolver]->ExplicitEuler_Iteration(geometry, solver_container, config);
        break;
      case (EULER_IMPLICIT):
        solver_container[MainSolver]->ImplicitEuler_Iteration(geometry, solver_container, config);
        break;
    }

   /*--- Structural time integration schemes ---*/
  
  }
  else if (KindSolver == FEM_ELASTICITY) {

    switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iRKStep);
      break;
    case (NEWMARK_IMPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Iteration(geometry, solver_container, config);
      break;
    case (GENERALIZED_ALPHA):
      solver_container[MainSolver]->ImplicitEuler_Iteration(geometry, solver_container, config);
      break;
    }
  }

}

void CIntegration::Time_Integration_FEM(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config,
                                    unsigned short RunTime_EqSystem) {

  unsigned short iMarker;

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Set the Jacobian according to the different time integration methods ---*/

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Iteration(geometry, solver_container, config);
      break;
    case (NEWMARK_IMPLICIT):
      solver_container[MainSolver]->ImplicitNewmark_Iteration(geometry, solver_container, config);
      break;
    case (GENERALIZED_ALPHA):
      solver_container[MainSolver]->GeneralizedAlpha_Iteration(geometry, solver_container, config);
      break;
    }

  /*--- Apply ESSENTIAL BOUNDARY CONDITIONS ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case CLAMPED_BOUNDARY:
          solver_container[MainSolver]->BC_Clamped(geometry, numerics[FEA_TERM], config, iMarker);
          break;
        case DISP_DIR_BOUNDARY:
          solver_container[MainSolver]->BC_DispDir(geometry, numerics[FEA_TERM], config, iMarker);
          break;
        case DISPLACEMENT_BOUNDARY:
          solver_container[MainSolver]->BC_Normal_Displacement(geometry, numerics[CONV_BOUND_TERM], config, iMarker);
          break;
      }

  /*--- Solver linearized system ---*/

    solver_container[MainSolver]->Solve_System(geometry, config);

  /*--- Update solution ---*/

    switch (config->GetKind_TimeIntScheme_FEA()) {
      case (CD_EXPLICIT):
        solver_container[MainSolver]->ImplicitNewmark_Update(geometry, solver_container, config);
        break;
      case (NEWMARK_IMPLICIT):
        solver_container[MainSolver]->ImplicitNewmark_Update(geometry, solver_container, config);
        break;
      case (GENERALIZED_ALPHA):
        solver_container[MainSolver]->GeneralizedAlpha_UpdateDisp(geometry, solver_container, config);
        break;
      }

  /*--- Reinforce ESSENTIAL BOUNDARY CONDITIONS: avoids accumulation of numerical error ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case CLAMPED_BOUNDARY:
      solver_container[MainSolver]->BC_Clamped_Post(geometry, numerics[FEA_TERM], config, iMarker);
      break;
//      case DISPLACEMENT_BOUNDARY:
//      solver_container[MainSolver]->BC_Normal_Displacement(geometry, numerics[CONV_BOUND_TERM], config, iMarker);
//      break;
    }

    /*--- Perform the MPI communication of the solution ---*/

    solver_container[MainSolver]->InitiateComms(geometry, config, SOLUTION_FEA);
    solver_container[MainSolver]->CompleteComms(geometry, config, SOLUTION_FEA);

}


void CIntegration::SetDualTime_Solver(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned short iMesh) {

  unsigned long iPoint;

  solver->GetNodes()->Set_Solution_time_n1();
  solver->GetNodes()->Set_Solution_time_n();

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    geometry->node[iPoint]->SetVolume_nM1();
    geometry->node[iPoint]->SetVolume_n();
    
    /*--- Store old coordinates in case there is grid movement ---*/
    
    if (config->GetGrid_Movement()) {
      geometry->node[iPoint]->SetCoord_n1();
      geometry->node[iPoint]->SetCoord_n();
    }
  }
  
  /*--- Store old aeroelastic solutions ---*/
  if (config->GetGrid_Movement() && config->GetAeroelastic_Simulation() && (iMesh == MESH_0)) {
    config->SetAeroelastic_n1();
    config->SetAeroelastic_n();
    
    /*--- Also communicate plunge and pitch to the master node. Needed for output in case of parallel run ---*/
    
#ifdef HAVE_MPI
    su2double plunge, pitch, *plunge_all = NULL, *pitch_all = NULL;
    unsigned short iMarker, iMarker_Monitoring;
    unsigned long iProcessor, owner, *owner_all = NULL;
    
    string Marker_Tag, Monitoring_Tag;
    int nProcessor = size;

    /*--- Only if master node allocate memory ---*/
    
    if (rank == MASTER_NODE) {
      plunge_all = new su2double[nProcessor];
      pitch_all  = new su2double[nProcessor];
      owner_all  = new unsigned long[nProcessor];
    }
    
    /*--- Find marker and give it's plunge and pitch coordinate to the master node ---*/
    
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag) { owner = 1; break;
        } else {
          owner = 0;
        }
        
      }
      plunge = config->GetAeroelastic_plunge(iMarker_Monitoring);
      pitch  = config->GetAeroelastic_pitch(iMarker_Monitoring);
      
      /*--- Gather the data on the master node. ---*/
      
      SU2_MPI::Gather(&plunge, 1, MPI_DOUBLE, plunge_all, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(&pitch, 1, MPI_DOUBLE, pitch_all, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(&owner, 1, MPI_UNSIGNED_LONG, owner_all, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      
      /*--- Set plunge and pitch on the master node ---*/
      
      if (rank == MASTER_NODE) {
        for (iProcessor = 0; iProcessor < (unsigned long)nProcessor; iProcessor++) {
          if (owner_all[iProcessor] == 1) {
            config->SetAeroelastic_plunge(iMarker_Monitoring, plunge_all[iProcessor]);
            config->SetAeroelastic_pitch(iMarker_Monitoring, pitch_all[iProcessor]);
            break;
          }
        }
      }
      
    }
    
    if (rank == MASTER_NODE) {
      delete [] plunge_all;
      delete [] pitch_all;
      delete [] owner_all;
    }
#endif
  }
  
}

void CIntegration::SetStructural_Solver(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned short iMesh) {

  solver->GetNodes()->Set_Solution_time_n();
  solver->GetNodes()->SetSolution_Vel_time_n();
  solver->GetNodes()->SetSolution_Accel_time_n();

  bool fsi = config->GetFSI_Simulation();

  /*--- If FSI problem, save the last Aitken relaxation parameter of the previous time step ---*/

  if (fsi) {

    su2double WAitk=0.0;

    WAitk = solver->GetWAitken_Dyn();
    solver->SetWAitken_Dyn_tn1(WAitk);

  }
}

void CIntegration::SetFEM_StructuralSolver(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  bool fsi = config->GetFSI_Simulation();
  
  /*--- Update the solution according to the integration scheme used ---*/
  
  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      break;
    case (NEWMARK_IMPLICIT):
      if (fsi) solver_container[FEA_SOL]->ImplicitNewmark_Relaxation(geometry, solver_container, config);
      break;
    case (GENERALIZED_ALPHA):
      //if (fsi)  solver_container[FEA_SOL]->Update_StructSolution(geometry, solver_container, config);
      solver_container[FEA_SOL]->GeneralizedAlpha_UpdateSolution(geometry, solver_container, config);
      solver_container[FEA_SOL]->GeneralizedAlpha_UpdateLoads(geometry, solver_container, config);
      break;
  }

  /*--- Store the solution at t+1 as solution at t, both for the local points and for the halo points ---*/

  solver_container[FEA_SOL]->GetNodes()->Set_Solution_time_n();
  solver_container[FEA_SOL]->GetNodes()->SetSolution_Vel_time_n();
  solver_container[FEA_SOL]->GetNodes()->SetSolution_Accel_time_n();

  /*--- If FSI problem, save the last Aitken relaxation parameter of the previous time step ---*/

  if (fsi) {

    su2double WAitk=0.0;

    WAitk = solver_container[FEA_SOL]->GetWAitken_Dyn();
    solver_container[FEA_SOL]->SetWAitken_Dyn_tn1(WAitk);

  }
}
