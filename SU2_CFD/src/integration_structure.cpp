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
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

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
        solver_container[MainSolver]->BC_Euler_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
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
      case DISPLACEMENT_BOUNDARY:
        solver_container[MainSolver]->BC_Normal_Displacement(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case LOAD_BOUNDARY:
        solver_container[MainSolver]->BC_Normal_Load(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case NEUMANN:
        solver_container[MainSolver]->BC_Neumann(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case LOAD_DIR_BOUNDARY:
        solver_container[MainSolver]->BC_Dir_Load(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case LOAD_SINE_BOUNDARY:
        solver_container[MainSolver]->BC_Sine_Load(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
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
      case CLAMPED_BOUNDARY:
        solver_container[MainSolver]->BC_Clamped(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
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
                                     unsigned short RunTime_EqSystem,
                                     unsigned long Iteration) {

    unsigned short iMarker;

    bool initial_calc = (config->GetExtIter() == 0);                  // Checks if it is the first calculation.
    bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool first_iter = (config->GetIntIter() == 0);                  // Checks if it is the first iteration
    bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.
    unsigned short IterativeScheme = config->GetKind_SpaceIteScheme_FEA();       // Iterative schemes: NEWTON_RAPHSON, MODIFIED_NEWTON_RAPHSON
    unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

    bool restart = config->GetRestart();                                  // Restart solution
    bool initial_calc_restart = (SU2_TYPE::Int(config->GetExtIter()) == config->GetDyn_RestartIter());  // Restart iteration

    /*--- Compute Mass Matrix ---*/
    /*--- The mass matrix is computed only once, at the beginning of the calculation, no matter whether the ---*/
    /*--- problem is linear or nonlinear. This is done in the preprocessing step. ---*/

    /*--- If the analysis is linear, only a the constitutive term of the stiffness matrix has to be computed ---*/
    /*--- This is done only once, at the beginning of the calculation. From then on, K is constant ---*/
    if ((linear_analysis && (initial_calc || dynamic)) ||
      (linear_analysis && restart && initial_calc_restart)) {
      solver_container[MainSolver]->Compute_StiffMatrix(geometry, solver_container, numerics, config);
    }
    else if (!linear_analysis) {
      /*--- If the analysis is nonlinear, also the stress terms need to be computed ---*/
      /*--- If the method is full Newton-Raphson, the stiffness matrix and the nodal term are updated every time ---*/
      /*--- They are calculated together to avoid looping twice over the elements ---*/
      if (IterativeScheme == NEWTON_RAPHSON) {
        /*--- The Jacobian is reinitialized every time in Preprocessing (before calling Space_Integration_FEM) */
        solver_container[MainSolver]->Compute_StiffMatrix_NodalStressRes(geometry, solver_container, numerics, config);
      }

      /*--- If the method is modified Newton-Raphson, the stiffness matrix is only computed once at the beginning of the time-step ---*/
      /*--- Nevertheless, the Nodal Stress Term has to be computed for each iteration ---*/
      else if (IterativeScheme == MODIFIED_NEWTON_RAPHSON) {

        if (first_iter) {
          solver_container[MainSolver]->Compute_StiffMatrix_NodalStressRes(geometry, solver_container, numerics, config);
        }

        else {
          solver_container[MainSolver]->Compute_NodalStressRes(geometry, solver_container, numerics, config);
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
        solver_container[MainSolver]->BC_Dir_Load(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
        break;
          case LOAD_SINE_BOUNDARY:
        solver_container[MainSolver]->BC_Sine_Load(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
        break;
        }
      }
    }

    /*--- Others are not, because they depend on the geometry ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case LOAD_BOUNDARY:
          solver_container[MainSolver]->BC_Normal_Load(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
          break;
        case DAMPER_BOUNDARY:
          solver_container[MainSolver]->BC_Damper(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
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
        SetRestricted_Gradient(RUNTIME_FLOW_SYS, solver_container[iZone][INST_0][iMGLevel][FLOW_SOL], solver_container[iZone][INST_0][iMGLevel+1][FLOW_SOL],
                               geometry[iZone][INST_0][iMGLevel], geometry[iZone][INST_0][iMGLevel+1], config[iZone]);
      }
      
    }
  }
  
}

void CIntegration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                    unsigned short RunTime_EqSystem, unsigned long Iteration) {
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
                                    unsigned short RunTime_EqSystem, unsigned long Iteration) {

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
          solver_container[MainSolver]->BC_Clamped(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
          break;
        case DISP_DIR_BOUNDARY:
          solver_container[MainSolver]->BC_DispDir(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
          break;
        case DISPLACEMENT_BOUNDARY:
          solver_container[MainSolver]->BC_Normal_Displacement(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
          break;
      }

  /*--- Solver linearized system ---*/

    solver_container[MainSolver]->Solve_System(geometry, solver_container, config);

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
      solver_container[MainSolver]->BC_Clamped_Post(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
      break;
//      case DISPLACEMENT_BOUNDARY:
//      solver_container[MainSolver]->BC_Normal_Displacement(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
//      break;
    }

    /*--- Perform the MPI communication of the solution ---*/

    solver_container[MainSolver]->InitiateComms(geometry, config, SOLUTION_FEA);
    solver_container[MainSolver]->CompleteComms(geometry, config, SOLUTION_FEA);

}

void CIntegration::Convergence_Monitoring(CGeometry *geometry, CConfig *config, unsigned long Iteration,
                                          su2double monitor, unsigned short iMesh) {
  
  unsigned short iCounter;

  /*--- Initialize some variables for controlling the output frequency. ---*/
  
  bool DualTime_Iteration = false;
  unsigned long iIntIter = config->GetIntIter();
  unsigned long iExtIter = config->GetExtIter();
  bool Unsteady = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                   (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config->GetWrt_Con_Freq() == 0));
  bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config->GetWrt_Con_Freq_DualTime() == 0));
  bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
  bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config->GetWrt_Con_Freq() == 0));
  bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config->GetWrt_Con_Freq() == 0));
  
  if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) {
    
    bool Already_Converged = Convergence;
    
    /*--- Cauchy based convergence criteria ---*/
    
    if (config->GetConvCriteria() == CAUCHY) {
      
      /*--- Initialize at the fist iteration ---*/
      
      if (Iteration  == 0) {
        Cauchy_Value = 0.0;
        Cauchy_Counter = 0;
        for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
          Cauchy_Serie[iCounter] = 0.0;
      }
      
      Old_Func = New_Func;
      New_Func = monitor;
      Cauchy_Func = fabs(New_Func - Old_Func);
      
      Cauchy_Serie[Cauchy_Counter] = Cauchy_Func;
      Cauchy_Counter++;
      
      if (Cauchy_Counter == config->GetCauchy_Elems()) Cauchy_Counter = 0;
      
      Cauchy_Value = 1;
      if (Iteration  >= config->GetCauchy_Elems()) {
        Cauchy_Value = 0;
        for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
          Cauchy_Value += Cauchy_Serie[iCounter];
      }
      
      if (Cauchy_Value >= config->GetCauchy_Eps()) { Convergence = false; Convergence_FullMG = false; }
      else { Convergence = true; Convergence_FullMG = true; }
      
    }
    
    /*--- Residual based convergence criteria ---*/
    
    if (config->GetConvCriteria() == RESIDUAL) {
      
      /*--- Compute the initial value ---*/
      
      if (Iteration == config->GetStartConv_Iter() ) InitResidual = monitor;
      if (monitor > InitResidual) InitResidual = monitor;
      
      /*--- Check the convergence ---*/
      
      if (((fabs(InitResidual - monitor) >= config->GetOrderMagResidual()) && (monitor < InitResidual))  ||
          (monitor <= config->GetMinLogResidual())) { Convergence = true; Convergence_FullMG = true; }
      else { Convergence = false; Convergence_FullMG = false; }
      
    }
    
    /*--- Do not apply any convergence criteria of the number
     of iterations is less than a particular value ---*/
    
    if (Iteration < config->GetStartConv_Iter()) {
      Convergence = false;
      Convergence_FullMG = false;
    }
    
    if (Already_Converged) { Convergence = true; Convergence_FullMG = true; }
    
    
    /*--- Apply the same convergence criteria to all the processors ---*/
    
#ifdef HAVE_MPI
    
    unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
    sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
    rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;
    
    /*--- Convergence criteria ---*/
    
    sbuf_conv[0] = Convergence;
    SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
    
    /*-- Compute global convergence criteria in the master node --*/
    
    sbuf_conv[0] = 0;
    if (rank == MASTER_NODE) {
      if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
      else sbuf_conv[0] = 0;
    }
    
    SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (sbuf_conv[0] == 1) { Convergence = true; Convergence_FullMG = true; }
    else { Convergence = false; Convergence_FullMG = false; }
    
    delete [] sbuf_conv;
    delete [] rbuf_conv;
    
#endif
    
    /*--- Stop the simulation in case a nan appears, do not save the solution ---*/
    
    if (monitor != monitor) {
      SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);
    }
    
    if (config->GetFinestMesh() != MESH_0 ) Convergence = false;
    
  }
  
}

void CIntegration::SetDualTime_Solver(CGeometry *geometry, CSolver *solver, CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    solver->node[iPoint]->Set_Solution_time_n1();
    solver->node[iPoint]->Set_Solution_time_n();
    
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
  
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    solver->node[iPoint]->SetSolution_time_n();
    solver->node[iPoint]->SetSolution_Vel_time_n();
    solver->node[iPoint]->SetSolution_Accel_time_n();
    
  }
  
  bool fsi = config->GetFSI_Simulation();
  
  /*--- If FSI problem, save the last Aitken relaxation parameter of the previous time step ---*/
  
  if (fsi) {
    
    su2double WAitk=0.0;
    
    WAitk = solver->GetWAitken_Dyn();
    solver->SetWAitken_Dyn_tn1(WAitk);
    
  }
  
  
}

void CIntegration::SetFEM_StructuralSolver(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  
  unsigned long iPoint;
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
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    solver_container[FEA_SOL]->node[iPoint]->SetSolution_time_n();
    solver_container[FEA_SOL]->node[iPoint]->SetSolution_Vel_time_n();
    solver_container[FEA_SOL]->node[iPoint]->SetSolution_Accel_time_n();
    
  }
  
  /*--- If FSI problem, save the last Aitken relaxation parameter of the previous time step ---*/
  
  if (fsi) {
    
    su2double WAitk=0.0;
    
    WAitk = solver_container[FEA_SOL]->GetWAitken_Dyn();
    solver_container[FEA_SOL]->SetWAitken_Dyn_tn1(WAitk);
    
  }
  
}

void CIntegration::Convergence_Monitoring_FEM(CGeometry *geometry, CConfig *config, CSolver *solver, unsigned long iOuterIter) {
  
  su2double Reference_UTOL, Reference_RTOL, Reference_ETOL;
  su2double Residual_UTOL, Residual_RTOL, Residual_ETOL;
  
  bool Already_Converged = Convergence;
  
  Reference_UTOL = config->GetResidual_FEM_UTOL();
  Reference_RTOL = config->GetResidual_FEM_RTOL();
  Reference_ETOL = config->GetResidual_FEM_ETOL();
  
  Residual_UTOL = log10(solver->GetRes_FEM(0));
  Residual_RTOL = log10(solver->GetRes_FEM(1));
  Residual_ETOL = log10(solver->GetRes_FEM(2));
  
  //  cout << "Reference - UTOL: " << Reference_UTOL << " ETOL: " << Reference_ETOL << " RTOL: " << Reference_RTOL << endl;
  //  cout << "Residual - UTOL: " << Residual_UTOL << " ETOL: " << Residual_ETOL << " RTOL: " << Residual_RTOL << endl;
  
  if ((Residual_UTOL <= Reference_UTOL) &&
      (Residual_ETOL <= Reference_ETOL) &&
      (Residual_RTOL <= Reference_RTOL)) {
    Convergence = true;
  }
  
  if (Already_Converged) Convergence = true;
  
  
  /*--- Apply the same convergence criteria to all the processors ---*/
  
#ifdef HAVE_MPI
  
  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;
  
  /*--- Convergence criteria ---*/
  
  sbuf_conv[0] = Convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  
  /*-- Compute global convergence criteria in the master node --*/
  
  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }
  
  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
  
  if (sbuf_conv[0] == 1) { Convergence = true; }
  else { Convergence = false; }
  
  delete [] sbuf_conv;
  delete [] rbuf_conv;
  
#endif
  
}

void CIntegration::Convergence_Monitoring_FEM_Adj(CGeometry *geometry, CConfig *config, CSolver *solver, unsigned long iOuterIter) {

  su2double val_I, Max_Val_I;

  bool Already_Converged = Convergence;

  Max_Val_I = config->GetCriteria_FEM_ADJ();
  val_I     = log10(solver->Get_val_I());

  //  cout << "Reference - UTOL: " << Reference_UTOL << " ETOL: " << Reference_ETOL << " RTOL: " << Reference_RTOL << endl;
  //  cout << "Residual - UTOL: " << Residual_UTOL << " ETOL: " << Residual_ETOL << " RTOL: " << Residual_RTOL << endl;

  if (val_I <= Max_Val_I){
    Convergence = true;
  }

  if (Already_Converged) Convergence = true;


  /*--- Apply the same convergence criteria to all the processors ---*/

#ifdef HAVE_MPI

  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

  /*--- Convergence criteria ---*/

  sbuf_conv[0] = Convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  /*-- Compute global convergence criteria in the master node --*/

  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }

  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

  if (sbuf_conv[0] == 1) { Convergence = true; }
  else { Convergence = false; }

  delete [] sbuf_conv;
  delete [] rbuf_conv;

#endif

}


void CIntegration::Convergence_Monitoring_FSI(CGeometry *fea_geometry, CConfig *fea_config, CSolver *fea_solver, unsigned long iOuterIter) {
  
  su2double FEA_check[2] = {0.0, 0.0};
  su2double magResidualFSI = 0.0, logResidualFSI_initial = 0.0, logResidualFSI = 0.0;
  su2double magResidualFSI_criteria, logResidualFSI_criteria;
  
  unsigned long iExtIter = fea_config->GetExtIter();
  
  unsigned long iPoint, iDim;
  unsigned long nPointDomain, nDim;
  su2double *dispPred, *dispPred_Old;
  su2double deltaU, deltaURad, deltaURes, deltaURes_recv = 0.0;
  
  magResidualFSI_criteria = -1*fea_config->GetOrderMagResidualFSI();
  logResidualFSI_criteria = fea_config->GetMinLogResidualFSI();
  
  deltaURes = 0.0;
  
  /*--- Only when there is movement it makes sense to check convergence (otherwise, it is always converged...) ---*/
  /*--- The same with the first iteration, if we are doing strongly coupled we need at least two. ---*/
  
  if (iOuterIter == 0) {
    /*--- Set the convergence values to 0.0 --*/
    fea_solver->SetFSI_ConvValue(0,0.0);
    fea_solver->SetFSI_ConvValue(1,0.0);
    
  }
  else if (iOuterIter > 0) {
    
    // We loop only over the points that belong to the processor
    nPointDomain = fea_geometry->GetnPointDomain();
    nDim = fea_geometry->GetnDim();
    
    for (iPoint=0; iPoint < nPointDomain; iPoint++) {
      
      deltaURad = 0.0;
      
      dispPred = fea_solver->node[iPoint]->GetSolution_Pred();
      dispPred_Old = fea_solver->node[iPoint]->GetSolution_Pred_Old();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        
        /*--- Compute the deltaU, and add deltaU2 to deltaURad ---*/
        deltaU = dispPred[iDim] - dispPred_Old[iDim];
        deltaURad += deltaU * deltaU;
        
      }
      
      /*--- The residual is the maximum of the values of sqrt(deltaURad) computed ---*/
      deltaURad = sqrt(deltaURad);
      deltaURes = max(deltaURes, deltaURad);
      
    }
    
    // We need to communicate the maximum residual throughout the different processors
    
#ifdef HAVE_MPI
    /*--- We sum the squares of the norms across the different processors ---*/
    SU2_MPI::Allreduce(&deltaURes, &deltaURes_recv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    deltaURes_recv         = deltaURes;
#endif

    /*--- Store the FSI residual ---*/
    fea_solver->SetFSI_Residual(deltaURes_recv);

    if (iOuterIter == 1) {
      fea_solver->SetFSI_ConvValue(0,deltaURes_recv);
      logResidualFSI_initial = log10(deltaURes_recv);
      
      if (logResidualFSI_initial < logResidualFSI_criteria) Convergence_FSI = true;

    }
    else {
      fea_solver->SetFSI_ConvValue(1,deltaURes_recv);
      FEA_check[0] = fea_solver->GetFSI_ConvValue(0);
      logResidualFSI_initial = log10(FEA_check[0]);
      logResidualFSI = log10(deltaURes_recv);
      
      magResidualFSI=logResidualFSI-logResidualFSI_initial;
      
      if ((logResidualFSI < logResidualFSI_criteria) || (magResidualFSI < magResidualFSI_criteria)) Convergence_FSI = true;
    }

  }
    
  /*--- Apply the same convergence criteria to all the processors ---*/
  
#ifdef HAVE_MPI
  
  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;
  
  /*--- Convergence criteria ---*/
  
  sbuf_conv[0] = Convergence_FSI;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  
  /*-- Compute global convergence criteria in the master node --*/
  
  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }
  
  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
  
  if (sbuf_conv[0] == 1) { Convergence_FSI = true; }
  else { Convergence_FSI = false; }
  
  delete [] sbuf_conv;
  delete [] rbuf_conv;
  
#endif
  
  if (rank == MASTER_NODE) {
    
    su2double WAitken;
    unsigned short RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();
    
    if (RelaxMethod_FSI == NO_RELAXATION) {
      WAitken = 1.0;
    }
    else if (RelaxMethod_FSI == FIXED_PARAMETER) {
      WAitken = fea_config->GetAitkenStatRelax();
    }
    else if (RelaxMethod_FSI == AITKEN_DYNAMIC) {
      WAitken = fea_solver->GetWAitken_Dyn();
    }
    else {
      WAitken = 1.0;
      cout << "No relaxation parameter used. " << endl;
    }
    
    /*--- Store the Relaxation coefficient residual ---*/
    fea_solver->SetRelaxCoeff(WAitken);

    cout.precision(6);
    if (iOuterIter == 0) cout << endl <<"BGS_Iter" << "        ExtIter" << "     Relaxation" <<  endl;
    else if (iOuterIter == 1) cout << endl <<"BGS_Iter" << "        ExtIter" << "     Relaxation" << "      Res[ATOL]"  <<  endl;
    else cout << endl <<"BGS_Iter" << "        ExtIter" << "     Relaxation" << "      Res[ATOL]"  << "      Res[OMAG]"<<  endl;
      
    cout.width(8); cout << iOuterIter;
    cout.width(15); cout << iExtIter;
    cout.width(15); cout << WAitken;
    cout.width(15);
    if (iOuterIter == 0) cout << " ";
    else if (iOuterIter == 1) cout << logResidualFSI_initial;
    else cout << logResidualFSI;
    cout.width(15);
    if (iOuterIter < 2) cout << " ";
    else cout << magResidualFSI;
    cout.setf(ios::fixed, ios::floatfield);
    cout << endl;
    cout << endl << "Simulation time: " << fea_config->GetCurrent_DynTime() << ". Time step: " << fea_config->GetDelta_DynTime() << ".";
    cout << endl;
    cout << endl << "------------------------------------------------------------------------- ";
    cout << endl;
  }
  
}
