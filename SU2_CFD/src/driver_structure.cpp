/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 4.0.1 "Cardinal"
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

#include "../include/driver_structure.hpp"

CDriver::CDriver(CIteration **iteration_container,
                 CSolver ****solver_container,
                 CGeometry ***geometry_container,
                 CIntegration ***integration_container,
                 CNumerics *****numerics_container,
                 CInterpolator ***interpolator_container,
                 CTransfer ***transfer_container,
                 CConfig **config_container,
                 unsigned short val_nZone,
                 unsigned short val_nDim) {
  
  unsigned short iMesh, iZone, jZone, iSol;
  unsigned short nZone, nDim;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  nZone = val_nZone;
  nDim = val_nDim;
  

  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Instantiate the type of physics iteration to be executed within each zone. For
     example, one can execute the same physics across multiple zones (mixing plane),
     different physics in different zones (fluid-structure interaction), or couple multiple
     systems tightly within a single zone by creating a new iteration class (e.g., RANS). ---*/
    if (rank == MASTER_NODE){
      cout << endl <<"------------------------ Iteration Preprocessing ------------------------" << endl;
    }
    Iteration_Preprocessing(iteration_container, config_container, iZone);

    /*--- Definition of the solver class: solver_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS].
     The solver classes are specific to a particular set of governing equations,
     and they contain the subroutines with instructions for computing each spatial
     term of the PDE, i.e. loops over the edges to compute convective and viscous
     fluxes, loops over the nodes to compute source terms, and routines for
     imposing various boundary condition type for the PDE. ---*/
    if (rank == MASTER_NODE)
      cout << endl <<"------------------------- Solver Preprocessing --------------------------" << endl;

    solver_container[iZone] = new CSolver** [config_container[iZone]->GetnMGLevels()+1];
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++)
      solver_container[iZone][iMesh] = NULL;
    
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
      solver_container[iZone][iMesh] = new CSolver* [MAX_SOLS];
      for (iSol = 0; iSol < MAX_SOLS; iSol++)
        solver_container[iZone][iMesh][iSol] = NULL;
    }
    Solver_Preprocessing(solver_container[iZone], geometry_container[iZone],
                         config_container[iZone]);
    
    
    if (rank == MASTER_NODE)
      cout << endl <<"----------------- Integration and Numerics Preprocessing ----------------" << endl;
    
    /*--- Definition of the integration class: integration_container[#ZONES][#EQ_SYSTEMS].
     The integration class orchestrates the execution of the spatial integration
     subroutines contained in the solver class (including multigrid) for computing
     the residual at each node, R(U) and then integrates the equations to a
     steady state or time-accurately. ---*/
    
    integration_container[iZone] = new CIntegration*[MAX_SOLS];
    Integration_Preprocessing(integration_container[iZone], geometry_container[iZone],
                              config_container[iZone]);
    
    
    if (rank == MASTER_NODE) cout << "Integration Preprocessing." << endl;
    
    /*--- Definition of the numerical method class:
     numerics_container[#ZONES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
     The numerics class contains the implementation of the numerical methods for
     evaluating convective or viscous fluxes between any two nodes in the edge-based
     data structure (centered, upwind, galerkin), as well as any source terms
     (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/
    
    numerics_container[iZone] = new CNumerics***[config_container[iZone]->GetnMGLevels()+1];
    Numerics_Preprocessing(numerics_container[iZone], solver_container[iZone],
                           geometry_container[iZone], config_container[iZone]);
    
    if (rank == MASTER_NODE) cout << "Numerics Preprocessing." << endl;
    
  }
  
	/*--- Definition of the interface and transfer conditions between different zones.
	 *--- The transfer container is defined for zones paired one to one.
	 *--- This only works for a multizone problem (nZone > 1).
	 *--- Also, at the moment this capability is limited to two zones (nZone < 3).
	 *--- This will change in the future. ---*/

	if (rank == MASTER_NODE)
		cout << endl <<"------------------- Multizone Interface Preprocessing -------------------" << endl;


	if ((nZone > 1) && (nZone < 3)) {

		for (iZone = 0; iZone < nZone; iZone++){
			transfer_container[iZone] = new CTransfer*[nZone];
			interpolator_container[iZone] = new CInterpolator*[nZone];
			for (jZone = 0; jZone < nZone; jZone++){
				transfer_container[iZone][jZone] = NULL;
				interpolator_container[iZone][jZone] = NULL;
			}
		}

		Interface_Preprocessing(transfer_container, interpolator_container, geometry_container,
				config_container, solver_container, nZone, nDim);

	}

}


void CDriver::Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry,
                                   CConfig *config) {
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  lin_euler, lin_ns,
  tne2_euler, tne2_ns,
  adj_tne2_euler, adj_tne2_ns,
  poisson, wave, fea, heat, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, machine_learning, transition,
  template_solver, disc_adj;
  
  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  lin_euler        = false;  lin_ns          = false;
  tne2_euler       = false;  tne2_ns         = false;
  adj_tne2_euler   = false;  adj_tne2_ns     = false;
  spalart_allmaras = false;  menter_sst      = false;   machine_learning = false;
  poisson          = false;  neg_spalart_allmaras = false;
  wave             = false;	 disc_adj        = false;
  fea              = false;  fem = false;
  heat             = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case TNE2_EULER : tne2_euler = true; break;
    case TNE2_NAVIER_STOKES: tne2_ns = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case LINEAR_ELASTICITY: fea = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case ADJ_TNE2_EULER : tne2_euler = true; adj_tne2_euler = true; break;
    case ADJ_TNE2_NAVIER_STOKES : tne2_ns = true; adj_tne2_ns = true; break;
    case LIN_EULER: euler = true; lin_euler = true; break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;           break;
      case ML:     machine_learning = true;     break;
        
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(EXIT_FAILURE); break;
    }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    
    /*--- Allocate solution for a template problem ---*/
    if (template_solver) {
      solver_container[iMGlevel][TEMPLATE_SOL] = new CTemplateSolver(geometry[iMGlevel], config);
    }
    
    /*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
    if (euler) {
      solver_container[iMGlevel][FLOW_SOL] = new CEulerSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (ns) {
      solver_container[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (tne2_euler) {
      solver_container[iMGlevel][TNE2_SOL] = new CTNE2EulerSolver(geometry[iMGlevel], config, iMGlevel);
      solver_container[iMGlevel][TNE2_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_TNE2_SYS, false);
    }
    if (tne2_ns) {
      solver_container[iMGlevel][TNE2_SOL] = new CTNE2NSSolver(geometry[iMGlevel], config, iMGlevel);
      solver_container[iMGlevel][TNE2_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_TNE2_SYS, false);
    }
    if (turbulent) {
      if (spalart_allmaras) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel, solver_container[iMGlevel][FLOW_SOL]->GetFluidModel() );
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      else if (neg_spalart_allmaras) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel, solver_container[iMGlevel][FLOW_SOL]->GetFluidModel() );
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      else if (menter_sst) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSSTSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      else if (machine_learning) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbMLSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      if (transition) {
        solver_container[iMGlevel][TRANS_SOL] = new CTransLMSolver(geometry[iMGlevel], config, iMGlevel);
      }
    }
    if (poisson) {
      solver_container[iMGlevel][POISSON_SOL] = new CPoissonSolver(geometry[iMGlevel], config);
    }
    if (wave) {
      solver_container[iMGlevel][WAVE_SOL] = new CWaveSolver(geometry[iMGlevel], config);
    }
    if (heat) {
      solver_container[iMGlevel][HEAT_SOL] = new CHeatSolver(geometry[iMGlevel], config);
    }
    if (fea) {
      solver_container[iMGlevel][FEA_SOL] = new CFEASolver(geometry[iMGlevel], config);
    }
    if (fem) {
      solver_container[iMGlevel][FEA_SOL] = new CFEM_ElasticitySolver(geometry[iMGlevel], config);
    }
    
    /*--- Allocate solution for adjoint problem ---*/
    if (adj_euler) {
      solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjEulerSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (adj_ns) {
      solver_container[iMGlevel][ADJFLOW_SOL] = new CAdjNSSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (adj_tne2_euler) {
      solver_container[iMGlevel][ADJTNE2_SOL] = new CAdjTNE2EulerSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (adj_tne2_ns) {
      solver_container[iMGlevel][ADJTNE2_SOL] = new CAdjTNE2NSSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (adj_turb) {
      solver_container[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolver(geometry[iMGlevel], config, iMGlevel);
    }
    
    /*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
    if (lin_euler) {
      solver_container[iMGlevel][LINFLOW_SOL] = new CLinEulerSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (lin_ns) {
      cout <<"Equation not implemented." << endl; exit(EXIT_FAILURE); break;
    }
    
    if (disc_adj) {
      solver_container[iMGlevel][ADJFLOW_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver_container[iMGlevel][FLOW_SOL], RUNTIME_FLOW_SYS, iMGlevel);
      if (turbulent)
        solver_container[iMGlevel][ADJTURB_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver_container[iMGlevel][TURB_SOL], RUNTIME_TURB_SYS, iMGlevel);
    }
  }
}

void CDriver::Integration_Preprocessing(CIntegration **integration_container,
                                        CGeometry **geometry, CConfig *config) {
  
  bool
  euler, adj_euler, lin_euler,
  ns, adj_ns, lin_ns,
  turbulent, adj_turb,
  tne2_euler, adj_tne2_euler,
  tne2_ns, adj_tne2_ns,
  poisson, wave, fea, fem, heat, template_solver, transition, disc_adj;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false; lin_euler         = false;
  ns               = false; adj_ns           = false; lin_ns            = false;
  turbulent        = false; adj_turb         = false;
  tne2_euler       = false; adj_tne2_euler   = false;
  tne2_ns          = false; adj_tne2_ns      = false;
  poisson          = false; disc_adj         = false;
  wave             = false;
  heat             = false;
  fea              = false; fem = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case TNE2_EULER : tne2_euler = true; break;
    case TNE2_NAVIER_STOKES: tne2_ns = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case LINEAR_ELASTICITY: fea = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_TNE2_EULER : tne2_euler = true; adj_tne2_euler = true; break;
    case ADJ_TNE2_NAVIER_STOKES : tne2_ns = true; adj_tne2_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case LIN_EULER: euler = true; lin_euler = true; break;
    case DISC_ADJ_EULER : euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS : ns = true; turbulent = true; disc_adj = true; break;
      
  }
  
  /*--- Allocate solution for a template problem ---*/
  if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);
  
  /*--- Allocate solution for direct problem ---*/
  if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (ns) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (tne2_euler) integration_container[TNE2_SOL] = new CMultiGridIntegration(config);
  if (tne2_ns) integration_container[TNE2_SOL] = new CMultiGridIntegration(config);
  if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
  if (transition) integration_container[TRANS_SOL] = new CSingleGridIntegration(config);
  if (poisson) integration_container[POISSON_SOL] = new CSingleGridIntegration(config);
  if (wave) integration_container[WAVE_SOL] = new CSingleGridIntegration(config);
  if (heat) integration_container[HEAT_SOL] = new CSingleGridIntegration(config);
  if (fea) integration_container[FEA_SOL] = new CStructuralIntegration(config);
  if (fem) integration_container[FEA_SOL] = new CStructuralIntegration(config);
  
  /*--- Allocate solution for adjoint problem ---*/
  if (adj_euler) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_ns) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_tne2_euler) integration_container[ADJTNE2_SOL] = new CMultiGridIntegration(config);
  if (adj_tne2_ns) integration_container[ADJTNE2_SOL] = new CMultiGridIntegration(config);
  if (adj_turb) integration_container[ADJTURB_SOL] = new CSingleGridIntegration(config);
  
  /*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
  if (lin_euler) integration_container[LINFLOW_SOL] = new CMultiGridIntegration(config);
  if (lin_ns) { cout <<"Equation not implemented." << endl; exit(EXIT_FAILURE); }
  
  if (disc_adj) integration_container[ADJFLOW_SOL] = new CIntegration(config);
  
}

void CDriver::Numerics_Preprocessing(CNumerics ****numerics_container,
                                     CSolver ***solver_container, CGeometry **geometry,
                                     CConfig *config) {
  
  unsigned short iMGlevel, iSol, nDim,
  
  nVar_Template         = 0,
  nVar_Flow             = 0,
  nVar_Trans            = 0,
  nVar_TNE2             = 0,
  nPrimVar_TNE2         = 0,
  nPrimVarGrad_TNE2     = 0,
  nVar_Turb             = 0,
  nVar_Adj_Flow         = 0,
  nVar_Adj_Turb         = 0,
  nVar_Adj_TNE2         = 0,
  nPrimVar_Adj_TNE2     = 0,
  nPrimVarGrad_Adj_TNE2 = 0,
  nVar_Poisson          = 0,
  nVar_FEA              = 0,
  nVar_FEM				= 0,
  nVar_Wave             = 0,
  nVar_Heat             = 0,
  nVar_Lin_Flow         = 0;
  
  su2double *constants = NULL;
  
  bool
  euler, adj_euler, lin_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  tne2_euler, adj_tne2_euler,
  tne2_ns, adj_tne2_ns,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, machine_learning,
  poisson,
  wave,
  fea, fem,
  heat,
  transition,
  template_solver;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;   turbulent        = false;
  poisson          = false;
  adj_euler        = false;   adj_ns           = false;   adj_turb         = false;
  wave             = false;   heat             = false;   fea              = false;   spalart_allmaras = false; neg_spalart_allmaras = false;
  tne2_euler       = false;   tne2_ns          = false;	  fem				= false;
  adj_tne2_euler   = false;	  adj_tne2_ns      = false;
  lin_euler        = false;   menter_sst       = false;    machine_learning = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case DISC_ADJ_EULER: euler = true; break;
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES: ns = true; break;
    case RANS : case DISC_ADJ_RANS:  ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case TNE2_EULER : tne2_euler = true; break;
    case TNE2_NAVIER_STOKES: tne2_ns = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case LINEAR_ELASTICITY: fea = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_TNE2_EULER : tne2_euler = true; adj_tne2_euler = true; break;
    case ADJ_TNE2_NAVIER_STOKES : tne2_ns = true; adj_tne2_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case LIN_EULER: euler = true; lin_euler = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case ML:     machine_learning = true;     break;
      case SST:    menter_sst = true; constants = solver_container[MESH_0][TURB_SOL]->GetConstants(); break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(EXIT_FAILURE); break;
    }
  
  /*--- Number of variables for the template ---*/
  
  if (template_solver) nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  
  /*--- Number of variables for direct problem ---*/
  
  if (euler)        nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (ns)           nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (turbulent)    nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
  if (transition)   nVar_Trans = solver_container[MESH_0][TRANS_SOL]->GetnVar();
  if ((tne2_euler) || (tne2_ns)) {
    nVar_TNE2         = solver_container[MESH_0][TNE2_SOL]->GetnVar();
    nPrimVar_TNE2     = solver_container[MESH_0][TNE2_SOL]->GetnPrimVar();
    nPrimVarGrad_TNE2 = solver_container[MESH_0][TNE2_SOL]->GetnPrimVarGrad();
  }
  if (poisson)      nVar_Poisson = solver_container[MESH_0][POISSON_SOL]->GetnVar();
  
  if (wave)				nVar_Wave = solver_container[MESH_0][WAVE_SOL]->GetnVar();
  if (fea)				nVar_FEA = solver_container[MESH_0][FEA_SOL]->GetnVar();
  if (fem)				nVar_FEM = solver_container[MESH_0][FEA_SOL]->GetnVar();
  if (heat)				nVar_Heat = solver_container[MESH_0][HEAT_SOL]->GetnVar();
  
  /*--- Number of variables for adjoint problem ---*/
  
  if (adj_euler)        nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_ns)           nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_turb)         nVar_Adj_Turb = solver_container[MESH_0][ADJTURB_SOL]->GetnVar();
  if ((adj_tne2_euler) || (adj_tne2_ns)) {
    nVar_Adj_TNE2         = solver_container[MESH_0][ADJTNE2_SOL]->GetnVar();
    nPrimVar_Adj_TNE2     = solver_container[MESH_0][ADJTNE2_SOL]->GetnPrimVar();
    nPrimVarGrad_Adj_TNE2 = solver_container[MESH_0][ADJTNE2_SOL]->GetnPrimVarGrad();
  }
  
  /*--- Number of variables for the linear problem ---*/
  
  if (lin_euler)  nVar_Lin_Flow = solver_container[MESH_0][LINFLOW_SOL]->GetnVar();
  
  /*--- Number of dimensions ---*/
  
  nDim = geometry[MESH_0]->GetnDim();
  
  /*--- Definition of the Class for the numerical method: numerics_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
  
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    numerics_container[iMGlevel] = new CNumerics** [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      numerics_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
  }
  
  /*--- Solver definition for the template problem ---*/
  if (template_solver) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][TEMPLATE_SOL][CONV_TERM] = new CConvective_Template(nDim, nVar_Template, config);
        break;
      default : cout << "Convective scheme not implemented (template_solver)." << endl; exit(EXIT_FAILURE); break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics_container[iMGlevel][TEMPLATE_SOL][VISC_TERM] = new CViscous_Template(nDim, nVar_Template, config);
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics_container[iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM] = new CSource_Template(nDim, nVar_Template, config);
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM] = new CConvective_Template(nDim, nVar_Template, config);
    }
    
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(EXIT_FAILURE);
        break;
        
      case SPACE_CENTERED :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config); break;
            case JST_KE : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_KE_Flow(nDim, nVar_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
          if (!config->GetLowFidelitySim()) {
            for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          }
          else {
            numerics_container[MESH_1][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config);
            for (iMGlevel = 2; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          }
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJSTArtComp_Flow(nDim, nVar_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
          
        }
        if (freesurface) {
          /*--- FreeSurface flow, use artificial compressibility method ---*/
          cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE);
        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              if (ideal_gas) {
                
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                  numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                }
              } else {
                
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                  numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                }
              }
              break;
              
            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case TURKEL:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case HLLC:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
        }
        if (freesurface) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwArtComp_FreeSurf_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_FreeSurf_Flow(nDim, nVar_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
        }
        
        break;
        
      default :
        cout << "Convective scheme not implemented (euler and ns)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible) {
      if (ideal_gas) {
        
        /*--- Compressible flow Ideal gas ---*/
        numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        
      } else{
        
        /*--- Compressible flow Realgas ---*/
        numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CGeneralAvgGradCorrected_Flow(nDim, nVar_Flow, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, config);
        
      }
    }
    if (incompressible) {
      /*--- Incompressible flow, use artificial compressibility method ---*/
      numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
      
      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
    }
    if (freesurface) {
      /*--- Freesurface flow, use artificial compressibility method ---*/
      numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
      
      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      if (config->GetRotating_Frame() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
      else if (config->GetAxisymmetric() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_Flow(nDim, nVar_Flow, config);
      else if (config->GetGravityForce() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceGravity(nDim, nVar_Flow, config);
      else if (config->GetWind_Gust() == YES)
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceWindGust(nDim, nVar_Flow, config);
      else
        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
      
      numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
    }
    
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((tne2_euler) || (tne2_ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_TNE2()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(EXIT_FAILURE);
        break;
        
      case SPACE_CENTERED :
        /*--- Compressible two-temperature flow ---*/
        switch (config->GetKind_Centered_TNE2()) {
          case NO_CENTERED : cout << "No centered scheme." << endl; break;
          case LAX :
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM]       = new CCentLax_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_TNE2(nDim, nVar_TNE2,  nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
            }
            break;
          default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
        }
        break;
        
      case SPACE_UPWIND :
        /*--- Compressible two-temperature flow ---*/
        switch (config->GetKind_Upwind_TNE2()) {
          case NO_UPWIND : cout << "No upwind scheme." << endl; break;
          case ROE:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwRoe_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_TNE2(nDim, nVar_TNE2,  nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
            }
            break;
            
          case MSW:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwMSW_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwMSW_TNE2(nDim, nVar_TNE2,  nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
            }
            break;
            
          case AUSM:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwAUSM_TNE2(nDim, nVar_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwAUSM_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          case AUSMPWPLUS:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwAUSMPWplus_TNE2(nDim, nVar_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwAUSMPWplus_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          case TURKEL:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwRoe_Turkel_TNE2(nDim, nVar_TNE2, config);
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_Turkel_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          case HLLC:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwHLLC_TNE2(nDim, nVar_TNE2, config);
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwHLLC_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
        }
        break;
        
      default :
        cout << "Convective scheme not implemented (TNE2)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    /*--- Compressible TNE2 ---*/
    numerics_container[MESH_0][TNE2_SOL][VISC_TERM] = new CAvgGradCorrected_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
    for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TNE2_SOL][VISC_TERM] = new CAvgGradCorrected_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
    }
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TNE2_SOL][VISC_BOUND_TERM] = new CAvgGradCorrected_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TNE2_SOL][SOURCE_FIRST_TERM] = new CSource_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
    }
    
  }
  
  /*--- Solver definition for the turbulent model problem ---*/
  
  if (turbulent) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          else if (neg_spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          else if (machine_learning) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbML(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        }
        break;
      default :
        cout << "Convective scheme not implemented (turbulent)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSA(nDim, nVar_Turb, config);
      else if (neg_spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSA_Neg(nDim, nVar_Turb, config);
      else if (machine_learning) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbML(nDim, nVar_Turb, config);
      else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSST(nDim, nVar_Turb, constants, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
      else if (neg_spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA_Neg(nDim, nVar_Turb, config);
      else if (machine_learning) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbML(nDim, nVar_Turb, config);
      else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, constants, config);
      numerics_container[iMGlevel][TURB_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Turb, config);
    }
    
    /*--- Definition of the boundary condition method ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
      }
      else if (neg_spalart_allmaras) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA_Neg(nDim, nVar_Turb, config);
      }
      else if (machine_learning) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbML(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbML(nDim, nVar_Turb, config);
      }
      else if (menter_sst) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, config);
      }
    }
  }
  
  /*--- Solver definition for the transition model problem ---*/
  if (transition) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][TRANS_SOL][CONV_TERM] = new CUpwSca_TransLM(nDim, nVar_Trans, config);
        }
        break;
      default :
        cout << "Convective scheme not implemented (transition)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TRANS_SOL][VISC_TERM] = new CAvgGradCorrected_TransLM(nDim, nVar_Trans, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TransLM(nDim, nVar_Trans, config);
      numerics_container[iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Trans, config);
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TRANS_SOL][CONV_BOUND_TERM] = new CUpwLin_TransLM(nDim, nVar_Trans, config);
    }
  }
  
  /*--- Solver definition for the poisson potential problem ---*/
  if (poisson) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][POISSON_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Poisson, config);
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][POISSON_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Poisson, config);
    numerics_container[MESH_0][POISSON_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Poisson, config);
    
  }
  
  /*--- Solver definition for the poisson potential problem ---*/
  if (heat) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][HEAT_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Heat, config);
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][HEAT_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Heat, config);
    numerics_container[MESH_0][HEAT_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Heat, config);
    
  }
  
  /*--- Solver definition for the flow adjoint problem ---*/
  
  if (adj_euler || adj_ns) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    
    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(EXIT_FAILURE);
        break;
        
      case SPACE_CENTERED :
        
        if (compressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        if (incompressible || freesurface) {
          
          /*--- Incompressible flow, use artificial compressibility method ---*/
          
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJSTArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
          
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        break;
        
      case SPACE_UPWIND :
        
        if (compressible) {
          
          /*--- Compressible flow ---*/
          
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
        }
        
        if (incompressible || freesurface) {
          
          /*--- Incompressible flow, use artificial compressibility method ---*/
          
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
          }
        }
        
        break;
        
      default :
        cout << "Convective scheme not implemented (adj_euler and adj_ns)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    
    if (compressible) {
      
      /*--- Compressible flow ---*/
      
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      }
      
    }
    
    if (incompressible || freesurface) {
      
      /*--- Incompressible flow, use artificial compressibility method ---*/
      
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
      numerics_container[MESH_0][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
      
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
        numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
      }
      
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      
      /*--- Note that RANS is incompatible with Axisymmetric or Rotational (Fix it!) ---*/
      
      if (compressible) {
        
        if (adj_ns) {
          
          numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceViscous_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          if (config->GetRotating_Frame() == YES)
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        
        else {
          
          if (config->GetRotating_Frame() == YES)
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else if (config->GetAxisymmetric() == YES)
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
          
          numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
          
        }
        
      }
      
      if (incompressible || freesurface) {
        
        numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
        numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Adj_Flow, config);
        
      }
      
    }
    
  }
  
  /*--- Solver definition for the flow adjoint problem ---*/
  if (adj_tne2_euler || adj_tne2_ns) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTNE2()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(EXIT_FAILURE);
        break;
      case SPACE_CENTERED :
        switch (config->GetKind_Centered_AdjTNE2()) {
          case NO_CENTERED : cout << "No centered scheme." << endl; break;
          case LAX : numerics_container[MESH_0][ADJTNE2_SOL][CONV_TERM] = new CCentLax_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config); break;
          case JST : numerics_container[MESH_0][ADJTNE2_SOL][CONV_TERM] = new CCentJST_AdjTNE2(nDim, nVar_Adj_TNE2, config); break;
          default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
        }
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][ADJTNE2_SOL][CONV_TERM] = new CCentLax_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][ADJTNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
        break;
      case SPACE_UPWIND :
        switch (config->GetKind_Upwind_AdjTNE2()) {
          case NO_UPWIND : cout << "No upwind scheme." << endl; break;
          case ROE:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_TERM] = new CUpwRoe_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
            }
            break;
          case SW:
            for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_TERM] = new CUpwSW_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_BOUND_TERM] = new CUpwSW_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
            }
            break;
          default : cout << "Upwind scheme not implemented." << endl; exit(EXIT_FAILURE); break;
        }
        break;
        
      default :
        cout << "Convective scheme not implemented (adj_tne2_euler and adj_tne2_ns)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][ADJTNE2_SOL][VISC_TERM] = new CAvgGrad_AdjTNE2(nDim, nVar_Adj_TNE2, config);
      numerics_container[iMGlevel][ADJTNE2_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjTNE2(nDim, nVar_Adj_TNE2, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][ADJTNE2_SOL][SOURCE_FIRST_TERM] = new CSource_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
      
      numerics_container[iMGlevel][ADJTNE2_SOL][SOURCE_SECOND_TERM] = new CSource_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
    }
    
  }
  
  /*--- Solver definition for the linearized flow problem ---*/
  if (lin_euler) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_LinFlow()) {
      case NONE :
        break;
      case SPACE_CENTERED :
        switch (config->GetKind_Centered_LinFlow()) {
          case LAX : numerics_container[MESH_0][LINFLOW_SOL][CONV_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config); break;
          case JST : numerics_container[MESH_0][LINFLOW_SOL][CONV_TERM] = new CCentJST_LinFlow(nDim, nVar_Lin_Flow, config); break;
          default : cout << "Centered scheme not implemented." << endl; exit(EXIT_FAILURE); break;
        }
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][LINFLOW_SOL][CONV_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
        break;
      default :
        cout << "Convective scheme not implemented (lin_euler)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics_container[iMGlevel][LINFLOW_SOL][CONV_BOUND_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
  }
  
  /*--- Solver definition for the turbulent adjoint problem ---*/
  if (adj_turb) {
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTurb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          if (spalart_allmaras) {
            numerics_container[iMGlevel][ADJTURB_SOL][CONV_TERM] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
          }
          else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
          else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
        break;
      default :
        cout << "Convective scheme not implemented (adj_turb)." << endl; exit(EXIT_FAILURE);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[iMGlevel][ADJTURB_SOL][VISC_TERM] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);
      }
      else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
      else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
        numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
      }
      else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
      else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);
      else if (neg_spalart_allmaras) {cout << "Adjoint Neg SA turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
      else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(EXIT_FAILURE);}
    }
    
  }
  
  /*--- Solver definition for the wave problem ---*/
  if (wave) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][WAVE_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Wave, config);
    
  }
  
  /*--- Solver definition for the FEA problem ---*/
  if (fea) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics_container[MESH_0][FEA_SOL][VISC_TERM] = new CGalerkin_FEA(nDim, nVar_FEA, config);
    
  }
  
  /*--- Solver definition for the FEM problem ---*/
  if (fem) {
	switch (config->GetGeometricConditions()) {
    	case SMALL_DEFORMATIONS :
    		switch (config->GetMaterialModel()) {
    			case LINEAR_ELASTIC: numerics_container[MESH_0][FEA_SOL][VISC_TERM] = new CFEM_LinearElasticity(nDim, nVar_FEM, config); break;
    			case NEO_HOOKEAN : cout << "Material model does not correspond to geometric conditions." << endl; exit(EXIT_FAILURE); break;
    			default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
    		}
    		break;
    	case LARGE_DEFORMATIONS :
    		switch (config->GetMaterialModel()) {
				case LINEAR_ELASTIC: cout << "Material model does not correspond to geometric conditions." << endl; exit(EXIT_FAILURE); break;
    			case NEO_HOOKEAN :
    				switch (config->GetMaterialCompressibility()) {
    					case COMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][VISC_TERM] = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config); break;
    					case INCOMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][VISC_TERM] = new CFEM_NeoHookean_Incomp(nDim, nVar_FEM, config); break;
    					default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
    				}
    				break;
    			default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
    		}
    		break;
    	default: cout << " Solver not implemented." << endl; exit(EXIT_FAILURE); break;
	}

  }

}

void CDriver::Iteration_Preprocessing(CIteration **iteration_container, CConfig **config, unsigned short iZone) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Initial print to console for this zone. ---*/
  
  if (rank == MASTER_NODE) cout << "Zone " << iZone+1;
  
  /*--- Loop over all zones and instantiate the physics iteration. ---*/

  switch (config[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
      if (rank == MASTER_NODE)
        cout << ": Euler/Navier-Stokes/RANS flow iteration." << endl;
      iteration_container[iZone] = new CMeanFlowIteration(config[iZone]);
      break;

    case TNE2_EULER: case TNE2_NAVIER_STOKES:
      if (rank == MASTER_NODE)
        cout << ": TNE2 iteration." << endl;
      iteration_container[iZone] = new CTNE2Iteration(config[iZone]);
      break;

    case WAVE_EQUATION:
      if (rank == MASTER_NODE)
        cout << ": wave iteration." << endl;
      iteration_container[iZone] = new CWaveIteration(config[iZone]);
      break;

    case HEAT_EQUATION:
      if (rank == MASTER_NODE)
        cout << ": heat iteration." << endl;
      iteration_container[iZone] = new CHeatIteration(config[iZone]);
      break;

    case POISSON_EQUATION:
      if (rank == MASTER_NODE)
        cout << ": poisson iteration." << endl;
      iteration_container[iZone] = new CPoissonIteration(config[iZone]);
      break;

    case LINEAR_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << ": FEA iteration." << endl;
      iteration_container[iZone] = new CFEAIteration(config[iZone]);
      break;

    case FEM_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << ": FEA iteration." << endl;
      iteration_container[iZone] = new CFEM_StructuralAnalysis(config[iZone]);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << ": adjoint Euler/Navier-Stokes/RANS flow iteration." << endl;
      iteration_container[iZone] = new CAdjMeanFlowIteration(config[iZone]);
      break;

    case ADJ_TNE2_EULER: case ADJ_TNE2_NAVIER_STOKES:
      if (rank == MASTER_NODE)
        cout << ": adjoint TNE2 iteration." << endl;
      iteration_container[iZone] = new CAdjTNE2Iteration(config[iZone]);
      break;

    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << ": discrete adjoint Euler/Navier-Stokes/RANS flow iteration." << endl;
      iteration_container[iZone] = new CDiscAdjMeanFlowIteration(config[iZone]);
      break;
  }

}


void CDriver::Interface_Preprocessing(CTransfer ***transfer_container, CInterpolator ***interpolator_container,
							 CGeometry ***geometry_container, CConfig **config_container,
							 CSolver ****solver_container, unsigned short nZone, unsigned short nDim) {

	int rank = MASTER_NODE;
	unsigned short donorZone, targetZone;
	unsigned short nVar, nVarTransfer;

	/*--- Initialize some useful booleans ---*/
	bool fluid_donor, structural_donor;
	bool fluid_target, structural_target;

	bool matching_mesh;

	fluid_donor  = false;  structural_donor  = false;
	fluid_target  = false;  structural_target  = false;


#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/*--- Coupling between zones (limited to two zones at the moment) ---*/
	for (donorZone = 0; donorZone < nZone; donorZone++){

		/*--- Initialize donor booleans ---*/
		fluid_donor  = false;  structural_donor  = false;
		matching_mesh = config_container[donorZone]->GetMatchingMesh();

		/*--- Set the donor boolean: as of now, only Fluid-Structure Interaction considered ---*/
		switch (config_container[donorZone]->GetKind_Solver()) {
			case EULER : case NAVIER_STOKES: case RANS: fluid_donor  = true; 		break;
			case FEM_ELASTICITY: 						structural_donor = true; 	break;
		}

		for (targetZone = 0; targetZone < nZone; targetZone++){

			/*--- Initialize donor booleans ---*/
			fluid_target  = false;  structural_target  = false;

			/*--- Set the target boolean: as of now, only Fluid-Structure Interaction considered ---*/
			switch (config_container[targetZone]->GetKind_Solver()) {
				case EULER : case NAVIER_STOKES: case RANS: fluid_target  = true; 		break;
				case FEM_ELASTICITY: 						structural_target = true; 	break;
			}

			/*--- Retrieve the number of conservative variables (for problems not involving structural analysis ---*/
			if (!structural_donor && !structural_target){
				nVar = solver_container[donorZone][MESH_0][FLOW_SOL]->GetnVar();
			}
			else{
				/*--- If at least one of the components is structural ---*/
				nVar = nDim;
			}

			/*--- Interface conditions are only defined between different zones ---*/
			if (donorZone != targetZone){

				if (rank == MASTER_NODE) cout << "From zone " << donorZone << " to zone " << targetZone << ": " << endl;

				/*--- Match Zones ---*/
				if (rank == MASTER_NODE) cout << "Setting coupling ";

				/*--- If the mesh is matching: match points ---*/
				if (matching_mesh){
					if (rank == MASTER_NODE) cout << "between matching meshes. " << endl;
					geometry_container[donorZone][MESH_0]->MatchZone(config_container[donorZone], geometry_container[targetZone][MESH_0],
							config_container[targetZone], donorZone, nZone);
				}
				/*--- Else: interpolate ---*/
				else {
					if (rank == MASTER_NODE) cout << "between non-matching meshes ";
					switch (config_container[donorZone]->GetKindInterpolation()){
						case NEAREST_NEIGHBOR:
							interpolator_container[donorZone][targetZone] = new CNearestNeighbor(geometry_container, config_container, donorZone, targetZone);
							if (rank == MASTER_NODE) cout << "using a nearest-neighbor approach." << endl;
							break;
						case ISOPARAMETRIC:
							interpolator_container[donorZone][targetZone] = new CIsoparametric(geometry_container, config_container, donorZone, targetZone);
							if (rank == MASTER_NODE) cout << "using an isoparametric approach." << endl;
							break;
					}
				}

				/*--- Initialize the appropriate transfer strategy ---*/
				if (rank == MASTER_NODE) cout << "Transferring ";

				if (fluid_donor && structural_target) {
					nVarTransfer = 2;
					transfer_container[donorZone][targetZone] = new CTransfer_FlowTraction(nVar, nVarTransfer, config_container[donorZone]);
					if (rank == MASTER_NODE) cout << "flow tractions. "<< endl;
				}
				else if (structural_donor && fluid_target){
					nVarTransfer = 0;
					transfer_container[donorZone][targetZone] = new CTransfer_StructuralDisplacements(nVar, nVarTransfer, config_container[donorZone]);
					if (rank == MASTER_NODE) cout << "structural displacements. "<< endl;
				}
				else {
					nVarTransfer = 0;
					transfer_container[donorZone][targetZone] = new CTransfer_ConservativeVars(nVar, nVarTransfer, config_container[donorZone]);
					if (rank == MASTER_NODE) cout << "generic conservative variables. " << endl;
				}

			}


		}

	}


}


CDriver::~CDriver(void) { }


CSingleZoneDriver::CSingleZoneDriver(CIteration **iteration_container,
                                     CSolver ****solver_container,
                                     CGeometry ***geometry_container,
                                     CIntegration ***integration_container,
                                     CNumerics *****numerics_container,
                                     CInterpolator ***interpolator_container,
                                     CTransfer ***transfer_container,
                                     CConfig **config_container,
                                     unsigned short val_nZone,
                                     unsigned short val_nDim) : CDriver(iteration_container,
                                                                         solver_container,
                                                                         geometry_container,
                                                                         integration_container,
                                                                         numerics_container,
                                                                         interpolator_container,
                                                                         transfer_container,
                                                                         config_container,
                                                                         val_nZone,
                                                                         val_nDim) { }

CSingleZoneDriver::~CSingleZoneDriver(void) { }

void CSingleZoneDriver::Run(CIteration **iteration_container,
                            COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            CInterpolator ***interpolator_container,
                            CTransfer ***transfer_container) {

  unsigned short iZone = ZONE_0;
  
  /*--- Run an iteration of the physics within this single zone.
   We assume that the zone of interest is in the ZONE_0 container position. ---*/
  
  iteration_container[ZONE_0]->Preprocess(); /*--- Does nothing for now. ---*/
  
  bool checkSubiter = true;
  
  if (!checkSubiter){
	  iteration_container[ZONE_0]->Iterate(output, integration_container, geometry_container,
	                                       solver_container, numerics_container, config_container,
	                                       surface_movement, grid_movement, FFDBox);
  }
  else {

	  cout << "YES! I'M SUBITERATING!!!!" << endl;
	  iteration_container[ZONE_0]->Subiterate(output, integration_container, geometry_container,
	                                       solver_container, numerics_container, config_container,
	                                       surface_movement, grid_movement, FFDBox, ZONE_0);

	  iteration_container[ZONE_0]->Update(output, integration_container, geometry_container,
	          solver_container, numerics_container, config_container,
	          surface_movement, grid_movement, FFDBox, ZONE_0);
  }

  iteration_container[ZONE_0]->Monitor(); /*--- Does nothing for now. ---*/
  
  iteration_container[ZONE_0]->Output(); /*--- Does nothing for now. ---*/
  
  iteration_container[ZONE_0]->Postprocess(); /*--- Does nothing for now. ---*/
  
}


CMultiZoneDriver::CMultiZoneDriver(CIteration **iteration_container,
                                   CSolver ****solver_container,
                                   CGeometry ***geometry_container,
                                   CIntegration ***integration_container,
                                   CNumerics *****numerics_container,
                                   CInterpolator ***interpolator_container,
                                   CTransfer ***transfer_container,
                                   CConfig **config_container,
                                   unsigned short val_nZone,
                                   unsigned short val_nDim) : CDriver(iteration_container,
                                                                       solver_container,
                                                                       geometry_container,
                                                                       integration_container,
                                                                       numerics_container,
                                                                       interpolator_container,
                                                                       transfer_container,
                                                                       config_container,
                                                                       val_nZone,
                                                                       val_nDim) { }


CMultiZoneDriver::~CMultiZoneDriver(void) { }

void CMultiZoneDriver::Run(CIteration **iteration_container,
                           COutput *output,
                           CIntegration ***integration_container,
                           CGeometry ***geometry_container,
                           CSolver ****solver_container,
                           CNumerics *****numerics_container,
                           CConfig **config_container,
                           CSurfaceMovement **surface_movement,
                           CVolumetricMovement **grid_movement,
                           CFreeFormDefBox*** FFDBox,
                           CInterpolator ***interpolator_container,
                           CTransfer ***transfer_container) {
  
  unsigned short iZone;
  
  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++) {
    
    iteration_container[iZone]->Preprocess();  /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox);
    
    iteration_container[iZone]->Update(output, integration_container, geometry_container,
            solver_container, numerics_container, config_container,
            surface_movement, grid_movement, FFDBox, iZone);      /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Monitor();     /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Output();      /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Postprocess(); /*--- Does nothing for now. ---*/
    
  }
  
}


CFSIDriver::CFSIDriver(CIteration **iteration_container,
                       CSolver ****solver_container,
                       CGeometry ***geometry_container,
                       CIntegration ***integration_container,
                       CNumerics *****numerics_container,
                       CInterpolator ***interpolator_container,
                       CTransfer ***transfer_container,
                       CConfig **config_container,
                       unsigned short val_nZone,
                       unsigned short val_nDim) : CDriver(iteration_container,
                                                           solver_container,
                                                           geometry_container,
                                                           integration_container,
                                                           numerics_container,
                                                           interpolator_container,
                                                           transfer_container,
                                                           config_container,
                                                           val_nZone,
                                                           val_nDim) { }

CFSIDriver::~CFSIDriver(void) { }

void CFSIDriver::Run(CIteration **iteration_container,
                     COutput *output,
                     CIntegration ***integration_container,
                     CGeometry ***geometry_container,
                     CSolver ****solver_container,
                     CNumerics *****numerics_container,
                     CConfig **config_container,
                     CSurfaceMovement **surface_movement,
                     CVolumetricMovement **grid_movement,
                     CFreeFormDefBox*** FFDBox,
                     CInterpolator ***interpolator_container,
                     CTransfer ***transfer_container) {

	/*--- As of now, we are coding it for just 2 zones. ---*/
	/*--- This will become more general, but we need to modify the configuration for that ---*/
	unsigned short ZONE_FLOW = 0, ZONE_STRUCT = 1;
	unsigned short iZone;

	unsigned long IntIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[IntIter]->SetIntIter(IntIter);
	unsigned long iFSIIter = 0;
	unsigned long nFSIIter = config_container[ZONE_FLOW]->GetnIterFSI();
	unsigned long ExtIter = config_container[ZONE_FLOW]->GetExtIter();

	bool fem_solver = (config_container[ZONE_1]->GetKind_Solver() == FEM_ELASTICITY);

	int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/*-----------------------------------------------------------------*/
	/*---------------- Predict structural displacements ---------------*/
	/*-----------------------------------------------------------------*/

//	FSI_Disp_Predictor(output, integration_container, geometry_container,
//                          solver_container, numerics_container, config_container,
//                          surface_movement, grid_movement, FFDBox);


	Predict_Displacements(output, integration_container, geometry_container,
            		      solver_container, numerics_container, config_container,
            		      surface_movement, grid_movement, FFDBox,
            		      ZONE_STRUCT, ZONE_FLOW);

	while (iFSIIter < nFSIIter){

		/*-----------------------------------------------------------------*/
		/*------------------------ Update mesh ----------------------------*/
		/*-----------------------------------------------------------------*/

//        FSI_Disp_Transfer(output, integration_container, geometry_container,
//                          solver_container, numerics_container, config_container,
//                          surface_movement, grid_movement, FFDBox, transfer_container);

		Transfer_Displacements(output, integration_container, geometry_container,
                solver_container, numerics_container, config_container,
                surface_movement, grid_movement, FFDBox, transfer_container,
                ZONE_STRUCT, ZONE_FLOW);

		/*-----------------------------------------------------------------*/
		/*-------------------- Fluid subiteration -------------------------*/
		/*-----------------------------------------------------------------*/

//        Flow_Subiteration(output, integration_container, geometry_container,
//                          solver_container, numerics_container, config_container,
//                          surface_movement, grid_movement, FFDBox);

		iteration_container[ZONE_FLOW]->Subiterate(output, integration_container, geometry_container,
		                                       solver_container, numerics_container, config_container,
		                                       surface_movement, grid_movement, FFDBox, ZONE_FLOW);

		/*--- Write the convergence history for the fluid (only screen output) ---*/

		output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_FLOW);

		/*-----------------------------------------------------------------*/
		/*------------------- Set FEA loads from fluid --------------------*/
		/*-----------------------------------------------------------------*/

//        FSI_Load_Transfer(output, integration_container, geometry_container,
//	                 	 solver_container, numerics_container, config_container,
//	                 	 surface_movement, grid_movement, FFDBox, transfer_container, ExtIter);

		Transfer_Tractions(output, integration_container, geometry_container,
                solver_container, numerics_container, config_container,
                surface_movement, grid_movement, FFDBox, transfer_container,
                ZONE_FLOW, ZONE_STRUCT);

		/*-----------------------------------------------------------------*/
		/*------------------ Structural subiteration ----------------------*/
		/*-----------------------------------------------------------------*/

//        if (fem_solver){
//    	    FEM_Subiteration(output, integration_container, geometry_container,
//    	                 	 solver_container, numerics_container, config_container,
//    	                 	 surface_movement, grid_movement, FFDBox);
//        }
//        else{
//    	    FEA_Subiteration(output, integration_container, geometry_container,
//    	                 	 solver_container, numerics_container, config_container,
//    	                 	 surface_movement, grid_movement, FFDBox);
//        }

		iteration_container[ZONE_STRUCT]->Subiterate(output, integration_container, geometry_container,
		                                       solver_container, numerics_container, config_container,
		                                       surface_movement, grid_movement, FFDBox, ZONE_STRUCT);

		/*--- Write the convergence history for the structure (only screen output) ---*/

		output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUCT);

		/*-----------------------------------------------------------------*/
		/*----------------- Displacements relaxation ----------------------*/
		/*-----------------------------------------------------------------*/

//	    FSI_Disp_Relaxation(output, geometry_container, solver_container, config_container, iFSIIter);

		Relaxation_Displacements(output, geometry_container, solver_container, config_container,
								 ZONE_STRUCT, ZONE_FLOW, iFSIIter);

		/*-----------------------------------------------------------------*/
		/*-------------------- Check convergence --------------------------*/
		/*-----------------------------------------------------------------*/

		integration_container[ZONE_STRUCT][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT],
																solver_container[ZONE_STRUCT][MESH_0][FEA_SOL], iFSIIter);

		if (integration_container[ZONE_STRUCT][FEA_SOL]->GetConvergence_FSI()) break;

		/*-----------------------------------------------------------------*/
		/*--------------------- Update iFSIIter ---------------------------*/
		/*-----------------------------------------------------------------*/

		iFSIIter++;

	}

	/*-----------------------------------------------------------------*/
  	/*-------------------- Update fluid solver ------------------------*/
	/*-----------------------------------------------------------------*/

//    Flow_Update(output, integration_container, geometry_container,
//                solver_container, numerics_container, config_container,
//                surface_movement, grid_movement, FFDBox, ExtIter);

	iteration_container[ZONE_FLOW]->Update(output, integration_container, geometry_container,
	          	  	  	  	  	  	  	  solver_container, numerics_container, config_container,
	          	  	  	  	  	  	  	  surface_movement, grid_movement, FFDBox, ZONE_FLOW);

	/*-----------------------------------------------------------------*/
  	/*----------------- Update structural solver ----------------------*/
	/*-----------------------------------------------------------------*/

//    if (fem_solver){
//    	FEM_Update(output, integration_container, geometry_container,
//                	solver_container, numerics_container, config_container, ExtIter);
//    }
//    else{
//    	FEA_Update(output, integration_container, geometry_container,
//                	solver_container, config_container, ExtIter);
//    }

	iteration_container[ZONE_STRUCT]->Update(output, integration_container, geometry_container,
	          	  	  	  	  	  	  	  solver_container, numerics_container, config_container,
	          	  	  	  	  	  	  	  surface_movement, grid_movement, FFDBox, ZONE_STRUCT);


	/*-----------------------------------------------------------------*/
	/*--------------- Update convergence parameter --------------------*/
	/*-----------------------------------------------------------------*/
	integration_container[ZONE_STRUCT][FEA_SOL]->SetConvergence_FSI(false);

  
}

void CFSIDriver::Predict_Displacements(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 unsigned short donorZone, unsigned short targetZone){

#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	solver_container[donorZone][MESH_0][FEA_SOL]->PredictStruct_Displacement(geometry_container[donorZone], config_container[donorZone],
			solver_container[donorZone]);

	/*--- For parallel simulations we need to communicate the predicted solution before updating the fluid mesh ---*/

	solver_container[donorZone][MESH_0][FEA_SOL]->Set_MPI_Solution_Pred(geometry_container[donorZone][MESH_0], config_container[donorZone]);


}

void CFSIDriver::Predict_Tractions(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 unsigned short donorZone, unsigned short targetZone){

}

void CFSIDriver::Transfer_Displacements(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short donorZone, unsigned short targetZone){

#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	bool MatchingMesh = config_container[targetZone]->GetMatchingMesh();

	/*--- Select the transfer method and the appropriate mesh properties (matching or nonmatching mesh) ---*/

	switch (config_container[targetZone]->GetKind_TransferMethod()) {
	case BROADCAST_DATA:
		if (MatchingMesh){
			transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
			/*--- Set the volume deformation for the fluid zone ---*/
			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);

		}
		else {
			transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
			/*--- Set the volume deformation for the fluid zone ---*/
			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);

		}
		break;
	case SCATTER_DATA:
		if (MatchingMesh){
			transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
			/*--- Set the volume deformation for the fluid zone ---*/
			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
		}
		else {
			cout << "Scatter method not implemented for non-matching meshes. Exiting..." << endl;
			exit(EXIT_FAILURE);
		}
		break;
	case ALLGATHER_DATA:
		if (MatchingMesh){
			cout << "Allgather method not yet implemented for matching meshes. Exiting..." << endl;
			exit(EXIT_FAILURE);
		}
		else {
			transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
			/*--- Set the volume deformation for the fluid zone ---*/
			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
		}
		break;
	case LEGACY_METHOD:
		if (MatchingMesh){
			solver_container[targetZone][MESH_0][FLOW_SOL]->SetFlow_Displacement(geometry_container[targetZone], grid_movement[targetZone],
					config_container[targetZone], config_container[donorZone],
					geometry_container[donorZone], solver_container[donorZone]);
		}
		else {
			solver_container[targetZone][MESH_0][FLOW_SOL]->SetFlow_Displacement_Int(geometry_container[targetZone], grid_movement[targetZone],
					config_container[targetZone], config_container[donorZone],
					geometry_container[donorZone], solver_container[donorZone]);
		}
		break;
	}

}

void CFSIDriver::Transfer_Tractions(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
		     CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short donorZone, unsigned short targetZone){

#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	bool MatchingMesh = config_container[donorZone]->GetMatchingMesh();

	/*--- Load transfer --  This will have to be modified for non-matching meshes ---*/

	unsigned short SolContainer_Position_fea = config_container[targetZone]->GetContainerPosition(RUNTIME_FEA_SYS);

	/*--- FEA equations -- Necessary as the SetFEA_Load routine is as of now contained in the structural solver ---*/
	unsigned long ExtIter = config_container[targetZone]->GetExtIter();
	config_container[targetZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

	/*--- Select the transfer method and the appropriate mesh properties (matching or nonmatching mesh) ---*/

	switch (config_container[donorZone]->GetKind_TransferMethod()) {
	case BROADCAST_DATA:
		if (MatchingMesh){
			transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Matching(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
		}
		else {
			transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
		}
		break;
	case SCATTER_DATA:
		if (MatchingMesh){
			transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
		}
		else {
			cout << "Scatter method not implemented for non-matching meshes. Exiting..." << endl;
			exit(EXIT_FAILURE);
		}
		break;
	case ALLGATHER_DATA:
		if (MatchingMesh){
			cout << "Allgather method not yet implemented for matching meshes. Exiting..." << endl;
			exit(EXIT_FAILURE);
		}
		else {
			transfer_container[donorZone][targetZone]->Allgather_InterfaceData(solver_container[donorZone][MESH_0][FLOW_SOL],solver_container[targetZone][MESH_0][FEA_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
		}
		break;
	case LEGACY_METHOD:
		if (MatchingMesh){
			solver_container[targetZone][MESH_0][FEA_SOL]->SetFEA_Load(solver_container[donorZone], geometry_container[targetZone], geometry_container[donorZone],
					config_container[targetZone], config_container[donorZone], numerics_container[targetZone][MESH_0][SolContainer_Position_fea][VISC_TERM]);
		}
		else {
			solver_container[targetZone][MESH_0][FEA_SOL]->SetFEA_Load_Int(solver_container[donorZone], geometry_container[targetZone], geometry_container[donorZone],
					config_container[targetZone], config_container[donorZone], numerics_container[targetZone][MESH_0][SolContainer_Position_fea][VISC_TERM]);
		}
		break;
	}

}

void CFSIDriver::Relaxation_Displacements(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter){

#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/*-------------------- Aitken's relaxation ------------------------*/

	/*------------------- Compute the coefficient ---------------------*/

	solver_container[donorZone][MESH_0][FEA_SOL]->ComputeAitken_Coefficient(geometry_container[donorZone], config_container[donorZone],
			solver_container[donorZone], iFSIIter);

	/*----------------- Set the relaxation parameter ------------------*/

	solver_container[donorZone][MESH_0][FEA_SOL]->SetAitken_Relaxation(geometry_container[donorZone], config_container[donorZone],
			solver_container[donorZone]);


	/*----------------- Communicate the predicted solution and the old one ------------------*/
	solver_container[donorZone][MESH_0][FEA_SOL]->Set_MPI_Solution_Pred_Old(geometry_container[donorZone][MESH_0], config_container[donorZone]);


}

void CFSIDriver::Relaxation_Tractions(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long iFSIIter){

}


