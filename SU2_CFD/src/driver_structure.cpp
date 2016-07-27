/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 4.2.0 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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
  unsigned short nDim;
  
  bool fsi = config_container[ZONE_0]->GetFSI_Simulation();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Store the number of zones for the problem (class data member). ---*/
  
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
	 *--- This only works for a multizone FSI problem (nZone > 1).
	 *--- Also, at the moment this capability is limited to two zones (nZone < 3).
	 *--- This will change in the future. ---*/

	if ((rank == MASTER_NODE) && (fsi))
		cout << endl <<"------------------- Multizone Interface Preprocessing -------------------" << endl;


	if (((nZone > 1) && (nZone < 3)) && (fsi)) {

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

void CDriver::Postprocessing(CIteration **iteration_container,
                 CSolver ****solver_container,
                 CGeometry ***geometry_container,
                 CIntegration ***integration_container,
                 CNumerics *****numerics_container,
                 CInterpolator ***interpolator_container,
                 CTransfer ***transfer_container,
                 CConfig **config_container,
                 unsigned short val_nZone){


  unsigned short iZone;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for (iZone = 0; iZone < nZone; iZone++) {
    /*
    if (rank == MASTER_NODE)
          cout << endl <<"----------------- Numerics Postprocessing ----------------" << endl;
    Numerics_Postprocessing(numerics_container[iZone], solver_container[iZone],
                                   geometry_container[iZone], config_container[iZone]);
    */
    if (rank == MASTER_NODE)
          cout << endl <<"----------------- Integration Postprocessing ----------------" << endl;

    Integration_Postprocessing(integration_container[iZone], geometry_container[iZone],
                                      config_container[iZone]);

    if (rank == MASTER_NODE)
          cout << endl <<"----------------- Solver Postprocessing ---------------------" << endl;


    Solver_Postprocessing(solver_container[iZone], geometry_container[iZone],
        config_container[iZone]);
    // TODO: implement Iteration Postprocessing
    /*
    if (rank == MASTER_NODE){
      cout << endl <<"------------------------ Iteration Postprocessing ------------------------" << endl;
    }
    Iteration_Postprocessing(iteration_container, config_container, iZone);
     */

    delete[] numerics_container[iZone];
    delete[] integration_container[iZone];
    delete[] solver_container[iZone];
  }
  delete[] numerics_container;
  delete[] integration_container;
  delete[] solver_container;

}

void CDriver::Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry,
                                   CConfig *config) {
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  poisson, wave, heat, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj;
  
  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;
  poisson          = false;  neg_spalart_allmaras = false;
  wave             = false;	 disc_adj        = false;
  fem = false;
  heat             = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;           break;
        
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
      solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    if (ns) {
      solver_container[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
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
    if (adj_turb) {
      solver_container[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolver(geometry[iMGlevel], config, iMGlevel);
    }
    
    if (disc_adj) {
      solver_container[iMGlevel][ADJFLOW_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver_container[iMGlevel][FLOW_SOL], RUNTIME_FLOW_SYS, iMGlevel);
      if (turbulent)
        solver_container[iMGlevel][ADJTURB_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver_container[iMGlevel][TURB_SOL], RUNTIME_TURB_SYS, iMGlevel);
    }
  }
}


void CDriver::Solver_Postprocessing(CSolver ***solver_container, CGeometry **geometry,
    CConfig *config) {
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  poisson, wave, heat, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj;

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;
  poisson          = false;  neg_spalart_allmaras = false;
  wave             = false;  disc_adj        = false;
  fem = false;
  heat             = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/

  switch (config->GetKind_Solver()) {
  case TEMPLATE_SOLVER: template_solver = true; break;
  case EULER : euler = true; break;
  case NAVIER_STOKES: ns = true; break;
  case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
  case POISSON_EQUATION: poisson = true; break;
  case WAVE_EQUATION: wave = true; break;
  case HEAT_EQUATION: heat = true; break;
  case FEM_ELASTICITY: fem = true; break;
  case ADJ_EULER : euler = true; adj_euler = true; break;
  case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
  case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
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
    }

  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
and potential are incompatible, they use the same position in sol container ---*/
  for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- DeAllocate solution for a template problem ---*/
    if (template_solver) {
      delete solver_container[iMGlevel][TEMPLATE_SOL];
    }

    /*--- DeAllocate solution for adjoint problem ---*/
    if (adj_euler || adj_ns || disc_adj) {
      delete solver_container[iMGlevel][ADJFLOW_SOL];
      if ((turbulent && disc_adj) || adj_turb){
        delete solver_container[iMGlevel][ADJTURB_SOL];
      }
    }

    /*--- DeAllocate solution for direct problem ---*/
    if (euler || ns) {
      delete solver_container[iMGlevel][FLOW_SOL];
    }

    if (turbulent) {
      if (spalart_allmaras || neg_spalart_allmaras || menter_sst ) {
        delete solver_container[iMGlevel][TURB_SOL];
      }
      if (transition) {
        delete solver_container[iMGlevel][TRANS_SOL];
      }
    }
    if (poisson) {
      delete solver_container[iMGlevel][POISSON_SOL];
    }
    if (wave) {
      delete solver_container[iMGlevel][WAVE_SOL];
    }
    if (heat) {
      delete solver_container[iMGlevel][HEAT_SOL];
    }
    if (fem) {
      delete solver_container[iMGlevel][FEA_SOL];
    }


    //delete[] solver_container[iMGlevel];
  }
}

void CDriver::Integration_Preprocessing(CIntegration **integration_container,
                                        CGeometry **geometry, CConfig *config) {
  
  bool
  euler, adj_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  poisson, wave, fem, heat, template_solver, transition, disc_adj;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  poisson          = false; disc_adj         = false;
  wave             = false;
  heat             = false;
  fem = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case DISC_ADJ_EULER : euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS : ns = true; turbulent = true; disc_adj = true; break;
      
  }
  
  /*--- Allocate solution for a template problem ---*/
  if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);
  
  /*--- Allocate solution for direct problem ---*/
  if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (ns) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
  if (transition) integration_container[TRANS_SOL] = new CSingleGridIntegration(config);
  if (poisson) integration_container[POISSON_SOL] = new CSingleGridIntegration(config);
  if (wave) integration_container[WAVE_SOL] = new CSingleGridIntegration(config);
  if (heat) integration_container[HEAT_SOL] = new CSingleGridIntegration(config);
  if (fem) integration_container[FEA_SOL] = new CStructuralIntegration(config);
  
  /*--- Allocate solution for adjoint problem ---*/
  if (adj_euler) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_ns) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_turb) integration_container[ADJTURB_SOL] = new CSingleGridIntegration(config);
  
  if (disc_adj) integration_container[ADJFLOW_SOL] = new CIntegration(config);
  
}

void CDriver::Integration_Postprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config){
  bool
  euler, adj_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  poisson, wave, fem, heat, template_solver, transition, disc_adj;

  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  poisson          = false; disc_adj         = false;
  wave             = false;
  heat             = false;
  fem = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case DISC_ADJ_EULER : euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS : ns = true; turbulent = true; disc_adj = true; adj_turb=true; break;

  }

  /*--- DeAllocate solution for a template problem ---*/
  if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);

  /*--- DeAllocate solution for direct problem ---*/
  if (euler || ns) delete integration_container[FLOW_SOL];
  if (turbulent) delete integration_container[TURB_SOL];
  if (transition) delete integration_container[TRANS_SOL];
  if (poisson) delete integration_container[POISSON_SOL];
  if (wave) delete integration_container[WAVE_SOL];
  if (heat) delete integration_container[HEAT_SOL];
  if (fem) delete integration_container[FEA_SOL];

  /*--- DeAllocate solution for adjoint problem ---*/
  if (adj_euler || adj_ns || disc_adj) delete integration_container[ADJFLOW_SOL];
  if (adj_turb) delete integration_container[ADJTURB_SOL];


}

void CDriver::Numerics_Preprocessing(CNumerics ****numerics_container,
                                     CSolver ***solver_container, CGeometry **geometry,
                                     CConfig *config) {
  
  unsigned short iMGlevel, iSol, nDim,
  
  nVar_Template         = 0,
  nVar_Flow             = 0,
  nVar_Trans            = 0,
  nVar_Turb             = 0,
  nVar_Adj_Flow         = 0,
  nVar_Adj_Turb         = 0,
  nVar_Poisson          = 0,
  nVar_FEM				= 0,
  nVar_Wave             = 0,
  nVar_Heat             = 0;
  
  su2double *constants = NULL;
  
  bool
  euler, adj_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  spalart_allmaras, neg_spalart_allmaras, menter_sst,
  poisson,
  wave,
  fem,
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
  wave             = false;   heat             = false;   fem				= false;
  spalart_allmaras = false; neg_spalart_allmaras = false;	menter_sst       = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case DISC_ADJ_EULER: euler = true; break;
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES: ns = true; break;
    case RANS : case DISC_ADJ_RANS:  ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
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
  if (poisson)      nVar_Poisson = solver_container[MESH_0][POISSON_SOL]->GetnVar();
  
  if (wave)				nVar_Wave = solver_container[MESH_0][WAVE_SOL]->GetnVar();
  if (fem)				nVar_FEM = solver_container[MESH_0][FEA_SOL]->GetnVar();
  if (heat)				nVar_Heat = solver_container[MESH_0][HEAT_SOL]->GetnVar();
  
  /*--- Number of variables for adjoint problem ---*/
  
  if (adj_euler)        nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_ns)           nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_turb)         nVar_Adj_Turb = solver_container[MESH_0][ADJTURB_SOL]->GetnVar();
  
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
		if (ideal_gas) {
		      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
		        numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
		        numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
		      }
		}
		else {
		      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
		        numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
		        numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
		      }
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
      else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSST(nDim, nVar_Turb, constants, config);
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
      else if (neg_spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA_Neg(nDim, nVar_Turb, config);
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

  /*--- Solver definition for the FEM problem ---*/
  if (fem) {
	switch (config->GetGeometricConditions()) {
    	case SMALL_DEFORMATIONS :
    		switch (config->GetMaterialModel()) {
    			case LINEAR_ELASTIC: numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_LinearElasticity(nDim, nVar_FEM, config); break;
    			case NEO_HOOKEAN : cout << "Material model does not correspond to geometric conditions." << endl; exit(EXIT_FAILURE); break;
    			default: cout << "Material model not implemented." << endl; exit(EXIT_FAILURE); break;
    		}
    		break;
    	case LARGE_DEFORMATIONS :
    		switch (config->GetMaterialModel()) {
				case LINEAR_ELASTIC: cout << "Material model does not correspond to geometric conditions." << endl; exit(EXIT_FAILURE); break;
    			case NEO_HOOKEAN :
    				switch (config->GetMaterialCompressibility()) {
    					case COMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config); break;
    					case INCOMPRESSIBLE_MAT : numerics_container[MESH_0][FEA_SOL][FEA_TERM] = new CFEM_NeoHookean_Incomp(nDim, nVar_FEM, config); break;
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


void CDriver::Numerics_Postprocessing(CNumerics ****numerics_container,
                            CSolver ***solver_container, CGeometry **geometry,
                            CConfig *config) {

  unsigned short iMGlevel, iSol;


  bool
  euler, adj_euler,
  ns, adj_ns,
  turbulent, adj_turb,
  spalart_allmaras, neg_spalart_allmaras, menter_sst,
  poisson,
  wave,
  fem,
  heat,
  transition,
  template_solver;

  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);

  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;   turbulent        = false;
  poisson          = false;
  adj_euler        = false;   adj_ns           = false;   adj_turb         = false;
  wave             = false;   heat             = false;   fem        = false;
  spalart_allmaras = false; neg_spalart_allmaras = false; menter_sst       = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case DISC_ADJ_EULER: euler = true; break;
    case NAVIER_STOKES: case DISC_ADJ_NAVIER_STOKES: ns = true; break;
    case RANS : case DISC_ADJ_RANS:  ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
  }

  /*--- Assign turbulence model booleans --- */

  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:     spalart_allmaras = true;     break;
      case SA_NEG: neg_spalart_allmaras = true; break;
      case SST:    menter_sst = true;  break;

    }

  /*--- Solver definition for the template problem ---*/
  if (template_solver) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          delete numerics_container[iMGlevel][TEMPLATE_SOL][CONV_TERM];
        break;
    }

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++){
      /*--- Definition of the viscous scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TEMPLATE_SOL][VISC_TERM];
      /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM];
      /*--- Definition of the boundary condition method ---*/
      delete numerics_container[iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM];
    }

  }

  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {

      case SPACE_CENTERED :
        if (compressible) {

          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case LAX : case JST :  case JST_KE : delete numerics_container[MESH_0][FLOW_SOL][CONV_TERM]; break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];

        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Centered_Flow()) {

            case LAX : case JST : delete numerics_container[MESH_0][FLOW_SOL][CONV_TERM]; break;

          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];

        }
        break;
      case SPACE_UPWIND :

        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE: case AUSM : case TURKEL: case HLLC: case MSW:  case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
              }

              break;
          }

        }
        if (incompressible || freesurface) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_TERM];
                delete numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM];
              }
              break;
          }
        }

        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible||incompressible||freesurface) {
      /*--- Compressible flow Ideal gas ---*/
      delete numerics_container[MESH_0][FLOW_SOL][VISC_TERM];
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        delete numerics_container[iMGlevel][FLOW_SOL][VISC_TERM];

      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        delete numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM];

    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      delete numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM];
      delete numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM];
    }

  }


  /*--- Solver definition for the turbulent model problem ---*/

  if (turbulent) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/

    switch (config->GetKind_ConvNumScheme_Turb()) {
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          if (spalart_allmaras || neg_spalart_allmaras ||menter_sst)
            delete numerics_container[iMGlevel][TURB_SOL][CONV_TERM];
        }
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (spalart_allmaras || neg_spalart_allmaras || menter_sst){
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          delete numerics_container[iMGlevel][TURB_SOL][VISC_TERM];
          delete numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM];
          delete numerics_container[iMGlevel][TURB_SOL][SOURCE_SECOND_TERM];
          /*--- Definition of the boundary condition method ---*/
          delete numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM];
          delete numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM];

      }
    }

  }

  /*--- Solver definition for the transition model problem ---*/
  if (transition) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          delete numerics_container[iMGlevel][TRANS_SOL][CONV_TERM];
        }
        break;
    }

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      /*--- Definition of the viscous scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TRANS_SOL][VISC_TERM];
      /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
      delete numerics_container[iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM];
      delete numerics_container[iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM];
      /*--- Definition of the boundary condition method ---*/
      delete numerics_container[iMGlevel][TRANS_SOL][CONV_BOUND_TERM];
    }
  }

  /*--- Solver definition for the poisson potential problem ---*/
  if (poisson || heat) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][POISSON_SOL][VISC_TERM];

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][POISSON_SOL][SOURCE_FIRST_TERM];
    delete numerics_container[MESH_0][POISSON_SOL][SOURCE_SECOND_TERM];

  }

  /*--- Solver definition for the flow adjoint problem ---*/

  if (adj_euler || adj_ns ) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/

    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case SPACE_CENTERED :

        if (compressible) {

          /*--- Compressible flow ---*/

          switch (config->GetKind_Centered_AdjFlow()) {
            case LAX : case JST:
              delete numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM];
              break;
          }

          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];

          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];

        }

        if (incompressible || freesurface) {

          /*--- Incompressible flow, use artificial compressibility method ---*/

          switch (config->GetKind_Centered_AdjFlow()) {
            case LAX : case JST:
              delete numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM]; break;
          }

          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];

          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];

        }

        break;

      case SPACE_UPWIND :

        if (compressible || incompressible || freesurface) {

          /*--- Compressible flow ---*/

          switch (config->GetKind_Upwind_AdjFlow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM];
                delete numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM];
              }
              break;
          }
        }

        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/

    if (compressible || incompressible || freesurface) {

      /*--- Compressible flow ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        delete numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM];
        delete numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM];
      }
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {


      if (compressible || incompressible || freesurface) {

        delete numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_FIRST_TERM];
        delete numerics_container[iMGlevel][ADJFLOW_SOL][SOURCE_SECOND_TERM];

      }
    }

  }


  /*--- Solver definition for the turbulent adjoint problem ---*/
  if (adj_turb) {
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTurb()) {

      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          if (spalart_allmaras) {
            delete numerics_container[iMGlevel][ADJTURB_SOL][CONV_TERM];
          }
        break;
    }


    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) {
        /*--- Definition of the viscous scheme for each equation and mesh level ---*/
        delete numerics_container[iMGlevel][ADJTURB_SOL][VISC_TERM];
        /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
        delete numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM];
        delete numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM];
        /*--- Definition of the boundary condition method ---*/
        delete numerics_container[iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM];
      }
    }
  }

  /*--- Solver definition for the wave problem ---*/
  if (wave) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][WAVE_SOL][VISC_TERM];

  }

  /*--- Solver definition for the FEA problem ---*/
  if (fem) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    delete numerics_container[MESH_0][FEA_SOL][VISC_TERM];

  }


  /*--- Definition of the Class for the numerical method: numerics_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    for (iSol = 0; iSol < MAX_SOLS; iSol++){
      delete numerics_container[iMGlevel][iSol];
    }
    delete[] numerics_container[iMGlevel];
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
      
    case FEM_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << ": FEM iteration." << endl;
      iteration_container[iZone] = new CFEM_StructuralAnalysis(config[iZone]);
      break;
    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << ": adjoint Euler/Navier-Stokes/RANS flow iteration." << endl;
      iteration_container[iZone] = new CAdjMeanFlowIteration(config[iZone]);
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
	for (targetZone = 0; targetZone < nZone; targetZone++){

    /*--- Initialize target booleans ---*/
    fluid_target  = false;  structural_target  = false;

    /*--- Set the target boolean: as of now, only Fluid-Structure Interaction considered ---*/
    switch (config_container[targetZone]->GetKind_Solver()) {
      case EULER : case NAVIER_STOKES: case RANS: fluid_target  = true;     break;
      case FEM_ELASTICITY:            structural_target = true;   break;
    }

		for (donorZone = 0; donorZone < nZone; donorZone++){
	    /*--- Initialize donor booleans ---*/
	    fluid_donor  = false;  structural_donor  = false;
	    matching_mesh = config_container[donorZone]->GetMatchingMesh();

	    /*--- Set the donor boolean: as of now, only Fluid-Structure Interaction considered ---*/
	    switch (config_container[donorZone]->GetKind_Solver()) {
	      case EULER : case NAVIER_STOKES: case RANS: fluid_donor  = true;    break;
	      case FEM_ELASTICITY:            structural_donor = true;  break;
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

				if (rank == MASTER_NODE) cout << "From zone " << donorZone << " to zone " << targetZone << ": ";

				/*--- Match Zones ---*/
				if (rank == MASTER_NODE) cout << "Setting coupling "<<endl;

				/*--- If the mesh is matching: match points ---*/
				if (matching_mesh){
					if (rank == MASTER_NODE) cout << "between matching meshes. " << endl;
					geometry_container[donorZone][MESH_0]->MatchZone(config_container[donorZone], geometry_container[targetZone][MESH_0],
							config_container[targetZone], donorZone, nZone);
				}
				/*--- Else: interpolate ---*/
				else {
          switch (config_container[donorZone]->GetKindInterpolation()){
            case NEAREST_NEIGHBOR:
              interpolator_container[donorZone][targetZone] = new CNearestNeighbor(geometry_container, config_container, donorZone, targetZone);
              if (rank == MASTER_NODE) cout << "using a nearest-neighbor approach." << endl;
              break;
            case ISOPARAMETRIC:
              interpolator_container[donorZone][targetZone] = new CIsoparametric(geometry_container, config_container, donorZone, targetZone);
              if (rank == MASTER_NODE) cout << "using an isoparametric approach." << endl;
              break;
            case CONSISTCONSERVE:
              if (targetZone>0 && structural_target){
                interpolator_container[donorZone][targetZone] = new CMirror(geometry_container, config_container, donorZone, targetZone);
                if (rank == MASTER_NODE) cout << "using a mirror approach: matching coefficients from opposite mesh." << endl;
              }
              else{
                interpolator_container[donorZone][targetZone] = new CIsoparametric(geometry_container, config_container, donorZone, targetZone);
                if (rank == MASTER_NODE) cout << "using an isoparametric approach." << endl;
              }
              if (targetZone==0 && structural_target){
                if (rank == MASTER_NODE) cout << "Consistent and conservative interpolation assumes the structure model mesh is evaluated second. Somehow this has not happened. The isoparametric coefficients will be calculated for both meshes, and are not guaranteed to be consistent." << endl;
              }
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


CDriver::~CDriver(void) {}


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

  /*--- Run an iteration of the physics within this single zone.
   We assume that the zone of interest is in the ZONE_0 container position. ---*/
  
  iteration_container[ZONE_0]->Preprocess(output, integration_container, geometry_container,
                                          solver_container, numerics_container, config_container,
                                          surface_movement, grid_movement, FFDBox, ZONE_0);
  
  iteration_container[ZONE_0]->Iterate(output, integration_container, geometry_container,
                                       solver_container, numerics_container, config_container,
                                       surface_movement, grid_movement, FFDBox, ZONE_0);
  
  iteration_container[ZONE_0]->Update(output, integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, ZONE_0);
  
  iteration_container[ZONE_0]->Monitor();     /*--- Does nothing for now. ---*/
  
  iteration_container[ZONE_0]->Output();      /*--- Does nothing for now. ---*/
  
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
    
    iteration_container[iZone]->Preprocess(output, integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, iZone);
    
    iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, iZone);
    
    iteration_container[iZone]->Update(output, integration_container, geometry_container,
                                       solver_container, numerics_container, config_container,
                                       surface_movement, grid_movement, FFDBox, iZone);
    
    iteration_container[iZone]->Monitor();     /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Output();      /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Postprocess(); /*--- Does nothing for now. ---*/
    
  }
  
}

CSpectralDriver::CSpectralDriver(CIteration **iteration_container,
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

CSpectralDriver::~CSpectralDriver(void) { }

void CSpectralDriver::Run(CIteration **iteration_container,
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
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- If this is the first iteration, set up the spectral operators,
   initialize the source terms, and compute any grid veocities, if necessary. ---*/
  
  if (ExtIter == 0) {
    SetTimeSpectral_Velocities(geometry_container, config_container, nZone);
    for (iZone = 0; iZone < nZone; iZone++)
      SetTimeSpectral(geometry_container, solver_container, config_container, nZone, (iZone+1)%nZone);
  }
  
  /*--- Run a single iteration of a spectral method problem. Preprocess all
   all zones before beginning the iteration. ---*/
  
  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone]->Preprocess(output, integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, iZone);
  
  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone]->Iterate(output, integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, iZone);
    
    
    /*--- Update the spectral source terms across all zones ---*/
    
    SetTimeSpectral(geometry_container, solver_container, config_container, nZone, (iZone+1)%nZone);
    
    iteration_container[iZone]->Update(output, integration_container, geometry_container,
                                       solver_container, numerics_container, config_container,
                                       surface_movement, grid_movement, FFDBox, iZone);
    
    iteration_container[iZone]->Monitor();     /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Output();      /*--- Does nothing for now. ---*/
    
    iteration_container[iZone]->Postprocess(); /*--- Does nothing for now. ---*/
    
  }
  
}

void CSpectralDriver::SetTimeSpectral(CGeometry ***geometry_container, CSolver ****solver_container,
                                      CConfig **config_container, unsigned short nZone, unsigned short iZone) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short iVar, jZone, kZone, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());
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
  su2double period = config_container[ZONE_0]->GetTimeSpectral_Period();
  
  /*--- allocate dynamic memory for D ---*/
  su2double **D = new su2double*[nZone];
  for (kZone = 0; kZone < nZone; kZone++) {
    D[kZone] = new su2double[nZone];
  }
  
  /*--- Build the time-spectral operator matrix ---*/
  ComputeTimeSpectral_Operator(D, period, nZone);
  //  for (kZone = 0; kZone < nZone; kZone++) {
  //    for (jZone = 0; jZone < nZone; jZone++) {
  //
  //      if (nZone%2 == 0) {
  //
  //        /*--- For an even number of time instances ---*/
  //        if (kZone == jZone) {
  //          D[kZone][jZone] = 0.0;
  //        }
  //        else {
  //          D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/tan(PI_NUMBER*(kZone-jZone)/nZone));
  //        }
  //      }
  //      else {
  //
  //        /*--- For an odd number of time instances ---*/
  //        if (kZone == jZone) {
  //          D[kZone][jZone] = 0.0;
  //        }
  //        else {
  //          D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/sin(PI_NUMBER*(kZone-jZone)/nZone));
  //        }
  //      }
  //
  //    }
  //  }
  
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
            solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar, Source[iVar]);
          }
          else {
            solver_container[iZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar, Source[iVar]);
          }
        }
        
      }
    }
  }
  
  //	/*--- Loop over all grid levels ---*/
  //	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {
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
  //					solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar, Source[iVar]);
  //			}
  //		}
  //	}
  
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
        solver_container[iZone][MESH_0][TURB_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar, Source_Turb[iVar]);
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

void CSpectralDriver::ComputeTimeSpectral_Operator(su2double **D, su2double period, unsigned short nZone) {
  
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

void CSpectralDriver::SetTimeSpectral_Velocities(CGeometry ***geometry_container,
                                                 CConfig **config_container, unsigned short nZone) {
  
  unsigned short iZone, jDegree, iDim, iMGlevel;
  unsigned short nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();
  
  su2double angular_interval = 2.0*PI_NUMBER/(su2double)(nZone);
  su2double *Coord;
  unsigned long iPoint;
  
  /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
  su2double period = config_container[ZONE_0]->GetTimeSpectral_Period();
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

	unsigned long IntIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetIntIter(IntIter);
	unsigned long FSIIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(FSIIter);
	unsigned long nFSIIter = config_container[ZONE_FLOW]->GetnIterFSI();

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	 /*--- If there is a restart, we need to get the old geometry from the fluid field ---*/
	 bool restart = (config_container[ZONE_FLOW]->GetRestart() || config_container[ZONE_FLOW]->GetRestart_Flow());
	 unsigned long ExtIter = config_container[ZONE_FLOW]->GetExtIter();

	 if (restart && (long)ExtIter == config_container[ZONE_FLOW]->GetUnst_RestartIter()){
		unsigned short ZONE_FLOW = 0;
		solver_container[ZONE_FLOW][MESH_0][FLOW_SOL]->Restart_OldGeometry(geometry_container[ZONE_FLOW][MESH_0],config_container[ZONE_FLOW]);
	 }

	/*-----------------------------------------------------------------*/
	/*---------------- Predict structural displacements ---------------*/
	/*-----------------------------------------------------------------*/

	Predict_Displacements(output, integration_container, geometry_container,
            		      solver_container, numerics_container, config_container,
            		      surface_movement, grid_movement, FFDBox,
            		      ZONE_STRUCT, ZONE_FLOW);

	while (FSIIter < nFSIIter){

		/*-----------------------------------------------------------------*/
		/*------------------- Transfer Displacements ----------------------*/
		/*-----------------------------------------------------------------*/

		Transfer_Displacements(output, integration_container, geometry_container,
                solver_container, numerics_container, config_container,
                surface_movement, grid_movement, FFDBox, transfer_container,
                ZONE_STRUCT, ZONE_FLOW);

		/*-----------------------------------------------------------------*/
		/*-------------------- Fluid subiteration -------------------------*/
		/*-----------------------------------------------------------------*/

		iteration_container[ZONE_FLOW]->Preprocess(output, integration_container, geometry_container,
		                                       solver_container, numerics_container, config_container,
		                                       surface_movement, grid_movement, FFDBox, ZONE_FLOW);

		iteration_container[ZONE_FLOW]->Iterate(output, integration_container, geometry_container,
		                                       solver_container, numerics_container, config_container,
		                                       surface_movement, grid_movement, FFDBox, ZONE_FLOW);

		/*--- Write the convergence history for the fluid (only screen output) ---*/

		output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_FLOW);

		/*--- Set the fluid convergence to false (to make sure FSI subiterations converge) ---*/

		integration_container[ZONE_FLOW][FLOW_SOL]->SetConvergence(false);

		/*-----------------------------------------------------------------*/
		/*------------------- Set FEA loads from fluid --------------------*/
		/*-----------------------------------------------------------------*/

		Transfer_Tractions(output, integration_container, geometry_container,
                solver_container, numerics_container, config_container,
                surface_movement, grid_movement, FFDBox, transfer_container,
                ZONE_FLOW, ZONE_STRUCT);

		/*-----------------------------------------------------------------*/
		/*------------------ Structural subiteration ----------------------*/
		/*-----------------------------------------------------------------*/

		iteration_container[ZONE_STRUCT]->Iterate(output, integration_container, geometry_container,
		                                       solver_container, numerics_container, config_container,
		                                       surface_movement, grid_movement, FFDBox, ZONE_STRUCT);

		/*--- Write the convergence history for the structure (only screen output) ---*/

		output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUCT);

		/*--- Set the fluid convergence to false (to make sure FSI subiterations converge) ---*/

		integration_container[ZONE_STRUCT][FEA_SOL]->SetConvergence(false);

		/*-----------------------------------------------------------------*/
		/*----------------- Displacements relaxation ----------------------*/
		/*-----------------------------------------------------------------*/

		Relaxation_Displacements(output, geometry_container, solver_container, config_container,
								 ZONE_STRUCT, ZONE_FLOW, FSIIter);

		/*-----------------------------------------------------------------*/
		/*-------------------- Check convergence --------------------------*/
		/*-----------------------------------------------------------------*/

		integration_container[ZONE_STRUCT][FEA_SOL]->Convergence_Monitoring_FSI(geometry_container[ZONE_STRUCT][MESH_0], config_container[ZONE_STRUCT],
																solver_container[ZONE_STRUCT][MESH_0][FEA_SOL], FSIIter);

		if (integration_container[ZONE_STRUCT][FEA_SOL]->GetConvergence_FSI()) break;

		/*-----------------------------------------------------------------*/
		/*--------------------- Update FSIIter ---------------------------*/
		/*-----------------------------------------------------------------*/

		FSIIter++; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetFSIIter(FSIIter);

	}

	/*-----------------------------------------------------------------*/
  	/*------------------ Update coupled solver ------------------------*/
	/*-----------------------------------------------------------------*/

	Update(output, integration_container, geometry_container,
           solver_container, numerics_container, config_container,
           surface_movement, grid_movement, FFDBox, transfer_container,
           ZONE_FLOW, ZONE_STRUCT);


	/*-----------------------------------------------------------------*/
  	/*-------------------- Update fluid solver ------------------------*/
	/*-----------------------------------------------------------------*/

	iteration_container[ZONE_FLOW]->Update(output, integration_container, geometry_container,
	          	  	  	  	  	  	  	  solver_container, numerics_container, config_container,
	          	  	  	  	  	  	  	  surface_movement, grid_movement, FFDBox, ZONE_FLOW);

	/*-----------------------------------------------------------------*/
  	/*----------------- Update structural solver ----------------------*/
	/*-----------------------------------------------------------------*/

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
//			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);

		}
		else {
			transfer_container[donorZone][targetZone]->Broadcast_InterfaceData_Interpolate(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
			/*--- Set the volume deformation for the fluid zone ---*/
//			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);

		}
		break;
	case SCATTER_DATA:
		if (MatchingMesh){
			transfer_container[donorZone][targetZone]->Scatter_InterfaceData(solver_container[donorZone][MESH_0][FEA_SOL],solver_container[targetZone][MESH_0][FLOW_SOL],
					geometry_container[donorZone][MESH_0],geometry_container[targetZone][MESH_0],
					config_container[donorZone], config_container[targetZone]);
			/*--- Set the volume deformation for the fluid zone ---*/
//			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
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
//			grid_movement[targetZone]->SetVolume_Deformation(geometry_container[targetZone][MESH_0], config_container[targetZone], true);
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
					config_container[targetZone], config_container[donorZone], numerics_container[targetZone][MESH_0][SolContainer_Position_fea][FEA_TERM]);
		}
		else {
			solver_container[targetZone][MESH_0][FEA_SOL]->SetFEA_Load_Int(solver_container[donorZone], geometry_container[targetZone], geometry_container[donorZone],
					config_container[targetZone], config_container[donorZone], numerics_container[targetZone][MESH_0][SolContainer_Position_fea][FEA_TERM]);
		}
		break;
	}

}

void CFSIDriver::Relaxation_Displacements(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long FSIIter){

#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/*-------------------- Aitken's relaxation ------------------------*/

	/*------------------- Compute the coefficient ---------------------*/

	solver_container[donorZone][MESH_0][FEA_SOL]->ComputeAitken_Coefficient(geometry_container[donorZone], config_container[donorZone],
			solver_container[donorZone], FSIIter);

	/*----------------- Set the relaxation parameter ------------------*/

	solver_container[donorZone][MESH_0][FEA_SOL]->SetAitken_Relaxation(geometry_container[donorZone], config_container[donorZone],
			solver_container[donorZone]);


	/*----------------- Communicate the predicted solution and the old one ------------------*/
	solver_container[donorZone][MESH_0][FEA_SOL]->Set_MPI_Solution_Pred_Old(geometry_container[donorZone][MESH_0], config_container[donorZone]);


}

void CFSIDriver::Relaxation_Tractions(COutput *output, CGeometry ***geometry_container, CSolver ****solver_container,
			CConfig **config_container, unsigned short donorZone, unsigned short targetZone, unsigned long FSIIter){

}

void CFSIDriver::Update(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
			 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
			 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox,
			 CTransfer ***transfer_container, unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT){

	unsigned long IntIter = 0; // This doesn't affect here but has to go into the function
	unsigned long ExtIter = config_container[ZONE_FLOW]->GetExtIter();


	/*-----------------------------------------------------------------*/
	/*--------------------- Enforce continuity ------------------------*/
	/*-----------------------------------------------------------------*/

	/*--- Enforces that the geometry of the flow corresponds to the converged, relaxed solution ---*/

	/*-------------------- Transfer the displacements --------------------*/

	Transfer_Displacements(output, integration_container, geometry_container,
            solver_container, numerics_container, config_container,
            surface_movement, grid_movement, FFDBox, transfer_container,
            ZONE_STRUCT, ZONE_FLOW);

	/*-------------------- Set the grid movement -------------------------*/

	SetGrid_Movement(geometry_container[ZONE_FLOW], surface_movement[ZONE_FLOW],
				grid_movement[ZONE_FLOW], FFDBox[ZONE_FLOW], solver_container[ZONE_FLOW], config_container[ZONE_FLOW],
				ZONE_FLOW, IntIter, ExtIter);

	/*----------- Store the solution_pred as solution_pred_old --------------*/

}


