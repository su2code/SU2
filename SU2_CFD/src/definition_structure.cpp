/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
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

#include "../include/definition_structure.hpp"


unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config) {
  string text_line, Marker_Tag;
  ifstream mesh_file;
  short nZone = 1; // Default value
  unsigned short iLine, nLine = 10;
  char cstr[200];
  string::size_type position;
  
  /*--- Search the mesh file for the 'NZONE' keyword. ---*/
  
  switch (val_format) {
    case SU2:
      
      /*--- Open grid file ---*/
      
      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      if (mesh_file.fail()) {
        cout << "cstr=" << cstr << endl;
        cout << "There is no geometry file (GetnZone))!" << endl;
        
#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
      }
      
      /*--- Read the SU2 mesh file ---*/
      
      for (iLine = 0; iLine < nLine ; iLine++) {
        
        getline (mesh_file, text_line);
        
        /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
        
        position = text_line.find ("NZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nZone = atoi(text_line.c_str());
        }
      }
      
      break;
      
  }
  
  /*--- For time spectral integration, nZones = nTimeInstances. ---*/
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    nZone = config->GetnTimeInstances();
  }
  
  return (unsigned short) nZone;
}

unsigned short GetnDim(string val_mesh_filename, unsigned short val_format) {
  
  string text_line, Marker_Tag;
  ifstream mesh_file;
  short nDim = 3;
  unsigned short iLine, nLine = 10;
  char cstr[200];
  string::size_type position;
  
  /*--- Open grid file ---*/
  
  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);
  
  switch (val_format) {
    case SU2:
      
      /*--- Read SU2 mesh file ---*/
      
      for (iLine = 0; iLine < nLine ; iLine++) {
        
        getline (mesh_file, text_line);
        
        /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
        
        position = text_line.find ("NDIME=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nDim = atoi(text_line.c_str());
        }
      }
      break;
      
    case CGNS:
      
#ifdef HAVE_CGNS
      
      /*--- Local variables which are needed when calling the CGNS mid-level API. ---*/
      
      int fn, nbases = 0, nzones = 0, file_type;
      int cell_dim = 0, phys_dim = 0;
      char basename[CGNS_STRING_SIZE];
      
      /*--- Check whether the supplied file is truly a CGNS file. ---*/
      
      if ( cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK ) {
        printf( "\n\n   !!! Error !!!\n" );
        printf( " %s is not a CGNS file.\n", val_mesh_filename.c_str());
        printf( " Now exiting...\n\n");
        exit(EXIT_FAILURE);
      }
      
      /*--- Open the CGNS file for reading. The value of fn returned
       is the specific index number for this file and will be
       repeatedly used in the function calls. ---*/
      
      if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn)) cg_error_exit();
      
      /*--- Get the number of databases. This is the highest node
       in the CGNS heirarchy. ---*/
      
      if (cg_nbases(fn, &nbases)) cg_error_exit();
      
      /*--- Check if there is more than one database. Throw an
       error if there is because this reader can currently
       only handle one database. ---*/
      
      if ( nbases > 1 ) {
        printf("\n\n   !!! Error !!!\n" );
        printf("CGNS reader currently incapable of handling more than 1 database.");
        printf("Now exiting...\n\n");
        exit(EXIT_FAILURE);
      }
      
      /*--- Read the databases. Note that the indexing starts at 1. ---*/
      
      for ( int i = 1; i <= nbases; i++ ) {
        
        if (cg_base_read(fn, i, basename, &cell_dim, &phys_dim)) cg_error_exit();
        
        /*--- Get the number of zones for this base. ---*/
        
        if (cg_nzones(fn, i, &nzones)) cg_error_exit();
        
      }
      
      /*--- Set the problem dimension as read from the CGNS file ---*/
      
      nDim = cell_dim;
      
#endif
      
      break;
      
  }
  
  mesh_file.close();
  
  return (unsigned short) nDim;
}

void Geometrical_Preprocessing(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {
  
  unsigned short iMGlevel, iZone;
  unsigned long iPoint; 
  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Compute elements surrounding points, points surrounding points ---*/
    
    if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
    geometry[iZone][MESH_0]->SetPoint_Connectivity();
    
    /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/
    
    if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
    geometry[iZone][MESH_0]->SetRCM_Ordering(config[iZone]);
    
    /*--- recompute elements surrounding points, points surrounding points ---*/
    
    if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
    geometry[iZone][MESH_0]->SetPoint_Connectivity();
    
    /*--- Compute elements surrounding elements ---*/
    
    if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
    geometry[iZone][MESH_0]->SetElement_Connectivity();
    
    /*--- Check the orientation before computing geometrical quantities ---*/
    
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry[iZone][MESH_0]->SetBoundVolume();
    geometry[iZone][MESH_0]->Check_IntElem_Orientation(config[iZone]);
    geometry[iZone][MESH_0]->Check_BoundElem_Orientation(config[iZone]);

    /*--- Create the edge structure ---*/
    
    if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
    geometry[iZone][MESH_0]->SetEdges();
    geometry[iZone][MESH_0]->SetVertex(config[iZone]);
    
    /*--- Compute cell center of gravity ---*/
    
    if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
    geometry[iZone][MESH_0]->SetCG();
    
    /*--- Create the control volume structures ---*/
    
    if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
    geometry[iZone][MESH_0]->SetControlVolume(config[iZone], ALLOCATE);
    geometry[iZone][MESH_0]->SetBoundControlVolume(config[iZone], ALLOCATE);
    
    /*--- Visualize a dual control volume if requested ---*/
    
    if ((config[iZone]->GetVisualize_CV() >= 0) &&
        (config[iZone]->GetVisualize_CV() < (long)geometry[iZone][MESH_0]->GetnPointDomain()))
      geometry[iZone][MESH_0]->VisualizeControlVolume(config[iZone], UPDATE);
    
    /*--- Identify closest normal neighbor ---*/
    
    if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
    geometry[iZone][MESH_0]->FindNormal_Neighbor(config[iZone]);
    
    /*--- Compute the surface curvature ---*/
    
    if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
    geometry[iZone][MESH_0]->ComputeSurf_Curvature(config[iZone]);
    
    if ((config[iZone]->GetnMGLevels() != 0) && (rank == MASTER_NODE))
      cout << "Setting the multigrid structure." << endl;
    
  }
  
  /*--- Loop over all the new grid ---*/
  
  for (iMGlevel = 1; iMGlevel <= config[ZONE_0]->GetnMGLevels(); iMGlevel++) {
    
    /*--- Loop over all zones at each grid level. ---*/
    
    for (iZone = 0; iZone < val_nZone; iZone++) {
      
      /*--- Create main agglomeration structure ---*/
      
      geometry[iZone][iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel, iZone);
      
      /*--- Compute points surrounding points. ---*/
      
      geometry[iZone][iMGlevel]->SetPoint_Connectivity(geometry[iZone][iMGlevel-1]);
      
      /*--- Create the edge structure ---*/
      
      geometry[iZone][iMGlevel]->SetEdges();
      geometry[iZone][iMGlevel]->SetVertex(geometry[iZone][iMGlevel-1], config[iZone]);
      
      /*--- Create the control volume structures ---*/
      
      geometry[iZone][iMGlevel]->SetControlVolume(config[iZone], geometry[iZone][iMGlevel-1], ALLOCATE);
      geometry[iZone][iMGlevel]->SetBoundControlVolume(config[iZone], geometry[iZone][iMGlevel-1], ALLOCATE);
      geometry[iZone][iMGlevel]->SetCoord(geometry[iZone][iMGlevel-1]);
      
      /*--- Find closest neighbor to a surface point ---*/
      
      geometry[iZone][iMGlevel]->FindNormal_Neighbor(config[iZone]);
      
    }
    
  }
  
  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    if (config[iZone]->GetUnsteady_Simulation() && config[iZone]->GetGrid_Movement()) {
      for (iMGlevel = 0; iMGlevel <= config[iZone]->GetnMGLevels(); iMGlevel++) {
        for (iPoint = 0; iPoint < geometry[iZone][iMGlevel]->GetnPoint(); iPoint++) {
          
          /*--- Update cell volume ---*/
          
          geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_n();
          geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_nM1();
          
          /*--- Update point coordinates ---*/
          geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n();
          geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n1();
          
        }
      }
    }
  }
  
}

void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry,
                          CConfig *config, unsigned short iZone) {
  
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
  adj_euler        = false;	 adj_ns          = false;	 adj_turb  = false;
  lin_euler        = false;	 lin_ns          = false;
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

void Integration_Preprocessing(CIntegration **integration_container,
                               CGeometry **geometry, CConfig *config,
                               unsigned short iZone) {
  
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

void Numerics_Preprocessing(CNumerics ****numerics_container,
                            CSolver ***solver_container, CGeometry **geometry,
                            CConfig *config, unsigned short iZone) {
  
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
  adj_euler        = false;	  adj_ns           = false;	 adj_turb         = false;
  wave             = false;   heat             = false;   fea              = false;   spalart_allmaras = false; neg_spalart_allmaras = false;
  tne2_euler       = false;   tne2_ns          = false;	fem				= false;
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
  
  if (euler)				nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (ns)	          nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (turbulent)		nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
  if (transition)		nVar_Trans = solver_container[MESH_0][TRANS_SOL]->GetnVar();
  if ((tne2_euler) || (tne2_ns)) {
    nVar_TNE2         = solver_container[MESH_0][TNE2_SOL]->GetnVar();
    nPrimVar_TNE2     = solver_container[MESH_0][TNE2_SOL]->GetnPrimVar();
    nPrimVarGrad_TNE2 = solver_container[MESH_0][TNE2_SOL]->GetnPrimVarGrad();
  }
  if (poisson)			nVar_Poisson = solver_container[MESH_0][POISSON_SOL]->GetnVar();
  
  if (wave)				nVar_Wave = solver_container[MESH_0][WAVE_SOL]->GetnVar();
  if (fea)				nVar_FEA = solver_container[MESH_0][FEA_SOL]->GetnVar();
  if (fem)				nVar_FEM = solver_container[MESH_0][FEA_SOL]->GetnVar();
  if (heat)				nVar_Heat = solver_container[MESH_0][HEAT_SOL]->GetnVar();
  
  /*--- Number of variables for adjoint problem ---*/
  
  if (adj_euler)    	  nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_ns)			      nVar_Adj_Flow = solver_container[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_turb)		      nVar_Adj_Turb = solver_container[MESH_0][ADJTURB_SOL]->GetnVar();
  if ((adj_tne2_euler) || (adj_tne2_ns)) {
    nVar_Adj_TNE2         = solver_container[MESH_0][ADJTNE2_SOL]->GetnVar();
    nPrimVar_Adj_TNE2     = solver_container[MESH_0][ADJTNE2_SOL]->GetnPrimVar();
    nPrimVarGrad_Adj_TNE2 = solver_container[MESH_0][ADJTNE2_SOL]->GetnPrimVarGrad();
  }
  
  /*--- Number of variables for the linear problem ---*/
  
  if (lin_euler)	nVar_Lin_Flow = solver_container[MESH_0][LINFLOW_SOL]->GetnVar();
  
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
    numerics_container[MESH_0][FEA_SOL][VISC_TERM] = new CGalerkin_FEA(nDim, nVar_Wave, config);
    
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
