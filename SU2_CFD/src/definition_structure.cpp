/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.0.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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
  short nZone = 1;
  bool isFound = false;
  char cstr[200];
  string::size_type position;
  int rank = MASTER_NODE;
  
#ifndef NO_MPI
  int size;
#ifdef WINDOWS
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
#else
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#endif
  if (size != 1) {
    unsigned short lastindex = val_mesh_filename.find_last_of(".");
    val_mesh_filename = val_mesh_filename.substr(0, lastindex);
    val_mesh_filename = val_mesh_filename + "_1.su2";
  }
#endif
  
  /*--- Search the mesh file for the 'NZONE' keyword. ---*/
  switch (val_format) {
    case SU2:
      
      /*--- Open grid file ---*/
      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      if (mesh_file.fail()) {
        cout << "cstr=" << cstr << endl;
        cout << "There is no geometry file (GetnZone))!" << endl;
#ifdef NO_MPI
        exit(1);
#else
#ifdef WINDOWS
		MPI_Abort(MPI_COMM_WORLD,1);
		MPI_Finalize();
#else
        MPI::COMM_WORLD.Abort(1);
        MPI::Finalize();
#endif
#endif
      }
      
      /*--- Open the SU2 mesh file ---*/
      while (getline (mesh_file,text_line)) {
        
        /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
        position = text_line.find ("NZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nZone = atoi(text_line.c_str()); isFound = true;
          if (rank == MASTER_NODE) {
            //					if (nZone == 1) cout << "SU2 mesh file format with a single zone." << endl;
            //					else if (nZone >  1) cout << "SU2 mesh file format with " << nZone << " zones." << endl;
            //					else
            if (nZone <= 0) {
              cout << "Error: Number of mesh zones is less than 1 !!!" << endl;
#ifdef NO_MPI
              exit(1);
#else
#ifdef WINDOWS
			  MPI_Abort(MPI_COMM_WORLD,1);
			  MPI_Finalize();
#else
              MPI::COMM_WORLD.Abort(1);
              MPI::Finalize();
#endif
#endif
            }
          }
        }
      }
      /*--- If the "NZONE" keyword was not found, assume this is an ordinary
       simulation on a single Zone ---*/
      if (!isFound) {
        nZone = 1;
        //			if (rank == MASTER_NODE) cout << "SU2 mesh file format with a single zone." << endl;
      }
      break;
      
    case CGNS:
      
      nZone = 1;
      //		if (rank == MASTER_NODE) cout << "CGNS mesh file format with a single zone." << endl;
      break;
      
    case NETCDF_ASCII:
      
      nZone = 1;
      //		if (rank == MASTER_NODE) cout << "NETCDF mesh file format with a single zone." << endl;
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
  bool isFound = false;
  char cstr[200];
  string::size_type position;
  
#ifndef NO_MPI
  int size;
#ifdef WINDOWS
  MPI_Comm_size(MPI_COMM_WORLD,&size);
#else
  size = MPI::COMM_WORLD.Get_size();
#endif
  if (size != 1) {
    unsigned short lastindex = val_mesh_filename.find_last_of(".");
    val_mesh_filename = val_mesh_filename.substr(0, lastindex);
    val_mesh_filename = val_mesh_filename + "_1.su2";
  }
#endif
  
  switch (val_format) {
    case SU2:
      
      /*--- Open grid file ---*/
      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      
      /*--- Read SU2 mesh file ---*/
      while (getline (mesh_file,text_line)) {
        /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
        position = text_line.find ("NDIME=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nDim = atoi(text_line.c_str()); isFound = true;
        }
      }
      break;
      
    case CGNS:
      nDim = 3;
      break;
      
    case NETCDF_ASCII:
      nDim = 3;
      break;
  }
  return (unsigned short) nDim;
}

void Geometrical_Preprocessing(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {
  
  unsigned short iMGlevel, iZone;
  unsigned long iPoint; 
  int rank = MASTER_NODE;

#ifndef NO_MPI
#ifdef WINDOWS
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
  rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Compute elements surrounding points, points surrounding points,
     and elements surrounding elements ---*/
    
    if (rank == MASTER_NODE) cout << "Setting local point and element connectivity." << endl;
    geometry[iZone][MESH_0]->SetEsuP();
    geometry[iZone][MESH_0]->SetPsuP();
    geometry[iZone][MESH_0]->SetEsuE();
    
    /*--- Check the orientation before computing geometrical quantities ---*/
    
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry[iZone][MESH_0]->SetBoundVolume();
    geometry[iZone][MESH_0]->Check_Orientation(config[iZone]);
    
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
    
    /*--- Visualize a control volume if requested ---*/
    
    if ((config[iZone]->GetVisualize_CV() >= 0) &&
        (config[iZone]->GetVisualize_CV() < geometry[iZone][MESH_0]->GetnPointDomain()))
      geometry[iZone][MESH_0]->VisualizeControlVolume(config[iZone], UPDATE);
    
    /*--- Identify closest normal neighbor ---*/
    
    if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
    geometry[iZone][MESH_0]->FindNormal_Neighbor(config[iZone]);
    
    /*--- Compute the surface curvature ---*/
    
    if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
    geometry[iZone][MESH_0]->ComputeSurf_Curvature(config[iZone]);
    
    if ((config[iZone]->GetMGLevels() != 0) && (rank == MASTER_NODE))
      cout << "Setting the multigrid structure." <<endl;
    
  }
  
#ifndef NO_MPI
  /*--- Synchronization point before the multigrid algorithm ---*/
#ifdef WINDOWS
  MPI_Barrier(MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Barrier();
#endif
#endif
  
  /*--- Loop over all the new grid ---*/
  
  for (iMGlevel = 1; iMGlevel <= config[ZONE_0]->GetMGLevels(); iMGlevel++) {
    
    /*--- Loop over all zones at each grid level. ---*/
    
    for (iZone = 0; iZone < val_nZone; iZone++) {
      
      /*--- Create main agglomeration structure ---*/
      
      geometry[iZone][iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel, iZone);
      
      /*--- Compute points surrounding points. ---*/
      
      geometry[iZone][iMGlevel]->SetPsuP(geometry[iZone][iMGlevel-1]);
      
      /*--- Create the edge structure ---*/
      
      geometry[iZone][iMGlevel]->SetEdges();
      geometry[iZone][iMGlevel]->SetVertex(geometry[iZone][iMGlevel-1], config[iZone]);
      
      /*--- Create the control volume structures ---*/
      
      geometry[iZone][iMGlevel]->SetControlVolume(config[iZone],geometry[iZone][iMGlevel-1], ALLOCATE);
      geometry[iZone][iMGlevel]->SetBoundControlVolume(config[iZone],geometry[iZone][iMGlevel-1], ALLOCATE);
      geometry[iZone][iMGlevel]->SetCoord(geometry[iZone][iMGlevel-1]);
      
      /*--- Find closest neighbor to a surface point ---*/
      
      geometry[iZone][iMGlevel]->FindNormal_Neighbor(config[iZone]);
      
    }
    
  }
  
  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    if (config[iZone]->GetUnsteady_Simulation() && config[iZone]->GetGrid_Movement()) {
      for (iMGlevel = 0; iMGlevel <= config[iZone]->GetMGLevels(); iMGlevel++) {
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
  lin_euler, lin_ns, lin_turb,
  tne2_euler, tne2_ns,
  adj_tne2_euler, adj_tne2_ns,
  poisson, wave, fea, heat,
  spalart_allmaras, menter_sst, machine_learning, transition,
  template_solver;
  
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;	 adj_ns          = false;	 adj_turb  = false;
  lin_euler        = false;	 lin_ns          = false;  lin_turb  = false;
  tne2_euler       = false;  tne2_ns         = false;
  adj_tne2_euler   = false;  adj_tne2_ns     = false;
  spalart_allmaras = false;  menter_sst      = false;   machine_learning = false;
  poisson          = false;
  wave             = false;
  fea              = false;
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
    case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
    case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
    case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case LINEAR_ELASTICITY: fea = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case ADJ_TNE2_EULER : tne2_euler = true; adj_tne2_euler = true; break;
    case ADJ_TNE2_NAVIER_STOKES : tne2_ns = true; adj_tne2_ns = true; break;
    case LIN_EULER: euler = true; lin_euler = true; break;
  }
  /*--- Assign turbulence model booleans --- */
  if (turbulent)
    switch (config->GetKind_Turb_Model()){
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; break;
      case ML: machine_learning = true;break;
        
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
    
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
        solver_container[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel);
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
      solver_container[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolver(geometry[iMGlevel], config);
    }
    
    /*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
    if (lin_euler) {
      solver_container[iMGlevel][LINFLOW_SOL] = new CLinEulerSolver(geometry[iMGlevel], config);
    }
    if (lin_ns) {
      cout <<"Equation not implemented." << endl; exit(1); break;
    }
    
  }
  
}

void Integration_Preprocessing(CIntegration **integration_container,
                               CGeometry **geometry, CConfig *config,
                               unsigned short iZone) {
  
  bool
  euler, adj_euler, lin_euler,
  ns, adj_ns, lin_ns,
  turbulent, adj_turb, lin_turb,
  spalart_allmaras, menter_sst,machine_learning,
  tne2_euler, adj_tne2_euler,
  tne2_ns, adj_tne2_ns,
  poisson, wave, fea, heat, template_solver, transition;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false; lin_euler         = false;
  ns               = false; adj_ns           = false; lin_ns            = false;
  turbulent        = false; adj_turb         = false; lin_turb          = false;
  spalart_allmaras = false; menter_sst       = false; machine_learning  = false;
  tne2_euler       = false; adj_tne2_euler   = false;
  tne2_ns          = false; adj_tne2_ns      = false;
  poisson          = false;
  wave             = false;
  heat             = false;
  fea              = false;
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
    case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
    case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
    case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case LINEAR_ELASTICITY: fea = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_TNE2_EULER : tne2_euler = true; adj_tne2_euler = true; break;
    case ADJ_TNE2_NAVIER_STOKES : tne2_ns = true; adj_tne2_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case LIN_EULER: euler = true; lin_euler = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  if (turbulent) {
    switch (config->GetKind_Turb_Model()) {
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; break;
      case ML: machine_learning = true; break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
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
  if (fea) integration_container[FEA_SOL] = new CSingleGridIntegration(config);
  
  /*--- Allocate solution for adjoint problem ---*/
  if (adj_euler) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_ns) integration_container[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_tne2_euler) integration_container[ADJTNE2_SOL] = new CMultiGridIntegration(config);
  if (adj_tne2_ns) integration_container[ADJTNE2_SOL] = new CMultiGridIntegration(config);
  if (adj_turb) integration_container[ADJTURB_SOL] = new CSingleGridIntegration(config);
  
  /*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
  if (lin_euler) integration_container[LINFLOW_SOL] = new CMultiGridIntegration(config);
  if (lin_ns) { cout <<"Equation not implemented." << endl; exit(1); }
  
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
  nVar_Wave             = 0,
  nVar_Heat             = 0,
  nVar_Lin_Flow         = 0;
  
  double *constants = NULL;
  
  bool
  euler, adj_euler, lin_euler,
  ns, adj_ns, lin_ns,
  turbulent, adj_turb, lin_turb,
  tne2_euler, adj_tne2_euler,
  tne2_ns, adj_tne2_ns,
  spalart_allmaras, menter_sst, machine_learning,
  poisson,
  wave,
  fea,
  heat,
  transition,
  template_solver;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;   turbulent        = false;
  poisson          = false;
  adj_euler        = false;	  adj_ns           = false;	 adj_turb         = false;
  wave             = false;   heat             = false;   fea              = false;   spalart_allmaras = false;
  tne2_euler       = false;   tne2_ns          = false;
  adj_tne2_euler   = false;	  adj_tne2_ns      = false;
  lin_euler        = false;   lin_ns           = false;   lin_turb         = false;	 menter_sst       = false;    machine_learning = false;
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
    case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
    case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
    case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
    case POISSON_EQUATION: poisson = true; break;
    case WAVE_EQUATION: wave = true; break;
    case HEAT_EQUATION: heat = true; break;
    case LINEAR_ELASTICITY: fea = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_TNE2_EULER : tne2_euler = true; adj_tne2_euler = true; break;
    case ADJ_TNE2_NAVIER_STOKES : tne2_ns = true; adj_tne2_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
    case LIN_EULER: euler = true; lin_euler = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  if (turbulent)
    switch (config->GetKind_Turb_Model()){
      case SA: spalart_allmaras = true; break;
      case ML: machine_learning = true; break;
      case SST: menter_sst = true; constants = solver_container[MESH_0][TURB_SOL]->GetConstants(); break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
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
  for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
    numerics_container[iMGlevel] = new CNumerics** [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      numerics_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
  }
  
  /*--- Solver definition for the template problem ---*/
  if (template_solver) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][TEMPLATE_SOL][CONV_TERM] = new CConvective_Template(nDim, nVar_Template, config);
        break;
      default : cout << "Convective scheme not implemented (template_solver)." << endl; exit(1); break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Template()) {
      case AVG_GRAD : case AVG_GRAD_CORRECTED : case GALERKIN :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][TEMPLATE_SOL][VISC_TERM] = new CViscous_Template(nDim, nVar_Template, config);
        break;
      default : cout << "Viscous scheme not implemented." << endl; exit(1); break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Template()) {
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM] = new CSource_Template(nDim, nVar_Template, config);
        break;
      default : cout << "Source term not implemented." << endl; exit(1); break;
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
      numerics_container[iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM] = new CConvective_Template(nDim, nVar_Template, config);
    }
    
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(1);
        break;
        
      case SPACE_CENTERED :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim,nVar_Flow, config); break;
            case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim,nVar_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(1); break;
          }
          
          if (!config->GetLowFidelitySim()) {
            for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          }
          else {
            numerics_container[MESH_1][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config);
            for (iMGlevel = 2; iMGlevel <= config->GetMGLevels(); iMGlevel++)
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          }
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJSTArtComp_Flow(nDim, nVar_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(1); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim,nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
          
        }
        if (freesurface) {
          /*--- FreeSurface flow, use artificial compressibility method ---*/
          cout << "Centered scheme not implemented." << endl; exit(1);
        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE_1ST : case ROE_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case AUSM_1ST : case AUSM_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case TURKEL_1ST : case TURKEL_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case HLLC_1ST : case HLLC_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case MSW_1ST : case MSW_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
          }
          
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE_1ST : case ROE_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_Flow(nDim, nVar_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
          }
        }
        if (freesurface) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE_1ST : case ROE_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwArtComp_FreeSurf_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwArtComp_FreeSurf_Flow(nDim, nVar_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
          }
        }
        
        break;
        
      default :
        cout << "Convective scheme not implemented (euler and ns)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Flow()) {
      case NONE :
        break;
      case AVG_GRAD :
        if (compressible) {
          /*--- Compressible flow ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
          }
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
          }
        }
        if (freesurface) {
          /*--- Freesurface flow, use artificial compressibility method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
          }
        }
        break;
      case AVG_GRAD_CORRECTED :
        if (compressible) {
          /*--- Compressible flow ---*/
          numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
        }
        if (freesurface) {
          /*--- Freesurface flow, use artificial compressibility method ---*/
          numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
        }
        break;
      case GALERKIN :
        cout << "Galerkin viscous scheme not implemented." << endl; exit(1); exit(1);
        break;
      default :
        cout << "Numerical viscous scheme not recognized." << endl; exit(1); exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Flow()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          
          if (config->GetRotating_Frame() == YES)
            numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
          else if (config->GetAxisymmetric() == YES)
            numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_Flow(nDim,nVar_Flow, config);
          else if (config->GetGravityForce() == YES)
            numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceGravity(nDim, nVar_Flow, config);
          else if (config->GetWind_Gust() == YES)
            numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceWindGust(nDim, nVar_Flow, config);
          else
            numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
          
          numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
        }
        
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
    
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((tne2_euler) || (tne2_ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_TNE2()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(1);
        break;
        
      case SPACE_CENTERED :
        /*--- Compressible two-temperature flow ---*/
        switch (config->GetKind_Centered_TNE2()) {
          case NO_CENTERED : cout << "No centered scheme." << endl; break;
          case LAX :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM]       = new CCentLax_TNE2(nDim,nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_TNE2(nDim, nVar_TNE2,  nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
            }
            break;
          default : cout << "Centered scheme not implemented." << endl; exit(1); break;
        }
        break;
        
      case SPACE_UPWIND :
        /*--- Compressible two-temperature flow ---*/
        switch (config->GetKind_Upwind_TNE2()) {
          case NO_UPWIND : cout << "No upwind scheme." << endl; break;
          case ROE_1ST : case ROE_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwRoe_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_TNE2(nDim, nVar_TNE2,  nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
            }
            break;
            
          case MSW_1ST : case MSW_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwMSW_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwMSW_TNE2(nDim, nVar_TNE2,  nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
            }
            break;
            
          case AUSM_1ST : case AUSM_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwAUSM_TNE2(nDim, nVar_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwAUSM_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          case AUSMPWPLUS_1ST : case AUSMPWPLUS_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwAUSMPWplus_TNE2(nDim, nVar_TNE2, config);
              numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwAUSMPWplus_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          case TURKEL_1ST : case TURKEL_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwRoe_Turkel_TNE2(nDim, nVar_TNE2, config);
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_Turkel_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          case HLLC_1ST : case HLLC_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_TERM] = new CUpwHLLC_TNE2(nDim, nVar_TNE2, config);
              //                numerics_container[iMGlevel][TNE2_SOL][CONV_BOUND_TERM] = new CUpwHLLC_TNE2(nDim, nVar_TNE2, config);
            }
            break;
            
          default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
        }
        break;
        
      default :
        cout << "Convective scheme not implemented (TNE2)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_TNE2()) {
      case NONE :
        break;
      case AVG_GRAD :
        /*--- Compressible TNE2 ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][TNE2_SOL][VISC_TERM]       = new CAvgGrad_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
          numerics_container[iMGlevel][TNE2_SOL][VISC_BOUND_TERM] = new CAvgGrad_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
        }
        break;
      case AVG_GRAD_CORRECTED :
        /*--- Compressible TNE2 ---*/
        numerics_container[MESH_0][TNE2_SOL][VISC_TERM] = new CAvgGradCorrected_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
        for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][TNE2_SOL][VISC_TERM] = new CAvgGradCorrected_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
        }
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][TNE2_SOL][VISC_BOUND_TERM] = new CAvgGradCorrected_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
        }
        break;
      case GALERKIN :
        cout << "Galerkin viscous scheme not implemented." << endl; exit(1); exit(1);
        break;
      default :
        cout << "Numerical viscous scheme not recognized." << endl; exit(1); exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_TNE2()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][TNE2_SOL][SOURCE_FIRST_TERM] = new CSource_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
        }
        
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
    
  }
  
  /*--- Solver definition for the turbulent model problem ---*/
  if (turbulent) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          else if (machine_learning) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbML(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        }
        break;
      default :
        cout << "Convective scheme not implemented (turbulent)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Turb()) {
      case NONE :
        break;
      case AVG_GRAD :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
          else if (machine_learning) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbML(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, config);
        }
        break;
      case AVG_GRAD_CORRECTED :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSA(nDim, nVar_Turb, config);
          else if (machine_learning) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbML(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSST(nDim, nVar_Turb, constants, config);
        }
        break;
      case GALERKIN :
        cout << "Viscous scheme not implemented." << endl;
        exit(1); break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Turb()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
          else if (machine_learning) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbML(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, constants, config);
          numerics_container[iMGlevel][TURB_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Turb, config);
        }
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
      if (spalart_allmaras) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
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
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          numerics_container[iMGlevel][TRANS_SOL][CONV_TERM] = new CUpwSca_TransLM(nDim, nVar_Trans, config);
        }
        break;
      default :
        cout << "Convective scheme not implemented (transition)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Turb()) {
      case NONE :
        break;
      case AVG_GRAD :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          numerics_container[iMGlevel][TRANS_SOL][VISC_TERM] = new CAvgGrad_TransLM(nDim, nVar_Trans, config);
        }
        break;
      case AVG_GRAD_CORRECTED :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          numerics_container[iMGlevel][TRANS_SOL][VISC_TERM] = new CAvgGradCorrected_TransLM(nDim, nVar_Trans, config);
        }
        break;
      case GALERKIN :
        cout << "Viscous scheme not implemented." << endl;
        exit(1); break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Turb()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][TRANS_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TransLM(nDim, nVar_Trans, config);
          numerics_container[iMGlevel][TRANS_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Trans, config);
        }
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
      numerics_container[iMGlevel][TRANS_SOL][CONV_BOUND_TERM] = new CUpwLin_TransLM(nDim, nVar_Trans, config);
    }
  }
  
  /*--- Solver definition for the poisson potential problem ---*/
  if (poisson) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Poisson()) {
      case GALERKIN :
        numerics_container[MESH_0][POISSON_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Poisson, config);
        break;
      default : cout << "Viscous scheme not implemented." << endl; exit(1); break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Poisson()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        numerics_container[MESH_0][POISSON_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Poisson, config);
        numerics_container[MESH_0][POISSON_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Poisson, config);
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
  }
  
  /*--- Solver definition for the poisson potential problem ---*/
  if (heat) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Heat()) {
      case GALERKIN :
        numerics_container[MESH_0][HEAT_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Heat, config);
        break;
      default : cout << "Viscous scheme not implemented." << endl; exit(1); break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Heat()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        numerics_container[MESH_0][HEAT_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Heat, config);
        numerics_container[MESH_0][HEAT_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Heat, config);
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
  }
  
  /*--- Solver definition for the flow adjoint problem ---*/
  if (adj_euler || adj_ns) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(1);
        break;
      case SPACE_CENTERED :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default : cout << "Centered scheme not implemented." << endl; exit(1); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Centered_AdjFlow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : cout << "Centered scheme not implemented." << endl; exit(1); break;
            default : cout << "Centered scheme not implemented." << endl; exit(1); break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CCentLaxArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          
        }
        if (freesurface) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          cout << "Centered scheme not implemented." << endl; exit(1);
        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE_1ST : case ROE_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
          }
        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE_1ST : case ROE_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
          }
        }
        if (freesurface) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_AdjFlow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE_1ST : case ROE_2ND :
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics_container[iMGlevel][ADJFLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
          }
        }
        break;
        
      default :
        cout << "Convective scheme not implemented (adj_euler and adj_ns)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_AdjFlow()) {
      case NONE :
        break;
      case AVG_GRAD :
        if (compressible) {
          /*--- Compressible flow ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
          }
        }
        if (incompressible || freesurface) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
        }
        
        break;
      case AVG_GRAD_CORRECTED :
        if (compressible) {
          /*--- Compressible flow ---*/
          numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
          numerics_container[MESH_0][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
          }
        }
        if (incompressible || freesurface) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          numerics_container[MESH_0][ADJFLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJFLOW_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_Flow, config);
        }
        
        break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_AdjFlow()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          
          /*--- Note that RANS is incompatible with Axisymmetric or Rotational (Fix it!) ---*/
          if ((adj_ns) && (!incompressible)) {
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
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
  }
  
  /*--- Solver definition for the flow adjoint problem ---*/
  if (adj_tne2_euler || adj_tne2_ns) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTNE2()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(1);
        break;
      case SPACE_CENTERED :
        switch (config->GetKind_Centered_AdjTNE2()) {
          case NO_CENTERED : cout << "No centered scheme." << endl; break;
          case LAX : numerics_container[MESH_0][ADJTNE2_SOL][CONV_TERM] = new CCentLax_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config); break;
          case JST : numerics_container[MESH_0][ADJTNE2_SOL][CONV_TERM] = new CCentJST_AdjTNE2(nDim, nVar_Adj_TNE2, config); break;
          default : cout << "Centered scheme not implemented." << endl; exit(1); break;
        }
        for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][ADJTNE2_SOL][CONV_TERM] = new CCentLax_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][ADJTNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
        break;
      case SPACE_UPWIND :
        switch (config->GetKind_Upwind_AdjTNE2()) {
          case NO_UPWIND : cout << "No upwind scheme." << endl; break;
          case ROE_1ST : case ROE_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_TERM] = new CUpwRoe_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_BOUND_TERM] = new CUpwRoe_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
            }
            break;
          case SW_1ST : case SW_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_TERM] = new CUpwSW_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
              numerics_container[iMGlevel][ADJTNE2_SOL][CONV_BOUND_TERM] = new CUpwSW_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
            }
            break;
          default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
        }
        break;
        
      default :
        cout << "Convective scheme not implemented (adj_tne2_euler and adj_tne2_ns)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_AdjTNE2()) {
      case NONE :
        break;
      case AVG_GRAD :
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJTNE2_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_TNE2, config);
        }
        else {
          /*--- Compressible flow ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][ADJTNE2_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_TNE2, config);
            numerics_container[iMGlevel][ADJTNE2_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_TNE2, config);
          }
        }
        break;
      case AVG_GRAD_CORRECTED :
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          numerics_container[MESH_0][ADJTNE2_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_AdjFlow(nDim, nVar_Adj_TNE2, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][ADJTNE2_SOL][VISC_TERM] = new CAvgGradArtComp_AdjFlow(nDim, nVar_Adj_TNE2, config);
        }
        else {
          /*--- Compressible flow ---*/
          numerics_container[MESH_0][ADJTNE2_SOL][VISC_TERM] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_TNE2, config);
          numerics_container[MESH_0][ADJTNE2_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_TNE2, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][ADJTNE2_SOL][VISC_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
            numerics_container[iMGlevel][ADJTNE2_SOL][VISC_BOUND_TERM] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
          }
          
        }
        break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_AdjTNE2()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][ADJTNE2_SOL][SOURCE_FIRST_TERM] = new CSource_TNE2(nDim, nVar_TNE2, nPrimVar_TNE2, nPrimVarGrad_TNE2, config);
          
          numerics_container[iMGlevel][ADJTNE2_SOL][SOURCE_SECOND_TERM] = new CSource_AdjTNE2(nDim, nVar_Adj_TNE2, nPrimVar_Adj_TNE2, nPrimVarGrad_Adj_TNE2, config);
        }
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
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
          default : cout << "Centered scheme not implemented." << endl; exit(1); break;
        }
        for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][LINFLOW_SOL][CONV_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
        break;
      default :
        cout << "Convective scheme not implemented (lin_euler)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
      numerics_container[iMGlevel][LINFLOW_SOL][CONV_BOUND_TERM] = new CCentLax_LinFlow(nDim, nVar_Lin_Flow, config);
  }
  
  /*--- Solver definition for the turbulent adjoint problem ---*/
  if (adj_turb) {
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTurb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          if (spalart_allmaras) {
            numerics_container[iMGlevel][ADJTURB_SOL][CONV_TERM] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
          }
          else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(1);}
        break;
      default :
        cout << "Convective scheme not implemented (adj_turb)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_AdjTurb()) {
      case NONE :
        break;
      case AVG_GRAD :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          if (spalart_allmaras){
            numerics_container[iMGlevel][ADJTURB_SOL][VISC_TERM] = new CAvgGrad_AdjTurb(nDim, nVar_Adj_Turb, config);
          }
          else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(1);}
        break;
      case AVG_GRAD_CORRECTED :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          if (spalart_allmaras){
            numerics_container[iMGlevel][ADJTURB_SOL][VISC_TERM] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);
          }
          else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(1);}
        break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_AdjTurb()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          if (spalart_allmaras) {
            numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
            numerics_container[iMGlevel][ADJTURB_SOL][SOURCE_SECOND_TERM] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
          }
          else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(1);}
        }
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics_container[iMGlevel][ADJTURB_SOL][CONV_BOUND_TERM] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);
      else if (menter_sst) {cout << "Adjoint SST turbulence model not implemented." << endl; exit(1);}
    }
    
  }
  
  /*--- Solver definition for the wave problem ---*/
  if (wave) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Wave()) {
      case NONE :
        break;
      case AVG_GRAD :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
      case AVG_GRAD_CORRECTED :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
      case GALERKIN :
        numerics_container[MESH_0][WAVE_SOL][VISC_TERM] = new CGalerkin_Flow(nDim, nVar_Wave, config);
        break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Wave()) {
      case NONE : break;
      case PIECEWISE_CONSTANT :
        break;
      default : break;
    }
  }
  
  /*--- Solver definition for the FEA problem ---*/
  if (fea) {
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_FEA()) {
      case NONE :
        break;
      case AVG_GRAD :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
      case AVG_GRAD_CORRECTED :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
      case GALERKIN :
        numerics_container[MESH_0][FEA_SOL][VISC_TERM] = new CGalerkin_FEA(nDim, nVar_Wave, config);
        break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Wave()) {
      case NONE : break;
      case PIECEWISE_CONSTANT :
        break;
      default : break;
    }
  }
  
}
