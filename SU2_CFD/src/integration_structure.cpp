/*!
 * \file integration_structure.cpp
 * \brief This subroutine includes the space and time integration structure
 * \author F. Palacios, T. Economon
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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

#include "/extra2/JKA/Cprograms/libm3l/Source/libm3l.h"
#include "/extra2/JKA/Cprograms/lsipdx/Source/lsipdx.h"
#include "/extra2/JKA/Cprograms/lsipdx/Source/socket_op.h"

void comm(void);

CIntegration::CIntegration(CConfig *config) {
	Cauchy_Value = 0;
	Cauchy_Func = 0;
	Old_Func = 0;
	New_Func = 0;
	Cauchy_Counter = 0;
	Convergence = false;
	Convergence_FullMG = false;
	Cauchy_Serie = new double [config->GetCauchy_Elems()+1];
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
  unsigned short iMarker;
  
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
  }
  
  
  /*--- Compute viscous residuals ---*/
  
  solver_container[MainSolver]->Viscous_Residual(geometry, solver_container, numerics[VISC_TERM], config, iMesh, iRKStep);
  
  
  /*--- Compute source term residuals ---*/
  
  solver_container[MainSolver]->Source_Residual(geometry, solver_container, numerics[SOURCE_FIRST_TERM], numerics[SOURCE_SECOND_TERM], config, iMesh);
  
  /*--- Add viscous and convective residuals, and compute the Dual Time Source term ---*/
  
  if (dual_time)
    solver_container[MainSolver]->SetResidual_DualTime(geometry, solver_container, config, iRKStep, iMesh, RunTime_EqSystem);
  
  /*--- Boundary conditions that depend on other boundaries (they require MPI sincronization)---*/
  
  solver_container[MainSolver]->BC_ActDisk_Boundary(geometry, solver_container, numerics[CONV_BOUND_TERM], config);
  
  solver_container[MainSolver]->BC_Interface_Boundary(geometry, solver_container, numerics[CONV_BOUND_TERM], config);

  solver_container[MainSolver]->BC_NearField_Boundary(geometry, solver_container, numerics[CONV_BOUND_TERM], config);

  
  /*--- Weak boundary conditions ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case EULER_WALL:
        solver_container[MainSolver]->BC_Euler_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case INLET_FLOW:
        solver_container[MainSolver]->BC_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
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
      case RIEMANN_BOUNDARY:
      	if (MainSolver == FLOW_SOL)
      		solver_container[MainSolver]->BC_Riemann(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
      	else if (MainSolver == TURB_SOL && config->GetKind_Data_Riemann(config->GetMarker_All_TagBound(iMarker)) == TOTAL_CONDITIONS_PT)
      		solver_container[MainSolver]->BC_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
      	else if (MainSolver == TURB_SOL && config->GetKind_Data_Riemann(config->GetMarker_All_TagBound(iMarker)) == STATIC_PRESSURE)
      		solver_container[MainSolver]->BC_Outlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
      	break;
      case FAR_FIELD:
        solver_container[MainSolver]->BC_Far_Field(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case SYMMETRY_PLANE:
        solver_container[MainSolver]->BC_Sym_Plane(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ENGINE_EXHAUST:
        solver_container[MainSolver]->BC_Engine_Exhaust(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ENGINE_INFLOW:
        solver_container[MainSolver]->BC_Engine_Inflow(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ENGINE_BLEED:
        solver_container[MainSolver]->BC_Engine_Bleed(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
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
      case FLOWLOAD_BOUNDARY:
        solver_container[MainSolver]->BC_Flow_Load(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case LOAD_BOUNDARY:
        solver_container[MainSolver]->BC_Normal_Load(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
      case PRESSURE_BOUNDARY:
        solver_container[MainSolver]->BC_Pressure(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
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
      case ISOTHERMAL_NONCATALYTIC:
        solver_container[MainSolver]->BC_IsothermalNonCatalytic_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case ISOTHERMAL_CATALYTIC:
        solver_container[MainSolver]->BC_IsothermalCatalytic_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case HEAT_FLUX:
        solver_container[MainSolver]->BC_HeatFlux_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case HEAT_FLUX_NONCATALYTIC:
        solver_container[MainSolver]->BC_HeatFluxNonCatalytic_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case HEAT_FLUX_CATALYTIC:
        solver_container[MainSolver]->BC_HeatFluxCatalytic_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
      case DIRICHLET:
        solver_container[MainSolver]->BC_Dirichlet(geometry, solver_container, config, iMarker);
        break;
      case CUSTOM_BOUNDARY:
        solver_container[MainSolver]->BC_Custom(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
        break;
    }
  
}

void CIntegration::Adjoint_Setup(CGeometry ***geometry, CSolver ****solver_container, CConfig **config,
                                 unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) {
  
	unsigned short iMGLevel;
  
	if ( ( ((RunTime_EqSystem == RUNTIME_ADJFLOW_SYS) ||
          (RunTime_EqSystem == RUNTIME_LINFLOW_SYS)) && (Iteration == 0) ) ) {
		for (iMGLevel = 0; iMGLevel <= config[iZone]->GetnMGLevels(); iMGLevel++) {
      
			/*--- Set the time step in all the MG levels ---*/
			solver_container[iZone][iMGLevel][FLOW_SOL]->SetTime_Step(geometry[iZone][iMGLevel], solver_container[iZone][iMGLevel], config[iZone], iMGLevel, Iteration);
      
			/*--- Set the force coefficients ---*/
			solver_container[iZone][iMGLevel][FLOW_SOL]->SetTotal_CDrag(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag());
			solver_container[iZone][iMGLevel][FLOW_SOL]->SetTotal_CLift(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CLift());
			solver_container[iZone][iMGLevel][FLOW_SOL]->SetTotal_CT(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CT());
			solver_container[iZone][iMGLevel][FLOW_SOL]->SetTotal_CQ(solver_container[iZone][MESH_0][FLOW_SOL]->GetTotal_CQ());
      
			/*--- Restrict solution and gradients to the coarse levels ---*/
			if (iMGLevel != config[iZone]->GetnMGLevels()) {
				SetRestricted_Solution(RUNTIME_FLOW_SYS, solver_container[iZone][iMGLevel][FLOW_SOL], solver_container[iZone][iMGLevel+1][FLOW_SOL],
                               geometry[iZone][iMGLevel], geometry[iZone][iMGLevel+1], config[iZone]);
				SetRestricted_Gradient(RUNTIME_FLOW_SYS, solver_container[iZone][iMGLevel][FLOW_SOL], solver_container[iZone][iMGLevel+1][FLOW_SOL],
                               geometry[iZone][iMGLevel], geometry[iZone][iMGLevel+1], config[iZone]);
			}
      
		}
  } else if ((RunTime_EqSystem == RUNTIME_ADJTNE2_SYS) && (Iteration == 0)) {
    for (iMGLevel = 0; iMGLevel <= config[iZone]->GetnMGLevels(); iMGLevel++) {
      
			/*--- Set the time step in all the MG levels ---*/
			solver_container[iZone][iMGLevel][TNE2_SOL]->SetTime_Step(geometry[iZone][iMGLevel],
                                                                solver_container[iZone][iMGLevel],
                                                                config[iZone], iMGLevel, Iteration);
      
			/*--- Set the force coefficients ---*/
			solver_container[iZone][iMGLevel][TNE2_SOL]->SetTotal_CDrag(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CDrag());
			solver_container[iZone][iMGLevel][TNE2_SOL]->SetTotal_CLift(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CLift());
			solver_container[iZone][iMGLevel][TNE2_SOL]->SetTotal_CT(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CT());
			solver_container[iZone][iMGLevel][TNE2_SOL]->SetTotal_CQ(solver_container[iZone][MESH_0][TNE2_SOL]->GetTotal_CQ());
      
			/*--- Restrict solution and gradients to the coarse levels ---*/
			if (iMGLevel != config[iZone]->GetnMGLevels()) {
				SetRestricted_Solution(RUNTIME_TNE2_SYS, solver_container[iZone][iMGLevel][TNE2_SOL], solver_container[iZone][iMGLevel+1][TNE2_SOL],
                               geometry[iZone][iMGLevel], geometry[iZone][iMGLevel+1], config[iZone]);
				SetRestricted_Gradient(RUNTIME_TNE2_SYS, solver_container[iZone][iMGLevel][TNE2_SOL], solver_container[iZone][iMGLevel+1][TNE2_SOL],
                               geometry[iZone][iMGLevel], geometry[iZone][iMGLevel+1], config[iZone]);
			}
      
		}
  }
}

void CIntegration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                    unsigned short RunTime_EqSystem, unsigned long Iteration) {
	unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
  
  /*--- Perform the time integration ---*/
  
  switch (config->GetKind_TimeIntScheme()) {
    case (RUNGE_KUTTA_EXPLICIT):
      solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iRKStep);
      break;
    case (EULER_EXPLICIT):
      solver_container[MainSolver]->ExplicitEuler_Iteration(geometry, solver_container, config);
      break;
    case (EULER_IMPLICIT):
      solver_container[MainSolver]->ImplicitEuler_Iteration(geometry, solver_container, config);
      break;
  }
  
}

void CIntegration::Convergence_Monitoring(CGeometry *geometry, CConfig *config, unsigned long Iteration,
                                          double monitor, unsigned short iMesh) {
  
  unsigned short iCounter;
  int rank = MASTER_NODE;
  
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
    
#ifdef HAVE_MPI
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    
    bool Already_Converged = Convergence;
    
    /*--- Cauchi based convergence criteria ---*/
    
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
    MPI_Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
    
    /*-- Compute global convergence criteria in the master node --*/
    
    sbuf_conv[0] = 0;
    if (rank == MASTER_NODE) {
      if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
      else sbuf_conv[0] = 0;
    }
    
    MPI_Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (sbuf_conv[0] == 1) { Convergence = true; Convergence_FullMG = true; }
    else { Convergence = false; Convergence_FullMG = false; }
    
    delete [] sbuf_conv;
    delete [] rbuf_conv;
    
#endif
    
    /*--- Stop the simulation in case a nan appears, do not save the solution ---*/
    
    if (monitor != monitor) {
      if (rank == MASTER_NODE)
      cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
#ifndef HAVE_MPI
      exit(EXIT_DIVERGENCE);
#else
      MPI_Abort(MPI_COMM_WORLD,1);
#endif
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
    double plunge, pitch, *plunge_all = NULL, *pitch_all = NULL;
    unsigned short iMarker, iMarker_Monitoring;
    unsigned long iProcessor, owner, *owner_all = NULL;
    
    string Marker_Tag, Monitoring_Tag;
	int rank, nProcessor;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

    /*--- Only if mater node allocate memory ---*/
    if (rank == MASTER_NODE) {
      plunge_all = new double[nProcessor];
      pitch_all  = new double[nProcessor];
      owner_all  = new unsigned long[nProcessor];
    }
    
    /*--- Find marker and give it's plunge and pitch coordinate to the master node ---*/
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
        Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag) { owner = 1; break;
        } else {
          owner = 0;
        }
        
      }
      plunge = config->GetAeroelastic_plunge(iMarker_Monitoring);
      pitch  = config->GetAeroelastic_pitch(iMarker_Monitoring);
      
      /*--- Gather the data on the master node. ---*/
      MPI_Gather(&plunge, 1, MPI_DOUBLE, plunge_all, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      MPI_Gather(&pitch, 1, MPI_DOUBLE, pitch_all, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      MPI_Gather(&owner, 1, MPI_UNSIGNED_LONG, owner_all, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
      
      /*--- Set plunge and pitch on the master node ---*/
      if (rank == MASTER_NODE) {
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
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



void comm(){
	node_t *Gnode=NULL, *TmpNode = NULL, *FoundNode = NULL;
	size_t dim[1], i, tot_dim;

	lmchar_t hostname[80], channel_name[80];
	lmint_t sockfd, portno;

	lmchar_t *name ="CFD2SIM";
	lmchar_t *name1="SIM2CFD";
	
	lmdouble_t Forces_moments[6], Angles[3], RotCenter[3], TransVec[3];
	
	lmchar_t name_i[80], name_o[80];

	lmdouble_t *tmpfloat;
	client_fce_struct_t InpPar, *PInpPar;
	opts_t *Popts_1, opts, opts_1, *Popts;
	find_t *SFounds;
	
	PInpPar = &InpPar;
	PInpPar->channel_name = name;
	PInpPar->SR_MODE = 'S';
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("socket_edge2simulink: wrong client mode");

	Popts   = &opts;
	Popts_1 = &opts_1;
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	m3l_set_Find(&Popts);
/*
 * create data structure which will be sent
 */
	if(  (Gnode = m3l_Mklist("CFD_2_SIM", "DIR", 0, 0, (node_t **)NULL, (const char *)NULL, (const char *)NULL, (char *)NULL)) == 0)
		Perror("socket_edge2simulink: m3l_Mklist");
	
	dim[0] = 6;
/*
 * store global forces moments
 */
	if(  (TmpNode = m3l_Mklist("ForcesMoments", "D", 1, dim, &Gnode, "/CFD_2_SIM", "./", "--no_malloc", (char *)NULL)) == 0)
		Error("socket_edge2simulink: m3l_Mklist");
	TmpNode->data.df = Forces_moments;
/*
 * add time
 */
// 	dim[0] = 1;
// 	if(  (TmpNode = m3l_Mklist("Time", "D", 1, dim, &Gnode, "/CFD_2_SIM", "./", (char *)NULL)) == 0)
// 		Error("socket_edge2simulink: m3l_Mklist");
// 	TmpNode->data.df[0] = *ttime;
/*
 * open socket
 */
	if( (sockfd = open_connection_to_server(hostname, portno, PInpPar, Popts_1)) < 1)
		Error("socket_edge2simulink: Error when opening socket");
/*
 * send data 
 */
	if ( client_sender(Gnode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL) !=1 )
		Error("socket_edge2simulink: client_sender()");
#pragma omp critical
{
	if( close(sockfd) == -1)
		Perror("socket_edge2simulink: close");
/*
 * free borrowed memory
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_edge2simulink: m3l_Umount");
/*
 * receive data 
 */

		
	PInpPar = &InpPar;
	PInpPar->channel_name = name1;
	PInpPar->SR_MODE = 'R';
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("socket_edge2simulink: wrong client mode");

	Popts   = &opts;
	Popts_1 = &opts_1;
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	m3l_set_Find(&Popts);

	if( (sockfd = open_connection_to_server(hostname, portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");
}
	if ( (Gnode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL)) == NULL)
		Error("socket_edge2simulink: client_receiver()");

#pragma omp critical
{
/*
 * close socket 
 */
	if( close(sockfd) == -1)
		Perror("socket_edge2simulink: close");
/*
 * find Angles - rotation matrix and copy the values to Edge allocated memory
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/Angles", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_edge2simulink: More then one Angles data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_edge2simulink: Did not find 1st data pointer");
// 		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 9)
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_edge2simulink: Wrong dimensions of Angles array");
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_edge2simulink: Did not find Angles data pointer");

		for (i=0; i<tot_dim; i++)
			Angles[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_edge2simulink: Angles not found\n");
	}
/*
 * find center of rotation
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/RotCenter", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_edge2simulink: More then one RotCenter data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_edge2simulink: Did not find 1st data pointer");
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_edge2simulink: Wrong dimensions of RotCenter array");
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_edge2simulink: Did not find RotCenter data pointer");

		for (i=0; i<tot_dim; i++)
			RotCenter[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_edge2simulink: RotCenter not found\n");
	}
/*
 * find center of translation
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/TransVec", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_edge2simulink: More then one TransVec data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_edge2simulink: Did not find 1st data pointer");
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_edge2simulink: Wrong dimensions of TransVec array");
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_edge2simulink: Did not find TransVec data pointer");

		for (i=0; i<tot_dim; i++)
			TransVec[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_edge2simulink: TransVec not found\n");
	}
/*
 * free borrowed memory
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_edge2simulink: m3l_Umount");
}
	
}
