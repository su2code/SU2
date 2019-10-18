/*!
 * \file integration_time.cpp
 * \brief Time dependent numerical methods
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

CMultiGridIntegration::CMultiGridIntegration(CConfig *config) : CIntegration(config) {}

CMultiGridIntegration::~CMultiGridIntegration(void) { }

void CMultiGridIntegration::MultiGrid_Iteration(CGeometry ****geometry,
                                                CSolver *****solver_container,
                                                CNumerics ******numerics_container,
                                                CConfig **config,
                                                unsigned short RunTime_EqSystem,
                                                unsigned long Iteration,
                                                unsigned short iZone,
                                                unsigned short iInst) {
  unsigned short FinestMesh, iMGLevel = 0;
  su2double monitor = 1.0;
  bool FullMG = false;
  
  const bool restart = (config[iZone]->GetRestart() || config[iZone]->GetRestart_Flow());
  const bool startup_multigrid = ((config[iZone]->GetRestart_Flow())     &&
                                  (RunTime_EqSystem == RUNTIME_FLOW_SYS) &&
                                  (Iteration == 0));
  const bool direct = ((config[iZone]->GetKind_Solver() == EULER)                         ||
                       (config[iZone]->GetKind_Solver() == NAVIER_STOKES)                 ||
                       (config[iZone]->GetKind_Solver() == RANS)                          ||
                       (config[iZone]->GetKind_Solver() == FEM_EULER)                     ||
                       (config[iZone]->GetKind_Solver() == FEM_NAVIER_STOKES)             ||
                       (config[iZone]->GetKind_Solver() == FEM_RANS)                      ||
                       (config[iZone]->GetKind_Solver() == FEM_LES)                       ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_EULER)                ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)        ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_FEM_EULER)            ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_FEM_NS)               ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_RANS));
  const bool turb = ((config[iZone]->GetKind_Solver() == RANS)                         ||
					 (config[iZone]->GetKind_Solver() == FEM_RANS)                     ||
					 (config[iZone]->GetKind_Solver() == DISC_ADJ_RANS));
  const unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  unsigned short RecursiveParam = config[iZone]->GetMGCycle();
  bool pressure_based = (config[iZone]->GetKind_Incomp_System() == PRESSURE_BASED);
  if (config[iZone]->GetMGCycle() == FULLMG_CYCLE) {
    RecursiveParam = V_CYCLE;
    FullMG = true;
  }
  
  /*--- If restart, update multigrid levels at the first multigrid iteration ---*/
	/*-- Since the restart takes care of this I dont think is required, but we should check after the new restart routines are added ---*/
  
  if ((restart && (Iteration == config[iZone]->GetnStartUpIter())) || startup_multigrid)
  {
    for (iMGLevel = 0; iMGLevel < config[iZone]->GetnMGLevels(); iMGLevel++) {
      
      SetRestricted_Solution(RunTime_EqSystem, solver_container[iZone][iInst][iMGLevel][SolContainer_Position],
                             solver_container[iZone][iInst][iMGLevel+1][SolContainer_Position],
                             geometry[iZone][iInst][iMGLevel], geometry[iZone][iInst][iMGLevel+1], config[iZone]);
      
      if (turb) SetRestricted_Solution(RUNTIME_TURB_SYS, solver_container[iZone][iInst][iMGLevel][SolContainer_Position], solver_container[iZone][iInst][iMGLevel+1][SolContainer_Position], geometry[iZone][iInst][iMGLevel], geometry[iZone][iInst][iMGLevel+1], config[iZone]);
      
      if (turb) SetRestricted_EddyVisc(RUNTIME_TURB_SYS, solver_container[iZone][iInst][iMGLevel][SolContainer_Position], solver_container[iZone][iInst][iMGLevel+1][SolContainer_Position], geometry[iZone][iInst][iMGLevel], geometry[iZone][iInst][iMGLevel+1], config[iZone]);
      
    }
  }
	
  /*--- Full multigrid strategy and start up with fine grid only works with the direct problem ---*/

  if (!config[iZone]->GetRestart() && FullMG && direct && ( Convergence_FullMG && (config[iZone]->GetFinestMesh() != MESH_0 ))) {
    SetProlongated_Solution(RunTime_EqSystem, solver_container[iZone][iInst][config[iZone]->GetFinestMesh()-1][SolContainer_Position],
                            solver_container[iZone][iInst][config[iZone]->GetFinestMesh()][SolContainer_Position],
                            geometry[iZone][iInst][config[iZone]->GetFinestMesh()-1], geometry[iZone][iInst][config[iZone]->GetFinestMesh()],
                            config[iZone]);
    config[iZone]->SubtractFinestMesh();
  }

  /*--- Set the current finest grid (full multigrid strategy) ---*/
  
  FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Perform the Full Approximation Scheme multigrid ---*/
  MultiGrid_Cycle(geometry, solver_container, numerics_container, config,
                  FinestMesh, RecursiveParam, RunTime_EqSystem,
                  Iteration, iZone, iInst);

  /*--- Computes primitive variables and gradients in the finest mesh (useful for the next solver (turbulence) and output ---*/

   solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][MESH_0],
                                                                         solver_container[iZone][iInst][MESH_0], config[iZone],
                                                                         MESH_0, NO_RK_ITER, RunTime_EqSystem, true);
  
  /*--- Compute non-dimensional parameters and the convergence monitor ---*/
  
  NonDimensional_Parameters(geometry[iZone][iInst], solver_container[iZone][iInst],
                            numerics_container[iZone][iInst], config[iZone],
                            FinestMesh, RunTime_EqSystem, Iteration, &monitor);
  
  /*--- Convergence strategy ---*/
  
  Convergence_Monitoring(geometry[iZone][iInst][FinestMesh], config[iZone], Iteration, monitor, FinestMesh);

}

void CMultiGridIntegration::MultiGrid_Cycle(CGeometry ****geometry,
                                            CSolver *****solver_container,
                                            CNumerics ******numerics_container,
                                            CConfig **config,
                                            unsigned short iMesh,
                                            unsigned short RecursiveParam,
                                            unsigned short RunTime_EqSystem,
                                            unsigned long Iteration,
                                            unsigned short iZone,
                                            unsigned short iInst) {
  
  unsigned short iPreSmooth, iPostSmooth, iRKStep, iRKLimit = 1;
  
  bool startup_multigrid = (config[iZone]->GetRestart_Flow() && (RunTime_EqSystem == RUNTIME_FLOW_SYS) && (Iteration == 0));
  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  
  /*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/
  
  for (iPreSmooth = 0; iPreSmooth < config[iZone]->GetMG_PreSmooth(iMesh); iPreSmooth++) {
    
    switch (config[iZone]->GetKind_TimeIntScheme()) {
      case RUNGE_KUTTA_EXPLICIT: iRKLimit = config[iZone]->GetnRKStep(); break;
      case CLASSICAL_RK4_EXPLICIT: iRKLimit = 4; break;
      case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break; }

    /*--- Time and space integration ---*/
    
    for (iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {
      
      /*--- Send-Receive boundary conditions, and preprocessing ---*/
      
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, iRKStep, RunTime_EqSystem, false);
      
      if (iRKStep == 0) {
        
        /*--- Set the old solution ---*/
        
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][iMesh]);

        if (config[iZone]->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT)
          solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_NewSolution(geometry[iZone][iInst][iMesh]);

        /*--- Compute time step, max eigenvalue, and integration scheme (steady and unsteady problems) ---*/
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, Iteration);
        
        /*--- Restrict the solution and gradient for the adjoint problem ---*/
        
        Adjoint_Setup(geometry, solver_container, config, RunTime_EqSystem, Iteration, iZone);
        
      }
      
      /*--- Space integration ---*/      
      Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], numerics_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone], iMesh, iRKStep, RunTime_EqSystem);
      
      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      Time_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iRKStep, RunTime_EqSystem, Iteration);
      
      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
      
    }
    
  }
  
  /*--- Compute Forcing Term $P_(k+1) = I^(k+1)_k(P_k+F_k(u_k))-F_(k+1)(I^(k+1)_k u_k)$ and update solution for multigrid ---*/
  
  if ( (iMesh < config[iZone]->GetnMGLevels() && ((Iteration >= config[iZone]->GetnStartUpIter()) || startup_multigrid)) ) {
    /*--- Compute $r_k = P_k + F_k(u_k)$ ---*/
    
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, NO_RK_ITER, RunTime_EqSystem, false);
    Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], numerics_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone], iMesh, NO_RK_ITER, RunTime_EqSystem);
    SetResidual_Term(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone]);
    
    /*--- Compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$ ---*/
    
    SetRestricted_Solution(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    solver_container[iZone][iInst][iMesh+1][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh+1], solver_container[iZone][iInst][iMesh+1], config[iZone], iMesh+1, NO_RK_ITER, RunTime_EqSystem, false);
    Space_Integration(geometry[iZone][iInst][iMesh+1], solver_container[iZone][iInst][iMesh+1], numerics_container[iZone][iInst][iMesh+1][SolContainer_Position], config[iZone], iMesh+1, NO_RK_ITER, RunTime_EqSystem);
    
    /*--- Compute $P_(k+1) = I^(k+1)_k(r_k) - r_(k+1) ---*/
    
    SetForcing_Term(solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone], iMesh+1);
    
    /*--- Recursive call to MultiGrid_Cycle ---*/
    
    for (unsigned short imu = 0; imu <= RecursiveParam; imu++) {
      if (iMesh == config[iZone]->GetnMGLevels()-2) MultiGrid_Cycle(geometry, solver_container, numerics_container, config, iMesh+1, 0, RunTime_EqSystem, Iteration, iZone, iInst);
      else MultiGrid_Cycle(geometry, solver_container, numerics_container, config, iMesh+1, RecursiveParam, RunTime_EqSystem, Iteration, iZone, iInst);
    }
    
    /*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +  Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/
    
    GetProlongated_Correction(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position],
                              geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    SmoothProlongated_Correction(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], geometry[iZone][iInst][iMesh],
                                 config[iZone]->GetMG_CorrecSmooth(iMesh), 1.25, config[iZone]);
    SetProlongated_Correction(solver_container[iZone][iInst][iMesh][SolContainer_Position], geometry[iZone][iInst][iMesh], config[iZone], iMesh);
    
    /*--- Solution postsmoothing in the prolongated grid ---*/
    
    for (iPostSmooth = 0; iPostSmooth < config[iZone]->GetMG_PostSmooth(iMesh); iPostSmooth++) {
      
      switch (config[iZone]->GetKind_TimeIntScheme()) {
        case RUNGE_KUTTA_EXPLICIT: iRKLimit = config[iZone]->GetnRKStep(); break;
        case CLASSICAL_RK4_EXPLICIT: iRKLimit = 4; break;
        case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break; }

      for (iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {
        
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, iRKStep, RunTime_EqSystem, false);
        
        if (iRKStep == 0) {
          solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][iMesh]);
          if (config[iZone]->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT)
            solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_NewSolution(geometry[iZone][iInst][iMesh]);
          solver_container[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, Iteration);
        }
        
        Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], numerics_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone], iMesh, iRKStep, RunTime_EqSystem);
        Time_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iRKStep, RunTime_EqSystem, Iteration);
        
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
        
      }
    }
  }
  
}

void CMultiGridIntegration::GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                                      CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short Boundary, iMarker, iChildren, iVar;
  su2double Area_Parent, Area_Children, *Solution_Fine, *Solution_Coarse;
  
  const unsigned short nVar = sol_coarse->GetnVar();
  
  su2double *Solution = new su2double[nVar];
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Solution_Fine = sol_fine->node[Point_Fine]->GetSolution();
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] -= Solution_Fine[iVar]*Area_Children/Area_Parent;
    }
    
    Solution_Coarse = sol_coarse->node[Point_Coarse]->GetSolution();
    
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] += Solution_Coarse[iVar];
    
    for (iVar = 0; iVar < nVar; iVar++)
      sol_coarse->node[Point_Coarse]->SetSolution_Old(Solution);
    
  }
  
  /*--- Remove any contributions from no-slip walls. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    if ((Boundary == HEAT_FLUX             ) ||
        (Boundary == ISOTHERMAL            ) ||
        (Boundary == CHT_WALL_INTERFACE    )) {
      
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        
        /*--- For dirichlet boundary condtions, set the correction to zero.
         Note that Solution_Old stores the correction not the actual value ---*/
        
        sol_coarse->node[Point_Coarse]->SetVelSolutionOldZero();
        
      }
      
    }
  }
  
  /*--- MPI the set solution old ---*/
  
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION_OLD);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION_OLD);

  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      sol_fine->LinSysRes.SetBlock(Point_Fine, sol_coarse->node[Point_Coarse]->GetSolution_Old());
    }
  }
  
  delete [] Solution;
  
}

void CMultiGridIntegration::SmoothProlongated_Correction (unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                                          unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {
  su2double *Residual_Old, *Residual_Sum, *Residual, *Residual_i, *Residual_j;
  unsigned short iVar, iSmooth, iMarker, nneigh;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  
  const unsigned short nVar = solver->GetnVar();
  
  if (val_nSmooth > 0) {
    
    Residual = new su2double [nVar];
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Residual_Old = solver->LinSysRes.GetBlock(iPoint);
      solver->node[iPoint]->SetResidual_Old(Residual_Old);
    }
    
    /*--- Jacobi iterations ---*/
    
    for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
        solver->node[iPoint]->SetResidualSumZero();
      
      /*--- Loop over Interior edges ---*/
      
      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
        iPoint = geometry->edge[iEdge]->GetNode(0);
        jPoint = geometry->edge[iEdge]->GetNode(1);
        
        Residual_i = solver->LinSysRes.GetBlock(iPoint);
        Residual_j = solver->LinSysRes.GetBlock(jPoint);
        
        /*--- Accumulate nearest neighbor Residual to Res_sum for each variable ---*/
        
        solver->node[iPoint]->AddResidual_Sum(Residual_j);
        solver->node[jPoint]->AddResidual_Sum(Residual_i);
      }
      
      /*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        nneigh = geometry->node[iPoint]->GetnPoint();
        Residual_Sum = solver->node[iPoint]->GetResidual_Sum();
        Residual_Old = solver->node[iPoint]->GetResidual_Old();
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] =(Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar])
          /(1.0 + val_smooth_coeff*su2double(nneigh));
        }
        solver->LinSysRes.SetBlock(iPoint, Residual);
      }
      
      /*--- Copy boundary values ---*/
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
        if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
            (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
          for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Residual_Old = solver->node[iPoint]->GetResidual_Old();
          solver->LinSysRes.SetBlock(iPoint, Residual_Old);
        }
        }
    }
    
    delete [] Residual;
    
  }
}

void CMultiGridIntegration::Smooth_Solution(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                            unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {
  su2double *Solution_Old, *Solution_Sum, *Solution, *Solution_i, *Solution_j;
  unsigned short iVar, iSmooth, iMarker, nneigh;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  
  const unsigned short nVar = solver->GetnVar();
  
  if (val_nSmooth > 0) {
    
    Solution = new su2double [nVar];
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Solution_Old = solver->node[iPoint]->GetSolution();
      solver->node[iPoint]->SetResidual_Old(Solution_Old);
    }
    
    /*--- Jacobi iterations ---*/
    
    for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
        solver->node[iPoint]->SetResidualSumZero();
      
      /*--- Loop over Interior edges ---*/
      
      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
        iPoint = geometry->edge[iEdge]->GetNode(0);
        jPoint = geometry->edge[iEdge]->GetNode(1);
        
        Solution_i = solver->node[iPoint]->GetSolution();
        Solution_j = solver->node[jPoint]->GetSolution();
        
        /*--- Accumulate nearest neighbor Residual to Res_sum for each variable ---*/
        
        solver->node[iPoint]->AddResidual_Sum(Solution_j);
        solver->node[jPoint]->AddResidual_Sum(Solution_i);
      }
      
      /*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        nneigh = geometry->node[iPoint]->GetnPoint();
        Solution_Sum = solver->node[iPoint]->GetResidual_Sum();
        Solution_Old = solver->node[iPoint]->GetResidual_Old();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] =(Solution_Old[iVar] + val_smooth_coeff*Solution_Sum[iVar])
          /(1.0 + val_smooth_coeff*su2double(nneigh));
        }
        solver->node[iPoint]->SetSolution(Solution);
      }
      
      /*--- Copy boundary values ---*/
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
        if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
            (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
          for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution_Old = solver->node[iPoint]->GetResidual_Old();
          solver->node[iPoint]->SetSolution(Solution_Old);
        }
        }
    }
    
    delete [] Solution;
    
  }
  
}

void CMultiGridIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh) {
  unsigned long Point_Fine;
  unsigned short iVar;
  su2double *Solution_Fine, *Residual_Fine;
  
  const unsigned short nVar = sol_fine->GetnVar();
  su2double factor = config->GetDamp_Correc_Prolong(); //pow(config->GetDamp_Correc_Prolong(), iMesh+1);
  
  su2double *Solution = new su2double [nVar];
  
  for (Point_Fine = 0; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
    Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
    Solution_Fine = sol_fine->node[Point_Fine]->GetSolution();
    for (iVar = 0; iVar < nVar; iVar++) {
      /*--- Prevent a fine grid divergence due to a coarse grid divergence ---*/
      if (Residual_Fine[iVar] != Residual_Fine[iVar]) Residual_Fine[iVar] = 0.0;
      Solution[iVar] = Solution_Fine[iVar]+factor*Residual_Fine[iVar];
    }
    sol_fine->node[Point_Fine]->SetSolution(Solution);
  }
  
  /*--- MPI the new interpolated solution ---*/
  
  sol_fine->InitiateComms(geo_fine, config, SOLUTION);
  sol_fine->CompleteComms(geo_fine, config, SOLUTION);
  
  delete [] Solution;
}


void CMultiGridIntegration::SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iChildren;
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      sol_fine->node[Point_Fine]->SetSolution(sol_coarse->node[Point_Coarse]->GetSolution());
    }
  }
}

void CMultiGridIntegration::SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config, unsigned short iMesh) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iMarker, iVar, iChildren;
  su2double *Residual_Fine;
  bool dirichlet_poisson = false;

  
  const unsigned short nVar = sol_coarse->GetnVar();
  su2double factor = config->GetDamp_Res_Restric(); //pow(config->GetDamp_Res_Restric(), iMesh);
  
  su2double *Residual = new su2double[nVar];
   
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->node[Point_Coarse]->SetRes_TruncErrorZero();
    
    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] += factor*Residual_Fine[iVar];
    }
    sol_coarse->node[Point_Coarse]->AddRes_TruncError(Residual);
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
	  
	  dirichlet_poisson = false;
	  if ((config->GetKind_Incomp_System() == PRESSURE_BASED) && 
	      ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) || (config->GetMarker_All_KindBC(iMarker) ==  FAR_FIELD ))) 
	      dirichlet_poisson = true;
        
    if (((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE    )) &&
		(!(config->GetKind_Incomp_System() == PRESSURE_BASED)))  {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->node[Point_Coarse]->SetVel_ResTruncError_Zero();
      }
      cout<<"Should not be here in SetForcing routine"<<endl;
    }
    if (dirichlet_poisson) {
		cout<<"Should be here in SetForcing routine"<<endl;
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->node[Point_Coarse]->SetRes_TruncErrorZero();
      }
    }    
  }
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->node[Point_Coarse]->SubtractRes_TruncError(sol_coarse->LinSysRes.GetBlock(Point_Coarse));
  }
  
  delete [] Residual;
}

void CMultiGridIntegration::SetResidual_Term(CGeometry *geometry, CSolver *solver, CConfig *config) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) 
    solver->LinSysRes.AddBlock(iPoint, solver->node[iPoint]->GetResTruncError());  
}

void CMultiGridIntegration::SetRestricted_Residual(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iVar, iChildren;
  su2double *Residual_Fine;
  bool dirichlet_poisson = false;
  
  const unsigned short nVar = sol_coarse->GetnVar();
  
  su2double *Residual = new su2double[nVar];
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->node[Point_Coarse]->SetRes_TruncErrorZero();
    
    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] += Residual_Fine[iVar];
    }
    sol_coarse->node[Point_Coarse]->AddRes_TruncError(Residual);
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)  {
	  
	  dirichlet_poisson = false;
	  if ((config->GetKind_Incomp_System() == PRESSURE_BASED) && 
	      ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) || (config->GetMarker_All_KindBC(iMarker) ==  FAR_FIELD ))) 
	      dirichlet_poisson = true;
        
    if (((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE    )) &&
		(!(config->GetKind_Incomp_System() == PRESSURE_BASED)))  {
			
      for (iVertex = 0; iVertex<geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->node[Point_Coarse]->SetVel_ResTruncError_Zero();
      }
      cout<<"Should not be here in SetRestRes routine"<<endl;
    }
    if (dirichlet_poisson) {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->node[Point_Coarse]->SetRes_TruncErrorZero();
      }
      cout<<"Should be here in SetRestSol routine"<<endl;
    }
  }
  
  delete [] Residual;
}

void CMultiGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iVar, iChildren, iDim;
  su2double Area_Parent, Area_Children, *Solution_Fine, *Grid_Vel, Vector[3];
  
  const unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
  const unsigned short nVar = sol_coarse->GetnVar();
  const unsigned short nDim = geo_fine->GetnDim();
  const bool grid_movement  = config->GetGrid_Movement();
  
  su2double *Solution = new su2double[nVar];
  
  /*--- Compute coarse solution from fine solution ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Solution_Fine = sol_fine->node[Point_Fine]->GetSolution();
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
      }
    }
    
    sol_coarse->node[Point_Coarse]->SetSolution(Solution);
    
  }
  
  /*--- Update the solution at the no-slip walls ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE    )) &&
        (SolContainer_Position != POISSON_SOL                             )) {
      
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        
        if (SolContainer_Position == FLOW_SOL) {
          
          /*--- At moving walls, set the solution based on the new density and wall velocity ---*/
          
          if (grid_movement) {
            Grid_Vel = geo_coarse->node[Point_Coarse]->GetGridVel();
            for (iDim = 0; iDim < nDim; iDim++)
              Vector[iDim] = sol_coarse->node[Point_Coarse]->GetSolution(0)*Grid_Vel[iDim];
            sol_coarse->node[Point_Coarse]->SetVelSolutionVector(Vector);
          } else {
            
            /*--- For stationary no-slip walls, set the velocity to zero. ---*/
            
            sol_coarse->node[Point_Coarse]->SetVelSolutionZero();
          }
          
        }
        cout<<"Should not be here in SetRestSol routine"<<endl;
        if (SolContainer_Position == ADJFLOW_SOL) {
          sol_coarse->node[Point_Coarse]->SetVelSolutionDVector();
        }
        
      }
    }
    if (((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW ) ||
        (config->GetMarker_All_KindBC(iMarker) == FAR_FIELD  ))    &&
        (SolContainer_Position == POISSON_SOL                   )   ) {
			for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
				Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
				sol_coarse->node[Point_Coarse]->SetSolutionZero();
			}
		}
  }
  
  /*--- MPI the new interpolated solution ---*/
  
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION);
  
  delete [] Solution;
  
}

void CMultiGridIntegration::SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                                   CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iVar, iDim, iChildren;
  su2double Area_Parent, Area_Children, **Gradient_fine;
  
  const unsigned short nDim = geo_coarse->GetnDim();
  const unsigned short nVar = sol_coarse->GetnVar();
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double [nDim];
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPoint(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Gradient[iVar][iDim] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Gradient_fine = sol_fine->node[Point_Fine]->GetGradient();
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (iDim = 0; iDim < nDim; iDim++)
          Gradient[iVar][iDim] += Gradient_fine[iVar][iDim]*Area_Children/Area_Parent;
    }
    sol_coarse->node[Point_Coarse]->SetGradient(Gradient);
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}





void CMultiGridIntegration::NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container,
                                                      CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem, unsigned long Iteration,
                                                      su2double *monitor) {
  
  const unsigned short nDim = geometry[FinestMesh]->GetnDim();
  
  switch (RunTime_EqSystem) {
      
    case RUNTIME_FLOW_SYS:
      
      /*--- Calculate the inviscid and viscous forces ---*/
      
      solver_container[FinestMesh][FLOW_SOL]->Pressure_Forces(geometry[FinestMesh], config);
      solver_container[FinestMesh][FLOW_SOL]->Momentum_Forces(geometry[FinestMesh], config);
      solver_container[FinestMesh][FLOW_SOL]->Friction_Forces(geometry[FinestMesh], config);
          
      /*--- Evaluate the buffet metric if requested ---*/
      
      if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){
          solver_container[FinestMesh][FLOW_SOL]->Buffet_Monitoring(geometry[FinestMesh], config);
      }
      
      /*--- Evaluate convergence monitor ---*/
      
      if (config->GetConvCriteria() == CAUCHY) {
        if (config->GetCauchy_Func_Flow() == DRAG_COEFFICIENT) (*monitor) = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CD();
        if (config->GetCauchy_Func_Flow() == LIFT_COEFFICIENT) (*monitor) = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CL();
        if (config->GetCauchy_Func_Flow() == NEARFIELD_PRESSURE) (*monitor) = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
      }
      
      if (config->GetConvCriteria() == RESIDUAL) {
        if (config->GetResidual_Func_Flow() == RHO_RESIDUAL) (*monitor) = log10(solver_container[FinestMesh][FLOW_SOL]->GetRes_RMS(0));
        else if (config->GetResidual_Func_Flow() == RHO_ENERGY_RESIDUAL) {
          if (nDim == 2) (*monitor) = log10(solver_container[FinestMesh][FLOW_SOL]->GetRes_RMS(3));
          else (*monitor) = log10(solver_container[FinestMesh][FLOW_SOL]->GetRes_RMS(4));
        }
      }
      
      break;
      
    case RUNTIME_ADJFLOW_SYS:
      
      /*--- Calculate the inviscid and viscous sensitivities ---*/
      
      solver_container[FinestMesh][ADJFLOW_SOL]->Inviscid_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh], numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
      solver_container[FinestMesh][ADJFLOW_SOL]->Viscous_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh], numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
      
      /*--- Smooth the inviscid and viscous sensitivities ---*/
      
      if (config->GetKind_SensSmooth() != NONE) solver_container[FinestMesh][ADJFLOW_SOL]->Smooth_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh], numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
      
      /*--- Evaluate convergence monitor ---*/
      
      if (config->GetConvCriteria() == CAUCHY) {
        if (config->GetCauchy_Func_AdjFlow() == SENS_GEOMETRY) (*monitor) = solver_container[FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
        if (config->GetCauchy_Func_AdjFlow() == SENS_MACH) (*monitor) = solver_container[FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
      }
      
      if (config->GetConvCriteria() == RESIDUAL) {
        if (config->GetResidual_Func_Flow() == RHO_RESIDUAL) (*monitor) = log10(solver_container[FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0));
        else if (config->GetResidual_Func_Flow() == RHO_ENERGY_RESIDUAL) {
          if (nDim == 2) (*monitor) = log10(solver_container[FinestMesh][ADJFLOW_SOL]->GetRes_RMS(3));
          else (*monitor) = log10(solver_container[FinestMesh][ADJFLOW_SOL]->GetRes_RMS(4));
        }
      }
      
      break;
          
  }
  
}

CSingleGridIntegration::CSingleGridIntegration(CConfig *config) : CIntegration(config) { }

CSingleGridIntegration::~CSingleGridIntegration(void) { }

void CSingleGridIntegration::SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone, unsigned short iInst) {
  unsigned short iMesh;
  su2double monitor = 0.0;
  
  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  unsigned short FinestMesh = config[iZone]->GetFinestMesh();
  
  /*--- Preprocessing ---*/
    
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh, 0, RunTime_EqSystem, false);
  
  /*--- Set the old solution ---*/
  
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][FinestMesh]);
  
  /*--- Time step evaluation ---*/
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh, 0);
  
  /*--- Space integration ---*/
  Space_Integration(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], numerics_container[iZone][iInst][FinestMesh][SolContainer_Position],
                    config[iZone], FinestMesh, NO_RK_ITER, RunTime_EqSystem);
  
  /*--- Time integration ---*/
  Time_Integration(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], NO_RK_ITER,
                   RunTime_EqSystem, Iteration);
  
  /*--- Postprocessing ---*/
  
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh);
  
  /*--- Compute adimensional parameters and the convergence monitor ---*/
  
  switch (RunTime_EqSystem) {
    case RUNTIME_FEA_SYS:     monitor = log10(solver_container[iZone][iInst][FinestMesh][FEA_SOL]->GetRes_RMS(0));     break;
    case RUNTIME_HEAT_SYS:    monitor = log10(solver_container[iZone][iInst][FinestMesh][HEAT_SOL]->GetRes_RMS(0));    break;
    case RUNTIME_POISSON_SYS: monitor = log10(solver_container[iZone][iInst][FinestMesh][POISSON_SOL]->GetRes_RMS(0)); break;
  }
  
  if (RunTime_EqSystem == RUNTIME_HEAT_SYS) {
    solver_container[iZone][iInst][FinestMesh][HEAT_SOL]->Heat_Fluxes(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone]);
  }
  
  /*--- Convergence strategy ---*/
  
  Convergence_Monitoring(geometry[iZone][iInst][FinestMesh], config[iZone], Iteration, monitor, FinestMesh);
  
  /*--- If turbulence model, copy the turbulence variables to the coarse levels ---*/
  
  if (RunTime_EqSystem == RUNTIME_TURB_SYS) {
    for (iMesh = FinestMesh; iMesh < config[iZone]->GetnMGLevels(); iMesh++) {
      SetRestricted_Solution(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
      SetRestricted_EddyVisc(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    }
  }
}

void CSingleGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iVar, iChildren;
  su2double Area_Parent, Area_Children, *Solution_Fine, *Solution;
  
  unsigned short nVar = sol_coarse->GetnVar();
  
  Solution = new su2double[nVar];
  
  /*--- Compute coarse solution from fine solution ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Solution_Fine = sol_fine->node[Point_Fine]->GetSolution();
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
    }
    
    sol_coarse->node[Point_Coarse]->SetSolution(Solution);
    
  }
  
  /*--- MPI the new interpolated solution ---*/
  
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION);
  
  delete [] Solution;
  
}

void CSingleGridIntegration::SetRestricted_EddyVisc(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iChildren;
  su2double Area_Parent, Area_Children, EddyVisc_Fine, EddyVisc;
  
  /*--- Compute coarse Eddy Viscosity from fine solution ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    EddyVisc = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      EddyVisc_Fine = sol_fine->node[Point_Fine]->GetmuT();
      EddyVisc += EddyVisc_Fine*Area_Children/Area_Parent;
    }
    
    sol_coarse->node[Point_Coarse]->SetmuT(EddyVisc);
    
  }
  
  /*--- Update solution at the no slip wall boundary, only the first
   variable (nu_tilde -in SA and SA_NEG- and k -in SST-), to guarantee that the eddy viscoisty
   is zero on the surface ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE     )) {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->node[Point_Coarse]->SetmuT(0);
      }
    }
  }

  /*--- MPI the new interpolated solution (this also includes the eddy viscosity) ---*/
    
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION_EDDY);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION_EDDY);
  
}



void CSingleGridIntegration::SetRestricted_PoissonSource(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iChildren;
  su2double Area_Parent, Area_Children, PoissonSource_Fine, PoissonSource;
  
  /*--- Compute coarse Source term from fine solution ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    PoissonSource = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      PoissonSource_Fine = sol_fine->node[Point_Fine]->GetSourceTerm();
      PoissonSource += PoissonSource_Fine*Area_Children/Area_Parent;
    }
    
    sol_coarse->node[Point_Coarse]->SetSourceTerm(PoissonSource);
    
  }
  
  

  /*--- MPI the new interpolated solution (this also includes the eddy viscosity) ---*/
  
  //sol_coarse->Set_MPI_Solution(geo_coarse, config);
  
}


CStructuralIntegration::CStructuralIntegration(CConfig *config) : CIntegration(config) { }

CStructuralIntegration::~CStructuralIntegration(void) { }

void CStructuralIntegration::Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config, unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone, unsigned short iInst) {

  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  /*--- Preprocessing ---*/

  solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0],
      config[iZone], numerics_container[iZone][iInst][MESH_0][SolContainer_Position], MESH_0, Iteration, RunTime_EqSystem, false);


  /*--- Space integration ---*/

  Space_Integration_FEM(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                    config[iZone], RunTime_EqSystem, Iteration);

  /*--- Time integration ---*/

  Time_Integration_FEM(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                config[iZone], RunTime_EqSystem, Iteration);

  /*--- Postprocessing ---*/

  solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0],
      config[iZone], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],  MESH_0);

  /*--- Convergence strategy ---*/

  switch (RunTime_EqSystem) {
    case RUNTIME_FEA_SYS: Convergence_Monitoring_FEM(geometry[iZone][iInst][MESH_0], config[iZone], solver_container[iZone][iInst][MESH_0][SolContainer_Position], Iteration); break;
    case RUNTIME_ADJFEA_SYS: Convergence_Monitoring_FEM_Adj(geometry[iZone][iInst][MESH_0], config[iZone], solver_container[iZone][iInst][MESH_0][SolContainer_Position], Iteration); break;
  }

}

CFEM_DG_Integration::CFEM_DG_Integration(CConfig *config) : CIntegration(config) { }

CFEM_DG_Integration::~CFEM_DG_Integration(void) { }

void CFEM_DG_Integration::SingleGrid_Iteration(CGeometry ****geometry,
                                               CSolver *****solver_container,
                                               CNumerics ******numerics_container,
                                               CConfig **config,
                                               unsigned short RunTime_EqSystem,
                                               unsigned long Iteration,
                                               unsigned short iZone,
                                               unsigned short iInst) {

  unsigned short iMesh, iStep, iLimit = 1;
  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /*--- For now, we assume no geometric multigrid. ---*/
  iMesh = FinestMesh;

  /*--- Check if only the Jacobian of the spatial discretization must
        be computed. If so, call the appropriate function and return. ---*/
  if (config[iZone]->GetJacobian_Spatial_Discretization_Only()) {
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->ComputeSpatialJacobian(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                         numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                                                                                         config[iZone], iMesh, RunTime_EqSystem);
    return;
  }

  /*--- Determine the number of stages in the time stepping algorithm.
        For the Runge-Kutta schemes this is the number of RK stages,
        while for ADER-DG this information is not used, because a more
        complicated algorithm must be used to facilitate time accurate
        local time stepping.  Note that we are currently hard-coding
        the classical RK4 scheme. ---*/
  bool useADER = false;
  switch (config[iZone]->GetKind_TimeIntScheme()) {
    case RUNGE_KUTTA_EXPLICIT: iLimit = config[iZone]->GetnRKStep(); break;
    case CLASSICAL_RK4_EXPLICIT: iLimit = 4; break;
    case ADER_DG: iLimit = 1; useADER = true; break;
    case EULER_EXPLICIT: case EULER_IMPLICIT: iLimit = 1; break; }

  /*--- In case an unsteady simulation is carried out, it is possible that a
        synchronization time step is specified. If so, set the boolean
        TimeSynSpecified to true, which leads to an outer loop in the
        algorithm below. ---*/
  bool TimeSyncSpecified   = false;
  const su2double TimeSync = config[iZone]->GetDelta_UnstTimeND();
  if(config[iZone]->GetUnsteady_Simulation() == TIME_STEPPING &&
     config[iZone]->GetUnst_CFL()            != 0.0           &&
     TimeSync                                != 0.0) TimeSyncSpecified = true;

  /*--- Outer loop, which is only active when a synchronization time has been
        specified for an unsteady simulation. ---*/
  bool syncTimeReached = false;
  su2double timeEvolved = 0.0;
  while( !syncTimeReached ) {

    /* Compute the time step for stability. */
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh],
                                                                               solver_container[iZone][iInst][iMesh],
                                                                               config[iZone], iMesh, Iteration);
    /* Possibly overrule the specified time step when a synchronization time was
       specified and determine whether or not the time loop must be continued.
       When TimeSyncSpecified is false, the loop is always terminated. */
    if( TimeSyncSpecified )
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->CheckTimeSynchronization(config[iZone],
                                                                                             TimeSync, timeEvolved,
                                                                                             syncTimeReached);
    else
      syncTimeReached = true;

    /*--- For ADER in combination with time accurate local time stepping, the
          space and time integration are tightly coupled and cannot be treated
          segregatedly. Therefore a different function is called for ADER to
          carry out the space and time integration. ---*/
    if( useADER ) {
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->ADER_SpaceTimeIntegration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                              numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                                                                                              config[iZone], iMesh, RunTime_EqSystem);
    }
    else {

      /*--- Time and space integration can be decoupled. ---*/
      for (iStep = 0; iStep < iLimit; iStep++) {

        /*--- Preprocessing ---*/
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                    config[iZone], iMesh, iStep, RunTime_EqSystem, false);

        /*--- Space integration ---*/
        Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                          numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                          config[iZone], iMesh, iStep, RunTime_EqSystem);

        /*--- Time integration, update solution using the old solution plus the solution increment ---*/
        Time_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                         config[iZone], iStep, RunTime_EqSystem, Iteration);

        /*--- Postprocessing ---*/
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                     config[iZone], iMesh);
      }
    }
  }

  /*--- Calculate the inviscid and viscous forces ---*/
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Pressure_Forces(geometry[iZone][iInst][iMesh], config[iZone]);

  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Friction_Forces(geometry[iZone][iInst][iMesh], config[iZone]);

  /*--- Convergence strategy ---*/

  //Convergence_Monitoring(geometry[iZone][iInst][FinestMesh], config[iZone], Iteration, monitor, FinestMesh);
}

void CFEM_DG_Integration::Space_Integration(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics **numerics,
                                            CConfig *config, unsigned short iMesh,
                                            unsigned short iStep,
                                            unsigned short RunTime_EqSystem) {

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Runge-Kutta type of time integration schemes. In the first step, i.e.
        if iStep == 0, set the old solution (working solution for the DG part),
        and if needed, the new solution. ---*/
  if (iStep == 0) {
    solver_container[MainSolver]->Set_OldSolution(geometry);

    if (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT) {
      solver_container[MainSolver]->Set_NewSolution(geometry);
    }
  }

  /*--- Compute the spatial residual by processing the task list. ---*/
  solver_container[MainSolver]->ProcessTaskList_DG(geometry, solver_container, numerics, config, iMesh);
}

void CFEM_DG_Integration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iStep,
                                    unsigned short RunTime_EqSystem, unsigned long Iteration) {

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Perform the time integration ---*/
  switch (config->GetKind_TimeIntScheme()) {
    case (RUNGE_KUTTA_EXPLICIT):
      solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iStep);
      break;
    case (CLASSICAL_RK4_EXPLICIT):
      solver_container[MainSolver]->ClassicalRK4_Iteration(geometry, solver_container, config, iStep);
      break;
    default:
      SU2_MPI::Error("Time integration scheme not implemented.", CURRENT_FUNCTION);
  }
}


void CMultiGridIntegration::MultiGrid_CyclePB(CGeometry ****geometry,
                                            CSolver *****solver_container,
                                            CNumerics ******numerics_container,
                                            CConfig **config,
                                            unsigned short iMesh,
                                            unsigned short RecursiveParam,
                                            unsigned short RunTime_EqSystem,
                                            unsigned long Iteration,
                                            unsigned short iZone,
                                            unsigned short iInst) {
  unsigned short iPreSmooth, iPostSmooth, iRKStep, iRKLimit = 1;
    
  bool startup_multigrid = (config[iZone]->GetRestart_Flow() && (RunTime_EqSystem == RUNTIME_FLOW_SYS) && (Iteration == 0));
  //bool piso = true;
  bool piso = (config[iZone]->GetKind_PBIter() == PISO);
  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  
  /*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/
  
  
  for (iPreSmooth = 0; iPreSmooth < config[iZone]->GetMG_PreSmooth(iMesh); iPreSmooth++) {
	  
	 /*--- Momentum solve (Predict velocities using existing pressure field)---*/
	 switch( config[iZone]->GetKind_Solver() ) {
		  case EULER: 
		  config[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, Iteration); break;
		  
		  case NAVIER_STOKES:
		  config[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, Iteration); break;
		  
		  case RANS:
		  config[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, Iteration); break;
     }
     
     config[iZone]->SetKind_TimeIntScheme(EULER_IMPLICIT);
     
     CurrentGridIteration(geometry, solver_container, numerics_container, config,
                  iMesh, RUNTIME_FLOW_SYS,
                  Iteration, iZone, iInst);
      
     /*--- Set source term for pressure correction equation based on current flow solution ---*/
     solver_container[iZone][iInst][iMesh][FLOW_SOL]->SetMomCoeff(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], false, iMesh);
     
     solver_container[iZone][iInst][iMesh][FLOW_SOL]->SetPoissonSourceTerm(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
	

     /*--- Solve the poisson equation (pressure correction) ---*/
     config[iZone]->SetGlobalParam(POISSON_EQUATION, RUNTIME_POISSON_SYS, Iteration);
     CurrentGridIteration(geometry, solver_container, numerics_container, config,
                  iMesh, RUNTIME_POISSON_SYS,
                  Iteration, iZone, iInst);     
     
     /*--- Correct pressure and velocities ---*/
      solver_container[iZone][iInst][iMesh][FLOW_SOL]->Flow_Correction(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone]);
      
      solver_container[iZone][iInst][iMesh][FLOW_SOL]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
      
     if (piso) {
		 config[iZone]->SetKind_TimeIntScheme(EULER_EXPLICIT);
		 
		 /*--- Momentum solve (Predict velocities using existing pressure field)---*/
		 switch( config[iZone]->GetKind_Solver() ) {
			 case EULER: 
			 config[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, Iteration); break;
			 
			 case NAVIER_STOKES:
			 config[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, Iteration); break;
			 
			 case RANS:
			 config[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, Iteration); break;
		}
		
		CurrentGridIteration(geometry, solver_container, numerics_container, config,
                  iMesh, RUNTIME_FLOW_SYS,
                  Iteration, iZone, iInst);
		
		config[iZone]->SetKind_TimeIntScheme(EULER_IMPLICIT);

		/*--- Set source term for pressure correction equation based on current flow solution ---*/
		solver_container[iZone][iInst][iMesh][FLOW_SOL]->SetMomCoeff(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], false, iMesh);
     
		solver_container[iZone][iInst][iMesh][FLOW_SOL]->SetPoissonSourceTerm(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
		 
		/*--- Solve the poisson equation (pressure correction) ---*/
		config[iZone]->SetGlobalParam(POISSON_EQUATION, RUNTIME_POISSON_SYS, Iteration);
		 
		CurrentGridIteration(geometry, solver_container, numerics_container, config,
                  iMesh, RUNTIME_POISSON_SYS,
                  Iteration, iZone, iInst);
         
        /*--- Correct pressure and velocities ---*/
        solver_container[iZone][iInst][iMesh][FLOW_SOL]->Flow_Correction(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone]); 
	}
    
  }
  /*--- Compute Forcing Term $P_(k+1) = I^(k+1)_k(P_k+F_k(u_k))-F_(k+1)(I^(k+1)_k u_k)$ and update solution for multigrid ---*/
  
  if ( (iMesh < config[iZone]->GetnMGLevels() && ((Iteration >= config[iZone]->GetnStartUpIter()) || startup_multigrid)) ) {
	  
	   //cout<<"Find flow res on "<<iMesh<<endl;
	  /* Find Residual for momentum equation - SpaceIntegration(imesh) with momentum eq, R_k*/
	    Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], numerics_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone], 
	                      iMesh, NO_RK_ITER, RunTime_EqSystem);
	    
	   //cout<<"Find massflux on "<<iMesh<<endl;
	   /* Find Mass imbalance - SetPoissonSourceterm(imesh) - Sets ResMassFlux m_k*/
	    solver_container[iZone][iInst][iMesh][FLOW_SOL]->SetPoissonSourceTerm(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
	    
	   //cout<<"Update R and M with trun errors on "<<iMesh<<endl;
	   /* Set Residuals - R_k = R_k + P_k^m and m_k = m_k + P_k^c*/
	    SetResidual_Term(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone]);
	     
	   //cout<<"Restrict Vars from "<<iMesh<<" to "<<iMesh+1<<endl;
	   /* RestrictVars - Restrict the primitive variables - Find U_k+1 and p_k+1 from U_k and p_k*/
	    SetRestricted_Variables(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][FLOW_SOL], solver_container[iZone][iInst][iMesh+1][FLOW_SOL]
	                           ,geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
	   
	   //cout<<"Find gradients on "<<iMesh+1<<endl;
	   /* Call Postprocessing routine to compute gradients of the restricted primitive varibales. Need pressure gradients during mass flux computation. */ 
	    solver_container[iZone][iInst][iMesh+1][FLOW_SOL]->Postprocessing(geometry[iZone][iInst][iMesh+1], solver_container[iZone][iInst][iMesh+1], config[iZone], iMesh+1);
	    
	   //cout<<"Find flow residuals on "<<iMesh+1<<endl;
	   /* Find Residual for momentum equation - SpaceIntegration(imesh+1) with momentum eq, R_k+1*/
	    Space_Integration(geometry[iZone][iInst][iMesh+1], solver_container[iZone][iInst][iMesh+1], numerics_container[iZone][iInst][iMesh+1][SolContainer_Position], config[iZone],
	                     iMesh+1, NO_RK_ITER, RunTime_EqSystem);
	    
	    //cout<<"Find massflux on "<<iMesh+1<<endl;
	   /* Find Mass imbalance - SetPoissonSourceterm(imesh+1) - Sets ResMassFlux m_k+1*/
	    solver_container[iZone][iInst][iMesh+1][FLOW_SOL]->SetPoissonSourceTerm(geometry[iZone][iInst][iMesh+1], solver_container[iZone][iInst][iMesh+1], config[iZone], iMesh+1);
	    
	   //cout<<"Find trunc errors based on restricted vars and restricted residuals on "<<iMesh+1<<endl;
	   /* SetForcingTerm(imesh+1,imesh) P_(k+1)^m = I^(k+1)_k R_k - R_(k+1), P_(k+1)^c = I^(k+1)_k m_k - m_k+1*/
	    SetPBForcing_Term(solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], 
	                      geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone], iMesh+1);
	    
	   //cout<<"Call MGCycle with "<<iMesh+1<<endl;
	   /* Multigrid_CyclePB(imesh+1) Find the solution at new iteration level (Only V cycle for now)*/
	    MultiGrid_CyclePB(geometry, solver_container, numerics_container, config, iMesh+1, 0, RunTime_EqSystem, Iteration, iZone, iInst);
        
        //cout<<"Get corrections from "<<iMesh+1<<" and interpolate back to "<<iMesh<<endl;
	   /* Get the prolongated correction to the variables - Find the update on the coarse grid level and prolongate to current grid (Have to prolongate primitive variables and not solutions)*/
	    GetProlongated_VarCorrection(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][FLOW_SOL], solver_container[iZone][iInst][iMesh+1][SolContainer_Position],
                              geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
	    
	    
	   /* Smooth the prolongated correction to the primitive variables (error) using weighted averages of the neighbors
	    SmoothProlongated_VarCorrection(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][FLOW_SOL], geometry[iZone][iInst][iMesh],
                                 config[iZone]->GetMG_CorrecSmooth(iMesh), 1.25, config[iZone]);*/
	    
	    //cout<<"Set corrections on "<<iMesh<<endl;
	   /* SetProlongated_VarCorrection - Apply the smoothed correction to the variables*/
	    SetProlongated_VarCorrection(solver_container[iZone][iInst][iMesh][FLOW_SOL], geometry[iZone][iInst][iMesh], config[iZone], iMesh);
	   
	   
	   /*--- Solution postsmoothing in the prolongated grid ---*/
    
       for (iPostSmooth = 0; iPostSmooth < config[iZone]->GetMG_PostSmooth(iMesh); iPostSmooth++) {
		   
		   /*--- Momentum solve (Predict velocities using existing pressure field)---*/
		   switch( config[iZone]->GetKind_Solver() ) {
			   case EULER: 
			   config[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, Iteration); break;
			   
			   case NAVIER_STOKES:
			   config[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, Iteration); break;
			   
			   case RANS:
			   config[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, Iteration); break;
			   
		    }
		    
		    CurrentGridIteration(geometry, solver_container, numerics_container, config,
                  iMesh, RUNTIME_FLOW_SYS,
                  Iteration, iZone, iInst);
      
            /*--- Set source term for pressure correction equation based on current flow solution ---*/
            solver_container[iZone][iInst][iMesh][FLOW_SOL]->SetMomCoeff(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], false, iMesh);
     
            solver_container[iZone][iInst][iMesh][FLOW_SOL]->SetPoissonSourceTerm(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone],  iMesh);
	

           /*--- Solve the poisson equation (pressure correction) ---*/
           config[iZone]->SetGlobalParam(POISSON_EQUATION, RUNTIME_POISSON_SYS, Iteration);
     
           CurrentGridIteration(geometry, solver_container, numerics_container, config,
                  iMesh, RUNTIME_POISSON_SYS,
                  Iteration, iZone, iInst);
     
     
           /*--- Correct pressure and velocities ---*/
           solver_container[iZone][iInst][iMesh][FLOW_SOL]->Flow_Correction(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone]);
       }

  }
  												


}

void CMultiGridIntegration::CurrentGridIteration(CGeometry ****geometry,
                                            CSolver *****solver,
                                            CNumerics ******numerics,
                                            CConfig **config,
                                            unsigned short iMesh,
                                            unsigned short RunTime_EqSystem,
                                            unsigned long Iteration,
                                            unsigned short iZone,
                                            unsigned short iInst) {

	unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
	su2double monitor = 0.0;
	unsigned short FinestMesh = config[iZone]->GetFinestMesh();
    
    
    /*--- Send-Receive boundary conditions, and preprocessing ---*/
      
    solver[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver[iZone][iInst][iMesh], config[iZone], iMesh, 0, RunTime_EqSystem, false);
             
    /*--- Set the old solution ---*/
        
    solver[iZone][iInst][iMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][iMesh]);

    /*--- Compute time step, max eigenvalue, and integration scheme (steady and unsteady problems) ---*/
    solver[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh], solver[iZone][iInst][iMesh], config[iZone], iMesh, Iteration);
      
    /*--- Space integration ---*/
    Space_Integration(geometry[iZone][iInst][iMesh], solver[iZone][iInst][iMesh], numerics[iZone][iInst][iMesh][SolContainer_Position], config[iZone], iMesh, NO_RK_ITER, RunTime_EqSystem);
      
    /*--- Time integration, update solution using the old solution plus the solution increment ---*/
    Time_Integration(geometry[iZone][iInst][iMesh], solver[iZone][iInst][iMesh], config[iZone], NO_RK_ITER, RunTime_EqSystem, Iteration);
      
    /*--- Send-Receive boundary conditions, and postprocessing ---*/
      
    solver[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver[iZone][iInst][iMesh], config[iZone], iMesh);

   /*--- Compute adimensional parameters and the convergence monitor ---*/
  
   switch (RunTime_EqSystem) {
    case RUNTIME_POISSON_SYS: monitor = log10(solver[iZone][iInst][iMesh][POISSON_SOL]->GetRes_RMS(0)); break;
    case RUNTIME_FLOW_SYS: monitor = log10(solver[iZone][iInst][iMesh][FLOW_SOL]->GetRes_RMS(0)); break;
  }
  
  /*--- Convergence strategy ---*/
  
  Convergence_Monitoring(geometry[iZone][iInst][FinestMesh], config[iZone], Iteration, monitor, FinestMesh);
  
  /*--- If turbulence model, copy the turbulence variables to the coarse levels ---*/
  
  if (RunTime_EqSystem == RUNTIME_TURB_SYS) {
    for (iMesh = FinestMesh; iMesh < config[iZone]->GetnMGLevels(); iMesh++) {
      SetRestricted_Solution(RunTime_EqSystem, solver[iZone][iInst][iMesh][SolContainer_Position], solver[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
      SetRestricted_EddyVisc(RunTime_EqSystem, solver[iZone][iInst][iMesh][SolContainer_Position], solver[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    }
  }
    
}


void CMultiGridIntegration::SetRestricted_Variables(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iVar, iChildren, iMarker, iDim;
  su2double Area_Parent, Area_Children, P_Outlet, *Var_Fine, *Var_Coarse, *Grid_Vel, Vector[3];
  su2double Flow_Dir[3], Flow_Dir_Mag, UnitFlowDir[3], Vel_Mag;
  
  const unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
  unsigned short nPrimVar = sol_coarse->GetnPrimVar();
  const unsigned short nDim = geo_fine->GetnDim();
  string Marker_Tag;
  const bool grid_movement  = config->GetGrid_Movement();
  ofstream Coarse_GridFile, Fine_GridFile;
  
  Coarse_GridFile.open("coarse_grid.txt", ios::out);
  Fine_GridFile.open("fine_grid.txt", ios::out);
  
  Var_Coarse = new su2double[nPrimVar];
  
  /*--- Compute coarse solution from fine solution ---*/
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nPrimVar; iVar++) Var_Coarse[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Var_Fine = sol_fine->node[Point_Fine]->GetPrimitive();
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        Var_Coarse[iVar] += Var_Fine[iVar]*Area_Children/Area_Parent;
	  }
	  
    }
    sol_coarse->node[Point_Coarse]->SetPrimitive(Var_Coarse);
  }
  
  
  
   /*--- Update the solution at the no-slip walls ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
	
	Marker_Tag  = config->GetMarker_All_TagBound(iMarker);
    
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             )) {
      
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        
        if (SolContainer_Position == FLOW_SOL) {
          
          /*--- At moving walls, set the solution based on the new density and wall velocity ---*/
          
          if (grid_movement) {
            Grid_Vel = geo_coarse->node[Point_Coarse]->GetGridVel();
            for (iDim = 0; iDim < nDim; iDim++) {
              Vector[iDim] = sol_coarse->node[Point_Coarse]->GetSolution(0)*Grid_Vel[iDim];
              sol_coarse->node[Point_Coarse]->SetPrimitive(iDim+1,Vector[iDim]);
		    }
          } else {
            
            /*--- For stationary no-slip walls, set the velocity to zero. ---*/
            for (iDim = 0; iDim < nDim; iDim++) 
            sol_coarse->node[Point_Coarse]->SetPrimitive(iDim+1,0.0);
          }
          
        }
      }
    }
     if ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) &&
         (config->GetKind_Inc_Outlet(Marker_Tag) == PRESSURE_OUTLET)) {
			 
		/*--- Retrieve the specified back pressure for this outlet. ---*/
        P_Outlet = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();
		 
		 for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
			 Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
			 sol_coarse->node[Point_Coarse]->SetPrimitive(0,P_Outlet);
		 }
	 }  
	 if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW)) {
			 
		 for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
			Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
			
			Vel_Mag  = sol_coarse->GetInlet_Ptotal(iMarker,iVertex)/config->GetVelocity_Ref();
			
			Flow_Dir_Mag = 0.0;
			for (iDim = 0; iDim < nDim; iDim++) {
				Flow_Dir[iDim] = sol_coarse->GetInlet_FlowDir(iMarker,iVertex,iDim);
				Flow_Dir_Mag += Flow_Dir[iDim]*Flow_Dir[iDim];
				
			}
			Flow_Dir_Mag = sqrt(Flow_Dir_Mag);
			
      
			/*--- Store the unit flow direction vector. ---*/
      
			for (iDim = 0; iDim < nDim; iDim++)
				UnitFlowDir[iDim] = Flow_Dir[iDim]/Flow_Dir_Mag;
			
			
			 
            for (iDim = 0; iDim < nDim; iDim++) {
              Vector[iDim] = Vel_Mag*UnitFlowDir[iDim];
              sol_coarse->node[Point_Coarse]->SetPrimitive(iDim+1,Vector[iDim]);
		    }
		 }
	 }  
  }
  
   /*--- Compute coarse solution from fine solution ---*/
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {     
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Var_Fine = sol_fine->node[Point_Fine]->GetPrimitive();

      Fine_GridFile<<geo_fine->node[Point_Fine]->GetCoord(0)<<"\t"<<geo_fine->node[Point_Fine]->GetCoord(1)<<"\t";
      Fine_GridFile<<Var_Fine[0]<<"\t"<<Var_Fine[1]<<"\t"<<Var_Fine[2]<<endl;
    }

    for (iVar = 0; iVar < nPrimVar; iVar++) Var_Coarse[iVar] = sol_coarse->node[Point_Coarse]->GetPrimitive(iVar);
    Coarse_GridFile<<geo_coarse->node[Point_Coarse]->GetCoord(0)<<"\t"<<geo_coarse->node[Point_Coarse]->GetCoord(1)<<"\t";
    Coarse_GridFile<<Var_Coarse[0]<<"\t"<<Var_Coarse[1]<<"\t"<<Var_Coarse[2]<<endl;
  }
  
  
  Coarse_GridFile.close();
  Fine_GridFile.close();
  delete [] Var_Coarse;
  
}

void CMultiGridIntegration::SetPBForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config, unsigned short iMesh) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iMarker, iVar, iChildren, iDim;
  su2double *Residual_Fine, Area_Children, Area_Parent;
  
  const unsigned short nVar = sol_coarse->GetnVar();
  const unsigned short nDim = geo_coarse->GetnDim();
  su2double factor = config->GetDamp_Res_Restric(); //pow(config->GetDamp_Res_Restric(), iMesh);
  su2double factor2 = 0.0;
  
  su2double *Residual = new su2double[nVar];
  su2double MassFlux, MassFlux_Fine;
  
  ofstream Coarse_GridFile, Fine_GridFile;
  
  Coarse_GridFile.open("coarse_mass.txt", ios::out);
  Fine_GridFile.open("fine_mass.txt", ios::out);

  /*--- Momentum equation forcing term. ---*/  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->node[Point_Coarse]->SetRes_TruncErrorZero();
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    /*--- Interpolation of fine grid residual (I^(k+1)_k R_k). ---*/
    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      
      /*for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] += factor*(Area_Children/Area_Parent)*Residual_Fine[iVar];*/
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] += factor*Residual_Fine[iVar];
    }
    sol_coarse->node[Point_Coarse]->AddRes_TruncError(Residual);
  }
  
  /*--- Find truncation error P_(k+1) = I^(k+1)_k R_k - R_(k+1) ---*/
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->node[Point_Coarse]->SubtractRes_TruncError(sol_coarse->LinSysRes.GetBlock(Point_Coarse));
  }
  
  /*--- Reset dirichlet conditions. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW             )) {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++) sol_coarse->node[Point_Coarse]->SetVal_ResTruncError_Zero(iDim);
      }
    }
  }

  /*--- Continuity equation forcing term. ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->node[Point_Coarse]->SetMass_TruncErrorZero();
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    MassFlux = 0.0;
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      MassFlux_Fine = sol_fine->node[Point_Fine]->GetMassFlux();
      //MassFlux += factor*(Area_Children/Area_Parent)*MassFlux_Fine;
      MassFlux += factor*MassFlux_Fine;
      Fine_GridFile<<Point_Fine<<"\t"<<geo_fine->node[Point_Fine]->GetCoord(0)<<"\t"<<geo_fine->node[Point_Fine]->GetCoord(1)<<"\t"<<MassFlux_Fine<<endl;
    }
    sol_coarse->node[Point_Coarse]->AddMass_TruncError(MassFlux);
    Coarse_GridFile<<Point_Coarse<<"\t"<<geo_coarse->node[Point_Coarse]->GetCoord(0)<<"\t"<<geo_coarse->node[Point_Coarse]->GetCoord(1)<<"\t";
    Coarse_GridFile<<MassFlux<<"\t"<<sol_coarse->node[Point_Coarse]->GetMassFlux()<<endl;
  }
  
  /*for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->node[Point_Coarse]->SubtractMass_TruncError(sol_coarse->node[Point_Coarse]->GetMassFlux());
    Coarse_GridFile<<Point_Coarse<<"\t"<<geo_coarse->node[Point_Coarse]->GetCoord(0)<<"\t"<<geo_coarse->node[Point_Coarse]->GetCoord(1)<<"\t";
    Coarse_GridFile<<sol_coarse->node[Point_Coarse]->GetMassTruncError()<<"\t"<<sol_coarse->node[Point_Coarse]->GetMassFlux()<<endl;
  }*/
  
  /*--- Reset dirichlet conditions. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->node[Point_Coarse]->SetMass_TruncErrorZero();
      }
    }
  }


  
  delete [] Residual;
  Coarse_GridFile.close();
  Fine_GridFile.close();
}



void CMultiGridIntegration::GetProlongated_VarCorrection(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                                      CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short Boundary, iMarker, iChildren, iVar, iDim;
  su2double Area_Parent, Area_Children, *PrimVar_Fine, *PrimVar_Coarse;
  
  const unsigned short nVar = sol_coarse->GetnPrimVar();
  const unsigned short nDim = geo_coarse->GetnDim();
  
  su2double *PrimVar_Temp = new su2double[nVar];
  
  ofstream Coarse_GridFile, Fine_GridFile;
  
  Coarse_GridFile.open("coarse_corr.txt", ios::out);
  Fine_GridFile.open("fine_corr.txt", ios::out);
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++) PrimVar_Temp[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      PrimVar_Fine = sol_fine->node[Point_Fine]->GetPrimitive();
      for (iVar = 0; iVar < nVar; iVar++)
        PrimVar_Temp[iVar] -= PrimVar_Fine[iVar]*Area_Children/Area_Parent;
    }
    
    PrimVar_Coarse = sol_coarse->node[Point_Coarse]->GetPrimitive();
    
    for (iVar = 0; iVar < nVar; iVar++)
      PrimVar_Temp[iVar] += PrimVar_Coarse[iVar];
    
    /*--- Here I have to set the updated primitive var back to coarse grid nodes. ---*/
    /*--- Note that the update primvar is actually I^(k+1)_k V_k -V_(k+1) i.e. the correction
     * on the primitive variables and not the actual value itself. This correction value is
     * stored in the primitive variable vector of the coarse mesh. ---*/
     
    sol_coarse->node[Point_Coarse]->SetPrimitiveMGCorr(PrimVar_Temp);
  }
  
  /*--- Reset updates to zero at dirichlet boundaries. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
	
    /*--- Velocity update is zero. ---*/
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX   ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL  ) ||
        (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW )) {
      
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        
       for (iDim = 0; iDim < nDim; iDim++) 
            sol_coarse->node[Point_Coarse]->SetPrimitiveMGCorr(iDim+1,0.0);

      }
    }
    /*--- Pressure update is zero. ---*/
    if ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW)) {

	   for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
		   Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
		   sol_coarse->node[Point_Coarse]->SetPrimitiveMGCorr(0,0.0);
	   }
	}  
  }
  
  /*--- Add MPI routine here. ---*/
  
  /*--- The update to primitive must be set to a central location that will be smoothed in the next routine. 
   * Here the code copied from the earlier routine saves the update in the residual vector. ---*/
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
	  Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Fine_GridFile<<Point_Fine<<"\t"<<geo_fine->node[Point_Fine]->GetCoord(0)<<"\t"<<geo_fine->node[Point_Fine]->GetCoord(1)<<"\t";
      for (iVar = 0; iVar < nVar; iVar++) {
		  PrimVar_Temp[iVar] = sol_coarse->node[Point_Coarse]->GetPrimitiveMGCorr(iVar);
		  //PrimVar_Temp[iVar] = (Area_Children/Area_Parent)*PrimVar_Temp[iVar];
		  Fine_GridFile<<PrimVar_Temp[iVar]<<"\t";
		  sol_fine->node[Point_Fine]->SetPrimitiveMGCorr(iVar,PrimVar_Temp[iVar]);
	  }
	  Fine_GridFile<<endl;
    }
    Coarse_GridFile<<Point_Coarse<<"\t"<<geo_coarse->node[Point_Coarse]->GetCoord(0)<<"\t"<<geo_coarse->node[Point_Coarse]->GetCoord(1)<<"\t";
    Coarse_GridFile<<sol_coarse->node[Point_Coarse]->GetPrimitiveMGCorr(0)<<"\t"<<sol_coarse->node[Point_Coarse]->GetPrimitiveMGCorr(1)<<endl;
  }
  
  delete [] PrimVar_Temp;
  Coarse_GridFile.close();
  Fine_GridFile.close();
}

void CMultiGridIntegration::SmoothProlongated_VarCorrection (unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                                          unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {
  su2double *Residual_Old, *Residual_Sum, *Residual, *Residual_i, *Residual_j;
  unsigned short iVar, iSmooth, iMarker, nneigh;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  
  const unsigned short nVar = solver->GetnPrimVar();
  if (val_nSmooth > 0) {
    
    Residual = new su2double [nVar];
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Residual_Old = solver->LinSysRes.GetBlock(iPoint);
      solver->node[iPoint]->SetResidual_Old(Residual_Old);
    }
    
    /*--- Jacobi iterations ---*/
    
    for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
        solver->node[iPoint]->SetResidualSumZero();
      
      /*--- Loop over Interior edges ---*/
      
      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
        iPoint = geometry->edge[iEdge]->GetNode(0);
        jPoint = geometry->edge[iEdge]->GetNode(1);
        
        Residual_i = solver->LinSysRes.GetBlock(iPoint);
        Residual_j = solver->LinSysRes.GetBlock(jPoint);
        
        /*--- Accumulate nearest neighbor Residual to Res_sum for each variable ---*/
        
        solver->node[iPoint]->AddResidual_Sum(Residual_j);
        solver->node[jPoint]->AddResidual_Sum(Residual_i);
      }
      
      /*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        nneigh = geometry->node[iPoint]->GetnPoint();
        Residual_Sum = solver->node[iPoint]->GetResidual_Sum();
        Residual_Old = solver->node[iPoint]->GetResidual_Old();
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] =(Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar])
          /(1.0 + val_smooth_coeff*su2double(nneigh));
        }
        solver->LinSysRes.SetBlock(iPoint, Residual);
      }
      
      /*--- Copy boundary values ---*/
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
        if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Residual_Old = solver->node[iPoint]->GetResidual_Old();
          solver->LinSysRes.SetBlock(iPoint, Residual_Old);
        }
    }
    
    delete [] Residual;
    
  }
}


void CMultiGridIntegration::SetProlongated_VarCorrection(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh) {
  unsigned long Point_Fine;
  unsigned short iVar;
  su2double *PrimVar_Fine, *MGCorrec_Fine;
  
  const unsigned short nPrimVar = sol_fine->GetnPrimVar();
  const unsigned short nVar = sol_fine->GetnVar();
  su2double factor = config->GetDamp_Correc_Prolong(); //pow(config->GetDamp_Correc_Prolong(), iMesh+1);
  //su2double factor = 0.0;
  
  su2double *PrimVar_New = new su2double [nPrimVar];
  su2double *Solution = new su2double [nVar];
  for (Point_Fine = 0; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
	  
	/*--- Update to primitive variables prolongated from coarse grid. ---*/  
    MGCorrec_Fine = sol_fine->node[Point_Fine]->GetPrimitiveMGCorr();
    
    /*--- Current primitive variables on the fine grid. ---*/
    PrimVar_Fine = sol_fine->node[Point_Fine]->GetPrimitive();
    
    for (iVar = 0; iVar < nPrimVar; iVar++) {
      
      /*--- Prevent a fine grid divergence due to a coarse grid divergence ---*/
      if (MGCorrec_Fine[iVar] != MGCorrec_Fine[iVar]) MGCorrec_Fine[iVar] = 0.0;
      
      PrimVar_New[iVar] = PrimVar_Fine[iVar] + factor*MGCorrec_Fine[iVar];
    }
    sol_fine->node[Point_Fine]->SetPrimitive(PrimVar_New);
    
    for (iVar = 0; iVar < nVar; iVar++) 
       Solution[iVar] = PrimVar_New[iVar+1]*PrimVar_New[nVar+1];
    
    sol_fine->node[Point_Fine]->SetSolution(Solution);
  }
  
  /*--- MPI the new interpolated solution ---*/
  //Add MPI routine here.
  
  delete [] Solution;
}

