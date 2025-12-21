/*!
 * \file CMultiGridIntegration.cpp
 * \brief Implementation of the multigrid integration class.
 * \author F. Palacios, T. Economon
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/integration/CMultiGridIntegration.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"

#include "ComputeLinSysResRMS.hpp"
#include <algorithm>

// Forward declaration of the adaptive damping helper so it can be used by
// MultiGrid_Cycle before the helper's definition further down the file.
static su2double Damp_Restric_Adapt(const unsigned short *lastPreSmoothIters,
                                   unsigned short lastPreSmoothCount,
                                   CConfig *config);


CMultiGridIntegration::CMultiGridIntegration() : CIntegration() { }

void CMultiGridIntegration::MultiGrid_Iteration(CGeometry ****geometry,
                                                CSolver *****solver_container,
                                                CNumerics ******numerics_container,
                                                CConfig **config,
                                                unsigned short RunTime_EqSystem,
                                                unsigned short iZone,
                                                unsigned short iInst) {

  bool direct;
  switch (config[iZone]->GetKind_Solver()) {
    case MAIN_SOLVER::EULER:
    case MAIN_SOLVER::NAVIER_STOKES:
    case MAIN_SOLVER::NEMO_EULER:
    case MAIN_SOLVER::NEMO_NAVIER_STOKES:
    case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::FEM_EULER:
    case MAIN_SOLVER::FEM_NAVIER_STOKES:
    case MAIN_SOLVER::FEM_RANS:
    case MAIN_SOLVER::FEM_LES:
    case MAIN_SOLVER::DISC_ADJ_EULER:
    case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER:
    case MAIN_SOLVER::DISC_ADJ_FEM_NS:
    case MAIN_SOLVER::DISC_ADJ_RANS:
      direct = true;
      break;
    default:
      direct = false;
      break;
  }

  const unsigned short Solver_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  /*--- Start an OpenMP parallel region covering the entire MG iteration, if the solver supports it. ---*/

  SU2_OMP_PARALLEL_(if(solver_container[iZone][iInst][MESH_0][Solver_Position]->GetHasHybridParallel()))
  {

  su2double monitor = 1.0;
  bool FullMG = false;

  unsigned short RecursiveParam = config[iZone]->GetMGCycle();

  if (config[iZone]->GetMGCycle() == FULLMG_CYCLE) {
    RecursiveParam = V_CYCLE;
    FullMG = true;
  }

  /*--- Full multigrid strategy and start up with fine grid only works with the direct problem ---*/

  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /// TODO: This was always false.
  const bool Convergence_FullMG = false;

  if (!config[iZone]->GetRestart() && FullMG && direct && ( Convergence_FullMG && (FinestMesh != MESH_0 ))) {

    SetProlongated_Solution(RunTime_EqSystem,
                            solver_container[iZone][iInst][FinestMesh-1][Solver_Position],
                            solver_container[iZone][iInst][FinestMesh][Solver_Position],
                            geometry[iZone][iInst][FinestMesh-1],
                            geometry[iZone][iInst][FinestMesh],
                            config[iZone]);

    SU2_OMP_SAFE_GLOBAL_ACCESS(config[iZone]->SubtractFinestMesh();)
  }

  /*--- Set the current finest grid (full multigrid strategy) ---*/

  FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Perform the Full Approximation Scheme multigrid ---*/

  // Allocate per-MG-run storage for the number of pre-smoothing iterations
  const unsigned short nMGLevels = config[iZone]->GetnMGLevels();
  unsigned short *lastPreSmoothIters = new unsigned short[nMGLevels + 1];
  for (unsigned short ii = 0; ii <= nMGLevels; ++ii) lastPreSmoothIters[ii] = 0;

  MultiGrid_Cycle(geometry, solver_container, numerics_container, config,
                  FinestMesh, RecursiveParam, RunTime_EqSystem, iZone, iInst,
                  lastPreSmoothIters);

  delete [] lastPreSmoothIters;


  /*--- Computes primitive variables and gradients in the finest mesh (useful for the next solver (turbulence) and output ---*/

  solver_container[iZone][iInst][MESH_0][Solver_Position]->Preprocessing(geometry[iZone][iInst][MESH_0],
                                                                         solver_container[iZone][iInst][MESH_0],
                                                                         config[iZone], MESH_0, NO_RK_ITER,
                                                                         RunTime_EqSystem, true);

  /*--- Compute non-dimensional parameters and the convergence monitor ---*/

  NonDimensional_Parameters(geometry[iZone][iInst], solver_container[iZone][iInst],
                            numerics_container[iZone][iInst], config[iZone],
                            FinestMesh, RunTime_EqSystem, &monitor);

  }
  END_SU2_OMP_PARALLEL

}

void CMultiGridIntegration::MultiGrid_Cycle(CGeometry ****geometry,
                                            CSolver *****solver_container,
                                            CNumerics ******numerics_container,
                                            CConfig **config_container,
                                            unsigned short iMesh,
                                            unsigned short RecursiveParam,
                                            unsigned short RunTime_EqSystem,
                                            unsigned short iZone,
                                            unsigned short iInst,
                                            unsigned short *lastPreSmoothIters) {

  CConfig* config = config_container[iZone];

  const unsigned short Solver_Position = config->GetContainerPosition(RunTime_EqSystem);
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Shorter names to refer to fine grid entities. ---*/

  CGeometry* geometry_fine = geometry[iZone][iInst][iMesh];
  CSolver** solver_container_fine = solver_container[iZone][iInst][iMesh];
  CSolver* solver_fine = solver_container_fine[Solver_Position];
  CNumerics** numerics_fine = numerics_container[iZone][iInst][iMesh][Solver_Position];

  /*--- Number of RK steps. ---*/

  unsigned short iRKLimit = config->GetnRKStep();

  // standard multigrid: pre-smoothing
  // start with solution on fine grid h
  // apply nonlinear smoother (e.g. Gauss-Seidel) to system to damp high-frequency errors.
  PreSmoothing(RunTime_EqSystem, geometry, solver_container, config_container, solver_fine, numerics_fine, geometry_fine, solver_container_fine, config, iMesh, iZone, iRKLimit, lastPreSmoothIters);


  /*--- Compute Forcing Term $P_(k+1) = I^(k+1)_k(P_k+F_k(u_k))-F_(k+1)(I^(k+1)_k u_k)$ and update solution for multigrid ---*/

  if ( iMesh < config->GetnMGLevels() ) {
    //cout << "imesh =" << iMesh << " " << config->GetnMGLevels() << endl;

    /*--- Shorter names to refer to coarse grid entities. ---*/

    CGeometry* geometry_coarse = geometry[iZone][iInst][iMesh+1];
    CSolver** solver_container_coarse = solver_container[iZone][iInst][iMesh+1];
    CSolver* solver_coarse = solver_container_coarse[Solver_Position];
    CNumerics** numerics_coarse = numerics_container[iZone][iInst][iMesh+1][Solver_Position];

    /*--- Temporarily disable implicit integration, for what follows we do not need the Jacobian. ---*/

    if (implicit) {
      SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetKind_TimeIntScheme(EULER_EXPLICIT);)
    }

    /*--- Compute $r_k = P_k + F_k(u_k)$ ---*/

    solver_fine->Preprocessing(geometry_fine, solver_container_fine, config, iMesh, NO_RK_ITER, RunTime_EqSystem, false);

    Space_Integration(geometry_fine, solver_container_fine, numerics_fine, config, iMesh, NO_RK_ITER, RunTime_EqSystem);

    SetResidual_Term(geometry_fine, solver_fine);

    /*--- Compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$ ---*/
    // restrict the fine-grid solution to the coarser grid. This must be FAS multigrid, because else we restrict the residual.
    SetRestricted_Solution(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);

    solver_coarse->Preprocessing(geometry_coarse, solver_container_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem, false);

    Space_Integration(geometry_coarse, solver_container_coarse, numerics_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem);

    /*--- Monitor coarse grid residual magnitude before forcing term ---*/
    if (SU2_MPI::GetRank() == 0) {
      su2double max_res_coarse = 0.0;
      for (unsigned long iPoint = 0; iPoint < geometry_coarse->GetnPointDomain(); iPoint++) {
        const su2double* res = solver_coarse->LinSysRes.GetBlock(iPoint);
        for (unsigned short iVar = 0; iVar < solver_coarse->GetnVar(); iVar++) {
          max_res_coarse = max(max_res_coarse, fabs(res[iVar]));
        }
      }
      cout << "  MG Level " << iMesh+1 << ": Max coarse residual before forcing = " << max_res_coarse << endl;
    }

    /*--- Compute $P_(k+1) = I^(k+1)_k(r_k) - r_(k+1) ---*/

    // Adapt damping based on recorded pre-smoothing iterations and apply to forcing term
    //cout << "calling damp_restric_Adapt" << endl;
    su2double adapted_factor = Damp_Restric_Adapt(lastPreSmoothIters, config->GetnMGLevels() + 1, config);
    SetForcing_Term(solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh+1, adapted_factor);

    /*--- Restore the time integration settings. ---*/

    if (implicit) {
      SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetKind_TimeIntScheme(EULER_IMPLICIT);)
    }

    /*--- Recursive call to MultiGrid_Cycle (this routine). ---*/

    for (unsigned short imu = 0; imu <= RecursiveParam; imu++) {

      unsigned short nextRecurseParam = RecursiveParam;
      if (iMesh == config->GetnMGLevels()-2)
        nextRecurseParam = 0;

      MultiGrid_Cycle(geometry, solver_container, numerics_container, config_container,
              iMesh+1, nextRecurseParam, RunTime_EqSystem, iZone, iInst, lastPreSmoothIters);
    }

    /*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +  Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/

    GetProlongated_Correction(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);
    SmoothProlongated_Correction(RunTime_EqSystem, solver_fine, geometry_fine, config->GetMG_CorrecSmooth(iMesh), config->GetMG_Smooth_Coeff(), config, iMesh);
    SetProlongated_Correction(solver_fine, geometry_fine, config, iMesh);


    /*--- Solution post-smoothing in the prolongated grid. ---*/
    PostSmoothing(RunTime_EqSystem, solver_fine, numerics_fine, geometry_fine, solver_container_fine, config, iMesh, iRKLimit);

  }

}


void CMultiGridIntegration::PreSmoothing(unsigned short RunTime_EqSystem,
CGeometry**** geometry,
CSolver***** solver_container,
CConfig **config_container,
CSolver* solver_fine,
CNumerics** numerics_fine,
CGeometry* geometry_fine,
CSolver** solver_container_fine,
CConfig *config,
unsigned short iMesh,
unsigned short iZone,
unsigned short iRKLimit,
unsigned short *lastPreSmoothIters) {
  const bool classical_rk4 = (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT);
  const unsigned short nPreSmooth = config->GetMG_PreSmooth(iMesh);
  const unsigned long timeIter = config->GetTimeIter();

  // Early exit settings from config
  const bool early_exit_enabled = config->GetMG_Smooth_EarlyExit();
  const su2double early_exit_threshold = config->GetMG_Smooth_Res_Threshold();
  const bool output_enabled = config->GetMG_Smooth_Output();

  su2double initial_rms = 0.0;
  if (early_exit_enabled || output_enabled) {
    initial_rms = ComputeLinSysResRMS(solver_fine, geometry_fine);
    if (output_enabled) {
      cout << "MG Pre-Smoothing Level " << iMesh << " Initial RMS: " << initial_rms << endl;
    }
  }

  /*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/
  unsigned short actual_iterations = nPreSmooth;
  for (unsigned short iPreSmooth = 0; iPreSmooth < nPreSmooth; iPreSmooth++) {
    /*--- Time and space integration ---*/
    for (unsigned short iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {
      /*--- Send-Receive boundary conditions, and preprocessing ---*/
      solver_fine->Preprocessing(geometry_fine, solver_container_fine, config, iMesh, iRKStep, RunTime_EqSystem, false);
      if (iRKStep == 0) {
        /*--- Set the old solution ---*/
        solver_fine->Set_OldSolution();
        if (classical_rk4) solver_fine->Set_NewSolution();
        solver_fine->SetTime_Step(geometry_fine, solver_container_fine, config, iMesh, timeIter);
        Adjoint_Setup(geometry, solver_container, config_container, RunTime_EqSystem, timeIter, iZone);
      }
      /*--- Space integration ---*/
      Space_Integration(geometry_fine, solver_container_fine, numerics_fine, config, iMesh, iRKStep, RunTime_EqSystem);
      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);
      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);
    }

    /*--- MPI sync after RK stage to ensure halos have updated solution for next smoothing iteration ---*/
    solver_fine->InitiateComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);
    solver_fine->CompleteComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);

    // Early exit check and output
    if (early_exit_enabled || output_enabled) {
      su2double current_rms = ComputeLinSysResRMS(solver_fine, geometry_fine);
      if (output_enabled) {
        cout << "MG Pre-Smoothing Level " << iMesh << " Iteration " << iPreSmooth + 1 << "/" << nPreSmooth << " RMS: " << current_rms << endl;
      }
      if (early_exit_enabled && current_rms < early_exit_threshold * initial_rms) {
        if (output_enabled) {
          cout << "MG Pre-Smoothing Level " << iMesh << " Early exit at iteration " << iPreSmooth + 1 << endl;
        }
        actual_iterations = iPreSmooth + 1;
        break;
      }
    }
  }

  // Record actual iterations performed for this MG level
  if (lastPreSmoothIters != nullptr) lastPreSmoothIters[iMesh] = actual_iterations;

  /*--- MPI communication after presmoothing to update solution at halo/boundary points ---*/
  solver_fine->InitiateComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);
  solver_fine->CompleteComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);
}


void CMultiGridIntegration::PostSmoothing(unsigned short RunTime_EqSystem, CSolver* solver_fine,CNumerics** numerics_fine, CGeometry* geometry_fine, CSolver** solver_container_fine, CConfig *config, unsigned short iMesh,unsigned short iRKLimit)
{
  const bool classical_rk4 = (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT);
  const unsigned short nPostSmooth = config->GetMG_PostSmooth(iMesh);
  const unsigned long timeIter = config->GetTimeIter();

  // Early exit settings from config
  const bool early_exit_enabled = config->GetMG_Smooth_EarlyExit();
  const su2double early_exit_threshold = config->GetMG_Smooth_Res_Threshold();
  const bool output_enabled = config->GetMG_Smooth_Output();

  su2double initial_rms = 0.0;
  if (early_exit_enabled || output_enabled) {
    initial_rms = ComputeLinSysResRMS(solver_fine, geometry_fine);
    if (output_enabled) {
      cout << "MG Post-Smoothing Level " << iMesh << " Initial RMS: " << initial_rms << endl;
    }
  }

  /*--- Do a postsmoothing on the grid iMesh after prolongation from the grid iMesh+1 ---*/
  for (unsigned short iPostSmooth = 0; iPostSmooth < nPostSmooth; iPostSmooth++) {
    for (unsigned short iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {
      solver_fine->Preprocessing(geometry_fine, solver_container_fine, config, iMesh, iRKStep, RunTime_EqSystem, false);
      if (iRKStep == 0) {
        /*--- Set the old solution ---*/
        solver_fine->Set_OldSolution();
        if (classical_rk4) solver_fine->Set_NewSolution();
        solver_fine->SetTime_Step(geometry_fine, solver_container_fine, config, iMesh,  timeIter);
      }
      /*--- Space integration ---*/
      Space_Integration(geometry_fine, solver_container_fine, numerics_fine, config, iMesh, iRKStep, RunTime_EqSystem);
      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);
      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);
    }

    /*--- MPI sync after RK stage to ensure halos have updated solution for next smoothing iteration ---*/
    solver_fine->InitiateComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);
    solver_fine->CompleteComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);

    // Early exit check and output
    if (early_exit_enabled || output_enabled) {
      su2double current_rms = ComputeLinSysResRMS(solver_fine, geometry_fine);
      if (output_enabled) {
        cout << "MG Post-Smoothing Level " << iMesh << " Iteration " << iPostSmooth + 1 << "/" << nPostSmooth << " RMS: " << current_rms << endl;
      }
      if (early_exit_enabled && current_rms < early_exit_threshold * initial_rms) {
        if (output_enabled) {
          cout << "MG Post-Smoothing Level " << iMesh << " Early exit at iteration " << iPostSmooth + 1 << endl;
        }
        break;
      }
    }
  }

  /*--- MPI communication after postsmoothing to update solution at halo/boundary points ---*/
  solver_fine->InitiateComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);
  solver_fine->CompleteComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION);
}


void CMultiGridIntegration::GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                      CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iMarker, iChildren, iVar;
  su2double Area_Parent, Area_Children;
  const su2double *Solution_Fine = nullptr, *Solution_Coarse = nullptr;

  const unsigned short nVar = sol_coarse->GetnVar();

  auto *Solution = new su2double[nVar];

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    Area_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);

    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;

    for (iChildren = 0; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
      Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] -= Solution_Fine[iVar]*Area_Children/Area_Parent;
    }

    Solution_Coarse = sol_coarse->GetNodes()->GetSolution(Point_Coarse);

    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] += Solution_Coarse[iVar];

    for (iVar = 0; iVar < nVar; iVar++)
      sol_coarse->GetNodes()->SetSolution_Old(Point_Coarse,Solution);
  }
  END_SU2_OMP_FOR

  delete [] Solution;

  /*--- Remove any contributions from no-slip walls. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {

      SU2_OMP_FOR_STAT(32)
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {

        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();

        /*--- For dirichlet boundary conditions, set the correction to zero.
         Note that Solution_Old stores the correction, not the actual value ---*/

        su2double zero[3] = {0.0};
        sol_coarse->GetNodes()->SetVelocity_Old(Point_Coarse, zero);

      }
      END_SU2_OMP_FOR
    }
  }

  /*--- MPI the set solution old ---*/

  sol_coarse->InitiateComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION_OLD);
  sol_coarse->CompleteComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION_OLD);

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (iChildren = 0; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      sol_fine->LinSysRes.SetBlock(Point_Fine, sol_coarse->GetNodes()->GetSolution_Old(Point_Coarse));
    }
  }
  END_SU2_OMP_FOR

}

void CMultiGridIntegration::SmoothProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                                         unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config, unsigned short iMesh) {

  /*--- Check if there is work to do. ---*/
  if (val_nSmooth == 0) return;

  const su2double *Residual_Old, *Residual_Sum, *Residual_j;
  unsigned short iVar, iSmooth, iMarker, iNeigh;
  unsigned long iPoint, jPoint, iVertex;

  const unsigned short nVar = solver->GetnVar();

  SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Residual_Old = solver->LinSysRes.GetBlock(iPoint);
    solver->GetNodes()->SetResidual_Old(iPoint,Residual_Old);
  }
  END_SU2_OMP_FOR

  /*--- Compute initial RMS for adaptive smoothing ---*/
  su2double initial_rms = 0.0;
  if (config->GetMG_Smooth_EarlyExit()) {
    initial_rms = ComputeLinSysResRMS(solver, geometry);
  }

  /*--- Output initial RMS if enabled ---*/
  if (config->GetMG_Smooth_Output()) {
    cout << "MG Correction-Smoothing Level " << iMesh << " Initial RMS: " << initial_rms << endl;
  }

  /*--- Jacobi iterations with adaptive early exit ---*/
  unsigned short actual_iterations = val_nSmooth;

  for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {

    /*--- Loop over domain points only (sum the residuals of direct neighbors). ---*/

    SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPointDomain(), omp_get_num_threads()))
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); ++iPoint) {

      solver->GetNodes()->SetResidualSumZero(iPoint);

      for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
        jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);

        /*--- Only include domain neighbors (skip halo points with stale LinSysRes) ---*/
        if (geometry->nodes->GetDomain(jPoint)) {
          Residual_j = solver->LinSysRes.GetBlock(jPoint);
          solver->GetNodes()->AddResidual_Sum(iPoint, Residual_j);
        }
      }

    }
    END_SU2_OMP_FOR

    /*--- Loop over domain points only (update residuals with the neighbor averages). ---*/

    SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPointDomain(), omp_get_num_threads()))
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); ++iPoint) {

      /*--- Count domain neighbors for proper averaging ---*/
      unsigned short nDomainNeighbors = 0;
      for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
        jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
        if (geometry->nodes->GetDomain(jPoint)) {
          nDomainNeighbors++;
        }
      }

      su2double factor = 1.0/(1.0+val_smooth_coeff*su2double(nDomainNeighbors));

      Residual_Sum = solver->GetNodes()->GetResidual_Sum(iPoint);
      Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);

      for (iVar = 0; iVar < nVar; iVar++)
        solver->LinSysRes(iPoint,iVar) = (Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar])*factor;
    }
    END_SU2_OMP_FOR

    /*--- Restore original residuals (without average) at boundary points. ---*/

    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

        SU2_OMP_FOR_STAT(32)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);
          solver->LinSysRes.SetBlock(iPoint, Residual_Old);
        }
        END_SU2_OMP_FOR
      }
    }

    /*--- Output RMS residual if enabled ---*/
    if (config->GetMG_Smooth_Output()) {
      const su2double RMS_Res = ComputeLinSysResRMS(solver, geometry);
      cout << "MG Correction-Smoothing Level " << iMesh << " Iteration " << iSmooth+1 << "/" << val_nSmooth << " RMS: " << RMS_Res << endl;
    }

    /*--- Adaptive early exit check after first iteration ---*/
    if (config->GetMG_Smooth_EarlyExit() && iSmooth == 0 && val_nSmooth > 1) {
      su2double current_rms = ComputeLinSysResRMS(solver, geometry);
      su2double reduction_ratio = current_rms / initial_rms;

      // If RMS reduction is sufficient (ratio <= threshold), additional iterations may not be necessary
      if (reduction_ratio <= config->GetMG_Smooth_Res_Threshold()) {
        if (config->GetMG_Smooth_Output()) {
          cout << "MG Correction-Smoothing Level " << iMesh << " Early exit: sufficient RMS reduction ("
               << reduction_ratio << " <= " << config->GetMG_Smooth_Res_Threshold() << ")" << endl;
        }
        actual_iterations = 1;  // Only do this one iteration
        break;
      }
      // If reduction is insufficient (ratio > threshold), continue with remaining iterations
    }

  }

  /*--- Log if iterations were skipped ---*/
  if (config->GetMG_Smooth_Output() && actual_iterations < val_nSmooth) {
    cout << "MG Correction-Smoothing Level " << iMesh << " completed " << actual_iterations
         << "/" << val_nSmooth << " iterations (adaptive early exit)" << endl;
  }
}

void CMultiGridIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine,
                                                      CConfig *config, unsigned short iMesh) {
  unsigned long Point_Fine;
  unsigned short iVar;
  su2double *Solution_Fine, *Residual_Fine;

  const unsigned short nVar = sol_fine->GetnVar();

  /*--- Apply level-dependent damping: more aggressive damping on deeper coarse levels ---*/
  const su2double base_damp = config->GetDamp_Correc_Prolong();
  const su2double level_factor = pow(0.75, iMesh);  // 0.75^iMesh reduces damping progressively
  const su2double factor = base_damp * level_factor;

  /*--- Track maximum correction magnitude for monitoring ---*/
  su2double max_correction_local = 0.0;
  su2double max_correction_global = 0.0;

  SU2_OMP_FOR_STAT(roundUpDiv(geo_fine->GetnPointDomain(), omp_get_num_threads()))
  for (Point_Fine = 0; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
    Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
    Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
    for (iVar = 0; iVar < nVar; iVar++) {
      /*--- Prevent a fine grid divergence due to a coarse grid divergence ---*/
      if (Residual_Fine[iVar] != Residual_Fine[iVar])
        Residual_Fine[iVar] = 0.0;

      su2double correction = factor * Residual_Fine[iVar];
      Solution_Fine[iVar] += correction;

      /*--- Track maximum correction ---*/
      SU2_OMP_CRITICAL
      max_correction_local = max(max_correction_local, fabs(correction));
    }
  }
  END_SU2_OMP_FOR

  /*--- Reduce maximum correction across all ranks ---*/
  SU2_MPI::Allreduce(&max_correction_local, &max_correction_global, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

  /*--- Output correction magnitude for monitoring ---*/
  if (SU2_MPI::GetRank() == 0) {
    cout << "  MG Level " << iMesh << ": Max correction = " << max_correction_global
         << ", Damping factor = " << factor << " (base=" << base_damp
         << " × level_factor=" << level_factor << ")" << endl;
  }

  /*--- MPI the new interpolated solution ---*/

  sol_fine->InitiateComms(geo_fine, config, MPI_QUANTITIES::SOLUTION);
  sol_fine->CompleteComms(geo_fine, config, MPI_QUANTITIES::SOLUTION);

}

void CMultiGridIntegration::SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                    CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iChildren;

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (iChildren = 0; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      sol_fine->GetNodes()->SetSolution(Point_Fine, sol_coarse->GetNodes()->GetSolution(Point_Coarse));
    }
  }
  END_SU2_OMP_FOR
}

void CMultiGridIntegration::SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                            CGeometry *geo_coarse, CConfig *config, unsigned short iMesh,
                                            su2double val_factor) {

  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iMarker, iVar, iChildren;
  const su2double *Residual_Fine;

  const unsigned short nVar = sol_coarse->GetnVar();
  su2double factor = config->GetDamp_Res_Restric();
  if (val_factor > 0.0) factor = val_factor;

  auto *Residual = new su2double[nVar];

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);

    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

    for (iChildren = 0; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] += factor*Residual_Fine[iVar];
    }
    sol_coarse->GetNodes()->AddRes_TruncError(Point_Coarse, Residual);
  }
  END_SU2_OMP_FOR

  delete [] Residual;

  /*--- Zero momentum components of truncation error on viscous walls ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {
      SU2_OMP_FOR_STAT(32)
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->GetNodes()->SetVel_ResTruncError_Zero(Point_Coarse);
      }
      END_SU2_OMP_FOR
    }
  }

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->GetNodes()->SubtractRes_TruncError(Point_Coarse, sol_coarse->LinSysRes.GetBlock(Point_Coarse));
  }
  END_SU2_OMP_FOR

}

/*
 * Damp_Restric_Adapt
 * ------------------
 * Compute an adapted restriction damping factor (Damp_Res_Restric) based on the
 * number of pre-smoothing iterations actually performed on each multigrid level.
 *
 * Arguments:
 *  - lastPreSmoothIters: array of length at least (nMGLevels+1) containing the
 *    recorded number of pre-smoothing iterations performed at each MG level.
 *    The convention used here is that index == level (0..nMGLevels).
 *  - lastPreSmoothCount: length of the array passed in (for safety checks).
 *  - config: the CConfig pointer used to read current settings (no write performed).
 *
 * Returns:
 *  - the adapted damping factor (may be equal to current factor if no change).
 *
 * Notes:
 *  - This routine performs MPI Allreduce operations to make the adapt decision
 *    globally consistent across ranks.
 *  - The function intentionally does NOT write back into CConfig; it returns the
 *    new factor so the caller can decide how to apply it (e.g., via a setter
 *    once available or by passing it into the next call).
 */
static su2double Damp_Restric_Adapt(const unsigned short *lastPreSmoothIters,
                                   unsigned short lastPreSmoothCount,
                                   CConfig *config) {
  //cout << "starting Damp_Restric_Adapt" << endl;
  if (config == nullptr) return 0.0;

  const unsigned short nMGLevels = config->GetnMGLevels();

  // Safety: require the provided array to have at least nMGLevels+1 entries
  if (lastPreSmoothIters == nullptr || lastPreSmoothCount < (nMGLevels + 1)) {
    // Not enough information to adapt; return current factor unchanged
    return config->GetDamp_Res_Restric();
  }

  int local_any_full = 0; // true if any inspected coarse level reached configured max
  int local_all_one = 1;  // true if all inspected coarse levels performed exactly 1 iter
  int local_inspected = 0; // number of coarse levels that have recorded a pre-smooth value

  // Inspect all coarse-grid levels (levels 1..nMGLevels). Ignore entries
  // that are still zero (not yet visited during this MG run).
  for (unsigned short lvl = 1; lvl <= nMGLevels; ++lvl) {
    const unsigned short performed = lastPreSmoothIters[lvl];
    const unsigned short configured = config->GetMG_PreSmooth(lvl);
    //cout << "  Level " << lvl << ": performed=" << performed
    //     << " configured=" << configured << endl;
    if (performed == 0) continue; // skip un-inspected level
    ++local_inspected;
    if (performed >= configured) local_any_full = 1;
    if (performed != 1) local_all_one = 0;
  }
  //cout << "local_inspected=" << local_inspected
  //     << " local_any_full=" << local_any_full
  //     << " local_all_one=" << local_all_one << endl;


    // Debug output: local decision stats
    //cout << "Damp_Restric_Adapt local: inspected=" << local_inspected
    //  << " any_full=" << local_any_full << " all_one=" << local_all_one << endl;

    // Make decision globally consistent across MPI ranks
  int global_any_full = 0;
  int global_all_one = 0;
  int global_inspected = 0;

  // Sum inspected counts across ranks; if nobody inspected any levels, skip adaptation
  SU2_MPI::Allreduce(&local_inspected, &global_inspected, 1, MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  if (global_inspected == 0) {
    // No ranks have inspected coarse levels yet — do not adapt
    return config->GetDamp_Res_Restric();
  }

    SU2_MPI::Allreduce(&local_any_full, &global_any_full, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&local_all_one, &global_all_one, 1, MPI_INT, MPI_MIN, SU2_MPI::GetComm());


  const su2double current = config->GetDamp_Res_Restric();
  //cout << "Current Damp_Res_Restric: " << current << endl;
  su2double new_factor = current;

  const su2double scale_down = 0.99;
  const su2double scale_up = 1.01;
  const su2double clamp_min = 0.1;
  // larger than this is bad for stability
  const su2double clamp_max = 0.95;

  if (global_any_full) {
    new_factor = current * scale_down;
  } else if (global_all_one) {
    new_factor = current * scale_up;
  } else {
    // no change
    return current;
  }

  // Clamp result
  new_factor = std::min(std::max(new_factor, clamp_min), clamp_max);

  // Optionally log the change on all ranks (config controls verbosity)
  if (config->GetMG_Smooth_Output()) {
    cout << "Adaptive Damp_Res_Restric: " << current << " -> " << new_factor << endl;
  }
  // Persist the adapted factor into the runtime config so subsequent cycles
  // observe the updated value. This is safe because the adapt decision was
  // already made through global MPI reductions, so calling the setter on all
  // ranks yields the same value everywhere.
  config->SetDamp_Res_Restric(new_factor);
  //cout << "ending Damp_Restric_Adapt, damping = " <<config->GetDamp_Res_Restric()  << endl;
  return new_factor;
}


void CMultiGridIntegration::SetResidual_Term(CGeometry *geometry, CSolver *solver) {

  AD::StartNoSharedReading();
  SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPointDomain(), omp_get_num_threads()))
  for (unsigned long iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    solver->LinSysRes.AddBlock(iPoint, solver->GetNodes()->GetResTruncError(iPoint));
  END_SU2_OMP_FOR
  AD::EndNoSharedReading();

}

void CMultiGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                   CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {

  const unsigned short Solver_Position = config->GetContainerPosition(RunTime_EqSystem);
  const bool grid_movement = config->GetGrid_Movement();

  /*--- Compute coarse solution from fine solution ---*/

  CSolver::MultigridRestriction(*geo_fine, sol_fine->GetNodes()->GetSolution(),
                                *geo_coarse, sol_coarse->GetNodes()->GetSolution());

  /*--- Update the solution at the no-slip walls ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {

      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {

        const auto Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();

        if (Solver_Position == FLOW_SOL) {

          /*--- At moving walls, set the solution based on the new density and wall velocity ---*/

          if (grid_movement) {
            const auto* Grid_Vel = geo_coarse->nodes->GetGridVel(Point_Coarse);
            sol_coarse->GetNodes()->SetVelSolutionVector(Point_Coarse, Grid_Vel);
          }
          else {
            // nijso asks: why only for the velocity?
            /*--- For stationary no-slip walls, set the velocity to zero. ---*/
            su2double zero[3] = {0.0};
            sol_coarse->GetNodes()->SetVelSolutionVector(Point_Coarse, zero);
          }

        }

        if (Solver_Position == ADJFLOW_SOL) {
          sol_coarse->GetNodes()->SetVelSolutionDVector(Point_Coarse);
        }

      }
      END_SU2_OMP_FOR
    }
  }

  /*--- MPI the new interpolated solution ---*/

  sol_coarse->InitiateComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION);

}

void CMultiGridIntegration::SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                   CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iVar, iDim, iChildren;
  su2double Area_Parent, Area_Children;

  const unsigned short nDim = geo_coarse->GetnDim();
  const unsigned short nVar = sol_coarse->GetnVar();

  auto **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double [nDim];

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPoint(), omp_get_num_threads()))
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPoint(); Point_Coarse++) {
    Area_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);

    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Gradient[iVar][iDim] = 0.0;

    for (iChildren = 0; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
      auto Gradient_fine = sol_fine->GetNodes()->GetGradient(Point_Fine);

      for (iVar = 0; iVar < nVar; iVar++)
        for (iDim = 0; iDim < nDim; iDim++)
          Gradient[iVar][iDim] += Gradient_fine[iVar][iDim]*Area_Children/Area_Parent;
    }
    sol_coarse->GetNodes()->SetGradient(Point_Coarse,Gradient);
  }
  END_SU2_OMP_FOR

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;

}

void CMultiGridIntegration::NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container,
                                                      CNumerics ****numerics_container, CConfig *config,
                                                      unsigned short FinestMesh, unsigned short RunTime_EqSystem,
                                                      su2double *monitor) {
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  switch (RunTime_EqSystem) {

    case RUNTIME_FLOW_SYS:

      /*--- Calculate the inviscid and viscous forces ---*/

      solver_container[FinestMesh][FLOW_SOL]->Pressure_Forces(geometry[FinestMesh], config);
      solver_container[FinestMesh][FLOW_SOL]->Momentum_Forces(geometry[FinestMesh], config);
      solver_container[FinestMesh][FLOW_SOL]->Friction_Forces(geometry[FinestMesh], config);

      /*--- Calculate the turbo performance ---*/
      if (config->GetBoolTurbomachinery()){

        /*--- Average quantities at the inflow and outflow boundaries ---*/

        solver_container[FinestMesh][FLOW_SOL]->TurboAverageProcess(solver_container[FinestMesh], geometry[FinestMesh],config,INFLOW);
        solver_container[FinestMesh][FLOW_SOL]->TurboAverageProcess(solver_container[FinestMesh], geometry[FinestMesh], config, OUTFLOW);

        /*--- Gather Inflow and Outflow quantities on the Master Node to compute performance ---*/

        solver_container[FinestMesh][FLOW_SOL]->GatherInOutAverageValues(config, geometry[FinestMesh]);

      }

      break;

    case RUNTIME_ADJFLOW_SYS:

      /*--- Calculate the inviscid and viscous sensitivities ---*/

      solver_container[FinestMesh][ADJFLOW_SOL]->Inviscid_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh],
                                                 numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);

      solver_container[FinestMesh][ADJFLOW_SOL]->Viscous_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh],
                                                 numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);

      /*--- Smooth the inviscid and viscous sensitivities ---*/

      if (config->GetKind_SensSmooth() != NONE)
        solver_container[FinestMesh][ADJFLOW_SOL]->Smooth_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh],
                                                   numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
      break;
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CMultiGridIntegration::Adjoint_Setup(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                                          unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) {

  if ((RunTime_EqSystem != RUNTIME_ADJFLOW_SYS) || (Iteration != 0)) return;

  for (unsigned short iMGLevel = 0; iMGLevel <= config[iZone]->GetnMGLevels(); iMGLevel++) {

    /*--- Set the time step in all the MG levels ---*/

    solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTime_Step(geometry[iZone][INST_0][iMGLevel],
                                                                      solver_container[iZone][INST_0][iMGLevel],
                                                                      config[iZone], iMGLevel, Iteration);

    /*--- Set the force coefficients ---*/

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CD(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CD());
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CL(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CL());
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CT(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CT());
      solver_container[iZone][INST_0][iMGLevel][FLOW_SOL]->SetTotal_CQ(solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->GetTotal_CQ());
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    /*--- Restrict solution and gradients to the coarse levels ---*/

    if (iMGLevel != config[iZone]->GetnMGLevels()) {
      SetRestricted_Solution(RUNTIME_FLOW_SYS,
                             solver_container[iZone][INST_0][iMGLevel][FLOW_SOL],
                             solver_container[iZone][INST_0][iMGLevel+1][FLOW_SOL],
                             geometry[iZone][INST_0][iMGLevel],
                             geometry[iZone][INST_0][iMGLevel+1],
                             config[iZone]);
//        ToDo: The flow solvers do not use the conservative variable gradients
//        SetRestricted_Gradient(RUNTIME_FLOW_SYS, solver_container[iZone][INST_0][iMGLevel][FLOW_SOL],
//                               solver_container[iZone][INST_0][iMGLevel+1][FLOW_SOL],
//                               geometry[iZone][INST_0][iMGLevel],
//                               geometry[iZone][INST_0][iMGLevel+1],
//                               config[iZone]);
    }

  }

}
