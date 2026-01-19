/*!
 * \file CMultiGridIntegration.cpp
 * \brief Implementation of the multigrid integration class.
 * \author F. Palacios, T. Economon
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

namespace {
/*!
 * \brief Helper function to enforce Euler wall BC by projecting momentum to tangent plane.
 * \param[in] geo_coarse - Coarse grid geometry.
 * \param[in] config - Problem configuration.
 * \param[in,out] sol_coarse - Coarse grid solver (to access and modify solution/correction).
 * \param[in] use_solution_old - If true, project Solution_Old (corrections); if false, project Solution.
 */
void ProjectEulerWallToTangentPlane(CGeometry* geo_coarse, const CConfig* config, CSolver* sol_coarse,
                                    bool use_solution_old) {
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {

      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        const auto Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();

        if (!geo_coarse->nodes->GetDomain(Point_Coarse)) continue;

        /*--- Get coarse grid normal ---*/
        su2double Normal[3] = {0.0};
        geo_coarse->vertex[iMarker][iVertex]->GetNormal(Normal);
        const auto nDim = geo_coarse->GetnDim();
        su2double Area = GeometryToolbox::Norm(nDim, Normal);

        if (Area < EPS) continue;

        /*--- Normalize normal vector ---*/
        su2double UnitNormal[3] = {0.0};
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim] / Area;
        }

        /*--- Get current solution or correction ---*/
        su2double* solution_coarse = use_solution_old ?
          sol_coarse->GetNodes()->GetSolution_Old(Point_Coarse) :
          sol_coarse->GetNodes()->GetSolution(Point_Coarse);

        /*--- Compute normal component of momentum (v·n or correction·n) ---*/
        su2double momentum_n = 0.0;
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          momentum_n += solution_coarse[iDim + 1] * UnitNormal[iDim];
        }

        /*--- Project to tangent plane: solution_coarse -= (solution_coarse·n)n ---*/
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          solution_coarse[iDim + 1] -= momentum_n * UnitNormal[iDim];
        }
      }
      END_SU2_OMP_FOR
    }
  }
}
}  // anonymous namespace

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

  MultiGrid_Cycle(geometry, solver_container, numerics_container, config,
                  FinestMesh, RecursiveParam, RunTime_EqSystem, iZone, iInst);


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
                                            unsigned short iInst) {

  CConfig* config = config_container[iZone];

  const unsigned short Solver_Position = config->GetContainerPosition(RunTime_EqSystem);
  const bool classical_rk4 = (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT);
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Shorter names to refer to fine grid entities. ---*/

  CGeometry* geometry_fine = geometry[iZone][iInst][iMesh];
  CSolver** solver_container_fine = solver_container[iZone][iInst][iMesh];
  CSolver* solver_fine = solver_container_fine[Solver_Position];
  CNumerics** numerics_fine = numerics_container[iZone][iInst][iMesh][Solver_Position];

  /*--- Number of RK steps. ---*/

  unsigned short iRKLimit = 1;

  switch (config->GetKind_TimeIntScheme()) {
    case RUNGE_KUTTA_EXPLICIT:
      iRKLimit = config->GetnRKStep();
      break;
    case CLASSICAL_RK4_EXPLICIT:
      iRKLimit = 4;
      break;
    case EULER_EXPLICIT:
    case EULER_IMPLICIT:
      iRKLimit = 1;
      break;
  }

  /*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/

  for (unsigned short iPreSmooth = 0; iPreSmooth < config->GetMG_PreSmooth(iMesh); iPreSmooth++) {

    /*--- Time and space integration ---*/

    for (unsigned short iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {

      /*--- Send-Receive boundary conditions, and preprocessing ---*/

      solver_fine->Preprocessing(geometry_fine, solver_container_fine, config, iMesh, iRKStep, RunTime_EqSystem, false);


      if (iRKStep == 0) {

        /*--- Set the old solution ---*/

        solver_fine->Set_OldSolution();

        if (classical_rk4) solver_fine->Set_NewSolution();

        /*--- Compute time step, max eigenvalue, and integration scheme (steady and unsteady problems) ---*/

        solver_fine->SetTime_Step(geometry_fine, solver_container_fine, config, iMesh, config->GetTimeIter());

        /*--- Restrict the solution and gradient for the adjoint problem ---*/

        Adjoint_Setup(geometry, solver_container, config_container, RunTime_EqSystem, config->GetTimeIter(), iZone);

      }

      /*--- Space integration ---*/

      Space_Integration(geometry_fine, solver_container_fine, numerics_fine, config, iMesh, iRKStep, RunTime_EqSystem);

      /*--- Time integration, update solution using the old solution plus the solution increment ---*/

      Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);

      /*--- Send-Receive boundary conditions, and postprocessing ---*/

      solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);

    }

  }

  /*--- Compute Forcing Term $P_(k+1) = I^(k+1)_k(P_k+F_k(u_k))-F_(k+1)(I^(k+1)_k u_k)$ and update solution for multigrid ---*/

  if ( iMesh < config->GetnMGLevels() ) {

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

    SetRestricted_Solution(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh);

    solver_coarse->Preprocessing(geometry_coarse, solver_container_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem, false);

    Space_Integration(geometry_coarse, solver_container_coarse, numerics_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem);

    /*--- Store initial coarse grid residual before MG cycle ---*/
    su2double max_res_coarse_before = 0.0;
    for (unsigned long iPoint = 0; iPoint < geometry_coarse->GetnPointDomain(); iPoint++) {
      const su2double* res = solver_coarse->LinSysRes.GetBlock(iPoint);
      for (unsigned short iVar = 0; iVar < solver_coarse->GetnVar(); iVar++) {
        max_res_coarse_before = max(max_res_coarse_before, fabs(res[iVar]));
      }
    }

#ifdef HAVE_MPI
    su2double global_max_res_coarse_before = 0.0;
    SU2_MPI::Allreduce(&max_res_coarse_before, &global_max_res_coarse_before, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    max_res_coarse_before = global_max_res_coarse_before;
#endif

    /*--- Compute $P_(k+1) = I^(k+1)_k(r_k) - r_(k+1) ---*/

    SetForcing_Term(solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh+1);

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
                      iMesh+1, nextRecurseParam, RunTime_EqSystem, iZone, iInst);
    }

    /*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +  Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/

    GetProlongated_Correction(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);

    /*--- Measure coarse grid residual after MG cycle to assess convergence rate ---*/
    su2double max_res_coarse_after = 0.0;
    for (unsigned long iPoint = 0; iPoint < geometry_coarse->GetnPointDomain(); iPoint++) {
      const su2double* res = solver_coarse->LinSysRes.GetBlock(iPoint);
      for (unsigned short iVar = 0; iVar < solver_coarse->GetnVar(); iVar++) {
        max_res_coarse_after = max(max_res_coarse_after, fabs(res[iVar]));
      }
    }

#ifdef HAVE_MPI
    su2double global_max_res_coarse_after = 0.0;
    SU2_MPI::Allreduce(&max_res_coarse_after, &global_max_res_coarse_after, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    max_res_coarse_after = global_max_res_coarse_after;
#endif

    /*--- Adaptive CFL based on coarse grid convergence rate ---*/
    if (max_res_coarse_before > 1e-16) {
      su2double convergence_rate = max_res_coarse_after / max_res_coarse_before;

      /*--- Get current CFL values ---*/
      su2double CFL_fine = config->GetCFL(iMesh);
      su2double CFL_coarse_current = config->GetCFL(iMesh+1);
      su2double current_coeff = CFL_fine / CFL_coarse_current;

      /*--- Adaptive CFL for coarse multigrid levels ---*/
      constexpr int MAX_MG_LEVELS = 10;
      static int consecutive_decreasing[MAX_MG_LEVELS] = {};
      static su2double prev_fine_res[MAX_MG_LEVELS] = {};
      static su2double min_fine_res[MAX_MG_LEVELS] = {};
      static int iterations_since_improvement[MAX_MG_LEVELS] = {};
      static bool initialized = false;

      const unsigned short nMGLevels = config->GetnMGLevels();
      const unsigned short max_lvl = min(nMGLevels, (unsigned short)MAX_MG_LEVELS);

      if (!initialized) {
        for (int lvl = 0; lvl < MAX_MG_LEVELS; lvl++) {
          consecutive_decreasing[lvl] = 0;
          prev_fine_res[lvl] = 0.0;
          min_fine_res[lvl] = 1e10;
          iterations_since_improvement[lvl] = 0;
        }
        initialized = true;
      }

      unsigned short lvl = min(iMesh, (unsigned short)(max_lvl - 1));

      /*--- Track fine grid RMS residual ---*/
      su2double rms_res_fine = solver_fine->GetRes_RMS(0);
      static su2double smoothed_fine_res[MAX_MG_LEVELS] = {};

      if (smoothed_fine_res[lvl] < 1e-30) {
        smoothed_fine_res[lvl] = rms_res_fine;
        prev_fine_res[lvl] = rms_res_fine;
        min_fine_res[lvl] = rms_res_fine;
      } else {
        smoothed_fine_res[lvl] = 0.9 * rms_res_fine + 0.1 * smoothed_fine_res[lvl];
      }

      /*--- Track consecutive decreasing iterations ---*/
      if (smoothed_fine_res[lvl] < prev_fine_res[lvl]) {
        consecutive_decreasing[lvl] = min(consecutive_decreasing[lvl] + 1, 10);
      } else {
        consecutive_decreasing[lvl] = 0;
      }
      prev_fine_res[lvl] = smoothed_fine_res[lvl];

      /*--- Track best residual for divergence detection only ---*/
      if (smoothed_fine_res[lvl] < min_fine_res[lvl]) {
        min_fine_res[lvl] = smoothed_fine_res[lvl];
        iterations_since_improvement[lvl] = 0;
      } else {
        iterations_since_improvement[lvl]++;
      }

      /*--- Detect actual divergence (not stall during slow convergence) ---*/
      bool fine_grid_diverging = (smoothed_fine_res[lvl] > 1.5 * min_fine_res[lvl] && min_fine_res[lvl] > 1e-30);

      /*--- Sustained convergence: residual decreasing for 5 consecutive iterations ---*/
      bool sustained_convergence = (consecutive_decreasing[lvl] >= 5);

      /*--- Adapt coefficient based on convergence state ---*/
      su2double new_coeff = current_coeff;

      if (fine_grid_diverging) {
        new_coeff = current_coeff * 1.1;
      } else if (convergence_rate > 1.0) {
        new_coeff = current_coeff * (1.0 + 0.10 * min(convergence_rate - 1.0, 1.0));
      } else if (sustained_convergence) {
        new_coeff = current_coeff * 0.99;
      }

      /*--- Clamp and enforce hierarchy ---*/
      new_coeff = min(3.0, max(1.0, new_coeff));
      su2double min_safe_coeff = CFL_fine / (0.95 * CFL_fine);
      new_coeff = max(new_coeff, min_safe_coeff);

      /*--- Update coarse grid CFL if changed significantly ---*/
      su2double CFL_coarse_new = CFL_fine / new_coeff;

      if (fabs(CFL_coarse_new - CFL_coarse_current) > 1e-8) {
        config->SetCFL(iMesh+1, CFL_coarse_new);

        /*--- Update LocalCFL at each coarse grid point ---*/
        SU2_OMP_FOR_STAT(roundUpDiv(geometry_coarse->GetnPoint(), omp_get_num_threads()))
        for (auto iPoint = 0ul; iPoint < geometry_coarse->GetnPoint(); iPoint++) {
          solver_coarse->GetNodes()->SetLocalCFL(iPoint, CFL_coarse_new);
        }
        END_SU2_OMP_FOR
      }

      /*--- Output monitoring information periodically ---*/
      if (SU2_MPI::GetRank() == 0 && config->GetTimeIter() % 50 == 0) {
        cout << "  MG Level " << iMesh+1 << ": conv_rate = " << convergence_rate
             << ", consecutive = " << consecutive_decreasing[lvl]
             << ", sustained = " << (sustained_convergence ? "YES" : "NO")
             << ", coeff = " << new_coeff << ", CFL = " << CFL_coarse_new << endl;
      }
    }

    SmoothProlongated_Correction(RunTime_EqSystem, solver_fine, geometry_fine, config->GetMG_CorrecSmooth(iMesh), 1.25, config);

    SetProlongated_Correction(solver_fine, geometry_fine, config, iMesh);


    /*--- Solution post-smoothing in the prolongated grid. ---*/

    for (unsigned short iPostSmooth = 0; iPostSmooth < config->GetMG_PostSmooth(iMesh); iPostSmooth++) {

      for (unsigned short iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {

        solver_fine->Preprocessing(geometry_fine, solver_container_fine, config, iMesh, iRKStep, RunTime_EqSystem, false);

        if (iRKStep == 0) {
          solver_fine->Set_OldSolution();

          if (classical_rk4) solver_fine->Set_NewSolution();

          solver_fine->SetTime_Step(geometry_fine, solver_container_fine, config, iMesh,  config->GetTimeIter());
        }

        Space_Integration(geometry_fine, solver_container_fine, numerics_fine, config, iMesh, iRKStep, RunTime_EqSystem);

        Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);

        solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);

      }
    }
  }

}

void CMultiGridIntegration::GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                      CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  const su2double *Solution_Fine = nullptr, *Solution_Coarse = nullptr;

  const unsigned short nVar = sol_coarse->GetnVar();

  auto *Solution = new su2double[nVar];

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    su2double Area_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);

    for (auto iVar = 0u; iVar < nVar; iVar++) Solution[iVar] = 0.0;

    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      su2double Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
      Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      for (auto iVar = 0u; iVar < nVar; iVar++)
        Solution[iVar] -= Solution_Fine[iVar]*Area_Children/Area_Parent;
    }

    Solution_Coarse = sol_coarse->GetNodes()->GetSolution(Point_Coarse);

    for (auto iVar = 0u; iVar < nVar; iVar++)
      Solution[iVar] += Solution_Coarse[iVar];

    for (auto iVar = 0u; iVar < nVar; iVar++)
      sol_coarse->GetNodes()->SetSolution_Old(Point_Coarse,Solution);
  }
  END_SU2_OMP_FOR

  delete [] Solution;

  /*--- Enforce Euler wall BC on corrections by projecting to tangent plane ---*/
  ProjectEulerWallToTangentPlane(geo_coarse, config, sol_coarse, true);

  /*--- Remove any contributions from no-slip walls. ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {

      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        auto Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();

        /*--- For dirichlet boundary condtions, set the correction to zero.
         Note that Solution_Old stores the correction not the actual value ---*/

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
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      sol_fine->LinSysRes.SetBlock(Point_Fine, sol_coarse->GetNodes()->GetSolution_Old(Point_Coarse));
    }
  }
  END_SU2_OMP_FOR

}

void CMultiGridIntegration::SmoothProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                                         unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {

  /*--- Check if there is work to do. ---*/
  if (val_nSmooth == 0) return;

  const su2double *Residual_Old, *Residual_Sum, *Residual_j;

  const unsigned short nVar = solver->GetnVar();

  SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); iPoint++) {
    Residual_Old = solver->LinSysRes.GetBlock(iPoint);
    solver->GetNodes()->SetResidual_Old(iPoint,Residual_Old);
  }
  END_SU2_OMP_FOR

  /*--- Jacobi iterations. ---*/

  for (auto iSmooth = 0u; iSmooth < val_nSmooth; iSmooth++) {

    /*--- Loop over all mesh points (sum the residuals of direct neighbors). ---*/

    SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {

      solver->GetNodes()->SetResidualSumZero(iPoint);

      for (auto iNeigh = 0u; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
        auto jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
        Residual_j = solver->LinSysRes.GetBlock(jPoint);
        solver->GetNodes()->AddResidual_Sum(iPoint, Residual_j);
      }

    }
    END_SU2_OMP_FOR

    /*--- Loop over all mesh points (update residuals with the neighbor averages). ---*/

    SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {

      su2double factor = 1.0/(1.0+val_smooth_coeff*su2double(geometry->nodes->GetnPoint(iPoint)));

      Residual_Sum = solver->GetNodes()->GetResidual_Sum(iPoint);
      Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);

      for (auto iVar = 0u; iVar < nVar; iVar++)
        solver->LinSysRes(iPoint,iVar) = (Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar])*factor;
    }
    END_SU2_OMP_FOR

    /*--- Restore original residuals (without average) at boundary points. ---*/

    for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

        SU2_OMP_FOR_STAT(32)
        for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);
          solver->LinSysRes.SetBlock(iPoint, Residual_Old);
        }
        END_SU2_OMP_FOR
      }
    }

  }

}

void CMultiGridIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine,
                                                      CConfig *config, unsigned short iMesh) {
  su2double *Solution_Fine, *Residual_Fine;

  const unsigned short nVar = sol_fine->GetnVar();
  //const su2double factor = config->GetDamp_Correc_Prolong();

  /*--- Apply level-dependent damping: more aggressive damping on deeper coarse levels ---*/
  const su2double base_damp = config->GetDamp_Correc_Prolong();
  const su2double level_factor = pow(0.75, iMesh);  // 0.75^iMesh reduces damping progressively
  const su2double factor = base_damp; // * level_factor;

  /*--- Track maximum correction magnitude for monitoring ---*/
  su2double max_correction_local = 0.0;
  su2double max_correction_global = 0.0;

  SU2_OMP_FOR_STAT(roundUpDiv(geo_fine->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Fine = 0ul; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
    Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
    Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
    for (auto iVar = 0u; iVar < nVar; iVar++) {
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

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      sol_fine->GetNodes()->SetSolution(Point_Fine, sol_coarse->GetNodes()->GetSolution(Point_Coarse));
    }
  }
  END_SU2_OMP_FOR
}

void CMultiGridIntegration::SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                            CGeometry *geo_coarse, CConfig *config, unsigned short iMesh) {

  const su2double *Residual_Fine;

  const unsigned short nVar = sol_coarse->GetnVar();
  su2double factor = config->GetDamp_Res_Restric();

  auto *Residual = new su2double[nVar];

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);

    for (auto iVar = 0u; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (auto iVar = 0u; iVar < nVar; iVar++)
        Residual[iVar] += factor*Residual_Fine[iVar];
    }
    sol_coarse->GetNodes()->AddRes_TruncError(Point_Coarse, Residual);
  }
  END_SU2_OMP_FOR

  delete [] Residual;

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {
      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        auto Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->GetNodes()->SetVel_ResTruncError_Zero(Point_Coarse);
      }
      END_SU2_OMP_FOR
    }
  }

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->GetNodes()->SubtractRes_TruncError(Point_Coarse, sol_coarse->LinSysRes.GetBlock(Point_Coarse));
  }
  END_SU2_OMP_FOR

}

void CMultiGridIntegration::SetResidual_Term(CGeometry *geometry, CSolver *solver) {

  AD::StartNoSharedReading();
  SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPointDomain(), omp_get_num_threads()))
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); iPoint++)
    solver->LinSysRes.AddBlock(iPoint, solver->GetNodes()->GetResTruncError(iPoint));
  END_SU2_OMP_FOR
  AD::EndNoSharedReading();

}

void CMultiGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                   CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config, unsigned short iMesh) {

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

  /*--- Enforce Euler wall BC by projecting velocity to tangent plane ---*/
  ProjectEulerWallToTangentPlane(geo_coarse, config, sol_coarse, false);

  /*--- MPI the new interpolated solution ---*/

  sol_coarse->InitiateComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION);

}

void CMultiGridIntegration::SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                   CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {

  const unsigned short nDim = geo_coarse->GetnDim();
  const unsigned short nVar = sol_coarse->GetnVar();

  auto **Gradient = new su2double* [nVar];
  for (auto iVar = 0u; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double [nDim];

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPoint(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPoint(); Point_Coarse++) {
    su2double Area_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);

    for (auto iVar = 0u; iVar < nVar; iVar++)
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Gradient[iVar][iDim] = 0.0;

    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      unsigned long Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      su2double Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
      auto Gradient_fine = sol_fine->GetNodes()->GetGradient(Point_Fine);

      for (auto iVar = 0u; iVar < nVar; iVar++)
        for (auto iDim = 0u; iDim < nDim; iDim++)
          Gradient[iVar][iDim] += Gradient_fine[iVar][iDim]*Area_Children/Area_Parent;
    }
    sol_coarse->GetNodes()->SetGradient(Point_Coarse,Gradient);
  }
  END_SU2_OMP_FOR

  for (auto iVar = 0u; iVar < nVar; iVar++)
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
