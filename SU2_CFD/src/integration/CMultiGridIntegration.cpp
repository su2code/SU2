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

#include "../../include/integration/ComputeLinSysResRMS.hpp"
#include "../../include/integration/CMultiGridIntegration.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"

/*!\cond PRIVATE Helper: shared logic for adapting a single MG damping factor.
 *  Inputs:
 *    performed[]  - actual iteration counts per level from this cycle
 *    progress[]   - whether residuals decreased per level
 *    getConfigured - returns the per-level configured maximum
 *    levelStart/End - level range to inspect
 *    getCurrent    - returns the current damping factor from config
 *    setPersist    - persists the updated factor back to config
 \endcond */
template <typename GetCfg, typename GetCur, typename SetPersist>
static void adaptMGDampingFactor(const unsigned short* performed,
                                  const bool* progress,
                                  GetCfg getConfigured,
                                  unsigned short levelStart, unsigned short levelEnd,
                                  GetCur getCurrent, SetPersist setPersist) {
  int local_any_stagnant = 0;  /*--- hit max iters AND residuals did not decrease: scale down. ---*/
  int local_all_early = 1;  /*--- all levels exited before max iters: scale up. ---*/
  int local_inspected = 0;

  for (unsigned short lvl = levelStart; lvl <= levelEnd; ++lvl) {
    const unsigned short configured = getConfigured(lvl);
    if (configured == 0) continue;
    ++local_inspected;
    const bool hit_max = (performed[lvl] >= configured);
    /*--- Scale-down signal: hit the cap AND residuals did not improve.
     *    Hitting the cap while still converging is not stagnation; no reduction needed. ---*/
    if (hit_max && !progress[lvl]) local_any_stagnant = 1;
    /*--- Scale-up signal: requires early exit on every level.
     *    Making progress at max iters is good, but the damping is already doing useful work;
     *    do not increase it further until the smoother actually exits early. ---*/
    if (hit_max) local_all_early = 0;
  }
  if (local_inspected == 0) return;

  /*--- performed[] and progress[] are derived from MPI-reduced ComputeLinSysResRMS values,
   *    so local_any_stagnant and local_all_early are already identical on every rank. ---*/
  const su2double SCALE_DOWN = 0.99;
  const su2double SCALE_UP = 1.01;
  const su2double CLAMP_MIN = 0.1;
  const su2double CLAMP_MAX = 0.95;

  su2double factor = getCurrent();
  if (local_any_stagnant) factor *= SCALE_DOWN;
  else if (local_all_early) factor *= SCALE_UP;
  /*--- else: hit max iters but still converging, or mixed — hold factor. ---*/
  factor = max(CLAMP_MIN, min(CLAMP_MAX, factor));
  setPersist(factor);
}

void CMultiGridIntegration::adaptRestrictionDamping(CConfig* config) {
  const auto& mgOpts = config->GetMGOptions();
  const unsigned short nMGLevels = config->GetnMGLevels();
  adaptMGDampingFactor(
    lastPreSmoothIters,
    lastPreSmoothProgress,
    [&mgOpts](unsigned short lvl){ return mgOpts.MG_PreSmooth[lvl]; },
    /*levelStart=*/1, nMGLevels,
    [config](){ return config->GetDamp_Res_Restric(); },
    [config](su2double v){ config->SetDamp_Res_Restric(v); });
}

void CMultiGridIntegration::adaptProlongationDamping(CConfig* config) {
  /*--- Post-smoothing directly measures whether the corrected fine-grid solution
   *    is well-behaved after prolongation.  If it exits early, increase damping.
   *    If it stagnates at max iters, decrease damping. ---*/
  const auto& mgOpts = config->GetMGOptions();
  const unsigned short nMGLevels = config->GetnMGLevels();
  if (nMGLevels == 0) return;
  adaptMGDampingFactor(
    lastPostSmoothIters,
    lastPostSmoothProgress,
    [&mgOpts](unsigned short lvl){ return mgOpts.MG_PostSmooth[lvl]; },
    /*levelStart=*/0, static_cast<unsigned short>(nMGLevels - 1),
    [config](){ return config->GetDamp_Correc_Prolong(); },
    [config](su2double v){ config->SetDamp_Correc_Prolong(v); });
}

passivedouble CMultiGridIntegration::computeMultigridCFL(CConfig* config, unsigned short iMesh,
                                                          passivedouble CFL_fine, passivedouble CFL_coarse_current,
                                                          passivedouble rms_res_coarse) {

  const bool wasActive = AD::BeginPassive();

  passivedouble current_coeff = CFL_coarse_current / CFL_fine;

  /*--- Adaptive CFL using Exponential Moving Average (EMA) ---*/
  constexpr int AVG_WINDOW = 5;

  passivedouble CFL_coarse_new = CFL_coarse_current; // Default: keep current value

  /*--- Get global iteration count first ---*/
  unsigned long current_iter;
  if (config->GetTime_Domain())
    current_iter = config->GetTimeIter();
  else
    current_iter = config->GetInnerIter();

  /*--- Reset state at the beginning of a new solve (iter 0 or 1) ---*/
  /*--- This ensures deterministic behavior across multiple runs ---*/
  if (current_iter <= 1 && last_reset_iter != current_iter) {
    for (int i = 0; i < MAX_MG_LEVELS; i++) {
      current_avg[i] = 0.0;
      prev_avg[i] = 0.0;
      last_res[i] = 0.0;
      last_was_increase[i] = false;
      oscillation_count[i] = 0;
      last_check_iter[i] = 0;
      last_update_iter[i] = 0;
    }
    last_reset_iter = current_iter;
  }

  unsigned short lvl = min(iMesh, (unsigned short)(MAX_MG_LEVELS - 1));
  unsigned long iter = current_iter;

  /*--- rms_res_coarse is passed in (from lastPreSmoothRMS, already MPI-reduced). ---*/

  /*--- Flip-flop detection: detect oscillating residuals (once per outer iteration) ---*/
  bool oscillation_detected = false;
  if (iter != last_check_iter[lvl]) {
    last_check_iter[lvl] = iter;

    if (last_res[lvl] > EPS) {
      bool current_is_increase = (rms_res_coarse > last_res[lvl]);
      if (current_is_increase != last_was_increase[lvl]) {
        /*--- Direction changed, increment oscillation counter ---*/
        oscillation_count[lvl]++;
        if (oscillation_count[lvl] >= 4) {
          /*--- Detected 4 consecutive direction changes = oscillation ---*/
          oscillation_detected = true;
          oscillation_count[lvl] = 0;  // Reset counter after detecting
        }
      } else {
        /*--- Same direction, reset counter ---*/
        oscillation_count[lvl] = 0;
      }
      last_was_increase[lvl] = current_is_increase;
    }
    last_res[lvl] = rms_res_coarse;
  }

  /*--- Update exponential moving average ---*/
  if (current_avg[lvl] < EPS) {
    current_avg[lvl] = rms_res_coarse;  // Initialize with first value
  } else {
    current_avg[lvl] = (current_avg[lvl] * (AVG_WINDOW - 1) + rms_res_coarse) / AVG_WINDOW;
  }

  /*--- Check if we should compare and adapt CFL ---*/
  passivedouble new_coeff = current_coeff;
  const passivedouble MIN_REDUCTION_FACTOR = 0.98;  // Require at least 2% reduction
  const int UPDATE_INTERVAL = 5;  // Update reference every N iterations

  /*--- Initialize prev_avg on first use ---*/
  if (prev_avg[lvl] < EPS) {
    prev_avg[lvl] = current_avg[lvl];
  }

  /*--- Periodically update prev_avg to allow ratio to reflect accumulated decrease ---*/
  bool should_update = (iter - last_update_iter[lvl] >= UPDATE_INTERVAL);

  /*--- Asymmetric adaptation for robustness ---*/
  if (prev_avg[lvl] > EPS) {
    passivedouble ratio = current_avg[lvl] / prev_avg[lvl];
    bool sufficient_decrease = (ratio < MIN_REDUCTION_FACTOR);
    bool increasing_trend = (ratio >= 1.0);

    if (increasing_trend) {
      /*--- Residual increasing: reduce CFL immediately for robustness ---*/
      new_coeff = current_coeff * 0.90;
      /*--- Update reference since we're reacting immediately ---*/
      prev_avg[lvl] = current_avg[lvl];
      last_update_iter[lvl] = iter;
    } else if (sufficient_decrease && should_update) {
      /*--- Residual decreasing sufficiently: increase CFL ---*/
      new_coeff = current_coeff * 1.05;
      /*--- Update reference only when we actually increase CFL ---*/
      prev_avg[lvl] = current_avg[lvl];
      last_update_iter[lvl] = iter;
    }
  }

  /*--- CFL reduction for oscillation detection ---*/
  if (oscillation_detected) {
    new_coeff = current_coeff * 0.75;
    /*--- Update reference after oscillation response ---*/
    prev_avg[lvl] = current_avg[lvl];
    last_update_iter[lvl] = iter;
  }

    /*--- Clamp coefficient between 0.5 and 1.0 ---*/
    new_coeff = max(0.5, min(1.0, new_coeff));

    /*--- Update coarse grid CFL ---*/
    CFL_coarse_new = max(0.5 * CFL_fine, min(CFL_fine, CFL_fine * new_coeff));

    config->SetCFL(iMesh+1, CFL_coarse_new);

  AD::EndPassive(wasActive);
  return CFL_coarse_new;
}

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

  unsigned short RecursiveParam = static_cast<unsigned short>(config[iZone]->GetMGCycle());

  if (config[iZone]->GetMGCycle() == MG_CYCLE::FULL) {
    RecursiveParam = static_cast<unsigned short>(MG_CYCLE::V);
    FullMG = true;
  }

  /*--- Full multigrid strategy and start up with fine grid only works with the direct problem ---*/

  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Initialize per-level smoothing iteration counters to the configured maximum.
   *    If early exit never fires, the output will show actual == max. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    const unsigned short nMGLevels = config[iZone]->GetnMGLevels();
    const auto& mgOpts = config[iZone]->GetMGOptions();
    for (unsigned short i = 0; i <= nMGLevels; ++i) {
      lastPreSmoothIters[i] = mgOpts.MG_PreSmooth[i];
      lastPostSmoothIters[i] = mgOpts.MG_PostSmooth[i];
      lastCorrecSmoothIters[i] = mgOpts.MG_CorrecSmooth[i];
      lastPreSmoothProgress[i] = false;
      lastPostSmoothProgress[i] = false;
      lastCorrecSmoothProgress[i] = false;
      lastPreSmoothRMS[i][0] = lastPreSmoothRMS[i][1] = 0.0;
      lastPostSmoothRMS[i][0] = lastPostSmoothRMS[i][1] = 0.0;
      lastCorrecSmoothRMS[i][0] = lastCorrecSmoothRMS[i][1] = 0.0;
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Full MG: advance to the next finer grid after a fixed number of
   *    outer iterations on the current coarsest active level.
   *    We use 100 iterations per level (nMGLevels levels total) ---*/
  const bool Convergence_FullMG =
      FullMG && (FinestMesh != MESH_0) &&
      (config[iZone]->GetInnerIter() % 100 == 99);

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

  /*--- Adapt coarse-grid CFL once per cycle using smoothing residuals gathered during the cycle.
   *    lastPreSmoothRMS[iMesh+1][1] is the final RMS after pre-smoothing at the coarse level ---*/
  const unsigned short nMGLevels = config[iZone]->GetnMGLevels();
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    for (unsigned short iMesh = FinestMesh; iMesh < nMGLevels; ++iMesh) {
      const passivedouble CFL_fine_p = SU2_TYPE::GetValue(config[iZone]->GetCFL(iMesh));
      const passivedouble CFL_coarse_p = SU2_TYPE::GetValue(config[iZone]->GetCFL(iMesh+1));
      computeMultigridCFL(config[iZone], iMesh, CFL_fine_p, CFL_coarse_p,
                          lastPreSmoothRMS[iMesh+1][1]);
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Propagate the updated coarse-grid CFL to every coarse-grid point (all threads). ---*/
  for (unsigned short iMesh = FinestMesh; iMesh < nMGLevels; ++iMesh) {
    const passivedouble CFL_coarse_new = SU2_TYPE::GetValue(config[iZone]->GetCFL(iMesh+1));
    CGeometry* geo_c = geometry[iZone][iInst][iMesh+1];
    CSolver* sol_c = solver_container[iZone][iInst][iMesh+1][Solver_Position];
    SU2_OMP_FOR_STAT(roundUpDiv(geo_c->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geo_c->GetnPoint(); iPoint++)
      sol_c->GetNodes()->SetLocalCFL(iPoint, CFL_coarse_new);
    END_SU2_OMP_FOR
  }

  /*--- Computes primitive variables and gradients in the finest mesh (useful for the next solver (turbulence) and output ---*/

  solver_container[iZone][iInst][MESH_0][Solver_Position]->Preprocessing(geometry[iZone][iInst][MESH_0],
                                                                         solver_container[iZone][iInst][MESH_0],
                                                                         config[iZone], MESH_0, NO_RK_ITER,
                                                                         RunTime_EqSystem, true);

  /*--- Compute non-dimensional parameters and the convergence monitor ---*/

  NonDimensional_Parameters(geometry[iZone][iInst], solver_container[iZone][iInst],
                            numerics_container[iZone][iInst], config[iZone],
                            FinestMesh, RunTime_EqSystem, &monitor);

  /*--- Adapt restriction damping based on coarse-level pre-smoothing workload from this cycle.
   *    Only effective when MG_SMOOTH_EARLY_EXIT= YES (otherwise all levels always run to completion
   *    and the signal would always point to "scale down"). ---*/
  const auto& mgOptsZone = config[iZone]->GetMGOptions();
  if (mgOptsZone.MG_Smooth_EarlyExit) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      adaptRestrictionDamping(config[iZone]);
      adaptProlongationDamping(config[iZone]);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Print compact smoothing summary when MG_SMOOTH_OUTPUT= YES. ---*/
  if (mgOptsZone.MG_Smooth_Output) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      if (SU2_MPI::GetRank() == MASTER_NODE) {

        /*--- Helper: format one cell as "act/max [init->final]". ---*/
        auto cellStr = [](unsigned short act, unsigned short mx,
                          su2double rms0, su2double rms1) -> std::string {
          std::ostringstream ss;
          ss << act << "/" << mx << " ["
             << std::scientific << std::setprecision(2)
             << rms0 << "->" << rms1 << "]";
          return ss.str();
        };

        PrintingToolbox::CTablePrinter table(&std::cout);
        table.AddColumn("Smoother", 13);
        for (unsigned short i = 0; i <= nMGLevels; ++i)
          table.AddColumn("Level " + std::to_string(i), 26);
        table.PrintHeader();

        /*--- Pre-smooth: defined on all levels 0..nMGLevels. ---*/
        table << "Pre-smooth";
        for (unsigned short i = 0; i <= nMGLevels; ++i)
          table << cellStr(lastPreSmoothIters[i], mgOptsZone.MG_PreSmooth[i],
                           lastPreSmoothRMS[i][0], lastPreSmoothRMS[i][1]);

        /*--- Post-smooth: defined on levels 0..nMGLevels-1; coarsest has none. ---*/
        table << "Post-smooth";
        for (unsigned short i = 0; i < nMGLevels; ++i)
          table << cellStr(lastPostSmoothIters[i], mgOptsZone.MG_PostSmooth[i],
                           lastPostSmoothRMS[i][0], lastPostSmoothRMS[i][1]);
        table << "-";

        /*--- Corr.-smooth: defined on levels 0..nMGLevels-1; coarsest has none. ---*/
        table << "Corr-smooth";
        for (unsigned short i = 0; i < nMGLevels; ++i)
          table << cellStr(lastCorrecSmoothIters[i], mgOptsZone.MG_CorrecSmooth[i],
                           lastCorrecSmoothRMS[i][0], lastCorrecSmoothRMS[i][1]);
        table << "-";

        table.PrintFooter();

        cout << std::fixed << std::setprecision(4)
             << "Damping [restrict | prolong] : " << config[iZone]->GetDamp_Res_Restric()
             << " | " << config[iZone]->GetDamp_Correc_Prolong() << "\n"
             << std::defaultfloat << std::setprecision(6);
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

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
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Shorter names to refer to fine grid entities. ---*/

  CGeometry* geometry_fine = geometry[iZone][iInst][iMesh];
  CSolver** solver_container_fine = solver_container[iZone][iInst][iMesh];
  CSolver* solver_fine = solver_container_fine[Solver_Position];
  CNumerics** numerics_fine = numerics_container[iZone][iInst][iMesh][Solver_Position];

  /*--- Number of RK steps. ---*/

  unsigned short iRKLimit = config->GetnRKStep();

  /*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/

  PreSmoothing(RunTime_EqSystem, geometry, solver_container, config_container, solver_fine, numerics_fine,
               geometry_fine, solver_container_fine, config, iMesh, iZone, iRKLimit);

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

    SetRestricted_Solution(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);

    solver_coarse->Preprocessing(geometry_coarse, solver_container_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem, false);

    Space_Integration(geometry_coarse, solver_container_coarse, numerics_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem);

    /*--- Compute $P_(k+1) = I^(k+1)_k(r_k) - r_(k+1) ---*/

    SetForcing_Term(solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh+1);

    /*--- Restore the time integration settings. ---*/

    if (implicit) {
      SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetKind_TimeIntScheme(EULER_IMPLICIT);)
    }

    /*--- Recursive call to MultiGrid_Cycle (this routine). ---*/
    /*--- Execute multigrid cycles sequentially to ensure deterministic recursion order ---*/
    /*--- This prevents accumulation of floating-point variations across recursive calls ---*/

    for (unsigned short imu = 0; imu <= RecursiveParam; imu++) {

      unsigned short nextRecurseParam = RecursiveParam;
      if (iMesh == config->GetnMGLevels()-2)
        nextRecurseParam = 0;

      MultiGrid_Cycle(geometry, solver_container, numerics_container, config_container,
                      iMesh+1, nextRecurseParam, RunTime_EqSystem, iZone, iInst);
    }

    /*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +  Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/

    GetProlongated_Correction(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);

    const auto& mgOpts = config->GetMGOptions();
    SmoothProlongated_Correction(RunTime_EqSystem, solver_fine, geometry_fine, mgOpts.MG_CorrecSmooth[iMesh], mgOpts.MG_Smooth_Coeff, config, iMesh);

    SetProlongated_Correction(solver_fine, geometry_fine, config, iMesh);

    /*--- Solution post-smoothing in the prolongated grid. ---*/

    PostSmoothing(RunTime_EqSystem, solver_fine, numerics_fine, geometry_fine, solver_container_fine,
                  config, iMesh, iRKLimit);
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
                                         unsigned short iRKLimit) {

  const auto& mgOpts = config->GetMGOptions();
  const bool classical_rk4 = (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT);
  const unsigned short nPreSmooth = mgOpts.MG_PreSmooth[iMesh];
  const unsigned long timeIter = config->GetTimeIter();
  const bool early_exit = mgOpts.MG_Smooth_EarlyExit && (nPreSmooth > 1);

  /*--- Reset the shared early-exit flag (master only). ---*/
  SU2_OMP_SAFE_GLOBAL_ACCESS(mg_early_exit_flag = false;)
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

      /*--- Capture initial RMS after the very first residual evaluation.
       *    This is the earliest point where LinSysRes = R(u_current) (not stale).
       *    ComputeLinSysResRMS must be called by all threads (uses parallel dot). ---*/
      if (iPreSmooth == 0 && iRKStep == 0) {
        const passivedouble initial_rms = ComputeLinSysResRMS(solver_fine);
        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
        {
          lastPreSmoothRMS[iMesh][0] = initial_rms;
          if (early_exit) mg_initial_smooth_rms = initial_rms;
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS
      }

      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);

      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);
    }

    /*--- Early exit: check if RMS has dropped sufficiently.
     *    ComputeLinSysResRMS must be called by all threads.
     *    only master uses the result inside the safe block. ---*/
    if (early_exit) {
      const passivedouble current_rms = ComputeLinSysResRMS(solver_fine);
      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      {
        mg_last_smooth_rms = current_rms;
        if (mg_last_smooth_rms < mgOpts.MG_Smooth_Res_Threshold * mg_initial_smooth_rms) {
          lastPreSmoothIters[iMesh] = iPreSmooth + 1;
          mg_early_exit_flag = true;
        }
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS
      if (mg_early_exit_flag) break;
    }
  }

  /*--- Record final RMS and progress flag.
   *    In the early-exit path mg_last_smooth_rms already holds the current value;
   *    in the normal path we compute it once here.
   *    The condition is the same for all threads so they all agree on whether to call. ---*/
  passivedouble final_pre_rms = mg_last_smooth_rms;
  if (!(early_exit && mg_early_exit_flag))
    final_pre_rms = ComputeLinSysResRMS(solver_fine);
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    mg_last_smooth_rms = final_pre_rms;
    lastPreSmoothRMS[iMesh][1] = final_pre_rms;
    lastPreSmoothProgress[iMesh] = mg_early_exit_flag ||
                                   (final_pre_rms < lastPreSmoothRMS[iMesh][0]);
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}


void CMultiGridIntegration::PostSmoothing(unsigned short RunTime_EqSystem,
                                          CSolver* solver_fine,
                                          CNumerics** numerics_fine,
                                          CGeometry* geometry_fine,
                                          CSolver** solver_container_fine,
                                          CConfig *config,
                                          unsigned short iMesh,
                                          unsigned short iRKLimit) {

  const auto& mgOpts = config->GetMGOptions();
  const bool classical_rk4 = (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT);
  const unsigned short nPostSmooth = mgOpts.MG_PostSmooth[iMesh];
  const unsigned long timeIter = config->GetTimeIter();
  const bool early_exit = mgOpts.MG_Smooth_EarlyExit && (nPostSmooth > 1);

  /*--- Reset the shared early-exit flag (master only). ---*/
  SU2_OMP_SAFE_GLOBAL_ACCESS(mg_early_exit_flag = false;)

  /*--- Do a postsmoothing on the grid iMesh after prolongation from the grid iMesh+1 ---*/
  for (unsigned short iPostSmooth = 0; iPostSmooth < nPostSmooth; iPostSmooth++) {

    for (unsigned short iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {
      solver_fine->Preprocessing(geometry_fine, solver_container_fine, config, iMesh, iRKStep, RunTime_EqSystem, false);
      if (iRKStep == 0) {

        /*--- Set the old solution ---*/
        solver_fine->Set_OldSolution();

        if (classical_rk4) solver_fine->Set_NewSolution();
        solver_fine->SetTime_Step(geometry_fine, solver_container_fine, config, iMesh, timeIter);
      }

      /*--- Space integration ---*/
      Space_Integration(geometry_fine, solver_container_fine, numerics_fine, config, iMesh, iRKStep, RunTime_EqSystem);

      /*--- Capture initial RMS after the very first residual evaluation.
       *    Before this point, LinSysRes held the smoothed correction (from SmoothProlongated_Correction),
       *    NOT the spatial residual R(u). This is the first valid R(u) after applying the correction.
       *    ComputeLinSysResRMS must be called by all threads (uses parallel dot). ---*/
      if (iPostSmooth == 0 && iRKStep == 0) {
        const passivedouble initial_rms = ComputeLinSysResRMS(solver_fine);
        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
        {
          lastPostSmoothRMS[iMesh][0] = initial_rms;
          if (early_exit) mg_initial_smooth_rms = initial_rms;
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS
      }

      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);

      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);

    }

    /*--- Early exit: check if RMS has dropped sufficiently.
     *    ComputeLinSysResRMS must be called by all threads (uses parallel dot). ---*/
    if (early_exit) {
      const passivedouble current_rms = ComputeLinSysResRMS(solver_fine);
      BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
      {
        mg_last_smooth_rms = current_rms;
        if (mg_last_smooth_rms < mgOpts.MG_Smooth_Res_Threshold * mg_initial_smooth_rms) {
          lastPostSmoothIters[iMesh] = iPostSmooth + 1;
          mg_early_exit_flag = true;
        }
      }
      END_SU2_OMP_SAFE_GLOBAL_ACCESS
      if (mg_early_exit_flag) break;
    }
  }

  /*--- Record final RMS after post-smoothing.
   *    In the early-exit path mg_last_smooth_rms is already current; otherwise compute once.
   *    The condition is the same for all threads so they all agree on whether to call. ---*/
  passivedouble final_post_rms = mg_last_smooth_rms;
  if (!(early_exit && mg_early_exit_flag))
    final_post_rms = ComputeLinSysResRMS(solver_fine);
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    mg_last_smooth_rms = final_post_rms;
    lastPostSmoothRMS[iMesh][1] = final_post_rms;
    lastPostSmoothProgress[iMesh] = mg_early_exit_flag ||
                                    (final_post_rms < lastPostSmoothRMS[iMesh][0]);
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
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

    /*--- Accumulate children contributions with stable ordering ---*/
    /*--- Process all children in sequential order to ensure deterministic FP summation ---*/
    auto nChildren = geo_coarse->nodes->GetnChildren_CV(Point_Coarse);
    for (auto iChildren = 0u; iChildren < nChildren; iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      su2double Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
      Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      su2double weight = Area_Children / Area_Parent;
      for (auto iVar = 0u; iVar < nVar; iVar++)
        Solution[iVar] -= Solution_Fine[iVar] * weight;
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
  sol_coarse->MultigridProjectEulerWall(geo_coarse, config, true);

  /*--- Remove any contributions from no-slip walls. ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {

      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        auto Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();

        /*--- For dirichlet boundary conditions, set the correction to zero.
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
                                                         unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config,
                                                         unsigned short iMesh) {

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

  /*--- Record initial correction norm for debugging output.
   *    ComputeLinSysResRMS must be called by all threads (uses parallel dot). ---*/
  {
    const passivedouble initial_corr_rms = ComputeLinSysResRMS(solver);
    SU2_OMP_SAFE_GLOBAL_ACCESS(lastCorrecSmoothRMS[iMesh][0] = initial_corr_rms;)
  }

  /*--- Jacobi iterations (no early exit — Jacobi targets high-frequency modes,
   *    so the global RMS is not a meaningful convergence indicator). ---*/

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

  /*--- Record final correction norm for debugging output.
   *    ComputeLinSysResRMS must be called by all threads (uses parallel dot). ---*/
    const passivedouble final_corr_rms = ComputeLinSysResRMS(solver);
    SU2_OMP_SAFE_GLOBAL_ACCESS(lastCorrecSmoothRMS[iMesh][1] = final_corr_rms;)

}

void CMultiGridIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine,
                                                      CConfig *config, unsigned short iMesh) {
  su2double *Solution_Fine, *Residual_Fine;

  const unsigned short nVar = sol_fine->GetnVar();

  /*--- Level-dependent damping: coarser prolongations produce noisier corrections
   *    due to larger cell-size jumps, so we reduce the factor progressively.
   *    iMesh=0: factor = base_damp * 1.0  (finest grid, full correction)
   *    iMesh=1: factor = base_damp * 0.75
   *    iMesh=2: factor = base_damp * 0.5625, etc. ---*/
  const su2double base_damp = config->GetDamp_Correc_Prolong();
  const su2double level_factor = pow(0.75, static_cast<su2double>(iMesh));
  const su2double factor = base_damp * level_factor;

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
    }
  }
  END_SU2_OMP_FOR

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
  const su2double factor = config->GetDamp_Res_Restric();

  auto *Residual = new su2double[nVar];

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);

    for (auto iVar = 0u; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (auto iVar = 0u; iVar < nVar; iVar++)
        Residual[iVar] += factor * Residual_Fine[iVar];
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
  sol_coarse->MultigridProjectEulerWall(geo_coarse, config, false);

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

      /*--- Calculate the turbo performance (only on the fine grid; turbo
       *    geometry data is only available on MESH_0). ---*/
      if (config->GetBoolTurbomachinery() && FinestMesh == MESH_0){

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
