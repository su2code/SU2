/*!
 * \file CMultiGridIntegration.cpp
 * \brief Implementation of the multigrid integration class.
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
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
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include <algorithm>


namespace {

/*!\cond PRIVATE
 *  Global-trend damping update.  Uses the cross-cycle EMA ratio to detect
 *  long-term divergence or convergence and adjusts both damping factors
 *  from a single, smooth signal.
 *
 *  crossCycleRatio < LO  : SCALE_UP   (residual below EMA trend)
 *  crossCycleRatio >= HI : SCALE_DOWN (residual above EMA trend)
 *  [LO, HI)              : no change  (neutral zone)
 \endcond */
static su2double applyGlobalTrend(su2double factor, passivedouble crossCycleRatio) {
  constexpr passivedouble SCALE_DOWN = 0.92;
  constexpr passivedouble SCALE_UP   = 1.02;
  constexpr passivedouble CLAMP_MIN  = 0.10;
  constexpr passivedouble CLAMP_MAX  = 0.90;
  constexpr passivedouble LO         = 0.95;  // ratio below this: converging, increase damping
  constexpr passivedouble HI         = 1.05;  // ratio above this: diverging, decrease damping

  if      (crossCycleRatio >= HI) factor *= SCALE_DOWN;
  else if (crossCycleRatio <  LO) factor *= SCALE_UP;
  return max(su2double{CLAMP_MIN}, min(su2double{CLAMP_MAX}, factor));
}

inline passivedouble ComputeLinSysResRMS(const CSolver* solver) {
  passivedouble result = 0;
  for (unsigned short iVar = 0; iVar < solver->GetnVar(); ++iVar) {
    result += pow(SU2_TYPE::GetValue(solver->GetRes_RMS(iVar)), 2);
  }
  return sqrt(result);
}

}  // anonymous namespace

void CMultiGridIntegration::adaptDampingFactors(CConfig* config, passivedouble crossCycleRatio) {
  SU2_ZONE_SCOPED
  if (config->GetnMGLevels() == 0) return;

  /*--- Both factors share the same global-trend signal.  The EMA already filters
   *    per-cycle noise; no per-level aggregation or floor counter needed. ---*/
  config->SetDamp_Res_Restric(
    applyGlobalTrend(config->GetDamp_Res_Restric(), crossCycleRatio));
  config->SetDamp_Correc_Prolong(
    applyGlobalTrend(config->GetDamp_Correc_Prolong(), crossCycleRatio));
}

CMultiGridIntegration::CMultiGridIntegration() : CIntegration() { }

void CMultiGridIntegration::MonitorFullMG_Startup(const vector<pair<string, passivedouble> >& convFields,
                                                 const CConfig *config) {

  const auto& mgOpts = config->GetMGOptions();

  /*--- No residual to promote on, MG_STARTUP_ITER and MG_STARTUP_STAGNATION are all that is left. ---*/

  if (convFields.empty()) return;

  const auto nFields = convFields.size();

  /*--- Fields are log10, so MG_STARTUP_CONVERGENCE is added to the level's starting value.
   *    All of them have to have dropped. ---*/

  if (mg_startup_conv_start.empty()) {
    for (const auto& field : convFields) mg_startup_conv_start.push_back(field.second);
  }
  else {
    const passivedouble drop = SU2_TYPE::GetValue(mgOpts.MG_Startup_Convergence);
    bool converged = true;
    for (auto iField = 0ul; iField < nFields; iField++)
      converged = converged && (convFields[iField].second <= mg_startup_conv_start[iField] + drop);

    if (converged && (mg_startup_promote_reason == MGStartupPromote::NONE))
      mg_startup_promote_reason = MGStartupPromote::CONVERGENCE;
  }

  /*--- Ratio of successive residuals, so a difference in log10. One slow iteration
   *    is not stagnation, and a field still coming down means the level is not stalled. ---*/

  const passivedouble stall_tol = SU2_TYPE::GetValue(mgOpts.MG_Startup_Stagnation);
  if (stall_tol > 0.0 && !mg_startup_conv_prev.empty()) {
    bool stalled = true;
    for (auto iField = 0ul; iField < nFields; iField++)
      stalled = stalled && (convFields[iField].second - mg_startup_conv_prev[iField] >= log10(stall_tol));

    if (stalled)
      mg_startup_stall_count++;
    else
      mg_startup_stall_count = 0;

    if ((mgOpts.MG_Startup_Stagnation_Iter > 0) &&
        (mg_startup_stall_count >= mgOpts.MG_Startup_Stagnation_Iter) &&
        (mg_startup_promote_reason == MGStartupPromote::NONE)) {
      mg_startup_promote_reason = MGStartupPromote::STAGNATION;
    }
  }

  mg_startup_conv_prev.clear();
  for (const auto& field : convFields) mg_startup_conv_prev.push_back(field.second);
}

void CMultiGridIntegration::SetCoarseGridCFL(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                                             unsigned short RunTime_EqSystem, unsigned short iZone,
                                             unsigned short iInst, unsigned short FinestMesh, bool FullMG,
                                             bool mesh0_ramp_window) {
  SU2_ZONE_SCOPED

  const unsigned short Solver_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  const unsigned short nMGLevels = config[iZone]->GetnMGLevels();
  const unsigned long startup_iter = config[iZone]->GetMGOptions().MG_Startup_Iter;
  CSolver* sol_f = solver_container[iZone][iInst][MESH_0][Solver_Position];

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    /*--- The active level is ramped, MESH_0 included, until it reaches its target. ---*/

    const bool ramping = FullMG && (FinestMesh <= nMGLevels) && (startup_iter > 0) &&
                         ((FinestMesh > MESH_0) || mesh0_ramp_window);
    const unsigned long iter_in_level = config[iZone]->GetInnerIter() - mg_ramp_level_start_iter;
    passivedouble ramp_progress = 1.0;
    if (ramping) {
      ramp_progress = min(passivedouble{1.0}, passivedouble(iter_in_level) / passivedouble(startup_iter));
    }

    mg_ramp_mesh0_active = FullMG_Mesh0Ramping(config[iZone], FinestMesh, FullMG);

    /*--- Coarse targets are derived from the level-0 CFL via MG_CFL_SCALING. ---*/
    passivedouble cfl_base = mesh0_ramp_window ? passivedouble{0.0} : SU2_TYPE::GetValue(sol_f->GetAvg_CFL_Local());
    if (cfl_base < EPS)
      cfl_base = SU2_TYPE::GetValue(config[iZone]->GetCFL(MESH_0));

    const auto& cflScaling = config[iZone]->GetMGOptions().MG_CflScaling;

    passivedouble CFL_target[MAX_MG_LEVELS+1];
    passivedouble CFL_scale[MAX_MG_LEVELS+1];
    CFL_target[0] = cfl_base;
    CFL_scale[0] = 1.0;
    for (unsigned short lvl = 1; lvl <= nMGLevels; ++lvl) {
      /*--- Entry lvl-1 is the transition lvl-1 -> lvl, clamped so a coarse level never runs
       *    at a higher CFL than the grid above it. ---*/
      const unsigned short iScale = lvl - 1;
      const passivedouble scale =
          max(passivedouble{1e-6}, min(passivedouble{1.0}, SU2_TYPE::GetValue(cflScaling[iScale])));
      CFL_scale[lvl] = scale;
      CFL_target[lvl] = CFL_target[lvl-1] * scale;
    }

    /*--- Ramp from the CFL the level below handed over at, which it may never have reached its
     *    target. The coarsest level has none and starts from its own target scaled down once more. ---*/

    const passivedouble CFL_handover =
        (mg_ramp_cfl_start > 0.0) ? mg_ramp_cfl_start : CFL_target[nMGLevels] * CFL_scale[nMGLevels];

    passivedouble CFL_mesh0 = CFL_target[MESH_0];

    for (unsigned short lvl = 0; lvl <= nMGLevels; ++lvl) {
      passivedouble CFL_local = CFL_target[lvl];
      if (ramping && (lvl == FinestMesh)) {
        CFL_local = (passivedouble(1.0) - ramp_progress) * CFL_handover + ramp_progress * CFL_target[lvl];
      }
      /*--- The configured level-0 CFL stays as the user wrote it, only the solver value is ramped. ---*/
      if (lvl > MESH_0) config[iZone]->SetCFL(lvl, CFL_local);
      else CFL_mesh0 = CFL_local;
    }
    /*--- Stashed in the solver's CFL stats for the point loop below to read back. ---*/
    if (mesh0_ramp_window) sol_f->SetCFL_Local_Stats(CFL_mesh0);
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Propagate the ramped CFL to every point of the finest grid. ---*/
  if (mesh0_ramp_window) {
    CGeometry* geo_f = geometry[iZone][iInst][MESH_0];
    const passivedouble cfl_mesh0 = SU2_TYPE::GetValue(sol_f->GetAvg_CFL_Local());
    SU2_OMP_FOR_STAT(roundUpDiv(geo_f->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geo_f->GetnPoint(); iPoint++)
      sol_f->GetNodes()->SetLocalCFL(iPoint, cfl_mesh0);
    END_SU2_OMP_FOR
  }

  /*--- Propagate the updated CFL to every coarse-grid point. ---*/
  for (unsigned short iMesh = 1; iMesh <= nMGLevels; ++iMesh) {
    const passivedouble CFL_coarse_new = SU2_TYPE::GetValue(config[iZone]->GetCFL(iMesh));
    CGeometry* geo_c = geometry[iZone][iInst][iMesh];
    CSolver* sol_c = solver_container[iZone][iInst][iMesh][Solver_Position];
    SU2_OMP_SAFE_GLOBAL_ACCESS(sol_c->SetCFL_Local_Stats(CFL_coarse_new);)
    SU2_OMP_FOR_STAT(roundUpDiv(geo_c->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geo_c->GetnPoint(); iPoint++)
      sol_c->GetNodes()->SetLocalCFL(iPoint, CFL_coarse_new);
    END_SU2_OMP_FOR
  }

  /*--- These solvers scale the flow's time step by CFL_scalar/CFL_flow, so keep the two in step
   *    on the active level for as long as the startup owns the flow's CFL. ---*/

  if ((Solver_Position != FLOW_SOL) || !FullMG || ((FinestMesh == MESH_0) && !mesh0_ramp_window)) return;

  CGeometry* geo_a = geometry[iZone][iInst][FinestMesh];
  const passivedouble cfl_flow =
      SU2_TYPE::GetValue(solver_container[iZone][iInst][FinestMesh][FLOW_SOL]->GetAvg_CFL_Local());

  for (const auto Scalar_Position : {TURB_SOL, TRANS_SOL, SPECIES_SOL}) {
    CSolver* sol_s = solver_container[iZone][iInst][FinestMesh][Scalar_Position];
    if (sol_s == nullptr) continue;

    const su2double cfl = cfl_flow * SU2_TYPE::GetValue((Scalar_Position == SPECIES_SOL)
                                                        ? config[iZone]->GetCFLRedCoeff_Species()
                                                        : config[iZone]->GetCFLRedCoeff_Turb());
    SU2_OMP_SAFE_GLOBAL_ACCESS(sol_s->SetCFL_Local_Stats(cfl);)
    SU2_OMP_FOR_STAT(roundUpDiv(geo_a->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geo_a->GetnPoint(); iPoint++)
      sol_s->GetNodes()->SetLocalCFL(iPoint, cfl);
    END_SU2_OMP_FOR
  }
}

void CMultiGridIntegration::MultiGrid_Iteration(CGeometry ****geometry,
                                                CSolver *****solver_container,
                                                CNumerics ******numerics_container,
                                                CConfig **config,
                                                unsigned short RunTime_EqSystem,
                                                unsigned short iZone,
                                                unsigned short iInst) {
  SU2_ZONE_SCOPED

  bool direct;
  switch (config[iZone]->GetKind_Solver()) {
    case MAIN_SOLVER::EULER:
    case MAIN_SOLVER::NAVIER_STOKES:
    case MAIN_SOLVER::INC_EULER:
    case MAIN_SOLVER::INC_NAVIER_STOKES:
    case MAIN_SOLVER::INC_RANS:
    case MAIN_SOLVER::DISC_ADJ_INC_EULER:
    case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_INC_RANS:
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

  auto RecursiveParam = static_cast<unsigned short>(config[iZone]->GetMGCycle());

  if (config[iZone]->GetMGCycle() == MG_CYCLE::FULL) {
    RecursiveParam = static_cast<unsigned short>(MG_CYCLE::V);
    FullMG = true;
  }

  /*--- Full multigrid strategy and start up with fine grid only works with the direct problem ---*/

  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Initialize per-level smoothing diagnostics for the current cycle. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    const unsigned short nMGLevels = config[iZone]->GetnMGLevels();
    const auto& mgOpts = config[iZone]->GetMGOptions();
    for (unsigned short i = 0; i <= nMGLevels; ++i) {
      lastPreSmoothIters[i] = 0;
      lastPostSmoothIters[i] = 0;
      lastCorrecSmoothIters[i] = mgOpts.MG_CorrecSmooth[i];
      lastCorrecSmoothRMS[i][0] = lastCorrecSmoothRMS[i][1] = 0.0;
      lastPreSmoothRMS[i][0] = lastPreSmoothRMS[i][1] = 0.0;
      lastPostSmoothRMS[i][0] = lastPostSmoothRMS[i][1] = 0.0;
      lastPreSmoothWorstStepRatio[i] = 0.0;
      lastPostSmoothWorstStepRatio[i] = 0.0;
      lastPreSmoothWorstStep[i] = 0;
      lastPostSmoothWorstStep[i] = 0;
      lastPreSmoothExitReason[i]  = ' ';
      lastPostSmoothExitReason[i] = ' ';
    }

    /*--- InnerIter restarts at every time step, the active level does not; re-anchor so the
     *    difference below cannot wrap around. ---*/

    if (config[iZone]->GetInnerIter() < mg_ramp_level_start_iter)
      mg_ramp_level_start_iter = config[iZone]->GetInnerIter();
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Promote to the next finer grid on MG_STARTUP_ITER, MG_STARTUP_CONVERGENCE or
   *    MG_STARTUP_STAGNATION, whichever comes first. ---*/
  const unsigned long startup_iter = config[iZone]->GetMGOptions().MG_Startup_Iter;
  const unsigned long iters_on_level = config[iZone]->GetInnerIter() - mg_ramp_level_start_iter;
  MGStartupPromote promote_reason = mg_startup_promote_reason;
  if ((promote_reason == MGStartupPromote::NONE) && (startup_iter > 0) && (iters_on_level >= startup_iter))
    promote_reason = MGStartupPromote::BUDGET;
  const bool Convergence_FullMG = FullMG && (FinestMesh != MESH_0) && (promote_reason != MGStartupPromote::NONE);

  if (Convergence_FullMG && direct && (RunTime_EqSystem == RUNTIME_FLOW_SYS)) {

    SetProlongated_Solution(RunTime_EqSystem,
                            solver_container[iZone][iInst][FinestMesh-1][Solver_Position],
                            solver_container[iZone][iInst][FinestMesh][Solver_Position],
                            geometry[iZone][iInst][FinestMesh-1],
                            geometry[iZone][iInst][FinestMesh],
                            config[iZone]);

    /*--- The scalar equations only restrict downward, so hand the new level their solution too,
     *    otherwise they restart from their initial condition at every promotion. ---*/

    for (const auto Scalar_Position : {TURB_SOL, TRANS_SOL, SPECIES_SOL, HEAT_SOL, RAD_SOL}) {

      CSolver* scalar_fine = solver_container[iZone][iInst][FinestMesh-1][Scalar_Position];
      CSolver* scalar_coarse = solver_container[iZone][iInst][FinestMesh][Scalar_Position];
      if ((scalar_fine == nullptr) || (scalar_coarse == nullptr)) continue;

      SetProlongated_ScalarSolution(scalar_fine, scalar_coarse,
                                    geometry[iZone][iInst][FinestMesh-1],
                                    geometry[iZone][iInst][FinestMesh],
                                    config[iZone], Scalar_Position == TURB_SOL);
    }

    /*--- Report the promotion before the startup state is reset below. ---*/
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    if (rank == MASTER_NODE) {
      cout << "Full-MG: mesh level " << FinestMesh << " -> " << FinestMesh - 1 << " after "
           << iters_on_level << " iteration(s) (";
      switch (promote_reason) {
        case MGStartupPromote::STAGNATION:
          cout << "residual stalled for " << config[iZone]->GetMGOptions().MG_Startup_Stagnation_Iter
               << " iteration(s)";
          break;
        case MGStartupPromote::CONVERGENCE:
          cout << "CONV_FIELD dropped "
               << fabs(SU2_TYPE::GetValue(config[iZone]->GetMGOptions().MG_Startup_Convergence))
               << " order(s) of magnitude";
          break;
        default:
          cout << "MG_STARTUP_ITER= " << startup_iter << " reached";
          break;
      }
      cout << ")." << endl;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    SU2_OMP_SAFE_GLOBAL_ACCESS(config[iZone]->SubtractFinestMesh();)
  }

  /*--- Set the current finest grid (full multigrid strategy) ---*/

  FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Rebuild the coarse-grid CFL before the cycle. ---*/
  const unsigned short nMGLevels = config[iZone]->GetnMGLevels();
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    /*--- Capture the configured damping factors, before adaptation has modified them. ---*/
    if (mg_damp_restric_initial < 0.0) {
      mg_damp_restric_initial = config[iZone]->GetDamp_Res_Restric();
      mg_damp_prolong_initial = config[iZone]->GetDamp_Correc_Prolong();
    }

    /*--- On a change of active level, restart the ramp window and the startup state. ---*/
    if (FinestMesh != mg_ramp_last_FinestMesh) {
      /*--- Only a promotion has a level below to take the handover CFL from. ---*/
      if (mg_ramp_last_FinestMesh == FinestMesh + 1)
        mg_ramp_cfl_start = SU2_TYPE::GetValue(config[iZone]->GetCFL(mg_ramp_last_FinestMesh));

      mg_ramp_level_start_iter = config[iZone]->GetInnerIter();
      mg_ramp_last_FinestMesh = FinestMesh;

      /*--- The EMA is measured on a different grid after a promotion. ---*/
      mg_fine_rms_ema = 0.0;
      config[iZone]->SetDamp_Res_Restric(mg_damp_restric_initial);
      config[iZone]->SetDamp_Correc_Prolong(mg_damp_prolong_initial);

      mg_startup_conv_start.clear();
      mg_startup_conv_prev.clear();
      mg_startup_stall_count = 0;
      mg_startup_promote_reason = MGStartupPromote::NONE;
    }

  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS


  /*--- While the startup ramps a level, its CFL has to be written before the cycle, not after
   *    it as in ordinary operation. ---*/

  const bool fmg_warmup = FullMG && (FinestMesh != MESH_0);
  const bool fmg_mesh0_ramp = FullMG_Mesh0RampWindow(config[iZone], FinestMesh, FullMG);

  if (fmg_warmup || fmg_mesh0_ramp)
    SetCoarseGridCFL(geometry, solver_container, config, RunTime_EqSystem, iZone, iInst, FinestMesh, FullMG,
                     fmg_mesh0_ramp);

  /*--- Perform the Full Approximation Scheme multigrid ---*/

  MultiGrid_Cycle(geometry, solver_container, numerics_container, config,
                  FinestMesh, RecursiveParam, RunTime_EqSystem, iZone, iInst);

  if (!fmg_warmup && !fmg_mesh0_ramp)
    SetCoarseGridCFL(geometry, solver_container, config, RunTime_EqSystem, iZone, iInst, FinestMesh, FullMG,
                     fmg_mesh0_ramp);

  /*--- Coarse-level residuals are per-rank partial sums; reduce once so the promotion criterion
   *    agrees on every rank. The smoothing early exit already reduces when it is on. ---*/

  if (fmg_warmup && !config[iZone]->GetMGOptions().MG_Smooth_EarlyExit) {
    solver_container[iZone][iInst][FinestMesh][Solver_Position]->SetResidual_RMS(
        geometry[iZone][iInst][FinestMesh], config[iZone], true);
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

  /*--- Adapt the damping factors, only meaningful when MG_SMOOTH_EARLY_EXIT= YES. ---*/
  const auto& mgOptsZone = config[iZone]->GetMGOptions();
  if (mgOptsZone.MG_Smooth_EarlyExit) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      constexpr passivedouble EMA_ALPHA = 0.02;
      const passivedouble fine_d0 = lastPreSmoothRMS[FinestMesh][0];
      const bool ema_ready = (mg_fine_rms_ema >= EPS);
      if (!ema_ready)
        mg_fine_rms_ema = fine_d0;
      else
        mg_fine_rms_ema = (1.0 - EMA_ALPHA) * mg_fine_rms_ema + EMA_ALPHA * fine_d0;
      const passivedouble crossCycleRatio = (mg_fine_rms_ema > EPS) ? fine_d0 / mg_fine_rms_ema : 1.0;

      /*--- Skip on the first cycle while the EMA is still being seeded. ---*/
      if (ema_ready) adaptDampingFactors(config[iZone], crossCycleRatio);
      last_crossCycleRatio = crossCycleRatio;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Print compact smoothing summary when MG_SMOOTH_OUTPUT= YES. ---*/
  if (mgOptsZone.MG_Smooth_Output) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    if (SU2_MPI::GetRank() == MASTER_NODE) {

      /*--- Helper: format one cell as "act/max [init->final]". ---*/
      auto cellStr = [](unsigned short act, unsigned short mx, char reason,
            su2double d0, su2double d1,
            passivedouble worstStepRatio,
            unsigned short worstStep) -> std::string {
        /*--- Steps taken/max, exit reason, initial defect and the d1/d0 ratio. ---*/
        std::ostringstream ss;
        ss << act << "/" << mx;
        if (act < mx) ss << reason;  /*--- only tag early exits ---*/
        ss << " [" << std::scientific << std::setprecision(2) << d0 << "]";
        if (d0 > 0.0) {
          ss << std::fixed << std::setprecision(3) << " r=" << d1 / d0;
        }
        if (worstStep > 0) {
          ss << " rw" << worstStep << "=" << std::fixed << std::setprecision(3) << worstStepRatio;
        }
        return ss.str();
      };

      PrintingToolbox::CTablePrinter table(&std::cout);
      table.AddColumn("Smoother", 13);
      for (unsigned short i = 0; i <= nMGLevels; ++i)
        table.AddColumn("Level " + std::to_string(i), 38);
      table.PrintHeader();

      /*--- Pre-smooth: defect [d0->d1] — what early exit and damping adaptation read. ---*/
      table << "Pre-smooth";
      for (unsigned short i = 0; i <= nMGLevels; ++i)
        table << cellStr(lastPreSmoothIters[i], mgOptsZone.MG_PreSmooth[i],
                          lastPreSmoothExitReason[i],
                          lastPreSmoothRMS[i][0], lastPreSmoothRMS[i][1],
                          lastPreSmoothWorstStepRatio[i], lastPreSmoothWorstStep[i]);

      /*--- Post-smooth: defect [d0->d1] — what early exit and prolongation-damping read. ---*/
      table << "Post-smooth";
      for (unsigned short i = 0; i < nMGLevels; ++i)
        table << cellStr(lastPostSmoothIters[i], mgOptsZone.MG_PostSmooth[i],
                          lastPostSmoothExitReason[i],
                          lastPostSmoothRMS[i][0], lastPostSmoothRMS[i][1],
                          lastPostSmoothWorstStepRatio[i], lastPostSmoothWorstStep[i]);
      table << "-";

      /*--- Corr.-smooth: defined on levels 0..nMGLevels-1; coarsest has none. ---*/
      table << "Corr-smooth";
      for (unsigned short i = 0; i < nMGLevels; ++i)
        table << cellStr(lastCorrecSmoothIters[i], mgOptsZone.MG_CorrecSmooth[i],
                          ' ', lastCorrecSmoothRMS[i][0], lastCorrecSmoothRMS[i][1],
                          0.0, 0);
      table << "-";

      table << "CFL";
      for (unsigned short i = 0; i <= nMGLevels; ++i) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4);
        if (i == MESH_0) {
          ss << SU2_TYPE::GetValue(solver_container[iZone][iInst][MESH_0][Solver_Position]->GetAvg_CFL_Local());
        } else {
          ss << SU2_TYPE::GetValue(config[iZone]->GetCFL(i));
        }
        table << ss.str();
      }

      table.PrintFooter();

      cout << std::fixed << std::setprecision(4)
            << "Damping [restrict | prolong] : " << config[iZone]->GetDamp_Res_Restric()
            << " | " << config[iZone]->GetDamp_Correc_Prolong()
            << "  cross-cycle: " << std::setprecision(3) << last_crossCycleRatio << "\n"
            << std::defaultfloat << std::setprecision(6);
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
  SU2_ZONE_SCOPED

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

    /*--- LinSysRes = R(u_N) here, before tau is added by SetResidual_Term. ---*/
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      lastPreSmoothRMS[iMesh][1] = ComputeLinSysResRMS(solver_fine);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

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

    /*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +
          Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/

    GetProlongated_Correction(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);

    const auto& mgOpts = config->GetMGOptions();
    SmoothProlongated_Correction(RunTime_EqSystem, solver_fine, geometry_fine, mgOpts.MG_CorrecSmooth[iMesh], mgOpts.MG_Smooth_Coeff, config, iMesh);

    SetProlongated_Correction(solver_fine, geometry_fine, config, iMesh);

    /*--- Solution post-smoothing in the prolongated grid. ---*/

    PostSmoothing(RunTime_EqSystem, solver_fine, numerics_fine, geometry_fine, solver_container_fine,
                  config, iMesh, iRKLimit);
  }

}

void CMultiGridIntegration::prePostEarlyExit(unsigned short iSmooth, unsigned short iMesh,
                                             passivedouble defect, const CMGOptions& mgOpts,
                                             passivedouble stag_tol, bool early_exit,
                                             passivedouble lastRMS[2], char& exitReason,
                                             passivedouble& worstStepRatio, unsigned short& worstStep) {
  if (iSmooth == 0) {
    mg_initial_smooth_rms = defect;
    lastRMS[0] = defect;
  } else {
    if (mg_prev_smooth_rms > EPS) {
      const passivedouble step_ratio = defect / mg_prev_smooth_rms;
      if (step_ratio > worstStepRatio) {
        worstStepRatio = step_ratio;
        worstStep = iSmooth + 1;
      }
    }

    if (early_exit) {
      if (defect < mgOpts.MG_Smooth_Res_Threshold * mg_initial_smooth_rms) {
        exitReason = 'T';
        mg_early_exit_flag = true;
      } else if (defect >= mg_prev_smooth_rms * stag_tol) {
        /*--- 'A' = amplifying-stagnation (defect grew vs previous step).
         *    'S' = clean stagnation (defect is not improving but also not growing). ---*/
        exitReason = (defect > mg_prev_smooth_rms) ? 'A' : 'S';
        mg_early_exit_flag = true;
      }
    }
  }
  mg_prev_smooth_rms = defect;
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
  SU2_ZONE_SCOPED

  const auto& mgOpts = config->GetMGOptions();
  const bool classical_rk4 = (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT);
  const unsigned short nPreSmooth = mgOpts.MG_PreSmooth[iMesh];
  const unsigned long timeIter = config->GetTimeIter();
  const bool early_exit = mgOpts.MG_Smooth_EarlyExit && (nPreSmooth > 1);
  const bool need_per_step_rms = early_exit || mgOpts.MG_Smooth_Output;
  /*--- Also capture initial RMS at MESH_0 for the cross-cycle EMA controller, even with nPreSmooth==1. ---*/
  const bool need_initial_rms = need_per_step_rms || (iMesh == MESH_0 && mgOpts.MG_Smooth_EarlyExit);
  const passivedouble stag_tol = (mgOpts.MG_Smooth_StagnationTol > 0.0)
                                 ? SU2_TYPE::GetValue(mgOpts.MG_Smooth_StagnationTol) : passivedouble(1.0);

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

      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);

      /*--- At iRKStep==0 LinSysRes = R(u_k). ---*/
      if (iRKStep == 0 && need_initial_rms) {
        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
          const passivedouble defect = ComputeLinSysResRMS(solver_fine);
          prePostEarlyExit(iPreSmooth, iMesh, defect, mgOpts, stag_tol, early_exit,
                           lastPreSmoothRMS[iMesh], lastPreSmoothExitReason[iMesh],
                           lastPreSmoothWorstStepRatio[iMesh], lastPreSmoothWorstStep[iMesh]);
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS
        if (mg_early_exit_flag) break;
      }


      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);
    }
    SU2_OMP_SAFE_GLOBAL_ACCESS(lastPreSmoothIters[iMesh] = iPreSmooth + 1;)
    if (mg_early_exit_flag) break;
  }

  /*--- Record the final pre-smooth defect; MultiGrid_Cycle overwrites it with the exact
   *    d_N on non-coarsest levels. ---*/
  if (nPreSmooth > 0) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      lastPreSmoothRMS[iMesh][1] = mg_prev_smooth_rms;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }
}


void CMultiGridIntegration::PostSmoothing(unsigned short RunTime_EqSystem,
                                          CSolver* solver_fine,
                                          CNumerics** numerics_fine,
                                          CGeometry* geometry_fine,
                                          CSolver** solver_container_fine,
                                          CConfig *config,
                                          unsigned short iMesh,
                                          unsigned short iRKLimit) {
  SU2_ZONE_SCOPED

  const auto& mgOpts = config->GetMGOptions();
  const bool classical_rk4 = (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT);
  const unsigned short nPostSmooth = mgOpts.MG_PostSmooth[iMesh];
  const unsigned long timeIter = config->GetTimeIter();
  const bool early_exit = mgOpts.MG_Smooth_EarlyExit && (nPostSmooth > 1);
  const bool need_per_step_rms = early_exit || mgOpts.MG_Smooth_Output;
  const passivedouble stag_tol = (mgOpts.MG_Smooth_StagnationTol > 0.0)
                                 ? SU2_TYPE::GetValue(mgOpts.MG_Smooth_StagnationTol) : passivedouble(1.0);

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

      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      Time_Integration(geometry_fine, solver_container_fine, config, iRKStep, RunTime_EqSystem);

      /*--- At iRKStep==0 LinSysRes = R(u_k) ---*/
      if (iRKStep == 0 && need_per_step_rms) {
        BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
          const passivedouble defect = ComputeLinSysResRMS(solver_fine);
          prePostEarlyExit(iPostSmooth, iMesh, defect, mgOpts, stag_tol, early_exit,
                           lastPostSmoothRMS[iMesh], lastPostSmoothExitReason[iMesh],
                           lastPostSmoothWorstStepRatio[iMesh], lastPostSmoothWorstStep[iMesh]);
        }
        END_SU2_OMP_SAFE_GLOBAL_ACCESS
        if (mg_early_exit_flag) break;
      }


      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      solver_fine->Postprocessing(geometry_fine, solver_container_fine, config, iMesh);

    }
    SU2_OMP_SAFE_GLOBAL_ACCESS(lastPostSmoothIters[iMesh] = iPostSmooth + 1;)
    if (mg_early_exit_flag) break;
  }

  /*--- Record the final post-smooth defect, for display only. ---*/
  if (nPostSmooth > 0) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      lastPostSmoothRMS[iMesh][1] = mg_prev_smooth_rms;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }
}


void CMultiGridIntegration::GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                      CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  SU2_ZONE_SCOPED

  const unsigned short nVar = sol_coarse->GetnVar();
  su2activevector Solution(nVar);

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    su2double Area_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);

    Solution = su2double(0);

    /*--- Accumulate children contributions with stable ordering ---*/
    /*--- Process all children in sequential order to ensure deterministic FP summation ---*/
    auto nChildren = geo_coarse->nodes->GetnChildren_CV(Point_Coarse);
    for (auto iChildren = 0u; iChildren < nChildren; iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      su2double Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
      const auto* Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      su2double weight = Area_Children / Area_Parent;
      for (auto iVar = 0u; iVar < nVar; iVar++)
        Solution[iVar] -= Solution_Fine[iVar] * weight;
    }

    const auto* Solution_Coarse = sol_coarse->GetNodes()->GetSolution(Point_Coarse);

    for (auto iVar = 0u; iVar < nVar; iVar++)
      Solution[iVar] += Solution_Coarse[iVar];

    for (auto iVar = 0u; iVar < nVar; iVar++)
      sol_coarse->GetNodes()->SetSolution_Old(Point_Coarse, Solution.data());
  }
  END_SU2_OMP_FOR

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
  SU2_ZONE_SCOPED

  /*--- Check if there is work to do. ---*/
  if (val_nSmooth == 0) return;

  const unsigned short nVar = solver->GetnVar();

  SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); iPoint++) {
    const auto* Residual_Old = solver->LinSysRes.GetBlock(iPoint);
    solver->GetNodes()->SetResidual_Old(iPoint, Residual_Old);
  }
  END_SU2_OMP_FOR

  /*--- Jacobi iterations (no early exit — Jacobi targets high-frequency modes,
   *    so the global RMS is not a meaningful convergence indicator). ---*/

  for (auto iSmooth = 0u; iSmooth < val_nSmooth; iSmooth++) {

    /*--- Loop over all mesh points (sum the residuals of direct neighbors). ---*/

    SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {

      solver->GetNodes()->SetResidualSumZero(iPoint);

      for (auto iNeigh = 0u; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
        auto jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
        const auto* Residual_j = solver->LinSysRes.GetBlock(jPoint);
        solver->GetNodes()->AddResidual_Sum(iPoint, Residual_j);
      }

    }
    END_SU2_OMP_FOR

    /*--- Loop over all mesh points (update residuals with the neighbor averages). ---*/

    SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {

      su2double factor = 1.0/(1.0+val_smooth_coeff*su2double(geometry->nodes->GetnPoint(iPoint)));

      const auto* Residual_Sum = solver->GetNodes()->GetResidual_Sum(iPoint);
      const auto* Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);

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
          const auto* Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);
          solver->LinSysRes.SetBlock(iPoint, Residual_Old);
        }
        END_SU2_OMP_FOR
      }
    }

  }

  /*--- Record final correction norm for debugging output. ---*/
  if (config->GetMGOptions().MG_Smooth_Output) {
    const su2double res = sqrt(solver->LinSysRes.squaredNorm() / (nVar * geometry->GetGlobal_nPointDomain()));
    SU2_OMP_SAFE_GLOBAL_ACCESS(lastCorrecSmoothRMS[iMesh][1] = SU2_TYPE::GetValue(res);)
  }
}

void CMultiGridIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine,
                                                      CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const unsigned short nVar = sol_fine->GetnVar();

  /*--- Use the adaptive damping factor uniformly across all prolongation levels. ---*/
  const su2double factor = config->GetDamp_Correc_Prolong();

  SU2_OMP_FOR_STAT(roundUpDiv(geo_fine->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Fine = 0ul; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
    auto* Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
    auto* Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
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
  SU2_ZONE_SCOPED

  const unsigned short Solver_Position = config->GetContainerPosition(RunTime_EqSystem);
  const bool grid_movement = config->GetGrid_Movement();

  /*--- Constant injection: every fine child of a coarse CV takes its parent's value. ---*/

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    const auto* Solution_Coarse = sol_coarse->GetNodes()->GetSolution(Point_Coarse);
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      const auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      sol_fine->GetNodes()->SetSolution(Point_Fine, Solution_Coarse);
    }
  }
  END_SU2_OMP_FOR

  /*--- The injected values do not satisfy the fine-grid wall conditions on their own. ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {

      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_fine->nVertex[iMarker]; iVertex++) {

        const auto Point_Fine = geo_fine->vertex[iMarker][iVertex]->GetNode();

        if (Solver_Position == FLOW_SOL) {

          /*--- At moving walls, set the solution based on the new density and wall velocity ---*/

          if (grid_movement) {
            const auto* Grid_Vel = geo_fine->nodes->GetGridVel(Point_Fine);
            sol_fine->GetNodes()->SetVelSolutionVector(Point_Fine, Grid_Vel);
          }
          else {
            /*--- For stationary no-slip walls, set the velocity to zero. ---*/
            su2double zero[3] = {0.0};
            sol_fine->GetNodes()->SetVelSolutionVector(Point_Fine, zero);
          }

        }

        if (Solver_Position == ADJFLOW_SOL) {
          sol_fine->GetNodes()->SetVelSolutionDVector(Point_Fine);
        }

      }
      END_SU2_OMP_FOR
    }
  }

  /*--- Project the velocity onto the fine-grid wall tangent plane. ---*/

  sol_fine->MultigridProjectEulerWall(geo_fine, config, false);

  /*--- MPI the new interpolated solution. ---*/

  sol_fine->InitiateComms(geo_fine, config, MPI_QUANTITIES::SOLUTION);
  sol_fine->CompleteComms(geo_fine, config, MPI_QUANTITIES::SOLUTION);

}

void CMultiGridIntegration::SetProlongated_ScalarSolution(CSolver *sol_fine, CSolver *sol_coarse,
                                                          CGeometry *geo_fine, CGeometry *geo_coarse,
                                                          CConfig *config, bool eddy_viscosity) {
  SU2_ZONE_SCOPED

  /*--- Constant injection: every fine child of a coarse CV takes its parent's value. ---*/

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    const auto* Solution_Coarse = sol_coarse->GetNodes()->GetSolution(Point_Coarse);
    const su2double muT_Coarse = eddy_viscosity ? sol_coarse->GetNodes()->GetmuT(Point_Coarse) : su2double(0.0);

    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      const auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      sol_fine->GetNodes()->SetSolution(Point_Fine, Solution_Coarse);
      if (eddy_viscosity) sol_fine->GetNodes()->SetmuT(Point_Fine, muT_Coarse);
    }
  }
  END_SU2_OMP_FOR

  /*--- The eddy viscosity is zero on a no-slip wall, the injected value is not. ---*/

  if (eddy_viscosity) {
    for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
      if (!config->GetViscous_Wall(iMarker)) continue;

      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_fine->nVertex[iMarker]; iVertex++) {
        const auto Point_Fine = geo_fine->vertex[iMarker][iVertex]->GetNode();
        sol_fine->GetNodes()->SetmuT(Point_Fine, 0.0);
      }
      END_SU2_OMP_FOR
    }
  }

  /*--- MPI the new interpolated solution, and the eddy viscosity with it. ---*/

  const auto commType = eddy_viscosity ? MPI_QUANTITIES::SOLUTION_EDDY : MPI_QUANTITIES::SOLUTION;
  sol_fine->InitiateComms(geo_fine, config, commType);
  sol_fine->CompleteComms(geo_fine, config, commType);

}

void CMultiGridIntegration::SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                            CGeometry *geo_coarse, CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const su2double *Residual_Fine;

  const unsigned short nVar = sol_coarse->GetnVar();
  const su2double factor = config->GetDamp_Res_Restric();

  su2activevector Residual(nVar);

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);

    Residual = su2double(0);
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (auto iVar = 0u; iVar < nVar; iVar++)
        Residual[iVar] += factor * Residual_Fine[iVar];
    }
    sol_coarse->GetNodes()->AddRes_TruncError(Point_Coarse, Residual.data());
  }
  END_SU2_OMP_FOR

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
  SU2_ZONE_SCOPED

  AD::StartNoSharedReading();
  SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPointDomain(), omp_get_num_threads()))
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); iPoint++)
    solver->LinSysRes.AddBlock(iPoint, solver->GetNodes()->GetResTruncError(iPoint));
  END_SU2_OMP_FOR
  AD::EndNoSharedReading();

}

void CMultiGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                   CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  SU2_ZONE_SCOPED

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
  SU2_ZONE_SCOPED

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
  SU2_ZONE_SCOPED

  if (RunTime_EqSystem == RUNTIME_FLOW_SYS) {
    /*--- Calculate the inviscid and viscous forces ---*/

    solver_container[FinestMesh][FLOW_SOL]->Pressure_Forces(geometry[FinestMesh], config);
    solver_container[FinestMesh][FLOW_SOL]->Momentum_Forces(geometry[FinestMesh], config);
    solver_container[FinestMesh][FLOW_SOL]->Friction_Forces(geometry[FinestMesh], config);
  }

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  if (RunTime_EqSystem == RUNTIME_ADJFLOW_SYS) {
    /*--- Calculate the inviscid and viscous sensitivities ---*/

    solver_container[FinestMesh][ADJFLOW_SOL]->Inviscid_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh],
                                                numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);

    solver_container[FinestMesh][ADJFLOW_SOL]->Viscous_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh],
                                                numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);

    /*--- Smooth the inviscid and viscous sensitivities ---*/

    if (config->GetKind_SensSmooth() != NONE)
      solver_container[FinestMesh][ADJFLOW_SOL]->Smooth_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh],
                                                  numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CMultiGridIntegration::Adjoint_Setup(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                                          unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone) {
  SU2_ZONE_SCOPED

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
