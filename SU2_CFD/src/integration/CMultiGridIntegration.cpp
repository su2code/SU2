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
#include "../../include/variables/CTurbVariable.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

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

static su2double GetMGLevelCorrectionScale(unsigned short iMesh) {
  switch (iMesh) {
    case 0: return 1.00;
    case 1: return 0.75;
    case 2: return 0.50;
    default: return 0.35;
  }
}

inline passivedouble ComputeLinSysResRMS(const CSolver* solver) {
  passivedouble result = 0;
  for (unsigned short iVar = 0; iVar < solver->GetnVar(); ++iVar) {
    result += pow(SU2_TYPE::GetValue(solver->GetRes_RMS(iVar)), 2);
  }
  return sqrt(result);
}

static void ApplyLineImplicitResidualSmoothing(CSolver* solver, CGeometry* geometry, unsigned short iMesh) {
  if (iMesh == 0) return;

  const auto nPoint = geometry->GetnPointDomain();
  if (nPoint < 3) return;

  const unsigned short nVar = solver->GetnVar();
  const unsigned short nDim = geometry->GetnDim();
  const su2double damping = 0.25;
  std::vector<bool> visited(nPoint, false);
  unsigned long nLines = 0;
  unsigned long totalLineLength = 0;
  su2double totalLineResidual = 0.0;

  for (auto iSeed = 0ul; iSeed < nPoint; ++iSeed) {
    if (visited[iSeed]) continue;

    std::vector<unsigned long> line;
    line.reserve(16);
    line.push_back(iSeed);
    visited[iSeed] = true;

    unsigned long current = iSeed;
    for (int step = 0; step < 8; ++step) {
      unsigned long next = std::numeric_limits<unsigned long>::max();
      const auto* coordCurrent = geometry->nodes->GetCoord(current);
      const unsigned short nNeigh = geometry->nodes->GetnPoint(current);
      su2double bestScore = -1e30;

      for (unsigned short iNeigh = 0; iNeigh < nNeigh; ++iNeigh) {
        const auto candidate = geometry->nodes->GetPoint(current, iNeigh);
        if (candidate >= nPoint || candidate == current || visited[candidate]) continue;

        const auto* coordCandidate = geometry->nodes->GetCoord(candidate);
        su2double score = 0.0;
        for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
          const su2double delta = fabs(coordCandidate[iDim] - coordCurrent[iDim]);
          score += delta;
        }
        if (score > bestScore) {
          bestScore = score;
          next = candidate;
        }
      }

      if (next == std::numeric_limits<unsigned long>::max()) break;
      line.push_back(next);
      visited[next] = true;
      current = next;
    }

    const auto nLine = static_cast<unsigned long>(line.size());
    if (nLine < 2) continue;

    ++nLines;
    totalLineLength += nLine;

    std::vector<su2double> lineAverage(nVar, 0.0);
    for (auto i = 0ul; i < nLine; ++i) {
      const auto* residual = solver->GetNodes()->GetResidual_Old(line[i]);
      if (residual == nullptr) continue;
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
        lineAverage[iVar] += residual[iVar];
        totalLineResidual += fabs(residual[iVar]);
      }
    }

    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      lineAverage[iVar] /= static_cast<su2double>(nLine);
    }

    for (auto i = 0ul; i < nLine; ++i) {
      const auto* oldResidual = solver->GetNodes()->GetResidual_Old(line[i]);
      std::vector<su2double> block(nVar, 0.0);
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
        block[iVar] = oldResidual[iVar] + damping * (lineAverage[iVar] - oldResidual[iVar]);
      }
      solver->LinSysRes.SetBlock(line[i], block.data());
    }
  }

  if (SU2_MPI::GetRank() == MASTER_NODE) {
    const su2double avgLineLength = (nLines > 0) ? static_cast<su2double>(totalLineLength) / static_cast<su2double>(nLines) : 0.0;
    const su2double avgLineResidual = (nLines > 0) ? totalLineResidual / static_cast<su2double>(nLines) : 0.0;
    cout << "[MG LINE] level=" << iMesh << " activated=" << (nLines > 0 ? "yes" : "no")
         << " lines=" << nLines << " avgLen=" << avgLineLength
         << " avgResidual=" << avgLineResidual << std::endl;
  }
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
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Full MG: advance to the next finer grid after a fixed number of
   *    outer iterations on the current coarsest active level.
   *    The number of iterations per level is controlled by MG_STARTUP_ITER config option. ---*/
  const unsigned long startup_iter = config[iZone]->GetMGOptions().MG_Startup_Iter;
  const bool Convergence_FullMG =
      FullMG && (FinestMesh != MESH_0) &&
      (config[iZone]->GetInnerIter() % startup_iter == startup_iter - 1);

  if (!config[iZone]->GetRestart() && FullMG && direct && ( Convergence_FullMG && (FinestMesh != MESH_0 )) &&
      RunTime_EqSystem == RUNTIME_FLOW_SYS) {

    SetProlongated_Solution(RunTime_EqSystem,
                            solver_container[iZone][iInst][FinestMesh-1][Solver_Position],
                            solver_container[iZone][iInst][FinestMesh][Solver_Position],
                            geometry[iZone][iInst][FinestMesh-1],
                            geometry[iZone][iInst][FinestMesh],
                            config[iZone]);

    /*--- Prolongate scalar solutions to the new finest mesh.
     *    All scalar solvers (turb, species, transition) run via SingleGrid_Iteration on
     *    GetFinestMesh().  Only turbulence additionally restricts its field downward to
     *    coarser meshes; no scalar ever propagates upward to finer meshes.  Consequently
     *    meshes finer than FinestMesh hold their iter-0 startup values for the entire
     *    warmup phase.  When FinestMesh is decremented these stale fields cause a large
     *    transient (e.g. +3 decade regression in rms[nu]).  Prolongating here mirrors
     *    what SetProlongated_Solution does for the flow and eliminates the regression. ---*/
    if (config[iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE) {
      SetProlongated_Solution(RUNTIME_TURB_SYS,
                              solver_container[iZone][iInst][FinestMesh-1][TURB_SOL],
                              solver_container[iZone][iInst][FinestMesh][TURB_SOL],
                              geometry[iZone][iInst][FinestMesh-1],
                              geometry[iZone][iInst][FinestMesh],
                              config[iZone]);
      /*--- Recompute mu_t on the new finest mesh from the prolongated nu_tilde/k/omega. ---*/
      solver_container[iZone][iInst][FinestMesh-1][TURB_SOL]->Postprocessing(
          geometry[iZone][iInst][FinestMesh-1],
          solver_container[iZone][iInst][FinestMesh-1],
          config[iZone], FinestMesh-1);
    }

    if (config[iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      SetProlongated_Solution(RUNTIME_TRANS_SYS,
                              solver_container[iZone][iInst][FinestMesh-1][TRANS_SOL],
                              solver_container[iZone][iInst][FinestMesh][TRANS_SOL],
                              geometry[iZone][iInst][FinestMesh-1],
                              geometry[iZone][iInst][FinestMesh],
                              config[iZone]);
    }

    if (config[iZone]->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
      SetProlongated_Solution(RUNTIME_SPECIES_SYS,
                              solver_container[iZone][iInst][FinestMesh-1][SPECIES_SOL],
                              solver_container[iZone][iInst][FinestMesh][SPECIES_SOL],
                              geometry[iZone][iInst][FinestMesh-1],
                              geometry[iZone][iInst][FinestMesh],
                              config[iZone]);
    }

    SU2_OMP_SAFE_GLOBAL_ACCESS(config[iZone]->SubtractFinestMesh();)
  }

  /*--- Set the current finest grid (full multigrid strategy) ---*/

  FinestMesh = config[iZone]->GetFinestMesh();

  /*--- For turbulence MG: before descending to coarse levels, ensure mu_t is computed
   *    at the finest level and restricted to all coarser levels. This prevents inf
   *    residuals from coarse-level turbulence solves using stale/uninitialized mu_t. ---*/
  if (RunTime_EqSystem == RUNTIME_TURB_SYS &&
      config[iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE) {

    solver_container[iZone][iInst][FinestMesh][TURB_SOL]->Postprocessing(
        geometry[iZone][iInst][FinestMesh],
        solver_container[iZone][iInst][FinestMesh],
        config[iZone], FinestMesh);

    RestrictTurbEddyViscToCoarseLevels(geometry[iZone][iInst],
                                       solver_container[iZone][iInst],
                                       config[iZone], FinestMesh,
                                       config[iZone]->GetnMGLevels());
  }

  /*--- Rebuild coarse-grid CFL before the cycle so the currently active FMG
   *    level uses the intended CFL in this iteration. ---*/
  const unsigned short nMGLevels = config[iZone]->GetnMGLevels();
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    /*--- Use the level-0 flow CFL as the base reference and derive all coarse
     *    levels from it via MG_CFL_SCALING[i] = CFL(i+1)/CFL(i). Fall back to
     *    config scalar when local level-0 CFL is unavailable. ---*/
    passivedouble cfl_base = SU2_TYPE::GetValue(
      solver_container[iZone][iInst][MESH_0][Solver_Position]->GetAvg_CFL_Local());
    if (cfl_base < EPS)
      cfl_base = SU2_TYPE::GetValue(config[iZone]->GetCFL(MESH_0));

    const auto& cflScaling = config[iZone]->GetMGOptions().MG_CflScaling;

    passivedouble CFL_local = cfl_base;
    for (unsigned short lvl = 1; lvl <= nMGLevels; ++lvl) {
      /*--- Index into cflScaling is (lvl-1): transition lvl-1 -> lvl. ---*/
      const unsigned short iScale = lvl - 1;
      const passivedouble scale = (iScale < cflScaling.size())
          ? max(passivedouble{1e-6}, SU2_TYPE::GetValue(cflScaling[iScale]))
          : passivedouble{0.25};
      CFL_local *= scale;
      config[iZone]->SetCFL(lvl, CFL_local);
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Propagate updated CFL to all coarse-grid points before the cycle. ---*/
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

  /*--- Perform the Full Approximation Scheme multigrid ---*/

  MultiGrid_Cycle(geometry, solver_container, numerics_container, config,
                  FinestMesh, RecursiveParam, RunTime_EqSystem, iZone, iInst);

  /*--- After a turb FAS V-cycle: recompute mu_t at the finest active level from the updated
   *    nu_hat/k/omega and restrict it to all coarser levels.  The flow FAS on the NEXT outer
   *    iteration uses these mu_t values at every coarse level for the eddy-viscosity coupling.
   *    (Postprocessing was already called on FinestMesh inside the last PreSmoothing step of
   *    MultiGrid_Cycle; we call it once more to be safe after the V-cycle correction is applied.) ---*/
  if (RunTime_EqSystem == RUNTIME_TURB_SYS &&
      config[iZone]->GetKind_Turb_Model() != TURB_MODEL::NONE) {
    solver_container[iZone][iInst][FinestMesh][TURB_SOL]->Postprocessing(
        geometry[iZone][iInst][FinestMesh],
        solver_container[iZone][iInst][FinestMesh],
        config[iZone], FinestMesh);
    RestrictTurbEddyViscToCoarseLevels(geometry[iZone][iInst],
                                       solver_container[iZone][iInst],
                                       config[iZone], FinestMesh,
                                       config[iZone]->GetnMGLevels());
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
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      constexpr passivedouble EMA_ALPHA = 0.02;
      const passivedouble fine_d0 = lastPreSmoothRMS[FinestMesh][0];
      const bool ema_ready = (mg_fine_rms_ema >= EPS);
      if (!ema_ready)
        mg_fine_rms_ema = fine_d0;
      else
        mg_fine_rms_ema = (1.0 - EMA_ALPHA) * mg_fine_rms_ema + EMA_ALPHA * fine_d0;
      const passivedouble crossCycleRatio = (mg_fine_rms_ema > EPS) ? fine_d0 / mg_fine_rms_ema : 1.0;

      /*--- Adapt both damping factors from the same global-trend signal.
       *    Skip on the first cycle while the EMA is still being seeded. ---*/
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
        /*--- Show: steps taken / max + exit reason + initial defect scale + d1/d0 ratio.
         *    r=d1/d0 < 1 means smoother reduced the defect (good).
         *    r > 1 means smoother grew the defect.
         *    Exit reason: T=threshold, S=clean stagnation, A=amplifying stagnation,
         *                 ' '=ran to completion. ---*/
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

      const string eqName = (RunTime_EqSystem == RUNTIME_FLOW_SYS)    ? "Flow"    :
                             (RunTime_EqSystem == RUNTIME_TURB_SYS)    ? "Turb"    :
                             (RunTime_EqSystem == RUNTIME_SPECIES_SYS) ? "Species" :
                             (RunTime_EqSystem == RUNTIME_TRANS_SYS)   ? "Trans"   : "Other";

      PrintingToolbox::CTablePrinter table(&std::cout);
      table.AddColumn("Smoother [" + eqName + "]", 13 + 7);
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
        if (i == FinestMesh) {
          ss << SU2_TYPE::GetValue(solver_container[iZone][iInst][FinestMesh][Solver_Position]->GetAvg_CFL_Local());
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

  /*--- Assemble the coarse-grid FAS defect term by restricting the fine-grid residual defect,
   *    solving the coarse-grid state, and prolongating only the state correction back to the fine grid. ---*/

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

    /*--- LinSysRes = R(u_N) here, before the fine-grid defect term is assembled. ---*/
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      lastPreSmoothRMS[iMesh][1] = ComputeLinSysResRMS(solver_fine);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    /*--- Assemble the fine-grid defect term that will be restricted to the coarse-grid FAS problem. ---*/
    SetResidual_Term(geometry_fine, solver_fine);

    /*--- Communicate frozen sources on finest level to halo cells before restriction aggregates from fine children. ---*/
    if (iMesh == MESH_0 && RunTime_EqSystem == RUNTIME_TURB_SYS && config->GetMGOptions().MG_Turb_Freeze_Source) {
      solver_fine->InitiateComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION_EDDY);
      solver_fine->CompleteComms(geometry_fine, config, MPI_QUANTITIES::SOLUTION_EDDY);
    }

    /*--- Restrict the fine-grid state to the coarse grid and initialize the coarse-grid state. ---*/

    SetRestricted_Solution(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);

    /*--- Restrict frozen source terms for turbulence multigrid if feature is enabled. ---*/
    if (RunTime_EqSystem == RUNTIME_TURB_SYS && config->GetMGOptions().MG_Turb_Freeze_Source) {
      SetRestricted_FrozenSource(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);
    }

    solver_coarse->Preprocessing(geometry_coarse, solver_container_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem, false);

    /*--- For turbulence: ensure flow primitives (density, laminar viscosity) are updated on the
     *    coarse level from the restricted conservative variables, THEN compute mu_t from the
     *    newly restricted turbulence variables. This ensures turbulence Postprocessing reads
     *    valid flow data and Space_Integration uses correct eddy viscosity. ---*/
    if (RunTime_EqSystem == RUNTIME_TURB_SYS && config->GetKind_Turb_Model() != TURB_MODEL::NONE) {
      solver_container_coarse[FLOW_SOL]->Preprocessing(geometry_coarse, solver_container_coarse, config, iMesh+1, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      solver_coarse->Postprocessing(geometry_coarse, solver_container_coarse, config, iMesh+1);
    }

    Space_Integration(geometry_coarse, solver_container_coarse, numerics_coarse, config, iMesh+1, NO_RK_ITER, RunTime_EqSystem);

    /*--- Restrict the fine-grid residual defect to the coarse-grid FAS forcing term. ---*/
    if (RunTime_EqSystem == RUNTIME_FLOW_SYS) {
      RestrictResidualToCoarseGrid(solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh+1);
    } else {
      SetForcing_Term(solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh+1);
    }

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

    /*--- Compute the coarse-grid state correction and prolongate it back to the fine grid. ---*/
    if (RunTime_EqSystem == RUNTIME_FLOW_SYS) {
      ProlongateCorrectionToFineGrid(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config, iMesh);
    } else {
      GetProlongated_Correction(RunTime_EqSystem, solver_fine, solver_coarse, geometry_fine, geometry_coarse, config);
    }

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

  /*--- Record d_{N-1} as the final pre-smooth defect (the last value captured inside the loop).
   *    For non-coarsest levels MultiGrid_Cycle overwrites this with the exact d_N at zero
   *    additional cost in the restriction block (Space_Integration already runs there).
   *    Skip when nPreSmooth==0: lastPreSmoothDefect[iMesh] stays {0,0} (initialized). ---*/
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

  /*--- Record d_{N-1} as the final post-smooth defect (display only).
   *    Skip when nPostSmooth==0: lastPostSmoothRMS[iMesh] stays {0,0} (initialized). ---*/
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

  /*--- Remove any contributions from no-slip walls (Dirichlet BC enforcement). ---*/

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {

      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        auto Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();

        /*--- For Dirichlet boundary conditions, set the correction to zero.
         *    Note that Solution_Old stores the correction, not the actual value. ---*/

        if (RunTime_EqSystem == RUNTIME_TURB_SYS) {
          /*--- Turbulence: explicitly zero all variables (scalar equations). ---*/
          for (auto iVar = 0u; iVar < sol_coarse->GetnVar(); iVar++) {
            sol_coarse->GetNodes()->SetSolution_Old(Point_Coarse, iVar, 0.0);
          }
        } else {
          /*--- Flow: zero velocity components only (momentum equations). ---*/
          su2double zero[3] = {0.0};
          sol_coarse->GetNodes()->SetVelocity_Old(Point_Coarse, zero);
        }

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
  const bool use_conservative_damping = (nVar <= 2);
  const su2double turbulence_base_damping = 0.50;

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

      for (auto iVar = 0u; iVar < nVar; iVar++) {
        su2double smoothed = (Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar])*factor;
        /*--- FIX 2: Do NOT apply turbulence_base_damping (0.50) to turbulence corrections.
         *    use_conservative_damping is true for nVar<=2 which includes SA (nVar=1), but the
         *    0.50 factor was designed for compressible flow density-velocity-pressure coupling
         *    near walls. Applying it to a scalar equation halves all corrections unnecessarily.
         *    Combined with SetProlongated_Correction's damping (~0.375 at Level 1) and the config
         *    damping factor (0.9), only 17% of the coarse correction was reaching the fine grid. ---*/
        if (use_conservative_damping && RunTime_EqSystem == RUNTIME_FLOW_SYS) smoothed *= turbulence_base_damping;
        solver->LinSysRes(iPoint,iVar) = smoothed;
      }
    }
    END_SU2_OMP_FOR

    if (iMesh > 0 && RunTime_EqSystem == RUNTIME_FLOW_SYS) {
      ApplyLineImplicitResidualSmoothing(solver, geometry, iMesh);
    }

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

    if (SU2_MPI::GetRank() == MASTER_NODE && use_conservative_damping) {
      cout << "[MG CORR-SMOOTH] turbulence nSmooth=" << val_nSmooth
           << " norm=" << res << "\n";
    }
  }
}

void CMultiGridIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine,
                                                      CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  const unsigned short nVar = sol_fine->GetnVar();
  /*--- Conservative damping prevents oscillations from density-velocity-pressure coupling
   *    in compressible flow (nVar > 2 equations). Scalar/turbulence equations (nVar <= 2)
   *    are single-variable and have no such coupling; applying flow-level damping (base ~0.375)
   *    reduces turbulence corrections to only 34% of their computed value, preventing the
   *    coarse-grid MG cycle from being effective. ---*/
  const bool use_conservative_damping = (nVar > 2);
  const su2double levelScale = GetMGLevelCorrectionScale(iMesh);
  const su2double base_damping = use_conservative_damping ? max(su2double{0.15}, 0.50 * levelScale) : 1.0;
  const su2double wall_damping = use_conservative_damping ? max(su2double{0.10}, 0.25 * levelScale) : 1.0;

  /*--- Use the adaptive damping factor uniformly across all prolongation levels. ---*/
  const su2double factor = config->GetDamp_Correc_Prolong();

  vector<bool> isWall(geo_fine->GetnPoint(), false);
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetViscous_Wall(iMarker))
      for (auto iVertex = 0ul; iVertex < geo_fine->nVertex[iMarker]; iVertex++)
        isWall[geo_fine->vertex[iMarker][iVertex]->GetNode()] = true;

  SU2_OMP_FOR_STAT(roundUpDiv(geo_fine->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Fine = 0ul; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
    auto* Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
    auto* Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);

    su2double residualMag = 0.0;
    su2double correctionMag = 0.0;
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      /*--- Prevent a fine grid divergence due to a coarse grid divergence ---*/
      if (Residual_Fine[iVar] != Residual_Fine[iVar]) {
        Residual_Fine[iVar] = 0.0;
      }

      const su2double corr = factor * Residual_Fine[iVar];
      residualMag = max(residualMag, fabs(Residual_Fine[iVar]));
      correctionMag = max(correctionMag, fabs(corr));
    }

    su2double correctionScale = 1.0;
    constexpr su2double maxAllowedRatio = 1.25;
    if (residualMag > 1e-30 && correctionMag > 1e-30) {
      const su2double ratio = correctionMag / residualMag;
      if (ratio > maxAllowedRatio) {
        correctionScale = maxAllowedRatio / ratio;
      }
    }

    const su2double localDamping = use_conservative_damping ? (isWall[Point_Fine] ? wall_damping : base_damping) : 1.0;
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      su2double correction = factor * Residual_Fine[iVar];
      correction *= localDamping;
      correction *= correctionScale;

      if (!std::isfinite(correction)) {
        correction = 0.0;
      }

      Solution_Fine[iVar] += correction;
    }
  }
  END_SU2_OMP_FOR

  /*--- DIAGNOSTIC: log the max applied correction (factor * LinSysRes) at fine-grid wall points
   *    vs interior. ---*/
  if (config->GetMGOptions().MG_Smooth_Output && SU2_MPI::GetRank() == MASTER_NODE) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      if (nVar > 2) {
        su2double maxWall0 = 0.0, maxWallN = 0.0, maxWallMom = 0.0;
        su2double maxInter0 = 0.0, maxInterN = 0.0, maxInterMom = 0.0;
        su2double maxWallApply0 = 0.0, maxWallApplyN = 0.0, maxWallApplyMom = 0.0;
        su2double maxInterApply0 = 0.0, maxInterApplyN = 0.0, maxInterApplyMom = 0.0;

        for (auto iPoint = 0ul; iPoint < geo_fine->GetnPointDomain(); iPoint++) {
          const auto* corr = sol_fine->LinSysRes.GetBlock(iPoint);
          const su2double localDamping = isWall[iPoint] ? wall_damping : base_damping;
          const su2double applied0 = fabs(localDamping * factor * corr[0]);
          const su2double appliedN = fabs(localDamping * factor * corr[nVar-1]);
          su2double appliedMom = 0.0;
          for (auto iVar = 1u; iVar < static_cast<unsigned short>(nVar-1); iVar++) {
            appliedMom = max(appliedMom, fabs(localDamping * factor * corr[iVar]));
          }

          if (isWall[iPoint]) {
            maxWall0 = max(maxWall0, fabs(factor * corr[0]));
            maxWallN = max(maxWallN, fabs(factor * corr[nVar-1]));
            maxWallMom = max(maxWallMom, fabs(factor * corr[0]));
            maxWallApply0 = max(maxWallApply0, applied0);
            maxWallApplyN = max(maxWallApplyN, appliedN);
            maxWallApplyMom = max(maxWallApplyMom, appliedMom);
          } else {
            maxInter0 = max(maxInter0, fabs(factor * corr[0]));
            maxInterN = max(maxInterN, fabs(factor * corr[nVar-1]));
            maxInterMom = max(maxInterMom, fabs(factor * corr[0]));
            maxInterApply0 = max(maxInterApply0, applied0);
            maxInterApplyN = max(maxInterApplyN, appliedN);
            maxInterApplyMom = max(maxInterApplyMom, appliedMom);
          }
        }
        auto ratio = [](su2double w, su2double i) { return (i > 1e-30) ? w/i : 0.0; };
        cout << "[MG APPL L" << iMesh << " wall/inter] rho=" << ratio(maxWallApply0, maxInterApply0)
             << "  mom=" << ratio(maxWallApplyMom, maxInterApplyMom)
             << "  E=" << ratio(maxWallApplyN, maxInterApplyN)
             << "  (raw wall/inter: rho=" << maxWall0 << "/" << maxInter0
             << ", E=" << maxWallN << "/" << maxInterN
             << "; applied wall/inter: rho=" << maxWallApply0 << "/" << maxInterApply0
             << ", E=" << maxWallApplyN << "/" << maxInterApplyN
             << "; levelScale=" << levelScale << ", damp(wall/inter)=" << wall_damping << "/" << base_damping << ")\n";
      } else {
        su2double maxWall = 0.0, maxInter = 0.0;
        su2double maxWallApply = 0.0, maxInterApply = 0.0;
        for (auto iPoint = 0ul; iPoint < geo_fine->GetnPointDomain(); iPoint++) {
          const auto* corr = sol_fine->LinSysRes.GetBlock(iPoint);
          const su2double localDamping = isWall[iPoint] ? wall_damping : base_damping;
          const su2double mag = fabs(factor * corr[0]);
          const su2double appliedMag = fabs(localDamping * factor * corr[0]);
          if (isWall[iPoint]) {
            maxWall = max(maxWall, mag);
            maxWallApply = max(maxWallApply, appliedMag);
          } else {
            maxInter = max(maxInter, mag);
            maxInterApply = max(maxInterApply, appliedMag);
          }
        }
        cout << "[MG TURB APPLY L" << iMesh << "] damp(wall/interior)= " << wall_damping << "/" << base_damping
             << "  raw max(wall/interior)= " << maxWall << "/" << maxInter
             << "  applied max(wall/interior)= " << maxWallApply << "/" << maxInterApply
             << "  levelScale=" << levelScale << "\n";
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- MPI the new interpolated solution ---*/

  sol_fine->InitiateComms(geo_fine, config, MPI_QUANTITIES::SOLUTION);
  sol_fine->CompleteComms(geo_fine, config, MPI_QUANTITIES::SOLUTION);

}

void CMultiGridIntegration::SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                    CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  SU2_ZONE_SCOPED

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
  SU2_ZONE_SCOPED

  const su2double *Residual_Fine;

  const unsigned short nVar = sol_coarse->GetnVar();
  const su2double factor = config->GetDamp_Res_Restric();

  su2activevector RestrictedDefect(nVar);

  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);

    RestrictedDefect = su2double(0);
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (auto iVar = 0u; iVar < nVar; iVar++)
        RestrictedDefect[iVar] += factor * Residual_Fine[iVar];
    }
    sol_coarse->GetNodes()->AddRes_TruncError(Point_Coarse, RestrictedDefect.data());
  }
  END_SU2_OMP_FOR

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetViscous_Wall(iMarker)) {
      SU2_OMP_FOR_STAT(32)
      for (auto iVertex = 0ul; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        auto Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        /*--- FIX 1: For turbulence (Dirichlet BC nu_tilde=0 at walls), zero the FULL truncation error.
         *    SetVel_ResTruncError_Zero is a no-op for CTurbVariable (only zeros velocity components
         *    in flow variables). Without this, the FAS forcing term injects spurious residuals at
         *    wall nodes where the coarse grid should simply enforce nu_tilde=0. For flow solvers,
         *    only the velocity components need zeroing (pressure and energy BCs are flux-based). ---*/
        if (sol_coarse->GetnVar() <= 2) {
          /*--- Turbulence/scalar: zero ALL truncation error components at the Dirichlet wall. ---*/
          sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);
        } else {
          /*--- Flow: zero only velocity components (pressure/energy BCs are flux-based). ---*/
          sol_coarse->GetNodes()->SetVel_ResTruncError_Zero(Point_Coarse);
        }
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

void CMultiGridIntegration::RestrictResidualToCoarseGrid(CSolver *sol_fine, CSolver *sol_coarse,
                                                        CGeometry *geo_fine, CGeometry *geo_coarse,
                                                        CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  /*--- This is the standard FAS restriction step: the fine-grid defect is passed to the
   *    coarse-grid problem as a forcing term. The existing SetForcing_Term routine already
   *    implements the conservative volume-weighted transfer and the damping factor in the
   *    same way the original MG cycle expects. ---*/
  SetForcing_Term(sol_fine, sol_coarse, geo_fine, geo_coarse, config, iMesh);
}

void CMultiGridIntegration::ProlongateCorrectionToFineGrid(unsigned short RunTime_EqSystem, CSolver *sol_fine,
                                                           CSolver *sol_coarse, CGeometry *geo_fine,
                                                           CGeometry *geo_coarse, CConfig *config,
                                                           unsigned short iMesh) {
  SU2_ZONE_SCOPED

  /*--- This is the standard FAS prolongation step: build the coarse-grid state correction,
   *    then transfer that correction to the fine-grid residual correction. The original
   *    GetProlongated_Correction routine already performs this transfer in the correct form;
   *    the additional scaling here would be equivalent to changing the correction operator.
   *    Keep the transfer operator unchanged and let the existing damping path control the size. ---*/
  GetProlongated_Correction(RunTime_EqSystem, sol_fine, sol_coarse, geo_fine, geo_coarse, config);
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

        if (Solver_Position == TURB_SOL) {
          /*--- CRITICAL FIX: Enforce Dirichlet BC for turbulence at walls.
           *    SA model requires nu_tilde = 0 at smooth walls. After restriction,
           *    wall values are averaged from fine grid children, violating the BC.
           *    This causes coarse grid solves with incorrect boundary conditions,
           *    producing corrupted corrections that stall convergence. ---*/
          for (auto iVar = 0u; iVar < sol_coarse->GetnVar(); iVar++) {
            sol_coarse->GetNodes()->SetSolution(Point_Coarse, iVar, 0.0);
          }
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

void CMultiGridIntegration::SetRestricted_FrozenSource(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                                        CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  SU2_ZONE_SCOPED

  /*--- Only applicable for turbulence equations ---*/
  if (RunTime_EqSystem != RUNTIME_TURB_SYS) return;

  auto* turbNodes_fine = su2staticcast_p<CTurbVariable*>(sol_fine->GetNodes());
  auto* turbNodes_coarse = su2staticcast_p<CTurbVariable*>(sol_coarse->GetNodes());

  /*--- Compute coarse frozen source density from fine grid using volume-weighted averaging.
   *    This is analogous to eddy viscosity restriction - both are intensive properties. ---*/
  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

    su2double Volume_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);
    su2double SourceDensity_Coarse = 0.0;
    su2double SourceJacobian_Coarse = 0.0;

    /*--- Volume-weighted average of source density and Jacobian from fine grid children ---*/
    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
      su2double Volume_Child = geo_fine->nodes->GetVolume(Point_Fine);
      su2double weight = Volume_Child / Volume_Parent;

      su2double SourceDensity_Fine = turbNodes_fine->GetFrozenSource(Point_Fine);
      SourceDensity_Coarse += SourceDensity_Fine * weight;

      /*--- Restrict Jacobian (intensive property, volume-weighted averaging) ---*/
      su2double SourceJac_Fine = turbNodes_fine->GetFrozenSourceJacobian(Point_Fine);
      SourceJacobian_Coarse += SourceJac_Fine * weight;
    }

    turbNodes_coarse->SetFrozenSource(Point_Coarse, SourceDensity_Coarse);
    turbNodes_coarse->SetFrozenSourceJacobian(Point_Coarse, SourceJacobian_Coarse);
  }
  END_SU2_OMP_FOR

  /*--- MPI communication of restricted frozen source to halo cells.
   *    Although frozen_source is not part of the solution vector, halo cells need correct values
   *    for proper parallel execution, especially when coarse points near partition boundaries
   *    aggregate from fine grid children. We reuse SOLUTION_EDDY communication infrastructure. ---*/

  sol_coarse->InitiateComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION_EDDY);
  sol_coarse->CompleteComms(geo_coarse, config, MPI_QUANTITIES::SOLUTION_EDDY);

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

void CMultiGridIntegration::RestrictTurbEddyViscToCoarseLevels(CGeometry** geometry, CSolver*** solver,
                                                               CConfig* config,
                                                               unsigned short FinestMesh,
                                                               unsigned short nMGLevels) {
  SU2_ZONE_SCOPED

  for (unsigned short iMesh = FinestMesh; iMesh < nMGLevels; iMesh++) {

    CGeometry* geo_fine   = geometry[iMesh];
    CGeometry* geo_coarse = geometry[iMesh + 1];
    CSolver*  sol_fine    = solver[iMesh][TURB_SOL];
    CSolver*  sol_coarse  = solver[iMesh + 1][TURB_SOL];

    /*--- Volume-weighted restriction of mu_t from fine to coarse. ---*/
    SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPointDomain(), omp_get_num_threads()))
    for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {

      const su2double Area_Parent = geo_coarse->nodes->GetVolume(Point_Coarse);
      su2double EddyVisc = 0.0;

      for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
        auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);
        su2double Area_Children = geo_fine->nodes->GetVolume(Point_Fine);
        su2double mu_t = sol_fine->GetNodes()->GetmuT(Point_Fine);
        EddyVisc += mu_t * Area_Children / Area_Parent;
      }

      sol_coarse->GetNodes()->SetmuT(Point_Coarse, EddyVisc);
    }
    END_SU2_OMP_FOR
  }
}
