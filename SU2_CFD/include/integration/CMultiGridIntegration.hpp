/*!
 * \file CMultiGridIntegration.hpp
 * \brief Declaration of class for time integration using a multigrid method.
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

#include "CIntegration.hpp"

/*!
 * \class CMultiGridIntegration
 * \ingroup Drivers
 * \brief Class for time integration using a multigrid method.
 * \author F. Palacios
 */
class CMultiGridIntegration final : public CIntegration {
public:
  /*!
   * \brief Constructor of the class.
   */
  CMultiGridIntegration();

  /*!
   * \brief This subroutine calls the MultiGrid_Cycle and also prepare the multigrid levels and the monitoring.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void MultiGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                           CNumerics ******numerics_container, CConfig **config,
                           unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) override;

  /*!
   * \brief Record CONV_FIELD on the active Full-MG level and decide whether the level is done,
   *        either because it has converged by MG_STARTUP_CONVERGENCE orders of magnitude or
   *        because it has stopped reducing the error (MG_STARTUP_STAGNATION).
   * \param[in] convFieldValue - Value of CONV_FIELD on the active level (log10 of the residual).
   * \param[in] config - Definition of the particular problem.
   */
  void MonitorFullMG_Startup(passivedouble convFieldValue, const CConfig *config) override;

private:
  /*!
   * \brief Perform a Full-Approximation Storage (FAS) Multigrid.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] mu - Variable for controlling the kind of multigrid algorithm.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void MultiGrid_Cycle(CGeometry ****geometry, CSolver *****solver_container,
                       CNumerics ******numerics_container, CConfig **config,
                       unsigned short iMesh, unsigned short mu, unsigned short RunTime_EqSystem,
                       unsigned short iZone, unsigned short iInst);

  /*!
   * \brief Compute the forcing term.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                       CGeometry *geo_coarse, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Add the truncation error to the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] flow - Flow solution.
   */
  void SetResidual_Term(CGeometry *geometry, CSolver *flow);

  /*!
   * \brief Set the value of the corrected fine grid solution.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Compute the gradient in coarse grid using the fine grid information.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                              CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*!
   * \brief Compute the non-dimensional parameters.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   * \param[in] monitor - value of the non-dimensional parameters for monitoring the convergence.
   */
  void NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container,
                                 CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem,
                                 su2double *monitor);

  /*!
   * \brief Compute the fine solution from a coarse solution.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                               CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);


  /*!
   * \brief Apply post-smoothing iterations on the fine grid after prolongation.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] solver_fine - Pointer to the solver on the fine grid.
   * \param[in] numerics_fine - Description of the numerical method on the fine grid.
   * \param[in] geometry_fine - Geometrical definition of the fine grid.
   * \param[in] solver_container_fine - Container with all solvers on the fine grid.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKLimit - Number of Runge-Kutta steps.
   */
  void PostSmoothing(unsigned short RunTime_EqSystem, CSolver* solver_fine, CNumerics** numerics_fine,
                     CGeometry* geometry_fine, CSolver** solver_container_fine, CConfig *config,
                     unsigned short iMesh, unsigned short iRKLimit);

  /*!
   * \brief Apply pre-smoothing iterations on the fine grid before restriction.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] geometry - Geometrical definition of the problem (all levels).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config_container - Definition of the particular problems.
   * \param[in] solver_fine - Pointer to the solver on the fine grid.
   * \param[in] numerics_fine - Description of the numerical method on the fine grid.
   * \param[in] geometry_fine - Geometrical definition of the fine grid.
   * \param[in] solver_container_fine - Container with all solvers on the fine grid.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iZone - Index of the zone.
   * \param[in] iRKLimit - Number of Runge-Kutta steps.
   */
  void PreSmoothing(unsigned short RunTime_EqSystem, CGeometry**** geometry, CSolver***** solver_container,
                    CConfig **config_container, CSolver* solver_fine, CNumerics** numerics_fine,
                    CGeometry* geometry_fine, CSolver** solver_container_fine, CConfig *config,
                    unsigned short iMesh, unsigned short iZone, unsigned short iRKLimit);

  /*!
   * \brief Compute the fine grid correction from the coarse solution.
   * \param[out] sol_fine - Pointer to the solution on the fine grid.
   * \param[in] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                                 CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*!
   * \brief Do an implicit smoothing of the prolongated correction.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] solution - Container vector with all the solutions on the finest grid.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_nSmooth - Number of smoothing iterations.
   * \param[in] val_smooth_coeff - Relaxation factor.
   * \param[in] config - Definition of the particular problem.
   */
  void SmoothProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                    unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config,
                                    unsigned short iMesh);

  /*!
   * \brief Restrict solution from fine grid to a coarse grid.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] sol_fine - Pointer to the solution on the fine grid.
   * \param[out] sol_coarse - Pointer to the solution on the coarse grid.
   * \param[in] geo_fine - Geometrical definition of the fine grid.
   * \param[in] geo_coarse - Geometrical definition of the coarse grid.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse,
                              CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config);

  /*!
   * \brief Initialize the adjoint solution using the primal problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Iteration - Current iteration.
   */
  void Adjoint_Setup(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                     unsigned short RunTime_EqSystem, unsigned long Iteration, unsigned short iZone);

  /*!
   * \brief Adapt both restriction and prolongation damping factors from the global-trend signal.
   *
   * Uses the cross-cycle EMA ratio (crossCycleRatio = fine_d0 / EMA(fine_d0)) to detect
   * long-term convergence or divergence, then adjusts both \c Damp_Res_Restric and
   * \c Damp_Correc_Prolong with a single shared signal.  The EMA filters per-cycle noise;
   * no per-level aggregation or floor counter is needed.
   *
   * \param[in,out] config          - Problem configuration.
   * \param[in]     crossCycleRatio - Current fine_d0 divided by the EMA of fine_d0.
   */
  void adaptDampingFactors(CConfig* config, passivedouble crossCycleRatio);

  /*!
   * \brief Set the CFL of every coarse multigrid level. level i starts with the final CFL of the coarser level i+1, and then linearly increases to the new final CFL for that level. Final CFL is determined from the CFL_NUMBER and MG_CFL_SCALING.
   *
   * \param[in]     geometry         - Geometrical definition of the problem.
   * \param[in,out] solver_container - Container vector with all the solutions.
   * \param[in,out] config           - Definition of the particular problem.
   * \param[in]     RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in]     iZone            - Zone index.
   * \param[in]     iInst            - Instance index.
   * \param[in]     FinestMesh       - Currently active finest mesh level.
   * \param[in]     FullMG           - Whether the Full-MG cycle is active.
   */
  void SetCoarseGridCFL(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
                        unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst,
                        unsigned short FinestMesh, bool FullMG);

  /*!
   * \brief Helper function for early-exit logic during pre/post-smoothing.
   * \param[in] iSmooth - Current smoothing iteration index.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] defect - Current RMS defect value.
   * \param[in] mgOpts - Reference to multigrid options.
   * \param[in] stag_tol - Stagnation tolerance value.
   * \param[in] early_exit - Whether early exit is enabled.
   * \param[out] lastRMS - Array to store RMS values [start, end].
   * \param[out] exitReason - Character for early exit reason ('T', 'S', 'A', or ' ').
   * \param[out] worstStepRatio - Worst step-to-step ratio seen.
   * \param[out] worstStep - Iteration number of worst step.
   */
  void prePostEarlyExit(unsigned short iSmooth, unsigned short iMesh,
                        passivedouble defect, const CMGOptions& mgOpts,
                        passivedouble stag_tol, bool early_exit,
                        passivedouble lastRMS[2], char& exitReason,
                        passivedouble& worstStepRatio, unsigned short& worstStep);

  static constexpr int MAX_MG_LEVELS = 10;

  /*--- Early-exit smoothing state (shared across OMP threads via master write + barrier). ---*/
  bool mg_early_exit_flag = false;              /*!< \brief Shared flag for early exit across OMP threads. */
  passivedouble mg_initial_smooth_rms = 0.0; /*!< \brief Initial RMS residual before current smoothing phase (FAS). */
  passivedouble mg_prev_smooth_rms = 0.0;    /*!< \brief RMS residual from previous smoothing step; used for stagnation detection. */
  passivedouble mg_fine_rms_ema = 0.0;      /*!< \brief EMA of fine-grid pre-smooth RMS across cycles; cross-cycle trend signal. */
  passivedouble last_crossCycleRatio = 1.0; /*!< \brief crossCycleRatio from the most recent cycle; stored for display only. */

  /*--- Actual iteration counts per MG level, filled each cycle for the compact output summary. ---*/
  unsigned short lastPreSmoothIters[MAX_MG_LEVELS+1] = {};
  unsigned short lastPostSmoothIters[MAX_MG_LEVELS+1] = {};
  unsigned short lastCorrecSmoothIters[MAX_MG_LEVELS+1] = {};
  /*--- Early-exit reason per level: 'T'=threshold, 'S'=stagnation, ' '=ran to completion. ---*/
  char lastPreSmoothExitReason[MAX_MG_LEVELS+1]  = {};
  char lastPostSmoothExitReason[MAX_MG_LEVELS+1] = {};

  /*--- Per-level start/end RMS residual for adaptive damping. ---*/
  passivedouble lastPreSmoothRMS[MAX_MG_LEVELS+1][2] = {};
  passivedouble lastPostSmoothRMS[MAX_MG_LEVELS+1][2] = {};
  passivedouble lastCorrecSmoothRMS[MAX_MG_LEVELS+1][2] = {};

  /*--- Per-level worst step-to-step amplification seen inside a smoothing phase.
   *    step==0 means no intra-smoother ratio was available (fewer than 2 sweeps). ---*/
  passivedouble lastPreSmoothWorstStepRatio[MAX_MG_LEVELS+1] = {};
  passivedouble lastPostSmoothWorstStepRatio[MAX_MG_LEVELS+1] = {};
  unsigned short lastPreSmoothWorstStep[MAX_MG_LEVELS+1] = {};
  unsigned short lastPostSmoothWorstStep[MAX_MG_LEVELS+1] = {};

  /*--- FMG startup CFL ramp bookkeeping: tracks the currently active FMG level
   *    and the InnerIter at which it became active, so its CFL can be ramped
   *    linearly towards the next (finer) level's target over MG_Startup_Iter
   *    iterations instead of jumping discontinuously at promotion. ---*/
  unsigned short mg_ramp_last_FinestMesh = MAX_MG_LEVELS + 1; /*!< \brief FinestMesh observed on the previous call; sentinel forces a reset on the first call. */
  unsigned long mg_ramp_level_start_iter = 0;                 /*!< \brief InnerIter at which the currently active FMG level became active. */

  /*--- User-configured damping factors, captured before any adaptation so they can be
   *    restored whenever the active FMG level changes (the cross-cycle EMA that drives
   *    the adaptation is only meaningful within a single level). ---*/
  bool mg_damp_initial_captured = false;  /*!< \brief Whether the configured damping factors have been stored yet. */
  su2double mg_damp_restric_initial = 0.0;  /*!< \brief MG_DAMP_RESTRICTION as configured. */
  su2double mg_damp_prolong_initial = 0.0;  /*!< \brief MG_DAMP_PROLONGATION as configured. */

  /*--- Full-MG startup promotion state, fed by MonitorFullMG_Startup once per iteration and
   *    consumed by the promotion check at the top of the next one. All of it tracks the currently
   *    active level and is reset when that level changes.
   *
   *    The values are CONV_FIELD as the output module reports it, i.e. log10 of the residual, so a
   *    drop of N orders of magnitude is a fall of N in these numbers. ---*/
  bool mg_startup_conv_captured = false;    /*!< \brief Whether the level's starting value has been recorded yet. */
  passivedouble mg_startup_conv_start = 0.0;/*!< \brief CONV_FIELD when the active level became active. */
  bool mg_startup_conv_prev_set = false;    /*!< \brief Whether a previous-iteration value is available. */
  passivedouble mg_startup_conv_prev = 0.0; /*!< \brief CONV_FIELD on the previous iteration. */
  unsigned long mg_startup_stall_count = 0; /*!< \brief Consecutive iterations without useful reduction. */
  bool mg_startup_promote = false;          /*!< \brief Set once a promotion criterion is met; consumed at the next promotion check. */
  bool mg_startup_promoted_on_stall = false;/*!< \brief Whether the pending promotion was triggered by stagnation (for reporting). */

};
