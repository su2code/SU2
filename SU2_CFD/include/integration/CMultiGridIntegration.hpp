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
   * \brief Compute adaptive CFL for multigrid coarse levels.
   * \param[in] config - Problem configuration.
   * \param[in] solver_coarse - Coarse grid solver.
   * \param[in] geometry_coarse - Coarse grid geometry.
   * \param[in] iMesh - Current multigrid level.
   * \param[in] CFL_fine - Fine grid CFL value (passive).
   * \param[in] CFL_coarse_current - Current coarse grid CFL value (passive).
   * \param[in] rms_res_coarse - Coarse-grid RMS residual (already MPI-reduced, from lastPreSmoothRMS).
   * \return New CFL value for the coarse grid.
   */
  passivedouble computeMultigridCFL(CConfig* config, unsigned short iMesh,
                                     passivedouble CFL_fine, passivedouble CFL_coarse_current,
                                     passivedouble rms_res_coarse);

  /*!
   * \brief Adapt the residual restriction damping factor.
   *
   * Uses \c lastPreSmoothIters[] (filled by the previous multigrid cycle) to assess
   * whether the pre-smoother is converging fast or slow on coarse levels, then adjusts
   * \c Damp_Res_Restric in \p config accordingly.
   *
   * Signal logic:
   *  - any coarse level ran its full configured iterations: reduce damping
   *  - all coarse levels exited early: increase damping
   *  - mixed (some full, some partial): no change
   *
   * \param[in,out] config - Problem configuration.
   */
  void adaptRestrictionDamping(CConfig* config);

  /*!
   * \brief Adapt the correction prolongation damping factor.
   *
   * Uses \c lastCorrecSmoothIters[] (filled by the previous multigrid cycle) to assess
   * whether the correction smoother is struggling or converging fast,
   * then adjusts \c Damp_Correc_Prolong in \p config accordingly.
   *
   * Signal logic:
   *  - any level ran its full correction-smooth iterations: reduce damping
   *  - all levels exited early: increase damping
   *  - mixed: no change
   *
   * \param[in,out] config - Problem configuration; \c SetDamp_Correc_Prolong is called to persist the result.
   */
  void adaptProlongationDamping(CConfig* config);

  /*--- CFL adaptation state variables.
   *    These must be passivedouble: AD::Reset() clears the tape between adjoint recordings,
   *    but class members survive. If these were su2double their stale AD indices would
   *    reference the cleared tape, causing invalid memory access during the backward pass. ---*/
  static constexpr int MAX_MG_LEVELS = 10;
  passivedouble current_avg[MAX_MG_LEVELS] = {};
  passivedouble prev_avg[MAX_MG_LEVELS] = {};
  passivedouble last_res[MAX_MG_LEVELS] = {};
  bool last_was_increase[MAX_MG_LEVELS] = {};
  int oscillation_count[MAX_MG_LEVELS] = {};
  unsigned long last_check_iter[MAX_MG_LEVELS] = {};
  unsigned long last_update_iter[MAX_MG_LEVELS] = {};
  unsigned long last_reset_iter = std::numeric_limits<unsigned long>::max();

  /*--- Early-exit smoothing state (shared across OMP threads via master write + barrier). ---*/
  bool mg_early_exit_flag = false;             /*!< \brief Shared flag for early exit across OMP threads. */
  passivedouble mg_initial_smooth_rms = 0.0;  /*!< \brief Initial RMS before current smoothing phase. */
  passivedouble mg_last_smooth_rms = 0.0;     /*!< \brief Last computed RMS; cached to avoid redundant Allreduce. */

  /*--- Actual iteration counts per MG level, filled each cycle for the compact output summary. ---*/
  unsigned short lastPreSmoothIters[MAX_MG_LEVELS+1] = {};
  unsigned short lastPostSmoothIters[MAX_MG_LEVELS+1] = {};
  unsigned short lastCorrecSmoothIters[MAX_MG_LEVELS+1] = {};

  /*--- Per-level residual progress flags: true if the final RMS after that phase was lower
   *    than the initial RMS.  Used by the adaptive damping routines to distinguish
   *    "hit max iters but still converging" from "hit max iters and stagnated". ---*/
  bool lastPreSmoothProgress[MAX_MG_LEVELS+1] = {};
  bool lastPostSmoothProgress[MAX_MG_LEVELS+1] = {};
  bool lastCorrecSmoothProgress[MAX_MG_LEVELS+1] = {};

  /*--- Per-level start/end RMS for the compact output summary.
   *    [0] = initial RMS before smoothing, [1] = final RMS after smoothing.
   *    Filled unconditionally (early-exit path and exhaustion path).
   *    Must be passivedouble: class members survive tape resets; su2double would
   *    carry stale AD indices referencing a cleared tape. ---*/
  passivedouble lastPreSmoothRMS[MAX_MG_LEVELS+1][2] = {};
  passivedouble lastPostSmoothRMS[MAX_MG_LEVELS+1][2] = {};
  passivedouble lastCorrecSmoothRMS[MAX_MG_LEVELS+1][2] = {};

};
