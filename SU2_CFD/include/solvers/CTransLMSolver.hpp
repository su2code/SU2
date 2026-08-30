/*!
 * \file CTransLMSolver.hpp
 * \brief Headers of the CTransLMSolver class
 * \author A. Aranake
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

#pragma once

#include "CTurbSolver.hpp"
#include "../numerics/turbulent/transition/trans_correlations.hpp"

/*!
 * \class CTransLMSolver
 * \brief Main class for defining the transition model solver.
 * \ingroup Turbulence_Model
 * \author A. Aranake, S. Kang.
 */

class CTransLMSolver final : public CTurbSolver {
private:

  LM_ParsedOptions options;
  TURB_FAMILY TurbFamily;

  TransLMCorrelations TransCorrelations;

  /*!
   * \brief Resolve the compile-time flow indices and dimension, and run the interior edge loop
   *        with the matching CScalarFlux_TransLM instantiation. nVar is fixed at 2 (gamma and
   *        Re_theta), so unlike SA's RunSA this dispatch has one axis fewer to resolve.
   */
  template <class Indices>
  void RunLM(CGeometry* geometry, CSolver** solver_container, CConfig* config, const ScalarFluxOptions& opt);

  template <class Indices, int nDim>
  void RunLM(CGeometry* geometry, CSolver** solver_container, CConfig* config, const ScalarFluxOptions& opt);

  /*!
   * \brief Same dispatch as RunLM, for a boundary's call into BoundaryFluxResidual.
   */
  template <class Indices>
  void RunLM_Boundary(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                      const ScalarFluxOptions& opt, unsigned short val_marker, bool implicit);

  template <class Indices, int nDim>
  void RunLM_Boundary(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                      const ScalarFluxOptions& opt, unsigned short val_marker, bool implicit);

public:
  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransLMSolver() = default;

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry,
                     CSolver **solver_container,
                     CConfig *config,
                     unsigned short iMesh,
                     unsigned short iRKStep,
                     unsigned short RunTime_EqSystem,
                     bool Output) override;

  /*!
   * \brief Computes the effective intermtittency.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry *geometry,
                      CSolver **solver_container,
                      CConfig *config,
                      unsigned short iMesh) override;

  /*!
   * \brief Compute the spatial integration using the CScalarFlux_TransLM edge kernel, which
   *        computes and writes both the convective and the diffusive term of every edge.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Unused, kept only for the boundary conditions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics_container,
                       CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Template(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Impose the Langtry Menter transition wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

 /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry *geometry,
                          CSolver **solver_container,
                          CNumerics *conv_numerics,
                          CNumerics *visc_numerics,
                          CConfig *config,
                          unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
                unsigned short val_marker) override;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry *geometry,
                 CSolver **solver_container,
                 CNumerics *conv_numerics,
                 CNumerics *visc_numerics,
                 CConfig *config,
                 unsigned short val_marker) override;

  /*!
   * \brief Impose the far-field boundary condition, via the CScalarFlux_TransLM edge kernel. Also
   *        used by BC_Outlet, matching this solver's pre-migration behavior of imposing the
   *        far-field state at both MARKER_FAR and MARKER_OUTLET.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] visc_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics *conv_numerics,
                    CNumerics *visc_numerics,
                    CConfig *config,
                    unsigned short val_marker) override;

  /*!
   * \brief Get the value of the intermittency.
   * \return Value of the turbulent kinetic energy.
   */
  inline su2double GetIntermittency_Inf(void) const override { return Solution_Inf[0]; }

  /*!
   * \brief Get the value of the intermittency.
   * \return Value of the turbulent kinetic energy.
   */
  inline su2double GetReThetaT_Inf(void) const override { return Solution_Inf[1]; }

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter, bool val_update_geo) final;

};
