/*!
 * \file CSpeciesSolver.hpp
 * \brief Headers of the CSpeciesSolver class
 * \author T. Kattmann.
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

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../variables/CSpeciesVariable.hpp"
#include "../variables/CEulerVariable.hpp"
#include "../variables/CIncEulerVariable.hpp"
#include "../variables/CNEMOEulerVariable.hpp"
#include "CScalarSolver.hpp"

/*!
 * \brief Main class for defining the species transport solver.
 * \author T. Kattmann.
 * \ingroup Scalar_Transport
 */
class CSpeciesSolver : public CScalarSolver<CSpeciesVariable> {
 protected:
  unsigned short Inlet_Position;             /*!< \brief Column index for scalar variables in inlet files. */
  vector<su2activematrix> Inlet_SpeciesVars; /*!< \brief Species variables at inlet profiles. */
  vector<su2activematrix> Wall_SpeciesVars; /*!< \brief Species variables at  profiles. */
  vector<su2matrix<su2double>> CustomBoundaryScalar;

  /*!
   * \brief Resolve the compile-time flow indices from the regime flag of config, and call f with
   *        a CIndicesTag of the result: f is a generic lambda, `[&](auto tag){ using Indices =
   *        typename decltype(tag)::type; ... }`. Unlike CTurbSolver::DispatchRegime, species
   *        transport is supported for NEMO, so this has a third branch.
   */
  template <class F>
  static void DispatchRegime(const CConfig* config, F&& f) {
    if (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) {
      f(CIndicesTag<CIncEulerVariable::CIndices<unsigned short>>{});
    } else if (config->GetNEMOProblem()) {
      f(CIndicesTag<CNEMOEulerVariable::CIndices<unsigned short>>{});
    } else {
      f(CIndicesTag<CEulerVariable::CIndices<unsigned short>>{});
    }
  }

  /*!
   * \brief Resolve the compile-time flow indices and dimension, and run the interior edge loop
   *        with the matching CScalarFlux_Species instantiation. The equation count is set at
   *        runtime (one per species), so unlike SST's/LM's RunXxx this dispatch never resolves it.
   */
  template <class Indices>
  void RunSpecies(CGeometry* geometry, CSolver** solver_container, CConfig* config, const ScalarFluxOptions& opt);

  template <class Indices, int nDim>
  void RunSpecies(CGeometry* geometry, CSolver** solver_container, CConfig* config, const ScalarFluxOptions& opt);

  /*!
   * \brief Same dispatch as RunSpecies, for a boundary's call into BoundaryFluxResidual.
   */
  template <class Indices>
  void RunSpecies_Boundary(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                           const ScalarFluxOptions& opt, unsigned short val_marker, bool implicit);

  template <class Indices, int nDim>
  void RunSpecies_Boundary(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                           const ScalarFluxOptions& opt, unsigned short val_marker, bool implicit);

  /*!
   * \brief Same dispatch as RunSpecies, for BC_Fluid_Interface's combined fill-and-flux donor loop.
   */
  template <class Indices>
  void RunSpecies_FluidInterface(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                 const ScalarFluxOptions& optConv, const ScalarFluxOptions& optVisc, bool implicit);

  template <class Indices, int nDim>
  void RunSpecies_FluidInterface(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                 const ScalarFluxOptions& optConv, const ScalarFluxOptions& optVisc, bool implicit);

 public:
  /*!
   * \brief Constructor of the class.
   */
  CSpeciesSolver();

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CSpeciesSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter, bool val_update_geo) final;

  void Initialize(CGeometry* geometry, CConfig* config, unsigned short iMesh, unsigned short nVar);

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
  void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh,
                     unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Compute the spatial integration using the CScalarFlux_Species edge kernel, which
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
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(const su2double* val_inlet, unsigned short iMarker, unsigned long iVertex) override;

  /*!
   * \brief Get the set of values imposed at an inlet.
   * \param[in] iMarker - Index of the surface marker.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in,out] val_inlet - vector returning the inlet values for the current vertex.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                             const CGeometry* geometry, su2double* val_inlet) const override;

  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(const CConfig* config, unsigned short iMarker) final;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                 CConfig* config, unsigned short val_marker) final;

  /*!
   * \brief Impose the far-field boundary condition, via the CScalarFlux_Species edge kernel.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] visc_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                    CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) final;

  /*!
   * \brief Impose the isothermal wall Dirichlet boundary condition (value).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container,
                          CNumerics* conv_numerics, CNumerics* visc_numerics,
                          CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the heat flux Neumann wall boundary condition (flux).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container,
                        CNumerics* conv_numerics, CNumerics* visc_numerics,
                        CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Generic wall boundary condition implementation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Wall_Generic(CGeometry* geometry, CSolver** solver_container,
                       CConfig* config, unsigned short val_marker);

  /*--- Note that BC_Sym_Plane, BC_HeatFlux_Wall, BC_Isothermal_Wall are all treated as zero-flux BC for the
   * mass-factions, which can be fulfilled by no additional residual contribution on these nodes.
   * If a specified mass-fractions flux (like BC_HeatFlux_Wall) or a constant mass-fraction on the boundary
   * (like BC_Isothermal_Wall) are implemented the respective CHeatSolver implementations can act as a  starting
   * point ---*/

  /*!
   * \brief Source term computation for axisymmetric flow.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Container for description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics_container, CConfig* config,
                       unsigned short iMesh) override;

  /*!
   * \brief Set the initial condition for the Species transport problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] TimeIter - Time iteration.
   */
  void SetInitialCondition(CGeometry **geometry,
                           CSolver ***solver_container,
                           CConfig *config,
                           unsigned long TimeIter) override;

  /*!
   * \brief Impose the fluid interface (sliding mesh) boundary condition, via the
   *        CScalarFlux_Species edge kernel. The convective term is a per-donor weighted average,
   *        computed in the same pass that fills the ghost row of each donor; the diffusive term is
   *        computed once per vertex, after the donor loop, from the interior point's own
   *        diffusivity mirrored into the ghost row (matching the pre-migration behavior, which
   *        likewise read the same point's diffusivity for both sides of the edge).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] visc_numerics - Unused, kept only for the boundary condition dispatch.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                          CNumerics *visc_numerics, CConfig *config) final;

  /*!
   * \brief Set custom boundary scalar values from Python.
   * \param[in] val_marker - Boundary marker index
   * \param[in] val_vertex - Boundary vertex index
   * \param[in] val_customBoundaryScalar - Vector of scalar values
   */
  inline void SetCustomBoundaryScalar(unsigned short val_marker, unsigned long val_vertex,
                                      vector<passivedouble> val_customBoundaryScalar) final {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      CustomBoundaryScalar[val_marker](val_vertex, iVar) = val_customBoundaryScalar[iVar];
    }
  }

};
