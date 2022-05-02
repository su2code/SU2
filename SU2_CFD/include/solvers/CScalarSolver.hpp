/*!
 * \file CScalarSolver.hpp
 * \brief Headers of the CScalarSolver class
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include <vector>

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../variables/CScalarVariable.hpp"
#include "CSolver.hpp"

/*!
 * \class CScalarSolver
 * \brief Main class for defining a scalar solver.
 * \tparam VariableType - Class of variable used by the solver inheriting from this template.
 */
template <class VariableType>
class CScalarSolver : public CSolver {
 protected:
  enum : size_t { MAXNDIM = 3 };      /*!< \brief Max number of space dimensions, used in some static arrays. */
  static constexpr size_t MAXNVAR = VariableType::MAXNVAR; /*!< \brief Max number of variables, for static arrays. */
  enum : size_t { MAXNVARFLOW = 12 }; /*!< \brief Max number of flow variables, used in some static arrays. */

  enum : size_t { OMP_MAX_SIZE = 512 }; /*!< \brief Max chunk size for light point loops. */
  enum : size_t { OMP_MIN_SIZE = 32 };  /*!< \brief Min chunk size for edge loops (max is color group size). */

  unsigned long omp_chunk_size; /*!< \brief Chunk size used in light point loops. */

  su2double lowerlimit[MAXNVAR]; /*!< \brief contains lower limits for turbulence variables. Note that ::min()
                                             returns the smallest positive value for floats. */
  su2double upperlimit[MAXNVAR]; /*!< \brief contains upper limits for turbulence variables. */

  su2double Solution_Inf[MAXNVAR]; /*!< \brief Far-field solution. */

  const bool Conservative; /*!< \brief Transported Variable is conservative. Solution has to be multiplied with rho. */

  /*--- Shallow copy of grid coloring for OpenMP parallelization. ---*/

#ifdef HAVE_OMP
  vector<GridColor<> > EdgeColoring; /*!< \brief Edge colors. */
  bool ReducerStrategy = false;      /*!< \brief If the reducer strategy is in use. */
#else
  array<DummyGridColor<>, 1> EdgeColoring;
  /*--- Never use the reducer strategy if compiling for MPI-only. ---*/
  static constexpr bool ReducerStrategy = false;
#endif

  /*--- Edge fluxes for reducer strategy (see the notes in CEulerSolver.hpp). ---*/
  CSysVector<su2double> EdgeFluxes; /*!< \brief Flux across each edge. */

  /*!
   * \brief The highest level in the variable hierarchy this solver can safely use.
   */
  VariableType* nodes = nullptr;

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }

  /*!
   * \brief Compute the viscous flux for the turbulent equation at a particular edge.
   * \tparam SolverSpecificNumericsFunc - lambda-function, that implements solver specific contributions to numerics.
   * \note The functor has to implement (iPoint, jPoint)
   * \param[in] iEdge - Edge for which we want to compute the flux
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  template <class SolverSpecificNumericsFunc>
  FORCEINLINE void Viscous_Residual_impl(const SolverSpecificNumericsFunc& SolverSpecificNumerics, unsigned long iEdge,
                                         CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                         CConfig* config) {
    const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
    CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

    /*--- Points in edge ---*/

    auto iPoint = geometry->edges->GetNode(iEdge, 0);
    auto jPoint = geometry->edges->GetNode(iEdge, 1);

    /*--- Points coordinates, and normal vector ---*/

    numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(jPoint));
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), flowNodes->GetPrimitive(jPoint));

    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

    numerics->SetScalarVar(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(jPoint));

    /*--- Call Numerics contribution which are Solver-Specifc. Implemented in the caller: Viscous_Residual.  ---*/

    SolverSpecificNumerics(iPoint, jPoint);

    /*--- Compute residual, and Jacobians ---*/

    auto residual = numerics->ComputeResidual(config);

    if (ReducerStrategy) {
      EdgeFluxes.SubtractBlock(iEdge, residual);
      if (implicit) Jacobian.UpdateBlocksSub(iEdge, residual.jacobian_i, residual.jacobian_j);
    } else {
      LinSysRes.SubtractBlock(iPoint, residual);
      LinSysRes.AddBlock(jPoint, residual);
      if (implicit) Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }
  }

  /*!
   * \brief Gradient and Limiter computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether the residual+jacobian should be zeroed.
   */
  void CommonPreprocessing(CGeometry *geometry, const CConfig *config, const bool Output);

 private:
  /*!
   * \brief Compute the viscous flux for the turbulent equation at a particular edge.
   * \param[in] iEdge - Edge for which we want to compute the flux
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \note Calls a generic implementation after defining a SolverSpecificNumerics object.
   */
  inline virtual void Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                       CNumerics* numerics, CConfig* config) {
    /*--- Define an empty object for solver specific numerics contribution. In case there are none, this default
     *--- implementation will be called ---*/
    auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {};

    /*--- Now instantiate the generic implementation with the functor above. ---*/

    Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
  }
  using CSolver::Viscous_Residual; /*--- Silence warning ---*/

  /*!
   * \brief Sum the edge fluxes for each cell to populate the residual vector, only used on coarse grids.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SumEdgeFluxes(CGeometry* geometry);

  /*!
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over
   * a nonlinear iteration for stability. Default value 1.0 set in ctor of CScalarVariable.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeUnderRelaxationFactor(const CConfig* config) {}

 public:
  /*!
   * \brief Destructor of the class.
   */
  ~CScalarSolver() override;

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CScalarSolver(CGeometry* geometry, CConfig* config, bool conservative);

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics_container, CConfig* config,
                       unsigned short iMesh) override;

  /*!
   * \brief Impose the Symmetry Plane boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline void BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                           CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override {
    /*--- Convective and viscous fluxes across symmetry plane are equal to zero. ---*/
  }

  /*!
   * \brief Impose via the residual the Euler wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  inline void BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                            CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) final {
    /*--- Convective fluxes across euler wall are equal to zero. ---*/
  }

  /*!
   * \brief Impose the supersonic inlet boundary condition (same as inlet, see BC_Inlet).
   */
  void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                           CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) final {
    BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
  }

  /*!
   * \brief Impose the supersonic outlet boundary condition (same as outlet, see BC_Outlet).
   */
  void BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                            CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) final {
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
  }

  /*!
   * \brief Impose a periodic boundary condition by summing contributions from the complete control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Periodic(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CConfig* config) final;

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  inline void SetFreeStream_Solution(const CConfig* config) final {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      nodes->SetSolution(iPoint, Solution_Inf);
    }
    END_SU2_OMP_FOR
  }

  /*!
   * \brief Prepare an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void PrepareImplicitIteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) final;

  /*!
   * \brief Complete an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void CompleteImplicitIteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) final;

  /*!
   * \brief Update the solution using the explicit Euler scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) final;

  /*!
   * \brief Update the solution using an implicit solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

  /*!
   * \brief Set the total residual adding the term that comes from the Dual Time-Stepping Strategy.
   * \param[in] geometry - Geometric definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
  void SetResidual_DualTime(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iRKStep,
                            unsigned short iMesh, unsigned short RunTime_EqSystem) final;

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                           bool val_update_geo) override = 0;

  /*!
   * \brief Scalar solvers support OpenMP+MPI.
   */
  inline bool GetHasHybridParallel() const override { return true; }
};
