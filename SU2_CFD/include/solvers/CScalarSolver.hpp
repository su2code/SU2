/*!
 * \file CScalarSolver.hpp
 * \brief Headers of the CScalarSolver class
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

#include <vector>

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../numerics/scalar/scalar_edge_flux.hpp"
#include "../variables/CScalarVariable.hpp"
#include "../variables/CEulerVariable.hpp"
#include "../variables/CFlowVariable.hpp"
#include "../variables/CGhostFlowVariable.hpp"
#include "../variables/CIncEulerVariable.hpp"
#include "../variables/CPrimitiveIndices.hpp"
#include "CSolver.hpp"

/*!
 * \brief Carries a type through a value, so a runtime branch can hand a compile-time type to a
 *        generic lambda (its parameter deduces as CTypeTag<T>, and the lambda recovers T as
 *        decltype(tag)::type). Standing in for a C++20 template lambda, which this project's
 *        C++17 baseline does not have.
 */
template <class T>
struct CTypeTag {
  using type = T;
};

/*!
 * \brief Main class for defining a scalar solver.
 * \tparam VariableType - Class of variable used by the solver inheriting from this template.
 * \ingroup Scalar_Transport
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

  su2double lowerlimit[MAXNVAR]; /*!< \brief contains lower limits for scalar variables. */
  su2double upperlimit[MAXNVAR]; /*!< \brief contains upper limits for scalar variables. */

  su2double Solution_Inf[MAXNVAR]; /*!< \brief Far-field solution. */

  const bool Conservative; /*!< \brief Transported Variable is conservative. Solution has to be multiplied with rho. */
  const bool BoundedScalar; /*!< \brief Whether the derived solver uses the bounded-scalar convective scheme. */

  const CPrimitiveIndices<unsigned short> prim_idx; /*!< \brief Indices of the primitive flow variables. */

  vector<su2matrix<su2double*> > SlidingState; // vector of matrix of pointers... inner dim alloc'd elsewhere (welcome, to the twilight zone)
  vector<vector<int> > SlidingStateNodes;

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
  CSysVector<su2double> EdgeFluxesDiff; /*!< \brief Flux difference between ij and ji for non-conservative discretisation. */

  /*--- Ghost states of the marker currently being processed by a boundary, indexed by vertex
   * and sized to the largest marker; same container types as the interior ones, so the flux
   * kernels read a boundary through the same accessors as an interior edge. Boundary loops run
   * one marker at a time, parallel over its vertices, so the buffers are written and consumed
   * before the next marker reaches them (see BoundaryFluxResidual). ---*/
  unique_ptr<VariableType> ghostNodes;         /*!< \brief Allocated by the derived solver, whose VariableType constructor it alone knows how to call. */
  unique_ptr<CGhostFlowVariable> ghostFlowNodes; /*!< \brief Sized from the flow solver, see the constructor. */
  su2activematrix ghostNormal; /*!< \brief Outward normals, sign flipped from the vertex normals. */
  su2activematrix ghostCoord;  /*!< \brief Reflected coordinates, read by the diffusion sites. */
  su2vector<uint8_t> ghostSkip; /*!< \brief Whether a vertex contributes no flux, set by the fill pass. */

  /*!
   * \brief The highest level in the variable hierarchy this solver can safely use.
   */
  VariableType* nodes = nullptr;

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }

  /*!
   * \brief Generic implementation of the fluid interface boundary condition for scalar solvers.
   * \tparam SolverSpecificNumericsFunc - lambda that implements solver specific contributions to viscous numerics.
   * \note The functor has to implement (iPoint)
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  template <class SolverSpecificNumericsFunc>
  void BC_Fluid_Interface_impl(const SolverSpecificNumericsFunc& SolverSpecificNumerics, CGeometry *geometry,
                               CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                               CConfig *config) {
    if (solver_container[FLOW_SOL] == nullptr) return;

    const auto nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();
    su2activevector PrimVar_j(nPrimVar);
    su2double solution_j[MAXNVAR] = {0.0};

    for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) != FLUID_INTERFACE) continue;

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (!geometry->nodes->GetDomain(iPoint)) continue;

        const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
        const auto nDonorVertex = GetnSlidingStates(iMarker,iVertex);

        su2double Normal[MAXNDIM] = {0.0};
        for (auto iDim = 0u; iDim < nDim; iDim++)
          Normal[iDim] = -geometry->vertex[iMarker][iVertex]->GetNormal()[iDim];

        su2double* PrimVar_i = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

        auto* Jacobian_i = Jacobian.GetBlock(iPoint, iPoint);

        /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

        for (auto jVertex = 0; jVertex < nDonorVertex; jVertex++) {

          for (auto iVar = 0u; iVar < nPrimVar; iVar++)
            PrimVar_j[iVar] = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, iVar, jVertex);

          /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

          const su2double weight = solver_container[FLOW_SOL]->GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

          /*--- Set primitive variables ---*/

          conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j.data() );

          /*--- Set the scalar variable states ---*/

          for (auto iVar = 0u; iVar < nVar; ++iVar)
            solution_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);

          conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);

          /*--- Set the normal vector ---*/

          conv_numerics->SetNormal(Normal);

          if (dynamic_grid)
            conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

          if (conv_numerics->GetBoundedScalar()) {
            const su2double* velocity = &PrimVar_j[prim_idx.Velocity()];
            const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
            conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, true, density, velocity, Normal));
          }

          auto residual = conv_numerics->ComputeResidual(config);

          /*--- Accumulate the residuals to compute the average ---*/

          for (auto iVar = 0u; iVar < nVar; iVar++) {
            LinSysRes(iPoint,iVar) += weight*residual[iVar];
            for (auto jVar = 0u; jVar < nVar; jVar++)
              Jacobian_i[iVar*nVar+jVar] += SU2_TYPE::GetValue(weight*residual.jacobian_i[iVar][jVar]);
          }
        }

        /*--- Set the normal vector and the coordinates ---*/

        visc_numerics->SetNormal(Normal);
        su2double Coord_Reflected[MAXNDIM];
        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                geometry->nodes->GetCoord(iPoint), Coord_Reflected);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

        /*--- Primitive variables ---*/

        visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j.data());

        /*--- Scalar variables and their gradients ---*/

        visc_numerics->SetScalarVar(nodes->GetSolution(iPoint), solution_j);
        visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), nodes->GetGradient(iPoint));

        /*--- Allow derived solvers to set more variables in numerics. ---*/

        SolverSpecificNumerics(iPoint);

        /*--- Compute and update residual ---*/

        auto residual = visc_numerics->ComputeResidual(config);

        LinSysRes.SubtractBlock(iPoint, residual);

        /*--- Jacobian contribution for implicit integration ---*/

        Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

      }
      END_SU2_OMP_FOR
    }
  }

  /*!
   * \brief Applies a convective flux correction to negate the effects of flow divergence at a BC node.
   * \note This function should be used for nodes that are part of a boundary marker, it computes a mass flux
   * from density and velocity at the node, and the outward-pointing normal (-1 * normal of vertex).
   * \return The mass flux.
   */
  inline su2double BoundedScalarBCFlux(unsigned long iPoint, bool implicit, const su2double& density,
                                       const su2double* velocity, const su2double* normal) {
    const su2double edgeMassFlux = density * GeometryToolbox::DotProduct(nDim, velocity, normal);
    LinSysRes.AddBlock(iPoint, nodes->GetSolution(iPoint), -edgeMassFlux);
    if (implicit) Jacobian.AddVal2Diag(iPoint, -edgeMassFlux);
    return edgeMassFlux;
  }

  /*!
   * \brief Move solution to previous time levels (for restarts).
   */
  void PushSolutionBackInTime(unsigned long TimeIter, bool restart, CSolver*** solver_container,
                              CGeometry** geometry, CConfig* config);

  /*!
   * \brief Gradient and Limiter computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether the residual+jacobian should be zeroed.
   */
  void CommonPreprocessing(CGeometry *geometry, const CConfig *config, const bool Output);

  /*!
   * \brief Sum the edge fluxes for each cell to populate the residual vector, only used on coarse grids.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SumEdgeFluxes(const CGeometry* geometry);

  /*!
   * \brief Generic interior edge loop for a scalar model expressed through the CUpwScalarBase
   *        CRTP chain (numerics/scalar/scalar_edge_flux.hpp), driving both the convective and
   *        the diffusive term of every edge with a single kernel, under either update strategy.
   * \tparam Scheme - A model instantiated from CUpwScalarBase, e.g. CScalarFlux_SA<...>.
   * \param[in] opt - Loop invariant flags built by the caller from the current CConfig state.
   */
  template <class Scheme>
  void EdgeFluxResidual(const CGeometry* geometry, CSolver** solver_container, const CConfig* config,
                        const ScalarFluxOptions& opt);

  /*!
   * \brief Write the four flow primitives the flux kernels read into one row of ghostFlowNodes.
   * \param[in] iVertex - Vertex of the marker currently being processed.
   * \param[in] V - Row of flow primitives to copy from (e.g. GetCharacPrimVar's or a sliding state's).
   */
  inline void SetGhostPrimitives(unsigned long iVertex, const su2double* V) {
    auto* ghostV = ghostFlowNodes->GetPrimitive(iVertex);
    ghostV[prim_idx.Density()] = V[prim_idx.Density()];
    for (auto iDim = 0u; iDim < nDim; ++iDim) ghostV[prim_idx.Velocity() + iDim] = V[prim_idx.Velocity() + iDim];
    ghostV[prim_idx.LaminarViscosity()] = V[prim_idx.LaminarViscosity()];
    ghostV[prim_idx.EddyViscosity()] = V[prim_idx.EddyViscosity()];
  }

  /*!
   * \brief Generic boundary flux pass, run after a boundary's fill pass has written the ghost
   *        row, the outward normal and (for the diffusion sites) the ghost gradient of every
   *        vertex of the marker. The ghost point has no row, so only the contribution to the
   *        interior point is assembled.
   * \tparam Scheme - Same model the interior loop uses, instantiated with muscl false.
   */
  template <class Scheme>
  void BoundaryFluxResidual(const CGeometry* geometry, CSolver** solver_container, const CConfig* config,
                            const ScalarFluxOptions& opt, unsigned short val_marker);

  /*!
   * \brief Generic fluid interface (sliding mesh) flux pass, shared by every model. The convective
   *        term is a per-donor weighted average, computed in the same pass that fills the ghost row
   *        of each donor; the diffusive term is computed once per vertex, after the donor loop,
   *        from the ghost state the last donor left behind. This does not fit the
   *        fill-pass-then-BoundaryFluxResidual shape the other boundaries use, so it drives the
   *        kernel directly.
   * \tparam Scheme - Same model the interior loop uses, instantiated with muscl false.
   * \param[in] fillGhostExtras - Functor (iVertex, iPoint) writing the auxiliary ghost fields the
   *            model's diffusion coefficients read, e.g. SST's blending function or the species
   *            mass diffusivities. Called once per vertex, before the diffusive flux.
   */
  template <class Scheme, class GhostFunc>
  void FluidInterfaceFluxResidual(const CGeometry* geometry, CSolver** solver_container, const CConfig* config,
                                  const ScalarFluxOptions& optConv, const ScalarFluxOptions& optVisc,
                                  const GhostFunc& fillGhostExtras);

  /*!
   * \brief Write the outward normal of one vertex into the ghost row and mark the vertex as
   *        contributing a flux.
   * \note Vertex normals point into the domain, the flux convention needs them outward.
   */
  inline void SetGhostGeometry(const CGeometry* geometry, unsigned short val_marker, unsigned long iVertex) {
    for (auto iDim = 0u; iDim < nDim; ++iDim)
      ghostNormal(iVertex, iDim) = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
    ghostSkip[iVertex] = false;
  }

  /*!
   * \brief Write what the diffusion term of a boundary reads beyond the ghost solution: the
   *        coordinate of the interior point reflected about the boundary, and the interior
   *        gradient mirrored into the ghost row.
   * \param[in] iPoint - Interior point of the vertex.
   * \param[in] jPoint - Point the interior one is reflected about, the vertex's normal neighbor.
   */
  inline void SetGhostDiffusionState(const CGeometry* geometry, unsigned long iVertex, unsigned long iPoint,
                                     unsigned long jPoint) {
    su2double Coord_Reflected[MAXNDIM];
    GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(jPoint), geometry->nodes->GetCoord(iPoint),
                                       Coord_Reflected);
    for (auto iDim = 0u; iDim < nDim; ++iDim) ghostCoord(iVertex, iDim) = Coord_Reflected[iDim];

    auto ghostGrad = ghostNodes->GetGradient(iVertex);
    const auto interiorGrad = nodes->GetGradient(iPoint);
    for (auto iVar = 0u; iVar < nVar; ++iVar)
      for (auto iDim = 0u; iDim < nDim; ++iDim) ghostGrad(iVar, iDim) = interiorGrad(iVar, iDim);
  }

  /*!
   * \brief Resolve the compile-time parameters of a scalar flux kernel, the flow indices, the
   *        dimension and the equation count, from the runtime state, and call f with a CTypeTag of
   *        the resulting scheme type: f is a generic lambda,
   *        `[&](auto tag){ using Scheme = typename decltype(tag)::type; ... }`.
   * \tparam Model - Model class template, e.g. CScalarFlux_SST, taking the four parameters of
   *         CUpwScalarBase: value type, flow indices, dimension and equation count.
   * \tparam nVarList - Equation counts to instantiate. Dynamic matches any count, a static one is
   *         taken when it equals the solver's nVar; the counts are tried in the order given.
   * \note NEMO is not one of the index branches: transported scalars are rejected for it at
   *       configuration.
   */
  template <template <class, class, int, size_t> class Model, size_t... nVarList, class F>
  void DispatchScheme(const CConfig* config, F&& f) {
    if (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) {
      DispatchDim<Model, CIncEulerVariable::CIndices<unsigned short>, nVarList...>(std::forward<F>(f));
    } else {
      DispatchDim<Model, CEulerVariable::CIndices<unsigned short>, nVarList...>(std::forward<F>(f));
    }
  }

 private:
  template <template <class, class, int, size_t> class Model, class Indices, size_t... nVarList, class F>
  void DispatchDim(F&& f) {
    if (nDim == 2) {
      DispatchEqnCount<Model, Indices, 2, nVarList...>(std::forward<F>(f));
    } else {
      DispatchEqnCount<Model, Indices, 3, nVarList...>(std::forward<F>(f));
    }
  }

  /*--- Recursion over the candidate equation counts, ending on the empty list. ---*/
  template <template <class, class, int, size_t> class Model, class Indices, int nDim_, class F>
  void DispatchEqnCount(F&&) {
    SU2_MPI::Error("The equation count of this model has no instantiation.", CURRENT_FUNCTION);
  }

  template <template <class, class, int, size_t> class Model, class Indices, int nDim_, size_t nVar_,
            size_t... others, class F>
  void DispatchEqnCount(F&& f) {
    if (nVar_ == Dynamic || nVar_ == nVar) {
      f(CTypeTag<Model<su2double, Indices, nDim_, nVar_>>{});
    } else {
      DispatchEqnCount<Model, Indices, nDim_, others...>(std::forward<F>(f));
    }
  }

  /*!
   * \brief Compute a suitable under-relaxation parameter to limit the change in the solution variables over
   * a nonlinear iteration for stability. Default value 1.0 set in ctor of CScalarVariable.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void ComputeUnderRelaxationFactor(CSolver** solver_container, const CConfig* config) {}

 public:
  /*!
   * \brief Destructor of the class.
   */
  ~CScalarSolver() override;

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] flow_solver - Flow solver the boundary ghost states mirror, null for a solver that
   *            has no flow to read (the heat equation on a solid zone).
   */
  CScalarSolver(CGeometry* geometry, CConfig* config, const CSolver* flow_solver, bool conservative,
                bool bounded_scalar);

  /*!
   * \brief Impose the Far Field boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) override;

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
   * \brief Impose the boundary condition using characteristic recostruction.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Riemann(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker) final;


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
   * \brief Impose the fluid interface boundary condition using transfer data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                          CNumerics *visc_numerics, CConfig *config) override {
    /*--- By default instantiate the generic implementation w/o extra variables, derived solvers can override. ---*/
    BC_Fluid_Interface_impl([](unsigned long){}, geometry, solver_container, conv_numerics, visc_numerics, config);
  }

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
   * \brief This base implementation simply copies the time step of the flow solver.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                    unsigned short iMesh, unsigned long Iteration) override;

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
  void ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) final;

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
                            unsigned short iMesh, unsigned short RunTime_EqSystem) override;

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

  /*!
  * \brief Get the outer state for fluid interface nodes.
  * \param[in] val_marker - marker index
  * \param[in] val_vertex - vertex index
  * \param[in] val_state  - requested state component
  * \param[in] donor_index- index of the donor node to get
  */
  inline su2double GetSlidingState(unsigned short val_marker,
                                   unsigned long val_vertex,
                                   unsigned short val_state,
                                   unsigned long donor_index) const final {
    return SlidingState[val_marker][val_vertex][val_state][donor_index];
  }

  /*!
   * \brief Allocates the final pointer of SlidingState depending on how many donor vertex donate to it. That number is stored in SlidingStateNodes[val_marker][val_vertex].
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   */
  inline void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex) final {
    int iVar;

    for( iVar = 0; iVar < nVar+1; iVar++){
      if( SlidingState[val_marker][val_vertex][iVar] != nullptr )
        delete [] SlidingState[val_marker][val_vertex][iVar];
    }

    for( iVar = 0; iVar < nVar+1; iVar++)
      SlidingState[val_marker][val_vertex][iVar] = new su2double[ GetnSlidingStates(val_marker, val_vertex) ];
  }

  /*!
   * \brief Set the outer state for fluid interface nodes.
   * \param[in] val_marker   - marker index
   * \param[in] val_vertex   - vertex index
   * \param[in] val_state    - requested state component
   * \param[in] donor_index  - index of the donor node to set
   * \param[in] component    - set value
   */
  inline void SetSlidingState(unsigned short val_marker,
                              unsigned long val_vertex,
                              unsigned short val_state,
                              unsigned long donor_index,
                              su2double component) final {
    SlidingState[val_marker][val_vertex][val_state][donor_index] = component;
  }

  /*!
   * \brief Set the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   * \param[in] value - number of outer states
   */
  inline void SetnSlidingStates(unsigned short val_marker,
                                unsigned long val_vertex,
                                int value) final { SlidingStateNodes[val_marker][val_vertex] = value; }

  /*!
   * \brief Get the number of outer state for fluid interface nodes.
   * \param[in] val_marker - marker index
   * \param[in] val_vertex - vertex index
   */
  inline int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex) const final {
    return SlidingStateNodes[val_marker][val_vertex];
  }

};
