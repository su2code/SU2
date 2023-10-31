/*!
 * \file CScalarSolver.hpp
 * \brief Headers of the CScalarSolver class
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
#include "../variables/CScalarVariable.hpp"
#include "../variables/CFlowVariable.hpp"
#include "../variables/CPrimitiveIndices.hpp"
#include "CSolver.hpp"

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

  /*!
   * \brief The highest level in the variable hierarchy this solver can safely use.
   */
  VariableType* nodes = nullptr;

  /*!
   * \brief Return nodes to allow CSolver::base_nodes to be set.
   */
  inline CVariable* GetBaseClassPointerToNodes() final { return nodes; }

  /*!
   * \brief Compute the viscous flux for the scalar equation at a particular edge.
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
    CFlowVariable* flowNodes = solver_container[FLOW_SOL] ?
        su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes()) : nullptr;

    /*--- Points in edge ---*/

    auto iPoint = geometry->edges->GetNode(iEdge, 0);
    auto jPoint = geometry->edges->GetNode(iEdge, 1);

    /*--- Points coordinates, and normal vector ---*/

    numerics->SetCoord(geometry->nodes->GetCoord(iPoint), geometry->nodes->GetCoord(jPoint));
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Conservative variables w/o reconstruction ---*/

    if (flowNodes) {
      numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), flowNodes->GetPrimitive(jPoint));
    }

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
  void SumEdgeFluxes(CGeometry* geometry);

 private:
  /*!
   * \brief Compute the viscous flux for the scalar equation at a particular edge.
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
  void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics_container,
                       CConfig* config, unsigned short iMesh) override;

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
                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) final;

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
