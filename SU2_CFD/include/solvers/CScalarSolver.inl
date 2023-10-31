/*!
 * \file CScalarSolver.inl
 * \brief Main subroutines of CScalarSolver class
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

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CScalarSolver.hpp"
#include "../../include/variables/CFlowVariable.hpp"

template <class VariableType>
CScalarSolver<VariableType>::CScalarSolver(CGeometry* geometry, CConfig* config, bool conservative)
    : CSolver(), Conservative(conservative),
      prim_idx(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE,
               config->GetNEMOProblem(), geometry->GetnDim(), config->GetnSpecies()) {
  nMarker = config->GetnMarker_All();

  /*--- Store the number of vertices on each marker for deallocation later ---*/
  nVertex.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = geometry->nVertex[iMarker];

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

#ifdef HAVE_OMP
  /*--- Get the edge coloring, see notes in CEulerSolver's constructor. ---*/
  su2double parallelEff = 1.0;
  const auto& coloring = geometry->GetEdgeColoring(&parallelEff);

  ReducerStrategy = parallelEff < COLORING_EFF_THRESH;

  if (ReducerStrategy && (coloring.getOuterSize() > 1)) geometry->SetNaturalEdgeColoring();

  if (!coloring.empty()) {
    auto groupSize = ReducerStrategy ? 1ul : geometry->GetEdgeColorGroupSize();
    auto nColor = coloring.getOuterSize();
    EdgeColoring.reserve(nColor);

    for (auto iColor = 0ul; iColor < nColor; ++iColor)
      EdgeColoring.emplace_back(coloring.innerIdx(iColor), coloring.getNumNonZeros(iColor), groupSize);
  }

  nPoint = geometry->GetnPoint();
  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), OMP_MAX_SIZE);
#else
  EdgeColoring[0] = DummyGridColor<>(geometry->GetnEdge());
#endif

  /*--- Initialize lower and upper limits for solution clipping. Solvers might overwrite these values. ---*/
  for (unsigned int iVar = 0; iVar < MAXNVAR; iVar++) {
    lowerlimit[iVar] = std::numeric_limits<su2double>::lowest();
    upperlimit[iVar] = std::numeric_limits<su2double>::max();
  }
}

template <class VariableType>
CScalarSolver<VariableType>::~CScalarSolver() {
  delete nodes;
}

template <class VariableType>
void CScalarSolver<VariableType>::CommonPreprocessing(CGeometry *geometry, const CConfig *config, const bool Output) {
  /*--- Define booleans that are solver specific through CConfig's GlobalParams which have to be set in CFluidIteration
   * before calling these solver functions. ---*/
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl = config->GetMUSCL();
  const bool limiter = (config->GetKind_SlopeLimit() != LIMITER::NONE) &&
                       (config->GetInnerIter() <= config->GetLimiterIter());

  /*--- Clear residual and system matrix, not needed for
   * reducer strategy as we write over the entire matrix. ---*/
  if (!ReducerStrategy && !Output) {
    LinSysRes.SetValZero();
    if (implicit) {
      Jacobian.SetValZero();
    } else {
      SU2_OMP_BARRIER
    }
  }

  /*--- Upwind second order reconstruction and gradients ---*/

  if (config->GetReconstructionGradientRequired()) {
    switch(config->GetKind_Gradient_Method_Recon()) {
      case GREEN_GAUSS: SetSolution_Gradient_GG(geometry, config, true); break;
      case LEAST_SQUARES: SetSolution_Gradient_LS(geometry, config, true); break;
      case WEIGHTED_LEAST_SQUARES: SetSolution_Gradient_LS(geometry, config, true); break;
    }
  }

  switch(config->GetKind_Gradient_Method()) {
    case GREEN_GAUSS: SetSolution_Gradient_GG(geometry, config); break;
    case WEIGHTED_LEAST_SQUARES: SetSolution_Gradient_LS(geometry, config); break;
  }

  if (limiter && muscl) SetSolution_Limiter(geometry, config);
}

template <class VariableType>
void CScalarSolver<VariableType>::Upwind_Residual(CGeometry* geometry, CSolver** solver_container,
                                                  CNumerics** numerics_container, CConfig* config,
                                                  unsigned short iMesh) {
  /*--- Define booleans that are solver specific through CConfig's GlobalParams which have to be set in CFluidIteration
   * before calling these solver functions. ---*/
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl = config->GetMUSCL();
  const bool limiter = (config->GetKind_SlopeLimit() != LIMITER::NONE) &&
                       (config->GetInnerIter() <= config->GetLimiterIter());

  /*--- Only reconstruct flow variables if MUSCL is on for flow (requires upwind) and turbulence. ---*/
  const bool musclFlow = config->GetMUSCL_Flow() && muscl && (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND);
  /*--- Only consider flow limiters for cell-based limiters, edge-based would need to be recomputed. ---*/
  const bool limiterFlow =
      (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  const auto& edgeMassFluxes = *(solver_container[FLOW_SOL]->GetEdgeMassFluxes());

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[CONV_TERM + omp_get_thread_num() * MAX_TERMS];

  /*--- Apply scalar advection correction terms for bounded scalar problems ---*/
  const bool bounded_scalar = numerics->GetBoundedScalar();

  /*--- Static arrays of MUSCL-reconstructed flow primitives and turbulence variables (thread safety). ---*/
  su2double solution_i[MAXNVAR] = {0.0}, flowPrimVar_i[MAXNVARFLOW] = {0.0};
  su2double solution_j[MAXNVAR] = {0.0}, flowPrimVar_j[MAXNVARFLOW] = {0.0};

  /*--- For hybrid parallel AD, pause preaccumulation if there is shared reading of
   * variables, otherwise switch to the faster adjoint evaluation mode. ---*/
  bool pausePreacc = false;
  if (ReducerStrategy)
    pausePreacc = AD::PausePreaccumulation();
  else
    AD::StartNoSharedReading();

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring) {
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];

      unsigned short iDim, iVar;

      /*--- Points in edge and normal vectors ---*/

      auto iPoint = geometry->edges->GetNode(iEdge, 0);
      auto jPoint = geometry->edges->GetNode(iEdge, 1);

      numerics->SetNormal(geometry->edges->GetNormal(iEdge));

      /*--- Primitive variables w/o reconstruction ---*/

      const auto V_i = flowNodes->GetPrimitive(iPoint);
      const auto V_j = flowNodes->GetPrimitive(jPoint);
      numerics->SetPrimitive(V_i, V_j);

      /*--- Scalar variables w/o reconstruction ---*/

      const auto Scalar_i = nodes->GetSolution(iPoint);
      const auto Scalar_j = nodes->GetSolution(jPoint);
      numerics->SetScalarVar(Scalar_i, Scalar_j);

      /*--- Grid Movement ---*/

      if (dynamic_grid) numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));

      if (muscl || musclFlow) {
        const su2double *Limiter_i = nullptr, *Limiter_j = nullptr;

        const auto Coord_i = geometry->nodes->GetCoord(iPoint);
        const auto Coord_j = geometry->nodes->GetCoord(jPoint);

        su2double Vector_ij[MAXNDIM] = {0.0};
        for (iDim = 0; iDim < nDim; iDim++) {
          Vector_ij[iDim] = 0.5 * (Coord_j[iDim] - Coord_i[iDim]);
        }

        if (musclFlow && !bounded_scalar) {
          /*--- Reconstruct mean flow primitive variables, note that in bounded scalar mode this is
           * not necessary because the edge mass flux is read directly from the flow solver, instead
           * of being computed from the primitive flow variables. ---*/

          auto Gradient_i = flowNodes->GetGradient_Reconstruction(iPoint);
          auto Gradient_j = flowNodes->GetGradient_Reconstruction(jPoint);

          if (limiterFlow) {
            Limiter_i = flowNodes->GetLimiter_Primitive(iPoint);
            Limiter_j = flowNodes->GetLimiter_Primitive(jPoint);
          }

          for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
            su2double Project_Grad_i = 0.0;
            su2double Project_Grad_j = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Project_Grad_i += Vector_ij[iDim] * Gradient_i[iVar][iDim];
              Project_Grad_j -= Vector_ij[iDim] * Gradient_j[iVar][iDim];
            }
            if (limiterFlow) {
              Project_Grad_i *= Limiter_i[iVar];
              Project_Grad_j *= Limiter_j[iVar];
            }
            flowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
            flowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
          }

          numerics->SetPrimitive(flowPrimVar_i, flowPrimVar_j);
        }

        if (muscl) {
          /*--- Reconstruct scalar variables. ---*/

          auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
          auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

          if (limiter) {
            Limiter_i = nodes->GetLimiter(iPoint);
            Limiter_j = nodes->GetLimiter(jPoint);
          }

          for (iVar = 0; iVar < nVar; iVar++) {
            su2double Project_Grad_i = 0.0, Project_Grad_j = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Project_Grad_i += Vector_ij[iDim] * Gradient_i[iVar][iDim];
              Project_Grad_j -= Vector_ij[iDim] * Gradient_j[iVar][iDim];
            }
            if (limiter) {
              Project_Grad_i *= Limiter_i[iVar];
              Project_Grad_j *= Limiter_j[iVar];
            }
            solution_i[iVar] = Scalar_i[iVar] + Project_Grad_i;
            solution_j[iVar] = Scalar_j[iVar] + Project_Grad_j;
          }

          numerics->SetScalarVar(solution_i, solution_j);
        }
      }

      /*--- Convective flux ---*/
      su2double EdgeMassFlux = 0.0;
      if (bounded_scalar) {
        EdgeMassFlux = edgeMassFluxes[iEdge];
        numerics->SetMassFlux(EdgeMassFlux);
      }

      /*--- Update convective residual value ---*/

      auto residual = numerics->ComputeResidual(config);

      if (ReducerStrategy) {
        EdgeFluxes.SetBlock(iEdge, residual);
        if (implicit) Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
      } else {
        LinSysRes.AddBlock(iPoint, residual);
        LinSysRes.SubtractBlock(jPoint, residual);
        if (implicit) Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
      }

      /*--- Apply convective flux correction to negate the effects of flow divergence in case of incompressible flow.
       * Note that for the bounded scalar model, we explicitly put div(v)=0.
       * If the ReducerStrategy is used, the corrections need to be applied in a loop over nodes
       * to avoid race conditions in accessing nodes shared by edges handled by different threads. ---*/

      if (bounded_scalar && !ReducerStrategy) {
        LinSysRes.AddBlock(iPoint, nodes->GetSolution(iPoint), -EdgeMassFlux);
        LinSysRes.AddBlock(jPoint, nodes->GetSolution(jPoint), EdgeMassFlux);

        if (implicit) {
          Jacobian.AddVal2Diag(iPoint, -EdgeMassFlux);
          Jacobian.AddVal2Diag(jPoint, EdgeMassFlux);
        }
      }

      /*--- Viscous contribution. ---*/

      Viscous_Residual(iEdge, geometry, solver_container,
                       numerics_container[VISC_TERM + omp_get_thread_num() * MAX_TERMS], config);
    }
    END_SU2_OMP_FOR
  }  // end color loop

  /*--- Restore preaccumulation and adjoint evaluation state. ---*/
  AD::ResumePreaccumulation(pausePreacc);
  if (!ReducerStrategy) AD::EndNoSharedReading();

  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    if (implicit) Jacobian.SetDiagonalAsColumnSum();

    /*--- Bounded scalar correction that cannot be applied in the edge loop when using the ReducerStrategy. ---*/
    if (bounded_scalar) {
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
        const auto* solution = nodes->GetSolution(iPoint);
        su2double divergence = 0;

        for (auto iEdge : geometry->nodes->GetEdges(iPoint)) {
          const auto sign = (iPoint == geometry->edges->GetNode(iEdge,0)) ? 1 : -1;
          const su2double EdgeMassFlux = sign * edgeMassFluxes[iEdge];
          divergence += EdgeMassFlux;
          LinSysRes.AddBlock(iPoint, solution, -EdgeMassFlux);
        }
        if (implicit) {
          Jacobian.AddVal2Diag(iPoint, -divergence);
        }
      }
      END_SU2_OMP_FOR
    }
  }
}

template <class VariableType>
void CScalarSolver<VariableType>::SumEdgeFluxes(CGeometry* geometry) {
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    LinSysRes.SetBlock_Zero(iPoint);

    for (auto iEdge : geometry->nodes->GetEdges(iPoint)) {
      if (iPoint == geometry->edges->GetNode(iEdge, 0))
        LinSysRes.AddBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
      else
        LinSysRes.SubtractBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
    }
  }
  END_SU2_OMP_FOR
}

template <class VariableType>
void CScalarSolver<VariableType>::BC_Periodic(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                              CConfig* config) {
  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated corectly during the communications. For implicit calculations
   the Jacobians and linear system are also correctly adjusted here. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic() / 2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }
}

template <class VariableType>
void CScalarSolver<VariableType>::BC_Far_Field(CGeometry* geometry, CSolver** solver_container,
                                               CNumerics* conv_numerics, CNumerics*, CConfig *config,
                                               unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the infinity ---*/

      auto V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Solution_Inf);

      /*--- Set Normal (it is necessary to change the sign) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      if (conv_numerics->GetBoundedScalar()) {
        const su2double* velocity = &V_infty[prim_idx.Velocity()];
        const su2double density = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        conv_numerics->SetMassFlux(BoundedScalarBCFlux(iPoint, implicit, density, velocity, Normal));
      }

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR
}

template <class VariableType>
void CScalarSolver<VariableType>::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                               unsigned short iMesh, unsigned long Iteration) {
  const auto flowNodes = solver_container[FLOW_SOL]->GetNodes();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    su2double dt = nodes->GetLocalCFL(iPoint) / flowNodes->GetLocalCFL(iPoint) * flowNodes->GetDelta_Time(iPoint);
    nodes->SetDelta_Time(iPoint, dt);
  }
  END_SU2_OMP_FOR
}

template <class VariableType>
void CScalarSolver<VariableType>::PrepareImplicitIteration(CGeometry* geometry, CSolver** solver_container,
                                                           CConfig* config) {
  /*--- Set shared residual variables to 0 and declare
   *    local ones for current thread to work on. ---*/

  SetResToZero();

  su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
  unsigned long idxMax[MAXNVAR] = {0};

  /*--- Build implicit system ---*/

  SU2_OMP_FOR_(schedule(static, omp_chunk_size) SU2_NOWAIT)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    /*--- Modify matrix diagonal to improve diagonal dominance. ---*/
    const su2double dt = nodes->GetDelta_Time(iPoint);

    if (dt != 0.0) {
      su2double Vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
      Jacobian.AddVal2Diag(iPoint, Vol / dt);
    } else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      LinSysRes.SetBlock_Zero(iPoint);
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      unsigned long total_index = iPoint * nVar + iVar;
      LinSysRes[total_index] = -LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;

      /*--- "Add" residual at (iPoint,iVar) to local residual variables. ---*/
      ResidualReductions_PerThread(iPoint, iVar, LinSysRes[total_index], resRMS, resMax, idxMax);
    }
  }
  END_SU2_OMP_FOR

  /*--- "Add" residuals from all threads to global residual variables. ---*/
  ResidualReductions_FromAllThreads(geometry, config, resRMS, resMax, idxMax);
}

template <class VariableType>
void CScalarSolver<VariableType>::CompleteImplicitIteration(CGeometry* geometry, CSolver** solver_container,
                                                            CConfig* config) {
  const bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  ComputeUnderRelaxationFactor(config);

  /*--- Update solution (system written in terms of increments) ---*/

  if (!adjoint) {
    /*--- Update the scalar solution. For transport equations, where Solution is not equivalent with the transported
     * quantity, multiply the respective factor.  ---*/
    if (Conservative) {
      const auto flowNodes = solver_container[FLOW_SOL]->GetNodes();

      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        /*--- Multiply the Solution var with density to get the conservative transported quantity, if necessary. ---*/
        const su2double density = flowNodes->GetDensity(iPoint);
        const su2double density_old = compressible ? flowNodes->GetSolution_Old(iPoint, 0) : density;

        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
          nodes->AddClippedSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint) * LinSysSol(iPoint, iVar),
                                    lowerlimit[iVar], upperlimit[iVar], density, density_old);
        }
      }
      END_SU2_OMP_FOR
    } else {
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
          nodes->AddClippedSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint) * LinSysSol(iPoint, iVar),
                                    lowerlimit[iVar], upperlimit[iVar]);
        }
      }
      END_SU2_OMP_FOR
    }
  }

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic() / 2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
}

template <class VariableType>
void CScalarSolver<VariableType>::ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container,
                                                          CConfig* config) {
  PrepareImplicitIteration(geometry, solver_container, config);

  /*--- Solve or smooth the linear system. ---*/

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned long iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }
  END_SU2_OMP_FOR

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    SetIterLinSolver(iter);
    SetResLinSolver(System.GetResidual());
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  CompleteImplicitIteration(geometry, solver_container, config);
}

template <class VariableType>
void CScalarSolver<VariableType>::ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container,
                                                          CConfig* config) {
  /*--- Local residual variables for current thread ---*/
  su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
  unsigned long idxMax[MAXNVAR] = {0};

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    const su2double dt = nodes->GetDelta_Time(iPoint);
    const su2double Vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);

    for (auto iVar = 0u; iVar < nVar; iVar++) {
      /*--- "Add" residual at (iPoint,iVar) to local residual variables. ---*/
      ResidualReductions_PerThread(iPoint, iVar, LinSysRes(iPoint, iVar), resRMS, resMax, idxMax);
      /*--- Explicit Euler step: ---*/
      LinSysSol(iPoint, iVar) = -dt / Vol * LinSysRes(iPoint, iVar);
    }
  }
  END_SU2_OMP_FOR

  /*--- "Add" residuals from all threads to global residual variables. ---*/
  ResidualReductions_FromAllThreads(geometry, config, resRMS, resMax, idxMax);

  /*--- Use LinSysSol for solution update. ---*/
  CompleteImplicitIteration(geometry, solver_container, config);
}

template <class VariableType>
void CScalarSolver<VariableType>::SetResidual_DualTime(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                                       unsigned short iRKStep, unsigned short iMesh,
                                                       unsigned short RunTime_EqSystem) {
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool first_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  const bool second_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);

  /*--- Flow solution, needed to get density. ---*/

  CVariable* flowNodes = solver_container[FLOW_SOL]->GetNodes();

  /*--- Store the physical time step ---*/

  const su2double TimeStep = config->GetDelta_UnstTimeND();

  /*--- Local variables ---*/

  unsigned short iVar, iMarker, iDim, iNeigh;
  unsigned long iPoint, jPoint, iVertex, iEdge;

  su2double *U_time_nM1 = nullptr, *U_time_n = nullptr, *U_time_nP1 = nullptr;
  /*--- For non-Conservative scalars (e.g. SA), multiply the Primitives with 1.0 instead of Density. ---*/
  su2double Density_nM1 = 1.0, Density_n = 1.0, Density_nP1 = 1.0;
  su2double Volume_nM1, Volume_nP1;
  const su2double *Normal = nullptr, *GridVel_i = nullptr, *GridVel_j = nullptr;
  su2double Residual_GCL;

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {
    /*--- Loop over all nodes (excluding halos) ---*/

    AD::StartNoSharedReading();

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (Conservative) {
        if (incompressible) {
          /*--- This is temporary and only valid for constant-density problems:
          density could also be temperature dependent, but as it is not a part
          of the solution vector it's neither stored for previous time steps
          nor updated with the solution at the end of each iteration. */
          Density_nM1 = flowNodes->GetDensity(iPoint);
          Density_n = flowNodes->GetDensity(iPoint);
          Density_nP1 = flowNodes->GetDensity(iPoint);
        } else {
          Density_nM1 = flowNodes->GetSolution_time_n1(iPoint)[0];
          Density_n = flowNodes->GetSolution_time_n(iPoint, 0);
          Density_nP1 = flowNodes->GetSolution(iPoint, 0);
        }
      }

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (first_order)
          LinSysRes(iPoint, iVar) +=
              (Density_nP1 * U_time_nP1[iVar] - Density_n * U_time_n[iVar]) * Volume_nP1 / TimeStep;
        if (second_order)
          LinSysRes(iPoint, iVar) += (3.0 * Density_nP1 * U_time_nP1[iVar] - 4.0 * Density_n * U_time_n[iVar] +
                                      1.0 * Density_nM1 * U_time_nM1[iVar]) *
                                     Volume_nP1 / (2.0 * TimeStep);
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1 / TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1 * 3.0) / (2.0 * TimeStep));
      }
    }
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();

  } else {
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      GridVel_i = geometry->nodes->GetGridVel(iPoint);
      U_time_n = nodes->GetSolution_time_n(iPoint);

      if (Conservative) {
        if (incompressible)
          Density_n = flowNodes->GetDensity(iPoint);  // Temporary fix
        else
          Density_n = flowNodes->GetSolution_time_n(iPoint, 0);
      }

      for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
        iEdge = geometry->nodes->GetEdge(iPoint, iNeigh);
        Normal = geometry->edges->GetNormal(iEdge);

        jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
        GridVel_j = geometry->nodes->GetGridVel(jPoint);

        /*--- Determine whether to consider the normal outward or inward. ---*/
        su2double dir = (iPoint < jPoint) ? 0.5 : -0.5;

        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Residual_GCL += dir * (GridVel_i[iDim] + GridVel_j[iDim]) * Normal[iDim];

        Residual_GCL *= Density_n;

        for (iVar = 0; iVar < nVar; iVar++) LinSysRes(iPoint, iVar) += U_time_n[iVar] * Residual_GCL;
      }
    }
    END_SU2_OMP_FOR

    /*--- Loop over the boundary edges ---*/

    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
        SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          /*--- Get the index for node i plus the boundary face normal ---*/

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          /*--- Grid velocities stored at boundary node i ---*/

          GridVel_i = geometry->nodes->GetGridVel(iPoint);

          /*--- Compute the GCL term by dotting the grid velocity with the face
           normal. The normal is negated to match the boundary convention. ---*/

          Residual_GCL = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) Residual_GCL -= 0.5 * (GridVel_i[iDim] + GridVel_i[iDim]) * Normal[iDim];

          /*--- Compute the GCL component of the source term for node i ---*/

          U_time_n = nodes->GetSolution_time_n(iPoint);

          /*--- Multiply by density at node i for the SST model ---*/

          if (Conservative) {
            if (incompressible)
              Density_n = flowNodes->GetDensity(iPoint);  // Temporary fix
            else
              Density_n = flowNodes->GetSolution_time_n(iPoint, 0);
          }

          for (iVar = 0; iVar < nVar; iVar++) LinSysRes(iPoint, iVar) += Density_n * U_time_n[iVar] * Residual_GCL;
        }
        END_SU2_OMP_FOR
      }
    }

    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/

    AD::EndNoSharedReading();

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/

      Volume_nM1 = geometry->nodes->GetVolume_nM1(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/

      if (Conservative) {
        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        if (incompressible) {
          /*--- This is temporary and only valid for constant-density problems:
          density could also be temperature dependent, but as it is not a part
          of the solution vector it's neither stored for previous time steps
          nor updated with the solution at the end of each iteration. */
          Density_nM1 = flowNodes->GetDensity(iPoint);
          Density_n = flowNodes->GetDensity(iPoint);
          Density_nP1 = flowNodes->GetDensity(iPoint);
        } else {
          Density_nM1 = flowNodes->GetSolution_time_n1(iPoint)[0];
          Density_n = flowNodes->GetSolution_time_n(iPoint, 0);
          Density_nP1 = flowNodes->GetSolution(iPoint, 0);
        }
      }

      for (iVar = 0; iVar < nVar; iVar++) {
        if (first_order)
          LinSysRes(iPoint, iVar) +=
              (Density_nP1 * U_time_nP1[iVar] - Density_n * U_time_n[iVar]) * (Volume_nP1 / TimeStep);
        if (second_order)
          LinSysRes(iPoint, iVar) +=
              (Density_nP1 * U_time_nP1[iVar] - Density_n * U_time_n[iVar]) * (3.0 * Volume_nP1 / (2.0 * TimeStep)) +
              (Density_nM1 * U_time_nM1[iVar] - Density_n * U_time_n[iVar]) * (Volume_nM1 / (2.0 * TimeStep));
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1 / TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1 * 3.0) / (2.0 * TimeStep));
      }
    }
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();

  }  // end dynamic grid
}
