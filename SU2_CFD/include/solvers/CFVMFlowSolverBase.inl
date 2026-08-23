/*!
 * \file CFVMFlowSolverBase.inl
 * \brief Base class template for all FVM flow solvers.
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

#include "../gradients/computeGradientsGreenGauss.hpp"
#include "../gradients/computeGradientsGreenGaussLimited.hpp"
#include "../gradients/computeGradientsLeastSquares.hpp"
#include "../limiters/computeLimiters.hpp"
#include "../numerics_simd/CNumericsSIMD.hpp"
#include "CFVMFlowSolverBase.hpp"

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::AeroCoeffsArray::allocate(int size) {
  _size = size;
  CD = new su2double[size];
  CL = new su2double[size];
  CSF = new su2double[size];
  CEff = new su2double[size];
  CFx = new su2double[size];
  CFy = new su2double[size];
  CFz = new su2double[size];
  CMx = new su2double[size];
  CMy = new su2double[size];
  CMz = new su2double[size];
  CoPx = new su2double[size];
  CoPy = new su2double[size];
  CoPz = new su2double[size];
  CT = new su2double[size];
  CQ = new su2double[size];
  CMerit = new su2double[size];
  setZero();
}

template <class V, ENUM_REGIME R>
CFVMFlowSolverBase<V, R>::AeroCoeffsArray::~AeroCoeffsArray() {
  delete[] CD;
  delete[] CL;
  delete[] CSF;
  delete[] CEff;
  delete[] CFx;
  delete[] CFy;
  delete[] CFz;
  delete[] CMx;
  delete[] CMy;
  delete[] CMz;
  delete[] CoPx;
  delete[] CoPy;
  delete[] CoPz;
  delete[] CT;
  delete[] CQ;
  delete[] CMerit;
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::Allocate(const CConfig& config) {
  SU2_ZONE_SCOPED

  /*--- Define some auxiliar vector related with the residual ---*/

  Residual_RMS.resize(nVar,0.0);
  Residual_Max.resize(nVar,0.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/

  if ((config.GetKind_ConvNumScheme_Flow() == SPACE_CENTERED && MGLevel == MESH_0) ||
      config.GetKind_Upwind_Flow() == UPWIND::MSW) {
    iPoint_UndLapl.resize(nPointDomain);
    jPoint_UndLapl.resize(nPointDomain);
  }

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- LinSysSol will always be init to 0. ---*/
  System.SetxIsZero(true);

  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/

  AllocVectorOfMatrices(nVertex, nPrimVar, CharacPrimVar);

  /*--- Store the value of the Total Pressure at the inlet BC ---*/

  AllocVectorOfVectors(nVertex, Inlet_Ttotal);

  /*--- Store the value of the Total Temperature at the inlet BC ---*/

  AllocVectorOfVectors(nVertex, Inlet_Ptotal);

  /*--- Store the value of the Flow direction at the inlet BC ---*/

  AllocVectorOfMatrices(nVertex, nDim, Inlet_FlowDir);

  /*--- Force definition and coefficient arrays for all of the markers ---*/

  AllocVectorOfVectors(nVertex, CPressure);
  AllocVectorOfVectors(nVertex, CPressureTarget);

  /*--- Non dimensional aerodynamic coefficients ---*/

  InvCoeff.allocate(nMarker);
  MntCoeff.allocate(nMarker);
  ViscCoeff.allocate(nMarker);
  SurfaceInvCoeff.allocate(config.GetnMarker_Monitoring());
  SurfaceMntCoeff.allocate(config.GetnMarker_Monitoring());
  SurfaceViscCoeff.allocate(config.GetnMarker_Monitoring());
  SurfaceCoeff.allocate(config.GetnMarker_Monitoring());

  /*--- Heat flux coefficients. ---*/

  HF_Visc.resize(nMarker,0.0);
  MaxHF_Visc.resize(nMarker,0.0);

  Surface_HF_Visc.resize(config.GetnMarker_Monitoring());
  Surface_MaxHF_Visc.resize(config.GetnMarker_Monitoring());

  /*--- Supersonic coefficients ---*/

  CNearFieldOF_Inv.resize(nMarker,0.0);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config.GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*--- Heat flux in all the markers ---*/

  AllocVectorOfVectors(nVertex, HeatFlux);
  AllocVectorOfVectors(nVertex, HeatFluxTarget);

  /*--- Y plus in all the markers ---*/

  AllocVectorOfVectors(nVertex, YPlus);

  /*--- U Tau in all the markers ---*/

  AllocVectorOfVectors(nVertex, UTau);

  /*--- wall eddy viscosity in all the markers ---*/

  AllocVectorOfVectors(nVertex, EddyViscWall);

  /*--- Skin friction in all the markers ---*/

  AllocVectorOfMatrices(nVertex, nDim, CSkinFriction);

  /*--- Wall Shear Stress in all the markers ---*/

  AllocVectorOfVectors(nVertex, WallShearStress);

  /*--- Store the values of the temperature and the heat flux density at the boundaries,
   used for coupling with a solid donor cell ---*/
  constexpr auto nHeatConjugateVar = 4u;
  AllocVectorOfMatrices(nVertex, nHeatConjugateVar, HeatConjugateVar, config.GetTemperature_FreeStreamND());

  if (MGLevel == MESH_0) {
    auto nSolidVertex = nVertex;
    for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
      if (!config.GetSolid_Wall(iMarker))
        nSolidVertex[iMarker] = 0;

    AllocVectorOfMatrices(nSolidVertex, nDim, VertexTraction);

    if (config.GetDiscrete_Adjoint()) AllocVectorOfMatrices(nSolidVertex, nDim, VertexTractionAdjoint);
  }

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (config.GetMultizone_Residual()) {
    Residual_BGS.resize(nVar,1.0);
    Residual_Max_BGS.resize(nVar,1.0);
    Point_Max_BGS.resize(nVar,0);
    Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
  }
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::AllocateTerribleLegacyTemporaryVariables() {
  SU2_ZONE_SCOPED

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual = new su2double[nVar]();
  Res_Conv = new su2double[nVar]();
  Res_Visc = new su2double[nVar]();
  Res_Sour = new su2double[nVar]();

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution = new su2double[nVar]();
  Solution_i = new su2double[nVar]();
  Solution_j = new su2double[nVar]();

  /*--- Define some auxiliary vectors related to the geometry ---*/

  Vector = new su2double[nDim]();
  Vector_i = new su2double[nDim]();
  Vector_j = new su2double[nDim]();

  /*--- Jacobian temporaries. ---*/

  Jacobian_i = new su2double*[nVar];
  Jacobian_j = new su2double*[nVar];
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double[nVar];
    Jacobian_j[iVar] = new su2double[nVar];
  }
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::CommunicateInitialState(CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  /*--- Define solver parameters needed for execution of destructor ---*/

  space_centered = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  euler_implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  least_squares = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);

  /*--- Communicate and store volume and the number of neighbors for
   any dual CVs that lie on on periodic markers. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic() / 2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
  }
  SetImplicitPeriodic(euler_implicit);
  if (MGLevel == MESH_0) SetRotatePeriodic(true);

  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  /*--- Store the initial CFL number for all grid points. ---*/

  const auto CFL = config->GetCFL(MGLevel);
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::HybridParallelInitialization(const CConfig& config, CGeometry& geometry) {
  SU2_ZONE_SCOPED

#ifdef HAVE_OMP
  /*--- Get the edge coloring. If the expected parallel efficiency becomes too low setup the
   *    reducer strategy. Where one loop is performed over edges followed by a point loop to
   *    sum the fluxes for each cell and set the diagonal of the system matrix. ---*/

  su2double parallelEff = 1.0;

#ifdef CODI_REVERSE_TYPE
  /*--- For the discrete adjoint, the reducer strategy is costly. Prefer coloring, possibly with reduced edge color
   *    group size. Find the maximum edge color group size that yields an efficient coloring. Also, allow larger numbers
   *    of colors. ---*/
  const bool relax =  config.GetEdgeColoringRelaxDiscAdj();
  const auto& coloring = geometry.GetEdgeColoring(&parallelEff, relax);
#else
  const auto& coloring = geometry.GetEdgeColoring(&parallelEff);
#endif

  /*--- The decision to use the strategy is local to each rank. ---*/
  ReducerStrategy = parallelEff < COLORING_EFF_THRESH;

  /*--- When using the reducer force a single color to reduce the color loop overhead. ---*/
  if (ReducerStrategy && (coloring.getOuterSize() > 1)) geometry.SetNaturalEdgeColoring();

  if (!coloring.empty()) {
    /*--- If the reducer strategy is used we are not constrained by group
     *    size as we have no other edge loops in the Euler/NS solvers. ---*/
    auto groupSize = static_cast<su2uint>(ReducerStrategy ? 1ul : geometry.GetEdgeColorGroupSize());
    auto nColor = coloring.getOuterSize();
    EdgeColoring.reserve(nColor);

    for (auto iColor = 0ul; iColor < nColor; ++iColor)
      EdgeColoring.emplace_back(coloring.innerIdx(iColor), coloring.getNumNonZeros(iColor), groupSize);
  }

  /*--- If the reducer strategy is not being forced (by EDGE_COLORING_GROUP_SIZE=0) print some messages. ---*/
  if (config.GetEdgeColoringGroupSize() != 1 << 30) {
    su2double minEff = 1.0;
    SU2_MPI::Reduce(&parallelEff, &minEff, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, SU2_MPI::GetComm());

    int tmp = ReducerStrategy, numRanksUsingReducer = 0;
    SU2_MPI::Reduce(&tmp, &numRanksUsingReducer, 1, MPI_INT, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

    if (minEff < COLORING_EFF_THRESH) {
      cout << "WARNING: On " << numRanksUsingReducer << " MPI ranks the coloring efficiency was less than "
           << COLORING_EFF_THRESH << " (min value was " << minEff << ").\n"
           << "         Those ranks will now use a fallback strategy, better performance may be possible\n"
           << "         with a different value of config option EDGE_COLORING_GROUP_SIZE (default 512)."
#ifdef HAVE_OPDI
           << "\n         The memory usage of the discrete adjoint solver is higher when using the fallback."
#endif
           << endl;
    } else {
      if (SU2_MPI::GetRank() == MASTER_NODE) {
        cout << "All ranks use edge coloring." << endl;
      }
    }

    const su2double coloredParallelEff = ReducerStrategy ? 1.0 : parallelEff;
    su2double minColoredParallelEff = 1.0;
    SU2_MPI::Reduce(&coloredParallelEff, &minColoredParallelEff, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, SU2_MPI::GetComm());

    const unsigned long coloredNumColors = ReducerStrategy ? 0 : coloring.getOuterSize();
    unsigned long maxColoredNumColors = 0;
    SU2_MPI::Reduce(&coloredNumColors, &maxColoredNumColors, 1, MPI_UNSIGNED_LONG, MPI_MAX, MASTER_NODE, SU2_MPI::GetComm());

    const unsigned long coloredEdgeColorGroupSize = ReducerStrategy ? 1 << 30 : geometry.GetEdgeColorGroupSize();
    unsigned long minColoredEdgeColorGroupSize = 1 << 30;
    SU2_MPI::Reduce(&coloredEdgeColorGroupSize, &minColoredEdgeColorGroupSize, 1, MPI_UNSIGNED_LONG, MPI_MIN, MASTER_NODE, SU2_MPI::GetComm());

    if (SU2_MPI::GetRank() == MASTER_NODE && numRanksUsingReducer != SU2_MPI::GetSize()) {
      cout << "Among the ranks that use edge coloring,\n"
           << "      the minimum efficiency is " << minColoredParallelEff << ",\n"
           << "      the maximum number of colors is " << maxColoredNumColors << ",\n"
           << "      the minimum edge color group size is " << minColoredEdgeColorGroupSize << "." << endl;
    }
  }

  if (ReducerStrategy) EdgeFluxes.Initialize(geometry.GetnEdge(), geometry.GetnEdge(), nVar, nullptr);

  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), OMP_MAX_SIZE);
#else
  EdgeColoring[0] = DummyGridColor<>(geometry.GetnEdge());
#endif
}

template <class V, ENUM_REGIME R>
CFVMFlowSolverBase<V, R>::~CFVMFlowSolverBase() {
  SU2_ZONE_SCOPED

  for (auto& mat : SlidingState) {
    for (auto ptr : mat) delete [] ptr;
  }

  delete nodes;
  delete edgeNumerics;
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetPrimitive_Gradient_GG(CGeometry* geometry, const CConfig* config,
                                                        bool reconstruction) {
  SU2_ZONE_SCOPED

  const auto& primitives = nodes->GetPrimitive();
  auto& gradient = reconstruction ? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
  const auto comm = reconstruction? MPI_QUANTITIES::PRIMITIVE_GRAD_REC : MPI_QUANTITIES::PRIMITIVE_GRADIENT;
  const auto commPer = reconstruction? PERIODIC_PRIM_GG_R : PERIODIC_PRIM_GG;

  computeGradientsGreenGauss(this, comm, commPer, *geometry, *config, primitives, 0, nPrimVarGrad, prim_idx.Velocity(), gradient);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetPrimitive_Gradient_GG_Limited(CGeometry* geometry, const CConfig* config,
                                                                bool reconstruction) {
  SU2_ZONE_SCOPED

  const auto& primitives = nodes->GetPrimitive();
  auto& gradient = reconstruction ? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
  const auto comm = reconstruction? MPI_QUANTITIES::PRIMITIVE_GRAD_REC : MPI_QUANTITIES::PRIMITIVE_GRADIENT;
  const auto commPer = reconstruction? PERIODIC_PRIM_GG_R : PERIODIC_PRIM_GG;

  computeGradientsGreenGaussLimited(config->GetKind_SlopeLimit_Flow(), this, comm, commPer, *geometry, *config,
                                    config->GetMUSCL_Kappa_Flow(), primitives, 0, nPrimVarGrad, prim_idx.Velocity(),
                                    gradient);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetPrimitive_Gradient(CGeometry* geometry, const CConfig* config,
                                                     bool reconstruction) {
  const auto kind = reconstruction ? config->GetKind_Gradient_Method_Recon() : config->GetKind_Gradient_Method();
  const auto limiterKind = config->GetKind_SlopeLimit_Flow();
  const bool limited = reconstruction ? config->GetLimitedGradientRecon(limiterKind)
                                      : config->GetLimitedGradient(limiterKind);
  switch (kind) {
    case GREEN_GAUSS_LIMITED:
      if (limited) { SetPrimitive_Gradient_GG_Limited(geometry, config, reconstruction); break; }
      /*--- Fall back to plain Green-Gauss, e.g. if no limiter is used. ---*/
      SetPrimitive_Gradient_GG(geometry, config, reconstruction); break;
    case GREEN_GAUSS:
      SetPrimitive_Gradient_GG(geometry, config, reconstruction); break;
    case LEAST_SQUARES:
    case WEIGHTED_LEAST_SQUARES:
      SetPrimitive_Gradient_LS(geometry, config, reconstruction); break;
    default: break;
  }
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetPrimitive_Gradient_LS(CGeometry* geometry, const CConfig* config,
                                                        bool reconstruction) {
  SU2_ZONE_SCOPED

  /*--- Set a flag for unweighted or weighted least-squares. ---*/
  bool weighted;
  PERIODIC_QUANTITIES commPer;

  if (reconstruction) {
    weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
    commPer = weighted? PERIODIC_PRIM_LS_R : PERIODIC_PRIM_ULS_R;
  }
  else {
    weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);
    commPer = weighted? PERIODIC_PRIM_LS : PERIODIC_PRIM_ULS;
  }

  const auto& primitives = nodes->GetPrimitive();
  auto& rmatrix = nodes->GetRmatrix();
  auto& gradient = reconstruction ? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
  const auto comm = reconstruction? MPI_QUANTITIES::PRIMITIVE_GRAD_REC : MPI_QUANTITIES::PRIMITIVE_GRADIENT;

  computeGradientsLeastSquares(this, comm, commPer, *geometry, *config, weighted,
                               primitives, 0, nPrimVarGrad, prim_idx.Velocity(), gradient, rmatrix);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetPrimitive_Limiter(CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  const auto kindLimiter = config->GetKind_SlopeLimit_Flow();
  const auto umusclKappa = config->GetMUSCL_Kappa_Flow();
  const auto& primitives = nodes->GetPrimitive();
  const auto& gradient = nodes->GetGradient_Reconstruction();
  auto& primMin = nodes->GetSolution_Min();
  auto& primMax = nodes->GetSolution_Max();
  auto& limiter = nodes->GetLimiter_Primitive();

  computeLimiters(kindLimiter, this, MPI_QUANTITIES::PRIMITIVE_LIMITER, PERIODIC_LIM_PRIM_1, PERIODIC_LIM_PRIM_2,
                  *geometry, *config, 0, nPrimVarGrad, umusclKappa, primitives, gradient, primMin, primMax, limiter);
}

template <class V, ENUM_REGIME R>
CNumerics::ResidualType<> CFVMFlowSolverBase<V, R>::Viscous_Residual_impl(unsigned long iEdge, CGeometry *geometry,
                                                                          CSolver **solver_container,
                                                                          CNumerics *numerics, CConfig *config) {
  SU2_ZONE_SCOPED

  const bool tkeNeeded = (config->GetKind_Turb_Model() == TURB_MODEL::SST);
  const bool backscatter = config->GetSBSParam().StochasticBackscatter;
  const bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR) ||
                         (config->GetKind_FluidModel() == IDEAL_GAS);

  CVariable* turbNodes = nullptr;
  if (tkeNeeded || backscatter) turbNodes = solver_container[TURB_SOL]->GetNodes();

  /*--- Points, coordinates and normal vector in edge ---*/

  auto iPoint = geometry->edges->GetNode(iEdge,0);
  auto jPoint = geometry->edges->GetNode(iEdge,1);

  numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                     geometry->nodes->GetCoord(jPoint));

  numerics->SetNormal(geometry->edges->GetNormal(iEdge));

  /*--- Primitive and secondary variables. ---*/

  numerics->SetPrimitive(nodes->GetPrimitive(iPoint),
                         nodes->GetPrimitive(jPoint));
  if (!ideal_gas) {
    numerics->SetSecondary(nodes->GetSecondary(iPoint),
                           nodes->GetSecondary(jPoint));
  }
  /*--- Gradients. ---*/

  numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                               nodes->GetGradient_Primitive(jPoint));

  /*--- Turbulent kinetic energy. ---*/

  if (tkeNeeded)
    numerics->SetTurbKineticEnergy(turbNodes->GetSolution(iPoint,0),
                                   turbNodes->GetSolution(jPoint,0));

  /*--- Stochastic variables from Langevin equations (Stochastic Backscatter Model). ---*/

  if (backscatter) {
    if (config->GetSBSParam().SBS_Ctau > 0.0) {
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        numerics->SetStochVar(iDim, turbNodes->GetSolution(iPoint, iDim+1),
                                    turbNodes->GetSolution(jPoint, iDim+1));
    } else {
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        numerics->SetStochVar(iDim, turbNodes->GetLangevinSourceTerms(iPoint, iDim),
                                    turbNodes->GetLangevinSourceTerms(jPoint, iDim));
    }
    su2double DES_length_i = max(turbNodes->GetDES_LengthScale(iPoint), 1e-10);
    su2double DES_length_j = max(turbNodes->GetDES_LengthScale(jPoint), 1e-10);
    su2double lesMode_i = turbNodes->GetLES_Mode(iPoint);
    su2double lesMode_j = turbNodes->GetLES_Mode(jPoint);
    numerics->SetDistance(DES_length_i, DES_length_j);
    numerics->SetLES_Mode(lesMode_i, lesMode_j);
  }

  /*--- Wall shear stress values (wall functions) ---*/

  numerics->SetTau_Wall(nodes->GetTau_Wall(iPoint),
                        nodes->GetTau_Wall(jPoint));

  /*--- Compute and update residual ---*/

  auto residual = numerics->ComputeResidual(config);

  if (ReducerStrategy) {
    EdgeFluxes.SubtractBlock(iEdge, residual);
  }
  else {
    LinSysRes.SubtractBlock(iPoint, residual);
    LinSysRes.AddBlock(jPoint, residual);
  }

  /*--- The Jacobians are applied by the caller, fused with the convective contribution. ---*/
  return residual;
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::ComputeVerificationError(CGeometry* geometry, CConfig* config) {
  SU2_ZONE_SCOPED

  /*--- The errors only need to be computed on the finest grid. ---*/
  if (MGLevel != MESH_0) return;

  /*--- If this is a verification case, we can compute the global
   error metrics by using the difference between the local error
   and the known solution at each DOF. This is then collected into
   RMS (L2) and maximum (Linf) global error norms. From these
   global measures, one can compute the order of accuracy. ---*/

  bool write_heads =
      ((((config->GetInnerIter() % (config->GetScreen_Wrt_Freq(2) * 40)) == 0) && (config->GetInnerIter() != 0)) ||
       (config->GetInnerIter() == 1));
  if (!write_heads) return;

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {

  /*--- Check if there actually is an exact solution for this
        verification case, if computed at all. ---*/
  if (VerificationSolution && VerificationSolution->ExactSolutionKnown()) {
    /*--- Get the physical time if necessary. ---*/
    su2double time = 0.0;
    if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

    /*--- Reset the global error measures to zero. ---*/
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      VerificationSolution->SetError_RMS(iVar, 0.0);
      VerificationSolution->SetError_Max(iVar, 0.0, 0);
    }

    /*--- Loop over all owned points. ---*/
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      /* Set the pointers to the coordinates and solution of this DOF. */
      const su2double* coor = geometry->nodes->GetCoord(iPoint);
      su2double* solDOF = nodes->GetSolution(iPoint);

      /* Get local error from the verification solution class. */
      vector<su2double> error(nVar, 0.0);
      VerificationSolution->GetLocalError(coor, time, solDOF, error.data());

      /* Increment the global error measures */
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        VerificationSolution->AddError_RMS(iVar, error[iVar] * error[iVar]);
        VerificationSolution->AddError_Max(iVar, fabs(error[iVar]), geometry->nodes->GetGlobalIndex(iPoint),
                                           geometry->nodes->GetCoord(iPoint));
      }
    }

    /* Finalize the calculation of the global error measures. */
    VerificationSolution->SetVerificationError(geometry->GetGlobal_nPointDomain(), config);

    /*--- Screen output of the error metrics. This can be improved
     once the new output classes are in place. ---*/

    PrintVerificationError(config);
  }

  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::CompleteImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {
  SU2_ZONE_SCOPED

  if constexpr (R == ENUM_REGIME::COMPRESSIBLE) ComputeUnderRelaxationFactor(config);

  /*--- Update solution with under-relaxation and communicate it. ---*/

  if (!config->GetContinuous_Adjoint()) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        nodes->AddSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint) * LinSysSol(iPoint, iVar));
      }
    }
    END_SU2_OMP_FOR
  }

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  /*--- For verification cases, compute the global error metrics. ---*/
  ComputeVerificationError(geometry, config);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::ImplicitEuler_Iteration(CGeometry *geometry, CSolver**, CConfig *config) {
  SU2_ZONE_SCOPED

  PrepareImplicitIteration(geometry, nullptr, config);

  /*--- Solve or smooth the linear system. ---*/

  SU2_OMP_FOR_(schedule(static,OMP_MIN_SIZE) SU2_NOWAIT)
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

  CompleteImplicitIteration(geometry, nullptr, config);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::ComputeVorticityAndStrainMag(const CConfig& config, const CGeometry *geometry, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  auto& StrainMag = nodes->GetStrainMag();

  ompMasterAssignBarrier(StrainMag_Max,0.0, Omega_Max,0.0);

  su2double strainMax = 0.0, omegaMax = 0.0;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    const auto VelocityGradient = nodes->GetVelocityGradient(iPoint);
    auto Vorticity = nodes->GetVorticity(iPoint);

    /*--- Vorticity ---*/

    Vorticity[0] = 0.0;
    Vorticity[1] = 0.0;
    Vorticity[2] = VelocityGradient(1,0)-VelocityGradient(0,1);

    if (nDim == 3) {
      Vorticity[0] = VelocityGradient(2,1)-VelocityGradient(1,2);
      Vorticity[1] = -(VelocityGradient(2,0)-VelocityGradient(0,2));
    }

    /*--- Strain Magnitude ---*/

    const su2double vy = nodes->GetVelocity(iPoint, 1);
    const su2double y = geometry->nodes->GetCoord(iPoint, 1);
    AD::StartPreacc();
    AD::SetPreaccIn(VelocityGradient, nDim, nDim);
    AD::SetPreaccIn(vy, y);

    StrainMag(iPoint) = 0.0;

    /*--- Add diagonal part ---*/

    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      StrainMag(iPoint) += pow(VelocityGradient(iDim, iDim), 2);
    }
    if (config.GetAxisymmetric() && y > EPS) {
      StrainMag(iPoint) += pow(vy / y, 2);
    }

    /*--- Add off diagonals ---*/

    StrainMag(iPoint) += 2.0*pow(0.5*(VelocityGradient(0,1) + VelocityGradient(1,0)), 2);

    if (nDim == 3) {
      StrainMag(iPoint) += 2.0*pow(0.5*(VelocityGradient(0,2) + VelocityGradient(2,0)), 2);
      StrainMag(iPoint) += 2.0*pow(0.5*(VelocityGradient(1,2) + VelocityGradient(2,1)), 2);
    }

    StrainMag(iPoint) = sqrt(2.0*StrainMag(iPoint));
    AD::SetPreaccOut(StrainMag(iPoint));
    AD::EndPreacc();

    /*--- The derivative with respect to strainMax and omegaMax is not required. ---*/
    bool wa = AD::PauseRecording();
    strainMax = max(strainMax, StrainMag(iPoint));
    omegaMax = max(omegaMax, GeometryToolbox::Norm(3, Vorticity));
    AD::ResumeRecording(wa);
  }
  END_SU2_OMP_FOR

  if ((iMesh == MESH_0) && (config.GetComm_Level() == COMM_FULL)) {
    atomicMax(strainMax, StrainMag_Max);
    atomicMax(omegaMax, Omega_Max);

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      su2double MyOmega_Max = Omega_Max;
      su2double MyStrainMag_Max = StrainMag_Max;

      SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetInletAtVertex(const su2double* val_inlet, unsigned short iMarker,
                                                unsigned long iVertex) {
  SU2_ZONE_SCOPED

  /*--- Alias positions within inlet file for readability ---*/

  unsigned short T_position = nDim;
  unsigned short P_position = nDim + 1;
  unsigned short FlowDir_position = nDim + 2;

  /*--- Note that it is not necessary anymore to use normalized normals for the inlet velocity ---*/


  /*--- Store the values in our inlet data structures. ---*/

  Inlet_Ttotal[iMarker][iVertex] = val_inlet[T_position];
  Inlet_Ptotal[iMarker][iVertex] = val_inlet[P_position];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Inlet_FlowDir[iMarker][iVertex][iDim] = val_inlet[FlowDir_position + iDim];
  }
}

template <class V, ENUM_REGIME R>
su2double CFVMFlowSolverBase<V, R>::GetInletAtVertex(unsigned short iMarker, unsigned long iVertex,
                                                     const CGeometry* geometry, su2double* val_inlet) const {
  SU2_ZONE_SCOPED

  const auto T_position = nDim;
  const auto P_position = nDim + 1;
  const auto FlowDir_position = nDim + 2;
  val_inlet[T_position] = Inlet_Ttotal[iMarker][iVertex];
  val_inlet[P_position] = Inlet_Ptotal[iMarker][iVertex];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    val_inlet[FlowDir_position + iDim] = Inlet_FlowDir[iMarker][iVertex][iDim];
  }

  /*--- Compute boundary face area for this vertex. ---*/

  su2double Normal[MAXNDIM] = {0.0};
  geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
  return GeometryToolbox::Norm(nDim, Normal);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetUniformInlet(const CConfig* config, unsigned short iMarker) {
  SU2_ZONE_SCOPED

  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
    const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    const su2double p_total = config->GetInletPtotal(Marker_Tag);
    const su2double t_total = config->GetInletTtotal(Marker_Tag);
    const su2double* flow_dir = config->GetInletFlowDir(Marker_Tag);

    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      Inlet_Ttotal[iMarker][iVertex] = t_total;
      Inlet_Ptotal[iMarker][iVertex] = p_total;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) Inlet_FlowDir[iMarker][iVertex][iDim] = flow_dir[iDim];
    }
  } else if (config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_INLET) {
    const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    const su2double p = config->GetInlet_Pressure(Marker_Tag);
    const su2double t = config->GetInlet_Temperature(Marker_Tag);
    const su2double* vel = config->GetInlet_Velocity(Marker_Tag);

    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      Inlet_Ttotal[iMarker][iVertex] = t;
      Inlet_Ptotal[iMarker][iVertex] = p;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) Inlet_FlowDir[iMarker][iVertex][iDim] = vel[iDim];
    }
  } else {
    /*--- For now, non-inlets just get set to zero. In the future, we
     can do more customization for other boundary types here. ---*/
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      Inlet_Ttotal[iMarker][iVertex] = 0.0;
      Inlet_Ptotal[iMarker][iVertex] = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) Inlet_FlowDir[iMarker][iVertex][iDim] = 0.0;
    }
  }
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::UpdateCustomBoundaryConditions(
    CGeometry** geometry_container, CSolver*** solver_container, CConfig *config) {
  SU2_ZONE_SCOPED

  struct {
    const CSolver* fine_solver{nullptr};
    CSolver* coarse_solver{nullptr};
    unsigned short marker{0};
    unsigned short var{0};

    su2double Get(unsigned long vertex) const {
      if (var == 0) {
        return fine_solver->GetInletTtotal(marker, vertex);
      } else if (var == 1) {
        return fine_solver->GetInletPtotal(marker, vertex);
      }
      return fine_solver->GetInletFlowDir(marker, vertex, var - 2);
    }

    void Set(unsigned long vertex, const su2double& val) const {
      if (var == 0) {
        coarse_solver->SetInletTtotal(marker, vertex, val);
      } else if (var == 1) {
        coarse_solver->SetInletPtotal(marker, vertex, val);
      }
      coarse_solver->SetInletFlowDir(marker, vertex, var - 2, val);
    }
  } inlet_values;

  for (auto mg_coarse = 1u; mg_coarse <= config->GetnMGLevels(); ++mg_coarse) {
    const auto mg_fine = mg_coarse - 1;
    inlet_values.fine_solver = solver_container[mg_fine][FLOW_SOL];
    inlet_values.coarse_solver = solver_container[mg_coarse][FLOW_SOL];

    for (auto marker = 0u; marker < config->GetnMarker_All(); ++marker) {
      if (config->GetMarker_All_KindBC(marker) != INLET_FLOW) continue;
      inlet_values.marker = marker;
      for (inlet_values.var = 0; inlet_values.var < 2 + nDim; ++inlet_values.var) {
        geometry_container[mg_coarse]->SetMultiGridMarkerQuantity(geometry_container[mg_fine], marker, inlet_values);
      }
    }
  }
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::LoadRestart_impl(CGeometry **geometry, CSolver ***solver, CConfig *config, int iter,
                                                bool update_geo, su2double* SolutionRestart,
                                                unsigned short nVar_Restart) {
  SU2_ZONE_SCOPED

  /*--- Restart the solution from file information ---*/

  string restart_filename = config->GetSolution_FileName();
  const bool static_fsi = ((config->GetTime_Marching() == TIME_MARCHING::STEADY) && config->GetFSI_Simulation());

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {

    if (nVar_Restart == 0) nVar_Restart = nVar;

    /*--- Skip coordinates ---*/

    unsigned short skipVars = nDim;

    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/
    if (config->GetRead_Binary_Restart()) {
      restart_filename = config->GetFilename(restart_filename, ".dat", iter);
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      restart_filename = config->GetFilename(restart_filename, ".csv", iter);
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    bool steady_restart = config->GetSteadyRestart();
    if (update_geo && dynamic_grid) {
      auto notFound = fields.end();
      if (find(fields.begin(), notFound, string("\"Grid_Velocity_x\"")) == notFound) {
        if (rank == MASTER_NODE)
          cout << "\nWARNING: The restart file does not contain grid velocities, these will be set to zero.\n" << endl;
        steady_restart = true;
      }
    }

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long counter = 0;
    for (auto iPoint_Global = 0ul; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {

      /*--- Retrieve local index. If this node from the restart file lives
      on the current processor, we will load and instantiate the vars. ---*/

      const auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {

        /*--- We need to store this point's data, so jump to the correct
        offset in the buffer of data from the restart file and load it. ---*/

        auto index = counter * Restart_Vars[1] + skipVars;

        if (SolutionRestart == nullptr) {
          for (auto iVar = 0u; iVar < nVar_Restart; iVar++)
            nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index+iVar]);
        }
        else {
          /*--- Used as buffer, allows defaults for nVar > nVar_Restart. ---*/
          for (auto iVar = 0u; iVar < nVar_Restart; iVar++)
            SolutionRestart[iVar] = Restart_Data[index + iVar];
          nodes->SetSolution(iPoint_Local, SolutionRestart);
        }

        /*--- For dynamic meshes, read in and store the
        grid coordinates and grid velocities for each node. ---*/

        if (dynamic_grid && update_geo) {

          /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
          /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
          /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter * Restart_Vars[1];
          const auto* Coord = &Restart_Data[index];

          su2double GridVel[MAXNDIM] = {0.0};
          if (!steady_restart) {
            /*--- Move the index forward to get the grid velocities. ---*/
            index += skipVars + nVar_Restart + config->GetnTurbVar();
            for (auto iDim = 0u; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
          }

          for (auto iDim = 0u; iDim < nDim; iDim++) {
            geometry[MESH_0]->nodes->SetCoord(iPoint_Local, iDim, Coord[iDim]);
            geometry[MESH_0]->nodes->SetGridVel(iPoint_Local, iDim, GridVel[iDim]);
          }
        }

        /*--- For static FSI problems, grid_movement is 0 but we need to read in and store the
        grid coordinates for each node (but not the grid velocities, as there are none). ---*/

        if (static_fsi && update_geo) {
        /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          const auto* Coord = &Restart_Data[index];

          for (auto iDim = 0u; iDim < nDim; iDim++) {
            geometry[MESH_0]->nodes->SetCoord(iPoint_Local, iDim, Coord[iDim]);
          }
        }

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }
    }

    /*--- Detect a wrong solution file ---*/

    if (counter != nPointDomain) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" does not match with the mesh file.\n") +
                     string("This can be caused by empty lines at the end of the file."), CURRENT_FUNCTION);
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- Update the geometry for flows on deforming meshes. ---*/

  if ((dynamic_grid || static_fsi) && update_geo) {

    CGeometry::UpdateGeometry(geometry, config);

    if (dynamic_grid) {
      for (auto iMesh = 0u; iMesh <= config->GetnMGLevels(); iMesh++) {

        /*--- Compute the grid velocities on the coarser levels. ---*/
        if (iMesh) geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMesh - 1]);
        else {
          geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, MPI_QUANTITIES::GRID_VELOCITY);
          geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, MPI_QUANTITIES::GRID_VELOCITY);
        }
      }
    }
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We also call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver[MESH_0][FLOW_SOL]->InitiateComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);
  solver[MESH_0][FLOW_SOL]->CompleteComms(geometry[MESH_0], config, MPI_QUANTITIES::SOLUTION);

  /*--- For turbulent/species simulations the flow preprocessing is done by the turbulence/species solver
   *    after it loads its variables (they are needed to compute flow primitives). In case turbulence and species, the
   *    species solver does all the Pre-/Postprocessing. ---*/
  if (config->GetKind_Turb_Model() == TURB_MODEL::NONE &&
      config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
    solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
  }

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {
    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][FLOW_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][FLOW_SOL]->GetNodes()->GetSolution());
    solver[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);
    solver[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, MPI_QUANTITIES::SOLUTION);

    if (config->GetKind_Turb_Model() == TURB_MODEL::NONE &&
        config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
      solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
    }
  }

  /*--- Update the old geometry (coordinates n and n-1) in dual time-stepping strategy. ---*/
  const bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                          (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
  if (dual_time && config->GetGrid_Movement() && !config->GetDeform_Mesh() &&
      (config->GetKind_GridMovement() != RIGID_MOTION)) {
    Restart_OldGeometry(geometry[MESH_0], config);
  }

  /*--- Go back to single threaded execution. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
  /*--- Delete the class memory that is used to load the restart. ---*/

    Restart_Vars = decltype(Restart_Vars){};
    Restart_Data = decltype(Restart_Data){};
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::LoadRestart(CGeometry **geometry, CSolver ***solver,
                                           CConfig *config, int iter, bool update_geo) {
  SU2_ZONE_SCOPED
  LoadRestart_impl(geometry, solver, config, iter, update_geo);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container,
                                                   CConfig *config, unsigned long TimeIter) {
  SU2_ZONE_SCOPED

  const bool restart = (config->GetRestart() || config->GetRestart_Flow());
  const bool rans = (config->GetKind_Turb_Model() != TURB_MODEL::NONE);
  const bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                          (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL {

  unsigned long iPoint;
  unsigned short iMesh;

  /*--- Check if a verification solution is to be computed. ---*/
  if ((VerificationSolution) && (TimeIter == 0) && !restart) {

    /*--- Loop over the multigrid levels. ---*/
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

      /*--- Loop over all grid points. ---*/
      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {

        /* Set the pointers to the coordinates and solution of this DOF. */
        const su2double *coor = geometry[iMesh]->nodes->GetCoord(iPoint);
        su2double *solDOF     = solver_container[iMesh][FLOW_SOL]->GetNodes()->GetSolution(iPoint);

        /* Set the solution in this DOF to the initial condition provided by
           the verification solution class. This can be the exact solution,
           but this is not necessary. */
        VerificationSolution->GetInitialCondition(coor, solDOF);
      }
      END_SU2_OMP_FOR
    }
  }

  /*--- The value of the solution for the first iteration of the dual time ---*/

  if (dual_time && TimeIter == config->GetRestart_Iter()) {
    PushSolutionBackInTime(TimeIter, restart, rans, solver_container, geometry, config);
  }

  }
  END_SU2_OMP_PARALLEL

}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::PushSolutionBackInTime(unsigned long TimeIter, bool restart, bool rans,
                                                      CSolver*** solver_container, CGeometry** geometry,
                                                      CConfig* config) {
  SU2_ZONE_SCOPED

  /*--- Push back the initial condition to previous solution containers
   for a 1st-order restart or when simply initializing to freestream. ---*/

  for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
    solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n1();
    if (rans) {
      solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
      solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1();
    }

    if (dynamic_grid) {
      geometry[iMesh]->nodes->SetVolume_n();
      geometry[iMesh]->nodes->SetVolume_nM1();
    }

    if (config->GetGrid_Movement()) {
      geometry[iMesh]->nodes->SetCoord_n();
      geometry[iMesh]->nodes->SetCoord_n1();
    }
  }

  if (restart && (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) {

    /*--- Load an additional restart file for a 2nd-order restart. ---*/

    solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, TimeIter-1, true);

    /*--- Load an additional restart file for the turbulence model. ---*/
    if (rans)
      solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, TimeIter-1, false);

    /*--- Push back this new solution to time level N. ---*/

    for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
      if (rans) solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();

      geometry[iMesh]->nodes->SetVolume_n();
      if (config->GetGrid_Movement()) geometry[iMesh]->nodes->SetCoord_n();
    }
  }
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                                     CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR) ||
                         (config->GetKind_FluidModel() == IDEAL_GAS);
  const auto iVel = prim_idx.Velocity();

  /*--- Blazek chapter 8.:
   * The components of the momentum residual normal to the symmetry plane are zeroed out.
   * The gradients have already been corrected acording to Eq. (8.40).
   * Contrary to Blazek we keep some scalar fluxes computed on the boundary to improve stability (see below). ---*/

  /*--- Loop over all the vertices on this boundary marker. ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Halo points do not need to be considered. ---*/
    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Get the normal of the current symmetry. This may be the original normal of the vertex
     * or a modified normal if there are intersecting symmetries. ---*/

    su2double Normal[MAXNDIM] = {}, UnitNormal[MAXNDIM] = {};
    geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
    const auto it = geometry->symmetryNormals[val_marker].find(iVertex);

    if (it != geometry->symmetryNormals[val_marker].end()) {
      for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = it->second[iDim];
    } else {
      const su2double Area = GeometryToolbox::Norm(nDim, Normal);
      for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = Normal[iDim] / Area;
    }

    /*--- Energy terms due to grid movement (aka work of pressure forces). ---*/
    if (dynamic_grid) {
      su2double* V_reflected = GetCharacPrimVar(val_marker, iVertex);

      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                geometry->nodes->GetGridVel(iPoint));

      /*--- Normal vector for this vertex (negate for outward convention). ---*/
      for (auto iDim = 0u; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      for (auto iVar = 0u; iVar < nPrimVar; iVar++)
        V_reflected[iVar] = nodes->GetPrimitive(iPoint, iVar);

      su2double ProjVelocity_i = nodes->GetProjVel(iPoint, UnitNormal);
      /*--- Adjustment to v.n due to grid movement. ---*/
      ProjVelocity_i -= GeometryToolbox::DotProduct(nDim, geometry->nodes->GetGridVel(iPoint), UnitNormal);

      for (auto iDim = 0u; iDim < nDim; iDim++)
        V_reflected[iDim + iVel] = nodes->GetVelocity(iPoint, iDim) - ProjVelocity_i * UnitNormal[iDim];

      /*--- Get current solution at this boundary node. ---*/
      const su2double* V_domain = nodes->GetPrimitive(iPoint);

      /*--- Set Primitive and Secondary for numerics class. ---*/
      conv_numerics->SetPrimitive(V_domain, V_reflected);
      if (!ideal_gas) {
        conv_numerics->SetSecondary(nodes->GetSecondary(iPoint), nodes->GetSecondary(iPoint));
      }

      /*--- Compute the residual using an upwind scheme. ---*/
      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Use just the energy fluxes to update the residual, adding the others would
       * increase numerical diffusion which we wish to avoid if possible. ---*/
      for (auto iVar = iVel + nDim; iVar < nVar; iVar++) {
        LinSysRes(iPoint, iVar) += residual.residual[iVar];
      }
      if (implicit) {
        auto* block = Jacobian.GetBlock(iPoint, iPoint);
        /*--- But in the Jacobian we also include the mass flux, this allows some cases with
         * motion to use larger CFL, for example pywrapper_translating_naca0012. ---*/
        for (auto iVar = 0u; iVar < nVar; iVar++) {
          if (iVar < iVel || iVar >= iVel + nDim) {
            for (auto jVar = 0u; jVar < nVar; jVar++) {
              block[iVar * nVar + jVar] += SU2_TYPE::GetValue(residual.jacobian_i[iVar][jVar]);
            }
          }
        }
      }
    }

    /*--- Explicitly set the velocity components normal to the symmetry plane to zero.
     * This is necessary because the modification of the residual leaves the problem
     * underconstrained (the normal residual is zero regardless of the normal velocity). ---*/

    su2double* solutionOld = nodes->GetSolution_Old(iPoint);

    su2double gridVel[MAXNDIM] = {};
    if (dynamic_grid) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        gridVel[iDim] = geometry->nodes->GetGridVel(iPoint)[iDim];
      }
      if (FlowRegime == ENUM_REGIME::COMPRESSIBLE) {
        for(auto iDim = 0u; iDim < nDim; iDim++) {
          /*--- Multiply by density since we are correcting conservative variables. ---*/
          gridVel[iDim] *= nodes->GetDensity(iPoint);
        }
      }
    }
    su2double vp = 0.0;
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      vp += (solutionOld[iVel + iDim] - gridVel[iDim]) * UnitNormal[iDim];
    }
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      solutionOld[iVel + iDim] -= vp * UnitNormal[iDim];
    }

    /*--- Keep only the tangential part of the momentum residuals. ---*/
    su2double normalRes = 0.0;
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      normalRes += LinSysRes(iPoint, iVel + iDim) * UnitNormal[iDim];
    }
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      LinSysRes(iPoint, iVel + iDim) -= normalRes * UnitNormal[iDim];
    }

    /*--- Jacobian contribution for implicit integration. ---*/
    if (implicit) {
      /*--- Modify the Jacobians according to the modification of the residual
       * J_new = (I - n * n^T) * J where n = {0, nx, ny, nz, 0, ...} ---*/
      su2double mat[MAXNVAR * MAXNVAR] = {};

      for (auto iVar = 0u; iVar < nVar; iVar++)
        mat[iVar * nVar + iVar] = 1;
      for (auto iDim = 0u; iDim < nDim; iDim++)
        for (auto jDim = 0u; jDim < nDim; jDim++)
          mat[(iDim + iVel) * nVar + jDim + iVel] -= UnitNormal[iDim] * UnitNormal[jDim];

      auto ModifyJacobian = [&](const unsigned long jPoint) {
        su2double jac[MAXNVAR * MAXNVAR], newJac[MAXNVAR * MAXNVAR];
        const auto view = Jacobian.GetBlockView(iPoint, jPoint);
        if (!view) return;
        for (auto iVar = 0u; iVar < nVar; iVar++)
          for (auto jVar = 0u; jVar < nVar; jVar++) jac[iVar * nVar + jVar] = view(iVar, jVar);

        CBlasStructure().gemm(nVar, nVar, nVar, mat, jac, newJac, config);

        Jacobian.SetBlock(iPoint, jPoint, newJac);
      };
      ModifyJacobian(iPoint);
      for (size_t iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
        ModifyJacobian(geometry->nodes->GetPoint(iPoint, iNeigh));
      }
    }

    /*--- Correction for multigrid. ---*/
    normalRes = 0.0;
    su2double* Res_TruncError = nodes->GetResTruncError(iPoint);
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      normalRes += Res_TruncError[iVel + iDim] * UnitNormal[iDim];
    }
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      Res_TruncError[iVel + iDim] -= normalRes * UnitNormal[iDim];
    }
  }
  END_SU2_OMP_FOR

}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::BC_Periodic(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                           CConfig* config) {
  SU2_ZONE_SCOPED

  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated correctly during the communications. For implicit calculations,
   the Jacobians and linear system are also correctly adjusted here. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic() / 2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::BC_Fluid_Interface(CGeometry* geometry, CSolver** solver_container,
                                                           CNumerics* conv_numerics, CNumerics* visc_numerics,
                                                           CConfig* config) {
  SU2_ZONE_SCOPED

  const bool ideal_gas = config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS;

  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, jVar, iMarker, nDonorVertex;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool viscous = config->GetViscous();

  su2double Normal[MAXNDIM] = {0.0};
  su2double PrimVar_i[MAXNVAR] = {0.0};
  su2double PrimVar_j[MAXNVAR] = {0.0};
  su2double Secondary_j[MAXNVAR] = {0.0};
  su2double Residual[MAXNVAR] = {0.0};
  su2double** Jacobian_i = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Jacobian_i[iVar] = new su2double[nVar];

  su2double weight;
  su2double P_static, rho_static;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->nodes->GetDomain(iPoint)) {
          nDonorVertex = GetnSlidingStates(iMarker, iVertex);

          /*--- Initialize Residual, this will serve to accumulate the average ---*/

          for (iVar = 0; iVar < nVar; iVar++) {
            Residual[iVar] = 0.0;
            for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          }

          /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

          for (jVertex = 0; jVertex < nDonorVertex; jVertex++) {
            Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = nodes->GetPrimitive(iPoint, iVar);
              PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }

            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

            weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set primitive variables ---*/

            conv_numerics->SetPrimitive(PrimVar_i, PrimVar_j);

            if (FlowRegime == ENUM_REGIME::COMPRESSIBLE && !ideal_gas) {
              P_static = PrimVar_j[nDim + 1];
              rho_static = PrimVar_j[nDim + 2];
              GetFluidModel()->SetTDState_Prho(P_static, rho_static);

              Secondary_j[0] = GetFluidModel()->GetdPdrho_e();
              Secondary_j[1] = GetFluidModel()->GetdPde_rho();

              conv_numerics->SetSecondary(nodes->GetSecondary(iPoint), Secondary_j);
            }

            /*--- Set the normal vector ---*/

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

            conv_numerics->SetNormal(Normal);

            if (dynamic_grid)
              conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

            /*--- Compute the convective residual using an upwind scheme ---*/

            auto residual = conv_numerics->ComputeResidual(config);

            /*--- Accumulate the residuals to compute the average ---*/

            for (iVar = 0; iVar < nVar; iVar++) {
              Residual[iVar] += weight * residual[iVar];
              for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] += weight * residual.jacobian_i[iVar][jVar];
            }
          }

          /*--- Add Residuals and Jacobians ---*/

          LinSysRes.AddBlock(iPoint, Residual);

          if (implicit) Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

          if (viscous) {
            /*--- Initialize Residual, this will serve to accumulate the average ---*/

            for (iVar = 0; iVar < nVar; iVar++) {
              Residual[iVar] = 0.0;
              for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
            }

            /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

            for (jVertex = 0; jVertex < nDonorVertex; jVertex++) {
              PrimVar_j[nDim + 5] = GetSlidingState(iMarker, iVertex, nDim + 5, jVertex);
              PrimVar_j[nDim + 6] = GetSlidingState(iMarker, iVertex, nDim + 6, jVertex);

              /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

              weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

              /*--- Set the normal vector and the coordinates ---*/

              visc_numerics->SetNormal(Normal);
              su2double Coord_Reflected[MAXNDIM];
              GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                        geometry->nodes->GetCoord(iPoint), Coord_Reflected);
              visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

              /*--- Primitive variables, and gradient ---*/

              visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);
              visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                                nodes->GetGradient_Primitive(iPoint));

              /*--- Turbulent kinetic energy ---*/

              if (config->GetKind_Turb_Model() == TURB_MODEL::SST)
                visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint, 0),
                                                    solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint, 0));

              /*--- Compute and update residual ---*/

              auto residual = visc_numerics->ComputeResidual(config);

              /*--- Accumulate the residuals to compute the average ---*/

              for (iVar = 0; iVar < nVar; iVar++) {
                Residual[iVar] += weight * residual[iVar];
                for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] += weight * residual.jacobian_i[iVar][jVar];
              }
            }

            LinSysRes.SubtractBlock(iPoint, Residual);

            /*--- Jacobian contribution for implicit integration ---*/

            if (implicit) Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
          }
        }
      }
      END_SU2_OMP_FOR
    }
  }

  for (iVar = 0; iVar < nVar; iVar++) delete[] Jacobian_i[iVar];
  delete[] Jacobian_i;
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::BC_Custom(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                         CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  /* Check for a verification solution. */

  if (VerificationSolution) {
    unsigned short iVar;
    unsigned long iVertex, iPoint;

    bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

    /*--- Get the physical time. ---*/

    su2double time = 0.0;
    if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

    /*--- Loop over all the vertices on this boundary marker ---*/

    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      /*--- Get the point index for the current node. ---*/

      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

      if (geometry->nodes->GetDomain(iPoint)) {
        /*--- Get the coordinates for the current node. ---*/

        const su2double* coor = geometry->nodes->GetCoord(iPoint);

        /*--- Get the conservative state from the verification solution. ---*/

        su2double Solution[MAXNVAR] = {0.0};
        VerificationSolution->GetBCState(coor, time, Solution);

        /*--- For verification cases, we will apply a strong Dirichlet
         condition by setting the solution values at the boundary nodes
         directly and setting the residual to zero at those nodes. ---*/

        nodes->SetSolution_Old(iPoint, Solution);
        nodes->SetSolution(iPoint, Solution);
        nodes->SetRes_TruncErrorZero(iPoint);
        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Adjust rows of the Jacobian (includes 1 in the diagonal) ---*/

        if (implicit) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Jacobian.DeleteValsRowi(iPoint, iVar);
          }
        }
      }
    }
    END_SU2_OMP_FOR

  } else {
    /* The user must specify the custom BC's here. */
    SU2_MPI::Error("Implement customized boundary conditions here.", CURRENT_FUNCTION);
  }
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::EdgeFluxResidual(const CGeometry *geometry,
                                                const CSolver* const* solvers,
                                                CConfig *config) {
  SU2_ZONE_SCOPED

  if (!edgeNumerics) {
    if (!ReducerStrategy && (omp_get_max_threads() > 1) &&
        (config->GetEdgeColoringGroupSize() % Double::Size != 0)) {
      SU2_MPI::Error("When using vectorization, the EDGE_COLORING_GROUP_SIZE must be divisible "
                     "by the SIMD length (2, 4, or 8).", CURRENT_FUNCTION);
    }
    InstantiateEdgeNumerics(solvers, config);

    /*--- The SIMD numerics do not use gradients of density and enthalpy. ---*/
    if (!config->GetContinuous_Adjoint()) {
      SU2_OMP_SAFE_GLOBAL_ACCESS(nPrimVarGrad = std::min<unsigned short>(nDim + 2, nPrimVarGrad);)
    }
  }

  /*--- Non-physical counter. ---*/
  unsigned long counterLocal = 0;
  SU2_OMP_MASTER
  ErrorCounter = 0;
  END_SU2_OMP_MASTER

  su2activevector* massFluxes = config->GetBounded_Scalar() ? &EdgeMassFluxes : nullptr;

  /*--- For hybrid parallel AD, pause preaccumulation if there is shared reading of
  * variables, otherwise switch to the faster adjoint evaluation mode. ---*/
  bool pausePreacc = false;
  if (ReducerStrategy) pausePreacc = AD::PausePreaccumulation();
  else AD::StartNoSharedReading();

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring) {
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for(auto k = 0ul; k < color.size; k += Double::Size) {
      Int iEdge;
      Double mask;
      for (auto j = 0ul; j < Double::Size; ++j) {
        bool in = (k+j < color.size);
        mask[j] = in;
        iEdge[j] = color.indices[k+j*in];
      }

      if (ReducerStrategy) {
        edgeNumerics->ComputeFlux(iEdge, *config, *geometry, *nodes, UpdateType::REDUCTION, mask, EdgeFluxes, Jacobian, massFluxes);
      } else {
        edgeNumerics->ComputeFlux(iEdge, *config, *geometry, *nodes, UpdateType::COLORING, mask, LinSysRes, Jacobian, massFluxes);
      }
      if (MGLevel == MESH_0) {
        for (auto j = 0ul; j < Double::Size; ++j)
          counterLocal += (nodes->NonPhysicalEdgeCounter[iEdge[j]] > 0);
      }
    }
    END_SU2_OMP_FOR
  }

  FinalizeResidualComputation(geometry, pausePreacc, counterLocal, config);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::SumEdgeFluxes(const CGeometry* geometry) {
  SU2_ZONE_SCOPED

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    LinSysRes.SetBlock_Zero(iPoint);

    for (auto iEdge : geometry->nodes->GetEdges(iPoint)) {
      if (iPoint == geometry->edges->GetNode(iEdge,0))
        LinSysRes.AddBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
      else
        LinSysRes.SubtractBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
    }
  }
  END_SU2_OMP_FOR
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container,
                                                             CConfig *config, unsigned short iRKStep, unsigned short iMesh,
                                                             unsigned short RunTime_EqSystem) {
  /*--- Local variables ---*/

  unsigned short iVar, iMarker, iDim, iNeigh;
  unsigned long iPoint, jPoint, iEdge, iVertex;

  const su2double *U_time_nM1 = nullptr, *U_time_n = nullptr, *U_time_nP1 = nullptr;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  const su2double *Normal = nullptr, *GridVel_i = nullptr, *GridVel_j = nullptr;
  su2double Residual_GCL;

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool first_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  const bool second_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);

  /*--- Store the physical time step ---*/

  TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {

    /*--- Loop over all nodes (excluding halos) ---*/

    AD::StartNoSharedReading();

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (first_order)
          LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (second_order)
          LinSysRes(iPoint,iVar) += ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                                     +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
      }
    }
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();

  }

  else {

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

      for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {

        iEdge = geometry->nodes->GetEdge(iPoint, iNeigh);
        Normal = geometry->edges->GetNormal(iEdge);

        jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
        GridVel_j = geometry->nodes->GetGridVel(jPoint);

        /*--- Determine whether to consider the normal outward or inward. ---*/
        su2double dir = (iPoint < jPoint)? 0.5 : -0.5;

        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL += dir*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
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
          for (iDim = 0; iDim < nDim; iDim++)
            Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];

          /*--- Compute the GCL component of the source term for node i ---*/

          U_time_n = nodes->GetSolution_time_n(iPoint);
          for (iVar = 0; iVar < nVar; iVar++)
            LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
        }
        END_SU2_OMP_FOR
      }
    }

    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/

    AD::StartNoSharedReading();

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/

      Volume_nM1 = geometry->nodes->GetVolume_nM1(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (first_order)
          LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (second_order)
          LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
                                     + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
      }
    }
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();
  }

}

/*--- Helpers shared by Pressure_Forces, Momentum_Forces and Friction_Forces. Free function
 *    templates rather than members: none need instance state, and the AeroCoeffs/AeroCoeffsArray
 *    types are simply deduced from the caller's arguments. ---*/
namespace {

/*!
 * \brief Whether a marker's boundary kind is one of the momentum (inlet/outlet/actuator-disk/
 *        engine) surfaces handled by Momentum_Forces, used both to gate accumulation into
 *        MntCoeff and to gate the later CEff/CMerit derivation from it.
 * \param[in] Boundary - Boundary kind of the marker (config->GetMarker_All_KindBC(iMarker)).
 */
inline bool IsMomentumBoundary(unsigned short Boundary) {
  return (Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) || (Boundary == ACTDISK_INLET) ||
         (Boundary == ACTDISK_OUTLET) || (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST);
}

/*!
 * \brief Find iMarker's index within the monitoring markers, and its reference origin if found.
 * \param[in] config - Problem definition.
 * \param[in] iMarker - Marker to look up.
 * \param[in] Monitoring - config->GetMarker_All_Monitoring(iMarker).
 * \param[in,out] Origin - Set to iMarker's reference origin if found, left unchanged otherwise.
 * \return Index within the monitoring markers, or -1 if iMarker is not monitored, or not found
 *         among the monitoring markers.
 */
inline int FindMonitoringIndex(const CConfig* config, unsigned long iMarker, unsigned short Monitoring,
                               std::array<su2double, 3>& Origin) {
  if (Monitoring != YES) return -1;
  const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
  const auto nMarker_Monitoring = static_cast<int>(config->GetnMarker_Monitoring());
  for (int iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
    if (Marker_Tag == config->GetMarker_Monitoring_TagBound(iMarker_Monitoring)) {
      Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      return iMarker_Monitoring;
    }
  }
  return -1;
}

/*!
 * \brief Fold one thread's partial force/moment coefficient contribution (accumulated over
 *        its share of a marker's vertices) into the per-marker, AllBound and per-surface
 *        aerodynamic coefficient totals, in a single critical section. Does not touch
 *        CEff/CMerit (nonlinear ratios), those are derived once from the fully-reduced totals
 *        after all threads have merged in.
 * \param[in] iMarker - Marker index the contribution belongs to.
 * \param[in] iMarker_Monitoring - If iMarker is monitored, index within the monitoring markers, -1 otherwise.
 * \param[in] partial - This thread's partial coefficients (CD, CL, ..., CQ only).
 * \param[in,out] coeffArray - Per-marker totals to update at iMarker.
 * \param[in,out] allBoundCoeff - Totals over all boundaries to update.
 * \param[in,out] surfaceCoeff - Per-monitoring-surface totals to update.
 */
template <class AeroCoeffsT, class AeroCoeffsArrayT>
void AddCoeffContribution(unsigned long iMarker, int iMarker_Monitoring, const AeroCoeffsT& partial,
                          AeroCoeffsArrayT& coeffArray, AeroCoeffsT& allBoundCoeff, AeroCoeffsArrayT& surfaceCoeff) {
  SU2_OMP_CRITICAL {
    coeffArray.CD[iMarker] += partial.CD;
    coeffArray.CL[iMarker] += partial.CL;
    coeffArray.CSF[iMarker] += partial.CSF;
    coeffArray.CFx[iMarker] += partial.CFx;
    coeffArray.CFy[iMarker] += partial.CFy;
    coeffArray.CFz[iMarker] += partial.CFz;
    coeffArray.CMx[iMarker] += partial.CMx;
    coeffArray.CMy[iMarker] += partial.CMy;
    coeffArray.CMz[iMarker] += partial.CMz;
    coeffArray.CoPx[iMarker] += partial.CoPx;
    coeffArray.CoPy[iMarker] += partial.CoPy;
    coeffArray.CoPz[iMarker] += partial.CoPz;
    coeffArray.CT[iMarker] += partial.CT;
    coeffArray.CQ[iMarker] += partial.CQ;

    allBoundCoeff.CD += partial.CD;
    allBoundCoeff.CL += partial.CL;
    allBoundCoeff.CSF += partial.CSF;
    allBoundCoeff.CFx += partial.CFx;
    allBoundCoeff.CFy += partial.CFy;
    allBoundCoeff.CFz += partial.CFz;
    allBoundCoeff.CMx += partial.CMx;
    allBoundCoeff.CMy += partial.CMy;
    allBoundCoeff.CMz += partial.CMz;
    allBoundCoeff.CoPx += partial.CoPx;
    allBoundCoeff.CoPy += partial.CoPy;
    allBoundCoeff.CoPz += partial.CoPz;
    allBoundCoeff.CT += partial.CT;
    allBoundCoeff.CQ += partial.CQ;

    /*--- Compute the coefficients per surface ---*/

    if (iMarker_Monitoring >= 0) {
      surfaceCoeff.CL[iMarker_Monitoring] += partial.CL;
      surfaceCoeff.CD[iMarker_Monitoring] += partial.CD;
      surfaceCoeff.CSF[iMarker_Monitoring] += partial.CSF;
      surfaceCoeff.CFx[iMarker_Monitoring] += partial.CFx;
      surfaceCoeff.CFy[iMarker_Monitoring] += partial.CFy;
      surfaceCoeff.CFz[iMarker_Monitoring] += partial.CFz;
      surfaceCoeff.CMx[iMarker_Monitoring] += partial.CMx;
      surfaceCoeff.CMy[iMarker_Monitoring] += partial.CMy;
      surfaceCoeff.CMz[iMarker_Monitoring] += partial.CMz;
    }
  }
  END_SU2_OMP_CRITICAL
}

/*!
 * \brief Project summed force/moment components onto the wind axes (Alpha, Beta) to get the
 *        standard aerodynamic coefficients. Identical formulas are used by Pressure_Forces,
 *        Momentum_Forces and Friction_Forces, applied respectively to their inviscid,
 *        momentum and viscous force/moment sums. Does not set CEff/CMerit (nonlinear ratios,
 *        derived later from fully-reduced totals) nor CSF/CMx/CMy/CFz/CoPz in 2D (n/a).
 * \param[in] nDim - Number of spatial dimensions (2 or 3).
 * \param[in] CosAlpha, SinAlpha, CosBeta, SinBeta - sin/cos of the angle of attack and sideslip,
 *            precomputed once by the caller (this is called once per monitored marker, redundantly
 *            by every thread, so recomputing the trig functions here would not be free).
 * \param[in] Force - Summed force components (size MAXNDIM).
 * \param[in] Moment - Summed moment components about Origin (size MAXNDIM).
 * \param[in] MomentX_Force, MomentY_Force, MomentZ_Force - Summed moment-of-force components
 *            about the coordinate axes, used for the center-of-pressure coordinates.
 * \return The wind-axis aerodynamic coefficients (CD, CL, CSF, CFx..CFz, CMx..CMz, CoPx..CoPz,
 *         CT, CQ). AeroCoeffsT is explicit at the call site (it can't be deduced, since it is
 *         only the return type): ComputeAeroCoeffsFromForceMoment<AeroCoeffs>(...).
 */
template <class AeroCoeffsT>
AeroCoeffsT ComputeAeroCoeffsFromForceMoment(unsigned short nDim, su2double CosAlpha, su2double SinAlpha,
                                              su2double CosBeta, su2double SinBeta, const su2double* Force,
                                              const su2double* Moment, const su2double* MomentX_Force,
                                              const su2double* MomentY_Force, const su2double* MomentZ_Force) {
  AeroCoeffsT c;

  if (nDim == 2) {
    c.CD = Force[0] * CosAlpha + Force[1] * SinAlpha;
    c.CL = -Force[0] * SinAlpha + Force[1] * CosAlpha;
    c.CMz = Moment[2];
    c.CoPx = MomentZ_Force[1];
    c.CoPy = -MomentZ_Force[0];
    c.CFx = Force[0];
    c.CFy = Force[1];
    c.CT = -c.CFx;
    c.CQ = -c.CMz;
  }
  if (nDim == 3) {
    c.CD = Force[0] * CosAlpha * CosBeta + Force[1] * SinBeta + Force[2] * SinAlpha * CosBeta;
    c.CL = -Force[0] * SinAlpha + Force[2] * CosAlpha;
    c.CSF = -Force[0] * SinBeta * CosAlpha + Force[1] * CosBeta - Force[2] * SinBeta * SinAlpha;
    c.CMx = Moment[0];
    c.CMy = Moment[1];
    c.CMz = Moment[2];
    c.CoPx = -MomentY_Force[0];
    c.CoPz = MomentY_Force[2];
    c.CFx = Force[0];
    c.CFy = Force[1];
    c.CFz = Force[2];
    c.CT = -c.CFz;
    c.CQ = -c.CMz;
  }

  return c;
}

/*!
 * \brief MPI-sum a single value across ranks (identity if not built with MPI).
 * \param[in] x - Value to reduce.
 * \return Sum of x over all ranks.
 */
inline su2double MPIReduceSum(su2double x) {
#ifdef HAVE_MPI
  su2double tmp = x;
  x = 0.0;
  SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
#endif
  return x;
}

/*!
 * \brief MPI-sum an array of per-monitoring-surface values across ranks, in place
 *        (no-op if not built with MPI).
 * \param[in,out] x - Array of size n to reduce in place.
 * \param[in] n - Number of entries in x.
 */
inline void MPIReduceSumInPlace(su2double* x, int n) {
#ifdef HAVE_MPI
  if (SU2_MPI::GetSize() == SINGLE_NODE) return;
  static vector<su2double> buffer;
  buffer.resize(n);
  SU2_MPI::Allreduce(x, buffer.data(), n, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  for (int i = 0; i < n; ++i) x[i] = buffer[i];
#endif
}

/*!
 * \brief MPI-reduce an AllBound/Surface aerodynamic coefficient pair across ranks (no-op if
 *        not built with MPI, or if the comm level does not require it).
 * \param[in] config - Definition of the particular problem.
 * \param[in,out] allBoundCoeff - Totals over all boundaries to reduce.
 * \param[in,out] surfaceCoeff - Per-monitoring-surface totals to reduce.
 */
template <class AeroCoeffsT, class AeroCoeffsArrayT>
void ReduceCoeffsMPI(const CConfig* config, AeroCoeffsT& allBoundCoeff, AeroCoeffsArrayT& surfaceCoeff) {
#ifdef HAVE_MPI
  if (config->GetComm_Level() != COMM_FULL) return;

  /*--- Add AllBound information using all the nodes ---*/

  allBoundCoeff.CD = MPIReduceSum(allBoundCoeff.CD);
  allBoundCoeff.CL = MPIReduceSum(allBoundCoeff.CL);
  allBoundCoeff.CSF = MPIReduceSum(allBoundCoeff.CSF);
  allBoundCoeff.CEff = allBoundCoeff.CL / (allBoundCoeff.CD + EPS);

  allBoundCoeff.CMx = MPIReduceSum(allBoundCoeff.CMx);
  allBoundCoeff.CMy = MPIReduceSum(allBoundCoeff.CMy);
  allBoundCoeff.CMz = MPIReduceSum(allBoundCoeff.CMz);

  allBoundCoeff.CoPx = MPIReduceSum(allBoundCoeff.CoPx);
  allBoundCoeff.CoPy = MPIReduceSum(allBoundCoeff.CoPy);
  allBoundCoeff.CoPz = MPIReduceSum(allBoundCoeff.CoPz);

  allBoundCoeff.CFx = MPIReduceSum(allBoundCoeff.CFx);
  allBoundCoeff.CFy = MPIReduceSum(allBoundCoeff.CFy);
  allBoundCoeff.CFz = MPIReduceSum(allBoundCoeff.CFz);

  allBoundCoeff.CT = MPIReduceSum(allBoundCoeff.CT);
  allBoundCoeff.CQ = MPIReduceSum(allBoundCoeff.CQ);
  allBoundCoeff.CMerit = allBoundCoeff.CT / (allBoundCoeff.CQ + EPS);

  /*--- Add the forces on the surfaces using all the nodes ---*/

  const int nMarkerMon = config->GetnMarker_Monitoring();

  MPIReduceSumInPlace(surfaceCoeff.CL, nMarkerMon);
  MPIReduceSumInPlace(surfaceCoeff.CD, nMarkerMon);
  MPIReduceSumInPlace(surfaceCoeff.CSF, nMarkerMon);

  for (int iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; iMarker_Monitoring++)
    surfaceCoeff.CEff[iMarker_Monitoring] =
        surfaceCoeff.CL[iMarker_Monitoring] / (surfaceCoeff.CD[iMarker_Monitoring] + EPS);

  MPIReduceSumInPlace(surfaceCoeff.CFx, nMarkerMon);
  MPIReduceSumInPlace(surfaceCoeff.CFy, nMarkerMon);
  MPIReduceSumInPlace(surfaceCoeff.CFz, nMarkerMon);

  MPIReduceSumInPlace(surfaceCoeff.CMx, nMarkerMon);
  MPIReduceSumInPlace(surfaceCoeff.CMy, nMarkerMon);
  MPIReduceSumInPlace(surfaceCoeff.CMz, nMarkerMon);
#endif
}

/*!
 * \brief Merge an AllBound/Surface aerodynamic coefficient pair into the running total/
 *        surfaceTotal grand totals.
 * \param[in] config - Definition of the particular problem.
 * \param[in] allBoundCoeff - Totals over all boundaries to merge in (already MPI-reduced).
 * \param[in] surfaceCoeff - Per-monitoring-surface totals to merge in (already MPI-reduced).
 * \param[in,out] total - Grand total to update (e.g. the solver's TotalCoeff).
 * \param[in,out] surfaceTotal - Per-surface grand total to update (e.g. the solver's SurfaceCoeff).
 * \param[in] overwrite - True to overwrite total/surfaceTotal (first contributor, i.e.
 *            Pressure_Forces, which also resets them to zero beforehand), false to add to them.
 */
template <class AeroCoeffsT, class AeroCoeffsArrayT>
void AccumulateTotalCoeffs(const CConfig* config, const AeroCoeffsT& allBoundCoeff,
                            const AeroCoeffsArrayT& surfaceCoeff, AeroCoeffsT& total, AeroCoeffsArrayT& surfaceTotal,
                            bool overwrite) {
  auto Update = [overwrite](su2double& dst, su2double src) {
    if (overwrite) dst = src;
    else dst += src;
  };

  Update(total.CD, allBoundCoeff.CD);
  Update(total.CL, allBoundCoeff.CL);
  Update(total.CSF, allBoundCoeff.CSF);
  total.CEff = total.CL / (total.CD + EPS);
  Update(total.CFx, allBoundCoeff.CFx);
  Update(total.CFy, allBoundCoeff.CFy);
  Update(total.CFz, allBoundCoeff.CFz);
  Update(total.CMx, allBoundCoeff.CMx);
  Update(total.CMy, allBoundCoeff.CMy);
  Update(total.CMz, allBoundCoeff.CMz);
  Update(total.CoPx, allBoundCoeff.CoPx);
  Update(total.CoPy, allBoundCoeff.CoPy);
  Update(total.CoPz, allBoundCoeff.CoPz);
  Update(total.CT, allBoundCoeff.CT);
  Update(total.CQ, allBoundCoeff.CQ);
  total.CMerit = total.CT / (total.CQ + EPS);

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring();
       iMarker_Monitoring++) {
    Update(surfaceTotal.CL[iMarker_Monitoring], surfaceCoeff.CL[iMarker_Monitoring]);
    Update(surfaceTotal.CD[iMarker_Monitoring], surfaceCoeff.CD[iMarker_Monitoring]);
    Update(surfaceTotal.CSF[iMarker_Monitoring], surfaceCoeff.CSF[iMarker_Monitoring]);
    surfaceTotal.CEff[iMarker_Monitoring] =
        surfaceTotal.CL[iMarker_Monitoring] / (surfaceTotal.CD[iMarker_Monitoring] + EPS);
    Update(surfaceTotal.CFx[iMarker_Monitoring], surfaceCoeff.CFx[iMarker_Monitoring]);
    Update(surfaceTotal.CFy[iMarker_Monitoring], surfaceCoeff.CFy[iMarker_Monitoring]);
    Update(surfaceTotal.CFz[iMarker_Monitoring], surfaceCoeff.CFz[iMarker_Monitoring]);
    Update(surfaceTotal.CMx[iMarker_Monitoring], surfaceCoeff.CMx[iMarker_Monitoring]);
    Update(surfaceTotal.CMy[iMarker_Monitoring], surfaceCoeff.CMy[iMarker_Monitoring]);
    Update(surfaceTotal.CMz[iMarker_Monitoring], surfaceCoeff.CMz[iMarker_Monitoring]);
  }
}

/*!
 * \brief Fold one vertex's Force/MomentDist/Coord into the running moment sums. Identical
 *        formulas are used by Pressure_Forces, Momentum_Forces and Friction_Forces, applied
 *        respectively to their inviscid, momentum and viscous force/moment sums.
 * \param[in] nDim - Number of spatial dimensions (2 or 3).
 * \param[in] RefLength - Reference length used to non-dimensionalize the moments.
 * \param[in] Force, MomentDist, Coord - This vertex's force, moment arm and position (size MAXNDIM).
 * \param[in,out] Moment - Running moment sum about Origin.
 * \param[in,out] MomentX_Force, MomentY_Force, MomentZ_Force - Running moment-of-force sums about
 *                the coordinate axes, used for the center-of-pressure coordinates.
 */
void AccumulateMoment(unsigned short nDim, su2double RefLength, const su2double* Force,
                      const su2double* MomentDist, const su2double* Coord, su2double* Moment,
                      su2double* MomentX_Force, su2double* MomentY_Force, su2double* MomentZ_Force) {
  if (nDim == 3) {
    Moment[0] += (Force[2] * MomentDist[1] - Force[1] * MomentDist[2]) / RefLength;
    MomentX_Force[1] += (-Force[1] * Coord[2]);
    MomentX_Force[2] += (Force[2] * Coord[1]);

    Moment[1] += (Force[0] * MomentDist[2] - Force[2] * MomentDist[0]) / RefLength;
    MomentY_Force[2] += (-Force[2] * Coord[0]);
    MomentY_Force[0] += (Force[0] * Coord[2]);
  }
  Moment[2] += (Force[1] * MomentDist[0] - Force[0] * MomentDist[1]) / RefLength;
  MomentZ_Force[0] += (-Force[0] * Coord[1]);
  MomentZ_Force[1] += (Force[1] * Coord[0]);
}

} // namespace

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::Pressure_Forces(const CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  const su2double Alpha = config->GetAoA() * PI_NUMBER / 180.0;
  const su2double Beta = config->GetAoS() * PI_NUMBER / 180.0;
  /*--- Precomputed once here since ComputeAeroCoeffsFromForceMoment is called once per monitored
   *    marker, redundantly by every thread (see below). ---*/
  const su2double CosAlpha = cos(Alpha), SinAlpha = sin(Alpha), CosBeta = cos(Beta), SinBeta = sin(Beta);
  const su2double RefArea = config->GetRefArea();
  const su2double RefLength = config->GetRefLength();
  auto Origin = config->GetRefOriginMoment(0);
  const bool axisymmetric = config->GetAxisymmetric();

  /*--- Variables initialization, and other writes to shared (possibly AD-active) state,
   *    are confined to the master thread, synchronized with a barrier so that the
   *    subsequent parallel loop over markers sees consistent, zeroed accumulators. ---*/

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    SetReferenceValues(*config);
    Total_CNearFieldOF = 0.0;
    Total_Heat = 0.0;
    Total_MaxHeat = 0.0;
    AllBound_CNearFieldOF_Inv = 0.0;
    /*--- AeroCoeffs::setZero is not parallel. ---*/
    AllBoundInvCoeff.setZero();
    TotalCoeff.setZero();
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  SurfaceInvCoeff.setZero();
  SurfaceCoeff.setZero();
  InvCoeff.setZero();

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) CNearFieldOF_Inv[iMarker] = 0.0;
  END_SU2_OMP_FOR

  const su2double factor = 1.0 / AeroCoeffForceRef;

  /*--- Reference pressure is always the far-field value. ---*/

  const su2double RefPressure = Pressure_Inf;

  /*--- Loop over the Euler and Navier-Stokes markers. Every thread runs this marker loop
   *    redundantly; only the per-vertex loop nested inside it is work-shared (SU2_OMP_FOR_STAT).
   *    Each thread accumulates its own partial force/moment sums and folds them into the shared
   *    per-marker/AllBound/Surface totals via AddCoeffContribution. CEff/CMerit are nonlinear
   *    ratios, so they are derived once after the loop, from the fully-reduced totals. ---*/

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    const auto Boundary = config->GetMarker_All_KindBC(iMarker);

    if (config->GetSolid_Wall(iMarker) || (Boundary == NEARFIELD_BOUNDARY) || IsMomentumBoundary(Boundary)) {

      /*--- Obtain the origin for the moment computation for a particular marker ---*/

      const auto Monitoring = config->GetMarker_All_Monitoring(iMarker);
      const int iMarker_Monitoring = FindMonitoringIndex(config, iMarker, Monitoring, Origin);

      su2double ForceInviscid[MAXNDIM] = {0.0}, MomentInviscid[MAXNDIM] = {0.0};
      su2double MomentX_Force[MAXNDIM] = {0.0}, MomentY_Force[MAXNDIM] = {0.0}, MomentZ_Force[MAXNDIM] = {0.0};

      su2double NFPressOF = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/

      SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
      for (unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        const su2double Pressure = nodes->GetPressure(iPoint);

        CPressure[iMarker][iVertex] = (Pressure - RefPressure) * factor * RefArea;

        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ((geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES)) {
          const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          const su2double* Coord = geometry->nodes->GetCoord(iPoint);

          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/

          NFPressOF += 0.5 * (Pressure - Pressure_Inf) * (Pressure - Pressure_Inf) * Normal[nDim - 1];

          su2double MomentDist[MAXNDIM] = {0.0};
          GeometryToolbox::Distance(nDim, Coord, Origin.data(), MomentDist);

          /*--- Axisymmetric simulations ---*/

          const su2double AxiFactor = axisymmetric ? su2double(2.0 * PI_NUMBER * geometry->nodes->GetCoord(iPoint, 1)) : su2double(1.0);

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/

          su2double Force[MAXNDIM] = {0.0};
          for (unsigned short iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf) * Normal[iDim] * factor * AxiFactor;
            ForceInviscid[iDim] += Force[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          AccumulateMoment(nDim, RefLength, Force, MomentDist, Coord, MomentInviscid, MomentX_Force, MomentY_Force,
                           MomentZ_Force);
        }
      }
      END_SU2_OMP_FOR

      if (Monitoring == YES) {
        if (Boundary != NEARFIELD_BOUNDARY) {
          const auto partial = ComputeAeroCoeffsFromForceMoment<AeroCoeffs>(
              nDim, CosAlpha, SinAlpha, CosBeta, SinBeta, ForceInviscid, MomentInviscid, MomentX_Force, MomentY_Force,
              MomentZ_Force);

          AddCoeffContribution(iMarker, iMarker_Monitoring, partial, InvCoeff, AllBoundInvCoeff, SurfaceInvCoeff);
        } else {
          /*--- At the Nearfield SU2 only cares about the pressure coeffient ---*/
          SU2_OMP_CRITICAL {
            CNearFieldOF_Inv[iMarker] += NFPressOF;
            AllBound_CNearFieldOF_Inv += NFPressOF;
          }
          END_SU2_OMP_CRITICAL
        }
      }
    }
  }
  /*--- For the SU2_NOWAIT in the vertex loop. ---*/
  SU2_OMP_BARRIER

  /*--- Derive the (nonlinear) per-marker, AllBound and Surface ratio coefficients from the
   *    now fully-reduced totals. This must happen once, after every thread has finished
   *    folding its partial contributions above (guaranteed by the barrier at the start of
   *    the safe-global-access section). ---*/

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    const auto Boundary = config->GetMarker_All_KindBC(iMarker);
    const auto Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if (Monitoring == YES && Boundary != NEARFIELD_BOUNDARY) {
      InvCoeff.CEff[iMarker] = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker] + EPS);
      InvCoeff.CMerit[iMarker] = InvCoeff.CT[iMarker] / (InvCoeff.CQ[iMarker] + EPS);
    }
  }
  END_SU2_OMP_FOR

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring();
       iMarker_Monitoring++) {
    SurfaceInvCoeff.CEff[iMarker_Monitoring] =
        SurfaceInvCoeff.CL[iMarker_Monitoring] / (SurfaceInvCoeff.CD[iMarker_Monitoring] + EPS);
  }
  END_SU2_OMP_FOR

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    AllBoundInvCoeff.CEff = AllBoundInvCoeff.CL / (AllBoundInvCoeff.CD + EPS);
    AllBoundInvCoeff.CMerit = AllBoundInvCoeff.CT / (AllBoundInvCoeff.CQ + EPS);

    ReduceCoeffsMPI(config, AllBoundInvCoeff, SurfaceInvCoeff);

    /*--- AllBound_CNearFieldOF_Inv, not covered by ReduceCoeffsMPI, is reduced separately. ---*/
    if (config->GetComm_Level() == COMM_FULL) {
      AllBound_CNearFieldOF_Inv = MPIReduceSum(AllBound_CNearFieldOF_Inv);
    }

    AccumulateTotalCoeffs(config, AllBoundInvCoeff, SurfaceInvCoeff, TotalCoeff, SurfaceCoeff, /*overwrite=*/true);
    Total_CNearFieldOF = AllBound_CNearFieldOF_Inv;
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::Momentum_Forces(const CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED

  const su2double Alpha = config->GetAoA() * PI_NUMBER / 180.0;
  const su2double Beta = config->GetAoS() * PI_NUMBER / 180.0;
  const su2double CosAlpha = cos(Alpha), SinAlpha = sin(Alpha), CosBeta = cos(Beta), SinBeta = sin(Beta);
  const su2double RefLength = config->GetRefLength();
  auto Origin = config->GetRefOriginMoment(0);
  const bool axisymmetric = config->GetAxisymmetric();

  const su2double factor = 1.0 / AeroCoeffForceRef;

  SU2_OMP_SAFE_GLOBAL_ACCESS(AllBoundMntCoeff.setZero();)
  SurfaceMntCoeff.setZero();
  MntCoeff.setZero();

  /*--- Loop over the Inlet-Outlet Markers (see Pressure_Forces for how the parallel
   *    reduction over threads is organized). ---*/

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    const auto Boundary = config->GetMarker_All_KindBC(iMarker);
    if (!IsMomentumBoundary(Boundary)) continue;

    const auto Monitoring = config->GetMarker_All_Monitoring(iMarker);

    /*--- Obtain the origin for the moment computation for a particular marker ---*/

    const int iMarker_Monitoring = FindMonitoringIndex(config, iMarker, Monitoring, Origin);

    su2double ForceMomentum[MAXNDIM] = {0.0}, MomentMomentum[MAXNDIM] = {0.0};
    su2double MomentX_Force[3] = {0.0}, MomentY_Force[3] = {0.0}, MomentZ_Force[3] = {0.0};

    /*--- Loop over the vertices to compute the forces (work-shared across threads, see
      *    Pressure_Forces for why the chunk size is computed and the barrier skipped). ---*/

    SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
    for (unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      /*--- Note that the pressure coefficient is computed at the
        halo cells (for visualization purposes), but not the forces ---*/

      if ((geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES)) {
        const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        const su2double* Coord = geometry->nodes->GetCoord(iPoint);
        const su2double Density = nodes->GetDensity(iPoint);
        su2double MassFlow = 0.0;
        su2double Velocity[MAXNDIM] = {0.0}, MomentDist[MAXNDIM] = {0.0};
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = nodes->GetVelocity(iPoint, iDim);
          MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          MassFlow -= Normal[iDim] * Velocity[iDim] * Density;
        }

        /*--- Axisymmetric simulations ---*/

        const su2double AxiFactor = axisymmetric ? su2double(2.0 * PI_NUMBER * geometry->nodes->GetCoord(iPoint, 1)) : su2double(1.0);

        /*--- Force computation, note the minus sign due to the
          orientation of the normal (outward) ---*/

        su2double Force[MAXNDIM] = {0.0};
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Force[iDim] = MassFlow * Velocity[iDim] * factor * AxiFactor;
          ForceMomentum[iDim] += Force[iDim];
        }

        /*--- Moment with respect to the reference axis ---*/

        AccumulateMoment(nDim, RefLength, Force, MomentDist, Coord, MomentMomentum, MomentX_Force, MomentY_Force,
                         MomentZ_Force);
      }
    }
    END_SU2_OMP_FOR

    if (Monitoring == YES) {
      const auto partial = ComputeAeroCoeffsFromForceMoment<AeroCoeffs>(
          nDim, CosAlpha, SinAlpha, CosBeta, SinBeta, ForceMomentum, MomentMomentum, MomentX_Force, MomentY_Force,
          MomentZ_Force);

      AddCoeffContribution(iMarker, iMarker_Monitoring, partial, MntCoeff, AllBoundMntCoeff, SurfaceMntCoeff);
    }
  }
  /*--- For the SU2_NOWAIT in the vertex loop. ---*/
  SU2_OMP_BARRIER

  /*--- Derive the ratio coefficients from the fully-reduced totals, once. ---*/

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    const auto Boundary = config->GetMarker_All_KindBC(iMarker);
    const auto Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if (Monitoring == YES && IsMomentumBoundary(Boundary)) {
      MntCoeff.CEff[iMarker] = MntCoeff.CL[iMarker] / (MntCoeff.CD[iMarker] + EPS);
      MntCoeff.CMerit[iMarker] = MntCoeff.CT[iMarker] / (MntCoeff.CQ[iMarker] + EPS);
    }
  }
  END_SU2_OMP_FOR

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring();
       iMarker_Monitoring++) {
    SurfaceMntCoeff.CEff[iMarker_Monitoring] =
        SurfaceMntCoeff.CL[iMarker_Monitoring] / (SurfaceMntCoeff.CD[iMarker_Monitoring] + EPS);
  }
  END_SU2_OMP_FOR

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    AllBoundMntCoeff.CEff = AllBoundMntCoeff.CL / (AllBoundMntCoeff.CD + EPS);
    AllBoundMntCoeff.CMerit = AllBoundMntCoeff.CT / (AllBoundMntCoeff.CQ + EPS);

    ReduceCoeffsMPI(config, AllBoundMntCoeff, SurfaceMntCoeff);
    AccumulateTotalCoeffs(config, AllBoundMntCoeff, SurfaceMntCoeff, TotalCoeff, SurfaceCoeff, /*overwrite=*/false);
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::Friction_Forces(const CGeometry* geometry, const CConfig* config) {
  SU2_ZONE_SCOPED
  if (!config->GetViscous()) return;

  constexpr int MaxNorm = 8;
  const su2double minYPlus = config->GetwallModel_MinYPlus();

  const su2double Alpha = config->GetAoA() * PI_NUMBER / 180.0;
  const su2double Beta = config->GetAoS() * PI_NUMBER / 180.0;
  const su2double CosAlpha = cos(Alpha), SinAlpha = sin(Alpha), CosBeta = cos(Beta), SinBeta = sin(Beta);
  const su2double RefLength = config->GetRefLength();
  const su2double RefHeatFlux = config->GetHeat_Flux_Ref();
  const su2double RefTemperature = config->GetTemperature_Ref();
  auto Origin = config->GetRefOriginMoment(0);

  const bool energy = config->GetEnergy_Equation();
  const bool QCR = config->GetSAParsedOptions().qcr2000;
  const bool axisymmetric = config->GetAxisymmetric();
  const bool roughwall = (config->GetnRoughWall() > 0);
  const bool nemo = config->GetNEMOProblem();

  const su2double factor = 1.0 / AeroCoeffForceRef;
  const su2double factorFric = config->GetRefArea() * factor;

  /*--- Variables initialization ---*/

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    AllBound_HF_Visc = 0.0;
    AllBound_MaxHF_Visc = 0.0;
    AllBoundViscCoeff.setZero();
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  SurfaceViscCoeff.setZero();

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring();
       iMarker_Monitoring++) {
    Surface_HF_Visc[iMarker_Monitoring] = 0.0;
    Surface_MaxHF_Visc[iMarker_Monitoring] = 0.0;
  }
  END_SU2_OMP_FOR

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (!config->GetViscous_Wall(iMarker)) continue;
    ViscCoeff.setZero(iMarker);
    HF_Visc[iMarker] = 0.0;
    MaxHF_Visc[iMarker] = 0.0;
  }
  END_SU2_OMP_FOR

  /*--- Loop over the Navier-Stokes markers (see Pressure_Forces for how the parallel
   *    reduction over threads is organized). The per-vertex loop below is the expensive
   *    part (stress-tensor and heat-flux evaluations) and is work-shared across threads. ---*/

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {

    if (!config->GetViscous_Wall(iMarker)) continue;
    const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    const auto Monitoring = config->GetMarker_All_Monitoring(iMarker);
    const bool py_custom = config->GetMarker_All_PyCustom(iMarker);

    /*--- Obtain the origin for the moment computation for a particular marker ---*/

    const int iMarker_Monitoring = FindMonitoringIndex(config, iMarker, Monitoring, Origin);

    su2double ForceViscous[MAXNDIM] = {0.0}, MomentViscous[MAXNDIM] = {0.0};
    su2double MomentX_Force[MAXNDIM] = {0.0}, MomentY_Force[MAXNDIM] = {0.0}, MomentZ_Force[MAXNDIM] = {0.0};
    su2double HF_Visc_Local = 0.0, MaxHF_Visc_Local = 0.0;

    /* --- check if wall functions are used --- */

    const bool wallfunctions = (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE);

    /*--- Marker-level lookups hoisted out of the per-vertex loop below: GetWallRoughnessProperties,
     *    GetWall_HeatFlux, GetIsothermal_Temperature and GetCatalytic_Wall all scan over marker
     *    lists (some by string comparison), so evaluating them once per vertex instead of once per
     *    marker was a real cost for markers with many vertices. ---*/

    WALL_TYPE WallType = WALL_TYPE::SMOOTH;
    if (roughwall) {
      su2double Roughness_Height;
      tie(WallType, Roughness_Height) = config->GetWallRoughnessProperties(Marker_Tag);
    }

    const auto KindBC = config->GetMarker_All_KindBC(iMarker);
    su2double Wall_HeatFlux_Value = 0.0, Twall = 0.0;
    if (!nemo && !py_custom) {
      if (KindBC == BC_TYPE::HEAT_FLUX) {
        Wall_HeatFlux_Value = -config->GetWall_HeatFlux(Marker_Tag);
        if (config->GetIntegrated_HeatFlux()) Wall_HeatFlux_Value /= geometry->GetSurfaceArea(config, iMarker);
      } else if (KindBC == BC_TYPE::ISOTHERMAL) {
        Twall = config->GetIsothermal_Temperature(Marker_Tag) / RefTemperature;
      }
    }
    const bool catalytic = nemo && config->GetCatalytic_Wall(iMarker);

    /*--- Loop over the vertices to compute the forces (work-shared across threads, see
     *    Pressure_Forces for why the chunk size is computed and the barrier skipped). ---*/

    SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
    for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      const su2double* Coord = geometry->nodes->GetCoord(iPoint);

      const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

      /*--- One view covering the whole velocity-gradient block, instead of nDim*nDim separate
       *    single-element lookups (ComputeStressTensor/AddQCR accept any [][]-indexable type). ---*/
      const auto Grad_Vel = nodes->GetVelocityGradient(iPoint);

      su2double Grad_Temp[3] = {0.0}, Grad_Temp_ve[3] = {0.0};
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        Grad_Temp[iDim] = nodes->GetGradient_Primitive(iPoint, prim_idx.Temperature(), iDim);
        if (nemo) Grad_Temp_ve[iDim] = nodes->GetGradient_Primitive(iPoint, prim_idx.Temperature_ve(), iDim);
      }

      su2double Viscosity = nodes->GetLaminarViscosity(iPoint);
      su2double EddyViscosity = 0.0;
      if (WallType == WALL_TYPE::ROUGH) {
        EddyViscosity = nodes->GetEddyViscosity(iPoint);
        Viscosity += EddyViscosity;
      }
      const su2double Density = nodes->GetDensity(iPoint);

      const su2double Area = GeometryToolbox::Norm(nDim, Normal);
      su2double UnitNormal[3] = {0.0};
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        UnitNormal[iDim] = Normal[iDim] / Area;
      }

      /*--- Evaluate Tau ---*/
      su2double Tau[3][3] = {{0.0}};
      CNumerics::ComputeStressTensor(nDim, Tau, Grad_Vel, Viscosity);

      /*--- If necessary evaluate the QCR contribution to Tau ---*/

      if (QCR) CNumerics::AddQCR(nDim, Grad_Vel, Tau, EddyViscosity / Viscosity);

      /*--- Project Tau in each surface element ---*/

      su2double TauElem[3] = {0.0};
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        for (unsigned short jDim = 0; jDim < nDim; jDim++) {
          TauElem[iDim] += Tau[iDim][jDim] * UnitNormal[jDim];
        }
      }

      /*--- Compute wall shear stress (using the stress tensor). Compute wall skin friction coefficient, and heat flux
       * on the wall ---*/

      su2double TauTangent[MAXNDIM] = {0.0};
      GeometryToolbox::TangentProjection(nDim, Tau, UnitNormal, TauTangent);

      WallShearStress[iMarker][iVertex] = GeometryToolbox::Norm(int(MAXNDIM), TauTangent);

      /*--- For wall functions, the wall stresses need to be scaled by the wallfunction stress Tau_Wall---*/
      if (wallfunctions && (YPlus[iMarker][iVertex] > minYPlus)){
        const su2double Tau_Wall = nodes->GetTau_Wall(iPoint);
        const su2double scale = Tau_Wall / WallShearStress[iMarker][iVertex];
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          TauTangent[iDim] *= scale;
          TauElem[iDim] *= scale;
        }

        WallShearStress[iMarker][iVertex] = Tau_Wall;
      }

      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        CSkinFriction[iMarker](iVertex,iDim) = TauTangent[iDim] * factorFric;
      }

      /*--- Compute non-dimensional velocity and y+ ---*/

      const su2double FrictionVel = sqrt(fabs(WallShearStress[iMarker][iVertex]) / Density);

      if (!wallfunctions && MGLevel == MESH_0 && geometry->nodes->GetDomain(iPoint)) {
        // for CMultiGridGeometry and halos, the nearest neighbor distance is not set
        const su2double WallDistMod = geometry->vertex[iMarker][iVertex]->GetNearestNeighborDistance();
        YPlus[iMarker][iVertex] = WallDistMod * FrictionVel * Density / Viscosity;
      }

      /*--- Compute total and maximum heat flux on the wall ---*/

      if (!nemo) {
        su2double thermal_conductivity = 0.0;
        if ((FlowRegime == ENUM_REGIME::COMPRESSIBLE) || (FlowRegime == ENUM_REGIME::INCOMPRESSIBLE)) {
          thermal_conductivity = nodes->GetThermalConductivity(iPoint);
        }

        if (KindBC == BC_TYPE::HEAT_FLUX) {
          HeatFlux[iMarker][iVertex] =
              py_custom ? -geometry->GetCustomBoundaryHeatFlux(iMarker, iVertex) : Wall_HeatFlux_Value;
        } else if (KindBC == BC_TYPE::ISOTHERMAL) {
          const su2double Twall_local =
              py_custom ? geometry->GetCustomBoundaryTemperature(iMarker, iVertex) / RefTemperature : Twall;
          const auto iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
          const su2double* Coord_Normal = geometry->nodes->GetCoord(iPointNormal);
          const su2double dist_ij = GeometryToolbox::NormalDistance(nDim, UnitNormal, Coord, Coord_Normal);
          const su2double There = nodes->GetTemperature(iPointNormal);
          HeatFlux[iMarker][iVertex] = thermal_conductivity * (There - Twall_local) / dist_ij * RefHeatFlux;
        } else {
          su2double dTdn = GeometryToolbox::DotProduct(nDim, Grad_Temp, UnitNormal);
          if (FlowRegime == ENUM_REGIME::INCOMPRESSIBLE && !energy) dTdn = 0.0;
          HeatFlux[iMarker][iVertex] = thermal_conductivity * dTdn * RefHeatFlux;
        }
      } else {

        const auto& thermal_conductivity_tr = nodes->GetThermalConductivity(iPoint);
        const auto& thermal_conductivity_ve = nodes->GetThermalConductivity_ve(iPoint);

        const su2double dTdn = -GeometryToolbox::DotProduct(nDim, Grad_Temp, UnitNormal);
        const su2double dTvedn = -GeometryToolbox::DotProduct(nDim, Grad_Temp_ve, UnitNormal);

        /*--- Surface energy balance: trans-rot heat flux, vib-el heat flux ---*/
        HeatFlux[iMarker][iVertex] = -(thermal_conductivity_tr*dTdn + thermal_conductivity_ve*dTvedn);

        /*--- Compute enthalpy transport to surface due to mass diffusion ---*/
        if (catalytic) {

          const auto nSpecies = config->GetnSpecies();
          const auto& Grad_PrimVar = nodes->GetGradient_Primitive(iPoint);
          const auto& PrimVar = nodes->GetPrimitive(iPoint);
          const auto& Ds = nodes->GetDiffusionCoeff(iPoint);
          const auto& hs = nodes->GetEnthalpys(iPoint);
          const su2double rho = PrimVar[prim_idx.Density()];

          su2double sumJhs = 0.0;
          for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++) {
            for (auto iDim = 0u; iDim < nDim; iDim++) {
              su2double dYdn = 1.0/rho*(Grad_PrimVar[iSpecies][iDim] - PrimVar[iSpecies]*Grad_PrimVar[prim_idx.Density()][iDim]/rho);
              sumJhs += rho*Ds[iSpecies]*hs[iSpecies]*dYdn*UnitNormal[iDim];
            }
          }
          /*--- Surface energy balance: mass diffusion ---*/
          HeatFlux[iMarker][iVertex] += sumJhs;

        }
      }

      /*--- Note that heat is computed at the
       halo cells (for visualization purposes), but not the forces ---*/

      if ((geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES)) {
        /*--- Axisymmetric simulations ---*/

        const su2double AxiFactor = axisymmetric ? su2double(2.0 * PI_NUMBER * geometry->nodes->GetCoord(iPoint, 1)) : su2double(1.0);

        /*--- Force computation ---*/

        su2double Force[MAXNDIM] = {0.0}, MomentDist[MAXNDIM] = {0.0};
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Force[iDim] = TauElem[iDim] * Area * factor * AxiFactor;
          ForceViscous[iDim] += Force[iDim];
          MomentDist[iDim] = Coord[iDim] - Origin[iDim];
        }

        /*--- Moment with respect to the reference axis ---*/

        AccumulateMoment(nDim, RefLength, Force, MomentDist, Coord, MomentViscous, MomentX_Force, MomentY_Force,
                         MomentZ_Force);

        HF_Visc_Local += HeatFlux[iMarker][iVertex] * Area;
        MaxHF_Visc_Local += pow(HeatFlux[iMarker][iVertex], MaxNorm);
      }
    }
    END_SU2_OMP_FOR

    /*--- MaxHF_Visc_Local (and the shared accumulators it feeds below) are left un-rooted,
     *    i.e. still a raw sum of HeatFlux^MaxNorm; pow(., 1/MaxNorm) is taken once at the end,
     *    which is equivalent since pow(x^(1/n), n) == x. ---*/

    if (Monitoring == YES) {
      const auto partial = ComputeAeroCoeffsFromForceMoment<AeroCoeffs>(
          nDim, CosAlpha, SinAlpha, CosBeta, SinBeta, ForceViscous, MomentViscous, MomentX_Force, MomentY_Force,
          MomentZ_Force);

      AddCoeffContribution(iMarker, iMarker_Monitoring, partial, ViscCoeff, AllBoundViscCoeff, SurfaceViscCoeff);

      /*--- Heat flux, not covered by AddCoeffContribution, is folded in its own critical section. ---*/

      SU2_OMP_CRITICAL {
        HF_Visc[iMarker] += HF_Visc_Local;
        AllBound_HF_Visc += HF_Visc_Local;
        MaxHF_Visc[iMarker] += MaxHF_Visc_Local;
        AllBound_MaxHF_Visc += MaxHF_Visc_Local;

        if (iMarker_Monitoring >= 0) {
          Surface_HF_Visc[iMarker_Monitoring] += HF_Visc_Local;
          Surface_MaxHF_Visc[iMarker_Monitoring] += MaxHF_Visc_Local;
        }
      }
      END_SU2_OMP_CRITICAL
    }
  }
  /*--- For the SU2_NOWAIT in the vertex loop. ---*/
  SU2_OMP_BARRIER

  /*--- Derive the ratio coefficients, and root the (still raw) per-marker maximum heat flux,
   *    from the fully-reduced totals, once. Surface_MaxHF_Visc and AllBound_MaxHF_Visc are
   *    rooted later below, after the MPI reduction. ---*/

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (!config->GetViscous_Wall(iMarker)) continue;
    if (config->GetMarker_All_Monitoring(iMarker) == YES) {
      ViscCoeff.CEff[iMarker] = ViscCoeff.CL[iMarker] / (ViscCoeff.CD[iMarker] + EPS);
      ViscCoeff.CMerit[iMarker] = ViscCoeff.CT[iMarker] / (ViscCoeff.CQ[iMarker] + EPS);
      MaxHF_Visc[iMarker] = pow(MaxHF_Visc[iMarker], 1.0 / MaxNorm);
    }
  }
  END_SU2_OMP_FOR

  SU2_OMP_FOR_(schedule(static, OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring();
       iMarker_Monitoring++) {
    SurfaceViscCoeff.CEff[iMarker_Monitoring] =
        SurfaceViscCoeff.CL[iMarker_Monitoring] / (SurfaceViscCoeff.CD[iMarker_Monitoring] + EPS);
  }
  END_SU2_OMP_FOR

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    ReduceCoeffsMPI(config, AllBoundViscCoeff, SurfaceViscCoeff);

    /*--- HF_Visc/MaxHF_Visc, not covered by ReduceCoeffsMPI, are reduced separately. ---*/
    if (config->GetComm_Level() == COMM_FULL) {
      AllBound_HF_Visc = MPIReduceSum(AllBound_HF_Visc);
      AllBound_MaxHF_Visc = MPIReduceSum(AllBound_MaxHF_Visc);

      const int nMarkerMon = config->GetnMarker_Monitoring();
      MPIReduceSumInPlace(Surface_HF_Visc.data(), nMarkerMon);
      MPIReduceSumInPlace(Surface_MaxHF_Visc.data(), nMarkerMon);
    }

    /*--- Complete the calculation of maximum heat flux. ---*/

    for (auto& hf : Surface_MaxHF_Visc) {
      hf = pow(hf, 1.0 / MaxNorm);
    }
    AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0 / MaxNorm);

    AccumulateTotalCoeffs(config, AllBoundViscCoeff, SurfaceViscCoeff, TotalCoeff, SurfaceCoeff, /*overwrite=*/false);
    Total_Heat = AllBound_HF_Visc;
    Total_MaxHeat = AllBound_MaxHF_Visc;

    /*--- Buffet_Monitoring is not thread-safe, hence confined to the master thread. ---*/

    Buffet_Monitoring(geometry, config);
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

template<class V, ENUM_REGIME R>
su2double CFVMFlowSolverBase<V,R>::EvaluateCommonObjFunc(const CConfig& config) const {

  su2double objFun = 0.0;

  /*--- Loop over all monitored markers, add to the 'combo' objective ---*/

  for (auto iMarker = 0u; iMarker < config.GetnMarker_Monitoring(); iMarker++) {

    const auto weight = config.GetWeight_ObjFunc(iMarker);

    switch (config.GetKind_ObjFunc(iMarker)) {
      case DRAG_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CD[iMarker];
        break;
      case LIFT_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CL[iMarker];
        break;
      case SIDEFORCE_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CSF[iMarker];
        break;
      case MOMENT_X_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CMx[iMarker];
        break;
      case MOMENT_Y_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CMy[iMarker];
        break;
      case MOMENT_Z_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CMz[iMarker];
        break;
      case FORCE_X_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CFx[iMarker];
        break;
      case FORCE_Y_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CFy[iMarker];
        break;
      case FORCE_Z_COEFFICIENT:
        objFun += weight * SurfaceCoeff.CFz[iMarker];
        break;
      case TOTAL_HEATFLUX:
        objFun += weight * Surface_HF_Visc[iMarker];
        break;
      case MAXIMUM_HEATFLUX:
        objFun += weight * Surface_MaxHF_Visc[iMarker];
        break;
      default:
        break;
    }
  }

  /*--- The following are not per-surface, and so to avoid that they are
   double-counted when multiple surfaces are specified, they have been
   placed outside of the loop above. In addition, multi-objective mode is
   also disabled for these objective functions (error thrown at start). ---*/

  const auto weight = config.GetWeight_ObjFunc(0);

  switch (config.GetKind_ObjFunc(0)) {
    case EFFICIENCY:
      objFun += weight * TotalCoeff.CEff;
      break;
    case INVERSE_DESIGN_PRESSURE:
      objFun += weight * Total_CpDiff;
      break;
    case INVERSE_DESIGN_HEATFLUX:
      objFun += weight * Total_HeatFluxDiff;
      break;
    case EQUIVALENT_AREA:
      objFun += weight*Total_CEquivArea;
      break;
    case THRUST_COEFFICIENT:
      objFun += weight * TotalCoeff.CT;
      break;
    case TORQUE_COEFFICIENT:
      objFun += weight * TotalCoeff.CQ;
      break;
    case FIGURE_OF_MERIT:
      objFun += weight * TotalCoeff.CMerit;
      break;
    case SURFACE_TOTAL_PRESSURE:
      objFun += weight * config.GetSurface_TotalPressure(0);
      break;
    case SURFACE_STATIC_PRESSURE:
      objFun += weight * config.GetSurface_Pressure(0);
      break;
    case SURFACE_STATIC_TEMPERATURE:
      objFun += weight * config.GetSurface_Temperature(0);
      break;
    case SURFACE_MASSFLOW:
      objFun += weight * config.GetSurface_MassFlow(0);
      break;
    case SURFACE_UNIFORMITY:
      objFun += weight * config.GetSurface_Uniformity(0);
      break;
    case SURFACE_SECONDARY:
      objFun += weight * config.GetSurface_SecondaryStrength(0);
      break;
    case SURFACE_MOM_DISTORTION:
      objFun += weight * config.GetSurface_MomentumDistortion(0);
      break;
    case SURFACE_SECOND_OVER_UNIFORM:
      objFun += weight * config.GetSurface_SecondOverUniform(0);
      break;
    case SURFACE_PRESSURE_DROP:
      objFun += weight * config.GetSurface_PressureDrop(0);
      break;
    case SURFACE_SPECIES_0:
      objFun += weight * config.GetSurface_Species_0(0);
      break;
    case SURFACE_SPECIES_VARIANCE:
      objFun += weight * config.GetSurface_Species_Variance(0);
      break;
    case CUSTOM_OBJFUNC:
      objFun += weight * Total_Custom_ObjFunc;
      break;
    default:
      break;
  }

  return objFun;
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::ComputeAxisymmetricAuxGradients(CGeometry *geometry, const CConfig* config) {

  /*--- Loop through all points to set the auxvargrad --*/
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    su2double yCoord          = geometry->nodes->GetCoord(iPoint, 1);
    su2double yVelocity       = nodes->GetVelocity(iPoint,1);
    su2double xVelocity       = nodes->GetVelocity(iPoint,0);
    su2double Total_Viscosity = nodes->GetLaminarViscosity(iPoint) + nodes->GetEddyViscosity(iPoint);

    if (yCoord > EPS){
      su2double nu_v_on_y = Total_Viscosity*yVelocity/yCoord;
      nodes->SetAuxVar(iPoint, 0, nu_v_on_y);
      nodes->SetAuxVar(iPoint, 1, nu_v_on_y*yVelocity);
      nodes->SetAuxVar(iPoint, 2, nu_v_on_y*xVelocity);
    }
  }
  END_SU2_OMP_FOR

  /*--- Compute the auxiliary variable gradient with GG or WLS. ---*/
  if (config->GreenGaussGradientMethod()) {
    SetAuxVar_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetAuxVar_Gradient_LS(geometry, config);
  }
}

template <class V, ENUM_REGIME FlowRegime>
void CFVMFlowSolverBase<V, FlowRegime>::MultigridProjectEulerWall(CGeometry* geometry, const CConfig* config,
                                                                   bool use_solution_old) {
  const auto iVel = prim_idx.Velocity();
  const auto nDim = geometry->GetnDim();

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != EULER_WALL) continue;

    SU2_OMP_FOR_STAT(32)
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      if (!geometry->nodes->GetDomain(iPoint)) continue;

      /*--- Use the Gram-Schmidt corrected normal for nodes on intersecting walls,
       *    consistent with BC_Sym_Plane.  Fall back to the raw marker normal otherwise. ---*/
      su2double UnitNormal[MAXNDIM] = {0.0};
      const auto it = geometry->symmetryNormals[iMarker].find(iVertex);
      if (it != geometry->symmetryNormals[iMarker].end()) {
        for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = it->second[iDim];
      } else {
        su2double Normal[MAXNDIM] = {0.0};
        geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
        const su2double Area = GeometryToolbox::Norm(nDim, Normal);
        if (Area < EPS) continue;
        for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = Normal[iDim] / Area;
      }

      su2double* sol = use_solution_old ? nodes->GetSolution_Old(iPoint) : nodes->GetSolution(iPoint);

      /*--- Compute normal component of the velocity / momentum vector.
       *    For dynamic grids subtract the grid velocity to enforce (v - v_grid).n = 0,
       *    multiplying by density for compressible flow (conservative variables). ---*/
      su2double gridVel[MAXNDIM] = {};
      if (dynamic_grid && !use_solution_old) {
        for (auto iDim = 0u; iDim < nDim; iDim++)
          gridVel[iDim] = geometry->nodes->GetGridVel(iPoint)[iDim];
        if constexpr (FlowRegime == ENUM_REGIME::COMPRESSIBLE) {
          for (auto iDim = 0u; iDim < nDim; iDim++)
            gridVel[iDim] *= nodes->GetDensity(iPoint);
        }
      }

      su2double momentum_n = 0.0;
      for (auto iDim = 0u; iDim < nDim; iDim++)
        momentum_n += (sol[iVel + iDim] - gridVel[iDim]) * UnitNormal[iDim];

      /*--- Project to tangent plane. ---*/
      for (auto iDim = 0u; iDim < nDim; iDim++) sol[iVel + iDim] -= momentum_n * UnitNormal[iDim];
    }
    END_SU2_OMP_FOR
  }
}
