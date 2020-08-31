/*!
 * \file CTurbSolver.cpp
 * \brief Main subrotuines of CTurbSolver class
 * \author F. Palacios, A. Bueno
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/solvers/CTurbSolver.hpp"
#include "../../include/solvers/CTurbSSTSolver.hpp"
#include "../../../Common/include/omp_structure.hpp"
#include "../../include/limiters/computeLimiters.hpp"


CTurbSolver::CTurbSolver(void) : CSolver() { }

CTurbSolver::CTurbSolver(CGeometry* geometry, CConfig *config) : CSolver() {

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  nMarker = config->GetnMarker_All();

  /*--- Store the number of vertices on each marker for deallocation later ---*/
  nVertex = new unsigned long[nMarker];
  for (auto iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Allocates a 2D array with variable "outer" sizes and init to 0. ---*/

  auto Alloc2D = [](unsigned long M, const unsigned long* N, su2double**& X) {
    X = new su2double* [M];
    for(unsigned long i = 0; i < M; ++i)
      X[i] = new su2double [N[i]] ();
  };

  /*--- Allocates a 3D array with variable "middle" sizes and init to 0. ---*/

  auto Alloc3D = [](unsigned long M, const unsigned long* N, unsigned long P, su2double***& X) {
    X = new su2double** [M];
    for(unsigned long i = 0; i < M; ++i) {
      X[i] = new su2double* [N[i]];
      for(unsigned long j = 0; j < N[i]; ++j)
        X[i][j] = new su2double [P] ();
    }
  };

  /*--- Store the value of the gradient basis function at the boundaries ---*/

  if (config->GetMUSCL_Flow() || config->GetMUSCL_Turb())
    Alloc3D(nMarker, nVertex, nVar, GradBasis_Aux);
  if (config->GetUse_Accurate_Turb_Jacobians())
    Alloc3D(nMarker, nVertex, nVar, GradBasis);

#ifdef HAVE_OMP
  /*--- Get the edge coloring, see notes in CEulerSolver's constructor. ---*/
  su2double parallelEff = 1.0;
  const auto& coloring = geometry->GetEdgeColoring(&parallelEff);

  ReducerStrategy = parallelEff < COLORING_EFF_THRESH;

  if (ReducerStrategy && (coloring.getOuterSize()>1))
    geometry->SetNaturalEdgeColoring();

  if (!coloring.empty()) {
    auto groupSize = ReducerStrategy? 1ul : geometry->GetEdgeColorGroupSize();
    auto nColor = coloring.getOuterSize();
    EdgeColoring.reserve(nColor);

    for(auto iColor = 0ul; iColor < nColor; ++iColor)
      EdgeColoring.emplace_back(coloring.innerIdx(iColor), coloring.getNumNonZeros(iColor), groupSize);
  }

  nPoint = geometry->GetnPoint();
  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), OMP_MAX_SIZE);
#else
  EdgeColoring[0] = DummyGridColor<>(geometry->GetnEdge());
#endif
}

CTurbSolver::~CTurbSolver(void) {

  if (Inlet_TurbVars != nullptr) {
    for (auto iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Inlet_TurbVars[iMarker] != nullptr) {
        for (auto iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          delete [] Inlet_TurbVars[iMarker][iVertex];
        }
        delete [] Inlet_TurbVars[iMarker];
      }
    }
    delete [] Inlet_TurbVars;
  }

  if (GradBasis != NULL) {
    for (auto iMarker = 0; iMarker < nMarker; iMarker++) {
      for (auto iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] GradBasis[iMarker][iVertex];
      delete [] GradBasis[iMarker];
    }
    delete [] GradBasis;
  }

  if (GradBasis_Aux != NULL) {
    for (auto iMarker = 0; iMarker < nMarker; iMarker++) {
      for (auto iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] GradBasis_Aux[iMarker][iVertex];
      delete [] GradBasis_Aux[iMarker];
    }
    delete [] GradBasis_Aux;
  }

  delete nodes;
}

void CTurbSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver,
                                  CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool muscl = config->GetMUSCL_Turb();
  const bool limiter = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER);
  const bool sst = ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST));

  /*--- Only reconstruct flow variables if MUSCL is on for flow (requires upwind) and turbulence. ---*/
  const bool musclFlow = (config->GetMUSCL_Flow()) && muscl && (config->GetKind_ConvNumScheme_Flow() == SPACE_UPWIND);
  /*--- Only consider flow limiters for cell-based limiters, edge-based would need to be recomputed. ---*/
  const bool limiterFlow = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER);

  CVariable* flowNodes = solver[FLOW_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Static arrays of MUSCL-reconstructed flow primitives and turbulence variables (thread safety). ---*/
  su2double solution_i[MAXNVAR] = {0.0}, flowPrimVar_i[MAXNVARFLOW] = {0.0};
  su2double solution_j[MAXNVAR] = {0.0}, flowPrimVar_j[MAXNVARFLOW] = {0.0};

  su2double ZeroVec[2] = {0.0}, OneVec[2]  = {0.0};

  su2double GradBasis_i[2] = {0.0}, GradBasis_j[2] = {0.0};

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring)
  {
  /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
  SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
  for(auto k = 0ul; k < color.size; ++k) {

    auto iEdge = color.indices[k];

    unsigned short iDim, iVar;

    /*--- Points in edge and normal vectors ---*/

    auto iPoint = geometry->edge[iEdge]->GetNode(0);
    auto jPoint = geometry->edge[iEdge]->GetNode(1);

    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Primitive variables w/o reconstruction ---*/

    const auto V_i = flowNodes->GetPrimitive(iPoint);
    const auto V_j = flowNodes->GetPrimitive(jPoint);
    numerics->SetPrimitive(V_i, V_j);

    /*--- Turbulent variables w/o reconstruction ---*/

    const auto T_i = (sst) ? nodes->GetPrimitive(iPoint) : nodes->GetSolution(iPoint);
    const auto T_j = (sst) ? nodes->GetPrimitive(jPoint) : nodes->GetSolution(jPoint);
    numerics->SetTurbVar(T_i, T_j);

    /*--- Grid Movement ---*/

    if (dynamic_grid)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                           geometry->node[jPoint]->GetGridVel());

    bool bad_i = false, bad_j = false;

    if (muscl) {
      su2double *Limiter_i = nullptr, *Limiter_j = nullptr;

      const auto Coord_i = geometry->node[iPoint]->GetCoord();
      const auto Coord_j = geometry->node[jPoint]->GetCoord();

      su2double Vector_ij[MAXNDIM] = {0.0};
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_ij[iDim] = (Coord_j[iDim] - Coord_i[iDim]);
      }

      const su2double Kappa = config->GetMUSCL_Kappa();

      if (muscl) {
        /*--- Reconstruct turbulence variables. ---*/

        auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
        auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

        if (limiter) {
          Limiter_i = nodes->GetLimiter(iPoint);
          Limiter_j = nodes->GetLimiter(jPoint);
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          
          const su2double T_ij = 0.5*(T_j[iVar] - T_i[iVar]);

          su2double Project_Grad_i = -T_ij;
          su2double Project_Grad_j = -T_ij;

          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_ij[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j += Vector_ij[iDim]*Gradient_j[iVar][iDim];
          }

          /*--- Blend upwind and centered differences ---*/

          Project_Grad_i = 0.5*((1.0-Kappa)*Project_Grad_i + (1.0+Kappa)*T_ij);
          Project_Grad_j = 0.5*((1.0-Kappa)*Project_Grad_j + (1.0+Kappa)*T_ij);

          /*--- Edge-based limiters ---*/

          if (limiter) {
            switch(config->GetKind_SlopeLimit_Turb()) {
              case VAN_ALBADA_EDGE:
                Limiter_i[iVar] = LimiterHelpers::vanAlbadaFunction(Project_Grad_i, T_ij);
                Limiter_j[iVar] = LimiterHelpers::vanAlbadaFunction(Project_Grad_j, T_ij);
                break;
              case PIPERNO:
                Limiter_i[iVar] = LimiterHelpers::pipernoFunction(Project_Grad_i, T_ij);
                Limiter_j[iVar] = LimiterHelpers::pipernoFunction(Project_Grad_j, T_ij);
                break;
            }

            /*--- Limit projection ---*/

            Project_Grad_i *= Limiter_i[iVar];
            Project_Grad_j *= Limiter_j[iVar];
          }

          solution_i[iVar] = T_i[iVar] + Project_Grad_i;
          solution_j[iVar] = T_j[iVar] - Project_Grad_j;

          bad_i = (solution_i[iVar] < 0.0) || (bad_i);
          bad_j = (solution_j[iVar] < 0.0) || (bad_j);

        }

      }
      else {
        for (iVar = 0; iVar < nVar; iVar++) {
          solution_i[iVar] = T_i[iVar];
          solution_j[iVar] = T_j[iVar];
        }
      }

      if (muscl) {
        /*--- Reconstruct mean flow primitive variables. ---*/

        auto Gradient_i = flowNodes->GetGradient_Reconstruction(iPoint);
        auto Gradient_j = flowNodes->GetGradient_Reconstruction(jPoint);

        if (limiterFlow) {
          Limiter_i = flowNodes->GetLimiter_Primitive(iPoint);
          Limiter_j = flowNodes->GetLimiter_Primitive(jPoint);
        }

        for (iVar = 0; iVar < solver[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {

          const su2double V_ij = 0.5*(V_j[iVar] - V_i[iVar]);

          su2double Project_Grad_i = -V_ij;
          su2double Project_Grad_j = -V_ij;

          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_ij[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j += Vector_ij[iDim]*Gradient_j[iVar][iDim];
          }

          /*--- Blend upwind and centered differences ---*/

          Project_Grad_i = 0.5*((1.0-Kappa)*Project_Grad_i + (1.0+Kappa)*V_ij);
          Project_Grad_j = 0.5*((1.0-Kappa)*Project_Grad_j + (1.0+Kappa)*V_ij);

          /*--- Edge-based limiters ---*/

          if (limiterFlow) {
            switch(config->GetKind_SlopeLimit_Flow()) {
              case VAN_ALBADA_EDGE:
                Limiter_i[iVar] = LimiterHelpers::vanAlbadaFunction(Project_Grad_i, V_ij);
                Limiter_j[iVar] = LimiterHelpers::vanAlbadaFunction(Project_Grad_j, V_ij);
                break;
              case PIPERNO:
                Limiter_i[iVar] = LimiterHelpers::pipernoFunction(Project_Grad_i, V_ij);
                Limiter_j[iVar] = LimiterHelpers::pipernoFunction(Project_Grad_j, V_ij);
                break;
            }

            /*--- Limit projection ---*/

            Project_Grad_i *= Limiter_i[iVar];
            Project_Grad_j *= Limiter_j[iVar];
          }

          flowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          flowPrimVar_j[iVar] = V_j[iVar] - Project_Grad_j;
        }

        const bool neg_pres_or_rho_i = (flowPrimVar_i[nDim+1] < 0.0) || (flowPrimVar_i[nDim+2] < 0.0);
        const bool neg_pres_or_rho_j = (flowPrimVar_j[nDim+1] < 0.0) || (flowPrimVar_j[nDim+2] < 0.0);

        su2double R = sqrt(fabs(flowPrimVar_j[nDim+2]/flowPrimVar_i[nDim+2]));
        su2double R_Plus_One = R+1.;
        su2double RoeSqVel = 0.0, SqVel_i = 0.0, SqVel_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          su2double RoeVelocity = (R*flowPrimVar_j[iDim+1]+flowPrimVar_i[iDim+1])/R_Plus_One;
          RoeSqVel += pow(RoeVelocity, 2);
          SqVel_i += pow(flowPrimVar_i[iDim+1],2);
          SqVel_j += pow(flowPrimVar_j[iDim+1],2);
        }
        su2double Energy_i = flowPrimVar_i[nDim+1]/(Gamma_Minus_One*flowPrimVar_i[nDim+2])+solution_i[0]+0.5*SqVel_i;
        su2double Energy_j = flowPrimVar_j[nDim+1]/(Gamma_Minus_One*flowPrimVar_j[nDim+2])+solution_j[0]+0.5*SqVel_j;
        su2double Enthalpy_i = Energy_i+flowPrimVar_i[nDim+1]/flowPrimVar_i[nDim+2];
        su2double Enthalpy_j = Energy_j+flowPrimVar_j[nDim+1]/flowPrimVar_j[nDim+2];
        su2double RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
        su2double RoeTke = (R*solution_j[0]+solution_i[0])/R_Plus_One;

        const bool bad_roe = (Gamma_Minus_One*(RoeEnthalpy-0.5*RoeSqVel-RoeTke) < 0.0);

        bad_i = neg_pres_or_rho_i || bad_roe || bad_i;
        bad_j = neg_pres_or_rho_j || bad_roe || bad_j;

      }
      else {
        for (iVar = 0; iVar < solver[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
          flowPrimVar_i[iVar] = V_i[iVar];
          flowPrimVar_j[iVar] = V_j[iVar];
        }
      }

      const bool bad_edge = bad_i || bad_j;

      /*--- Store extrapolated state ---*/

      numerics->SetPrimitive(bad_edge ? V_i : flowPrimVar_i, 
                             bad_edge ? V_j : flowPrimVar_j);
      numerics->SetTurbVar(bad_edge ? T_i : solution_i, 
                           bad_edge ? T_j : solution_j);

      /*--- Store values for limiter, even if limiter isn't being used ---*/

      if (muscl) {
        SetGradBasis(solver, geometry, config,
                   GradBasis_i, GradBasis_j,
                   limiter ? nodes->GetLimiter(iPoint) : OneVec, 
                   limiter ? nodes->GetLimiter(jPoint) : OneVec,
                   iPoint, jPoint);
        numerics->SetLimiter(bad_edge ? ZeroVec : GradBasis_i, 
                             bad_edge ? ZeroVec : GradBasis_j);
        
        // if (limiter) {
        //   numerics->SetLimiter(bad_edge ? ZeroVec : nodes->GetLimiter(iPoint), 
        //                        bad_edge ? ZeroVec : nodes->GetLimiter(jPoint));
        // }
        // else {
        //   numerics->SetLimiter(bad_edge ? ZeroVec : OneVec, 
        //                        bad_edge ? ZeroVec : OneVec);
        // }
      }

      /*--- Store nodal values ---*/

      numerics->SetNodalPrimitive(V_i, V_j);
      numerics->SetNodalTurbVar(T_i, T_j);
    }

    /*--- Update convective residual value ---*/

    auto residual = numerics->ComputeResidual(config);

    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, residual);
      Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
    }
    else {
      LinSysRes.AddBlock(iPoint, residual);
      LinSysRes.SubtractBlock(jPoint, residual);
      Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }

    /*--- Viscous contribution. ---*/

    Viscous_Residual(iEdge, geometry, solver,
                     numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);
  }
  } // end color loop

  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    Jacobian.SetDiagonalAsColumnSum();
  }
}

void CTurbSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver,
                                   CNumerics *numerics, CConfig *config) {

  const bool sst = ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST));

  CVariable* flowNodes = solver[FLOW_SOL]->GetNodes();

  /*--- Points in edge ---*/

  auto iPoint = geometry->edge[iEdge]->GetNode(0);
  auto jPoint = geometry->edge[iEdge]->GetNode(1);

  /*--- Points coordinates, and normal vector ---*/

  numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());
  numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

  /*--- Conservative variables w/o reconstruction ---*/

  numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint),
                         flowNodes->GetPrimitive(jPoint));
  
  numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint),
                               flowNodes->GetGradient_Primitive(jPoint));

  /*--- Turbulent variables w/o reconstruction, and its gradients ---*/

  numerics->SetTurbVar(sst ? nodes->GetPrimitive(iPoint) : nodes->GetSolution(iPoint),
                       sst ? nodes->GetPrimitive(jPoint) : nodes->GetSolution(jPoint));

  numerics->SetTurbVarGradient(nodes->GetGradient(iPoint),
                               nodes->GetGradient(jPoint));

  /*--- Menter's first blending function (only SST)---*/
  if (sst) {
    numerics->SetF1blending(nodes->GetF1blending(iPoint),
                            nodes->GetF1blending(jPoint));
    numerics->SetF2blending(nodes->GetF2blending(iPoint),
                            nodes->GetF2blending(jPoint));
    numerics->SetVorticity(flowNodes->GetVorticity(iPoint),
                           flowNodes->GetVorticity(jPoint));
    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint),
                           flowNodes->GetStrainMag(jPoint));
  }

  /*--- Compute residual, and Jacobians ---*/

  auto residual = numerics->ComputeResidual(config);

  if (ReducerStrategy) {
    EdgeFluxes.SubtractBlock(iEdge, residual);
    Jacobian.UpdateBlocksSub(iEdge, residual.jacobian_i, residual.jacobian_j);
  }
  else {
    LinSysRes.SubtractBlock(iPoint, residual);
    LinSysRes.AddBlock(jPoint, residual);
    Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);

    if (config->GetUse_Accurate_Turb_Jacobians()) {
      CorrectJacobian(geometry, solver, config, iPoint, jPoint, residual.jacobian_ic);
      CorrectJacobian(geometry, solver, config, jPoint, iPoint, residual.jacobian_jc);
    }
  }
  
}

void CTurbSolver::SumEdgeFluxes(CGeometry* geometry) {

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0; iPoint < nPoint; ++iPoint) {

    LinSysRes.SetBlock_Zero(iPoint);

    for (auto iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {

      auto iEdge = geometry->node[iPoint]->GetEdge(iNeigh);

      if (iPoint == geometry->edge[iEdge]->GetNode(0))
        LinSysRes.AddBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
      else
        LinSysRes.SubtractBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
    }
  }

}

void CTurbSolver::SetGradBasis(CSolver** solver, CGeometry *geometry, CConfig *config,
                               su2double* gradBasis_i, su2double* gradBasis_j,
                               su2double* limiter_i, su2double* limiter_j,
                               unsigned long iPoint, unsigned long jPoint) {
  const bool reconRequired = config->GetReconstructionGradientRequired();
  const unsigned short kindRecon = reconRequired ? config->GetKind_Gradient_Method_Recon()
                                                 : config->GetKind_Gradient_Method();

  const su2double kappa = config->GetMUSCL_Kappa();

  for (auto iVar = 0; iVar < nVar; iVar++) {
    gradBasis_i[iVar] = 0.5*kappa*limiter_i[iVar];
    gradBasis_j[iVar] = 0.5*kappa*limiter_j[iVar];
  }

  auto node_i = geometry->node[iPoint], node_j = geometry->node[jPoint];
  switch (kindRecon) {
    case GREEN_GAUSS:
      break;
    case WEIGHTED_LEAST_SQUARES:
    case LEAST_SQUARES:
      su2double weight = 1.0;
      su2double dist_ij[MAXNDIM] = {0.0};
      for (auto iDim = 0; iDim < nDim; iDim++)
        dist_ij[iDim] = node_j->GetCoord(iDim) - node_i->GetCoord(iDim);
      if (kindRecon == WEIGHTED_LEAST_SQUARES) {
        weight = 0.0;
        for (auto iDim = 0; iDim < nDim; iDim++)
          weight += pow(dist_ij[iDim],2); 
      }
      weight *= 0.5*(1.-kappa);

      auto var = solver[TURB_SOL]->GetNodes();
      auto S_i = reconRequired ? var->GetSmatrix_Aux(iPoint)
                               : var->GetSmatrix(iPoint);
      auto S_j = reconRequired ? var->GetSmatrix_Aux(jPoint)
                               : var->GetSmatrix(jPoint);
      for (auto iVar = 0; iVar < nVar; iVar++) {
        for (auto iDim = 0; iDim < nDim; iDim++) {
          for (auto jDim = 0; jDim < nDim; jDim++) {
            gradBasis_i[iVar] += weight*S_i[iDim][jDim]*dist_ij[iDim]*dist_ij[jDim]*limiter_i[iVar];
            gradBasis_j[iVar] += weight*S_j[iDim][jDim]*dist_ij[iDim]*dist_ij[jDim]*limiter_j[iVar];
          }
        }
      }
      break;
  }
}

void CTurbSolver::CorrectJacobian(CGeometry           *geometry,
                                  CSolver             **solver,
                                  CConfig             *config,
                                  const unsigned long iPoint,
                                  const unsigned long jPoint,
                                  const su2double     *const *const *const Jacobian_ic) {

  const bool wasActive = AD::BeginPassive();

  su2double Weight = 0.0;
  const bool gg  = config->GetKind_Gradient_Method() == GREEN_GAUSS;
  const bool ls  = config->GetKind_Gradient_Method() == LEAST_SQUARES;
  const bool wls = config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES;

  /*--- First we compute contributions of first neighbors to the Jacobian.
        In Green-Gauss, this contribution is scaled by 0.5*Sum(n_v)/r = 0 for
        volume nodes and (0.5*Sum(n_v)+n_s)/r for surface nodes. So only add to
        the Jacobian if iPoint is on a physical boundary. Least squares gradients
        do not have a surface term. ---*/

  /*--- Common factors for all Jacobian terms --*/
  const su2double HalfOnVol = 0.5/geometry->node[iPoint]->GetVolume();
  const su2double sign      = 1.0 - 2.0*(iPoint > jPoint);
    
  CVariable *nodesFlo = solver[FLOW_SOL]->GetNodes();

  if ((gg) && (geometry->node[iPoint]->GetPhysicalBoundary())) {

    for (auto iVar = 0; iVar < nVar; iVar++)
      for (auto jVar = 0; jVar < nVar; jVar++)
        Jacobian_i[iVar][jVar] = 0.0;

    Weight = -HalfOnVol*sign;
    
    /*--- Influence of boundary i on R(i,j) ---*/
    for (auto iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      const long iVertex = geometry->node[iPoint]->GetVertex(iMarker);
      if (iVertex != -1) {
        const su2double *Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (auto iDim = 0; iDim < nDim; iDim++)
          for (auto iVar = 0; iVar < nVar; iVar++)
            Jacobian_i[iVar][iVar] += Weight*Jacobian_ic[iDim][iVar][iVar]*Normal[iDim];
      }// iVertex
    }// iMarker

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);

  }// physical boundary

  /*--- Next we compute contributions of neighbor nodes to the Jacobian.
        To reduce extra communication overhead, we only consider first neighbors on
        the current rank. Note that Jacobian_ic is already weighted by 0.5 ---*/

  for (auto iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {

    for (auto iVar = 0; iVar < nVar; iVar++)
      for (auto jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }

    auto kPoint = geometry->node[iPoint]->GetPoint(iNeigh);

    auto kEdge = geometry->node[iPoint]->GetEdge(iNeigh);

    su2double Basis[MAXNDIM] = {0.0};
    if (gg) {
      const su2double signk = 1.0 - 2.0*(iPoint > kPoint);
      const su2double denom = nodesFlo->GetDensity(kPoint)/nodesFlo->GetDensity(iPoint);
      Weight = HalfOnVol*sign*signk/denom;
      for (auto iDim = 0; iDim < nDim; iDim++)
        Basis[iDim] = Weight*geometry->edge[kEdge]->GetNormal()[iDim];

      for (auto iDim = 0; iDim < nDim; iDim++)
        for (auto iVar = 0; iVar < nVar; iVar++)
          Jacobian_i[iVar][iVar] += Jacobian_ic[iDim][iVar][iVar]*Basis[iDim];
    }
    else {
      const su2double signk = 1.0 - 2.0*(iPoint > kPoint);
      const su2double denom = nodesFlo->GetDensity(kPoint)/nodesFlo->GetDensity(iPoint);

      su2double dist_ik[MAXNDIM] = {0.0};
      for (auto iDim = 0; iDim < nDim; iDim++)
        dist_ik[iDim] = geometry->node[kPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim);

      if (ls)
        Weight = 1.0;
      else {
        Weight = 0.0;
        for (auto iDim = 0; iDim < nDim; iDim++)
          Weight += pow(dist_ik[iDim],2);
      }

      const auto Smat = nodes->GetSmatrix(iPoint);
      for (auto iDim = 0; iDim < nDim; iDim++)
        for (auto jDim = 0; jDim < nDim; jDim++)
          Basis[iDim] += Weight*Smat[iDim][jDim]*dist_ik[jDim];

      for (auto iDim = 0; iDim < nDim; iDim++)
        for (auto jDim = 0; jDim < nDim; jDim++)
          for (auto iVar = 0; iVar < nVar; iVar++) {
            Jacobian_i[iVar][iVar] += Jacobian_ic[iDim][iVar][iVar]*Basis[iDim]/denom;
            Jacobian_j[iVar][iVar] -= Jacobian_ic[iDim][iVar][iVar]*Basis[iDim];
          }
    }

    Jacobian.SubtractBlock(iPoint, kPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, kPoint, Jacobian_i);

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_j);

  }// iNeigh

  AD::EndPassive(wasActive);
  
}


void CTurbSolver::BC_Sym_Plane(CGeometry      *geometry,
                               CSolver        **solver,
                               CNumerics      *conv_numerics,
                               CNumerics      *visc_numerics,
                               CConfig        *config,
                               unsigned short val_marker) {

  /*--- Convective and viscous fluxes across symmetry plane are equal to zero. ---*/

}

void CTurbSolver::BC_Euler_Wall(CGeometry      *geometry,
                                CSolver        **solver,
                                CNumerics      *conv_numerics,
                                CNumerics      *visc_numerics,
                                CConfig        *config,
                                unsigned short val_marker) {

  /*--- Convective fluxes across euler wall are equal to zero. ---*/

}

void CTurbSolver::BC_Riemann(CGeometry *geometry, CSolver **solver, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}


void CTurbSolver::BC_Giles(CGeometry *geometry, CSolver **solver, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Giles(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT:case TOTAL_CONDITIONS_PT_1D: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    break;
  case MIXING_IN:
    if (config->GetBoolTurbMixingPlane()){
      BC_Inlet_MixingPlane(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    }
    else{
      BC_Inlet_Turbo(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    }
    break;

  case STATIC_PRESSURE: case MIXING_OUT: case STATIC_PRESSURE_1D: case RADIAL_EQUILIBRIUM:
    BC_Outlet(geometry, solver, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::BC_Periodic(CGeometry *geometry, CSolver **solver,
                                  CNumerics *numerics, CConfig *config) {

  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated corectly during the communications. For implicit calculations
   the Jacobians and linear system are also correctly adjusted here. ---*/

  for (auto iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }

}

void CTurbSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver, CConfig *config) {

  const bool adjoint = config->GetContinuous_Adjoint() || (config->GetDiscrete_Adjoint() && config->GetFrozen_Visc_Disc());
  const bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  const bool sst = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);

  CVariable* flowNodes = solver[FLOW_SOL]->GetNodes();

  /*--- Set shared residual variables to 0 and declare
   *    local ones for current thread to work on. ---*/

  SU2_OMP_MASTER
  for (auto iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  SU2_OMP_BARRIER

  su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
  const su2double* coordMax[MAXNVAR] = {nullptr};
  unsigned long idxMax[MAXNVAR] = {0};

  /*--- Build implicit system ---*/

  SU2_OMP(for schedule(static,omp_chunk_size) nowait)
  for (auto iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the volume ---*/

    su2double Vol = (geometry->node[iPoint]->GetVolume() + geometry->node[iPoint]->GetPeriodicVolume());

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/

    su2double Delta = (sst) ? su2double(Vol / nodes->GetDelta_Time(iPoint))
                            : su2double(Vol / ((nodes->GetLocalCFL(iPoint)/flowNodes->GetLocalCFL(iPoint))*flowNodes->GetDelta_Time(iPoint)));
    // su2double Delta = Vol / ((nodes->GetLocalCFL(iPoint)/flowNodes->GetLocalCFL(iPoint))*flowNodes->GetDelta_Time(iPoint));
    Jacobian.AddVal2Diag(iPoint, Delta);

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (auto iVar = 0; iVar < nVar; iVar++) {
      unsigned long total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = -LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;

      su2double Res = fabs(LinSysRes[total_index]);
      resRMS[iVar] += Res*Res;
      if (Res > resMax[iVar]) {
        resMax[iVar] = Res;
        idxMax[iVar] = iPoint;
        coordMax[iVar] = geometry->node[iPoint]->GetCoord();
      }
    }
  }
  SU2_OMP_CRITICAL
  for (auto iVar = 0; iVar < nVar; iVar++) {
    AddRes_RMS(iVar, resRMS[iVar]);
    AddRes_Max(iVar, resMax[iVar], geometry->node[idxMax[iVar]]->GetGlobalIndex(), coordMax[iVar]);
  }

  /*--- Initialize residual and solution at the ghost points ---*/

  SU2_OMP(sections)
  {
    SU2_OMP(section)
    for (auto iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      LinSysRes.SetBlock_Zero(iPoint);

    SU2_OMP(section)
    for (auto iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      LinSysSol.SetBlock_Zero(iPoint);
  }

  /*--- Solve or smooth the linear system ---*/

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  SU2_OMP_MASTER
  {
    SetIterLinSolver(iter);
    SetResLinSolver(System.GetResidual());
  }
  SU2_OMP_BARRIER


  ComputeUnderRelaxationFactor(solver, config);

  /*--- Update solution (system written in terms of increments) ---*/

  if (!adjoint) {

    /*--- Update the turbulent solution. Only SST variants are clipped. ---*/

    switch (config->GetKind_Turb_Model()) {

      case SA: case SA_E: case SA_COMP: case SA_E_COMP: case SA_NEG:

        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (auto iPoint = 0; iPoint < nPointDomain; iPoint++) {
          nodes->AddSolution(iPoint, 0, nodes->GetUnderRelaxation(iPoint)*LinSysSol[iPoint]);
        }
        break;

      case SST: case SST_SUST:

        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (auto iPoint = 0; iPoint < nPointDomain; iPoint++) {

          for (auto iVar = 0; iVar < nVar; iVar++)
            nodes->AddSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint)*LinSysSol[iPoint*nVar+iVar]);
        }
        break;

    }
  }

  SU2_OMP_MASTER
  {
    for (auto iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
      InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
      CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    }

    /*--- MPI solution ---*/

    InitiateComms(geometry, config, SOLUTION);
    CompleteComms(geometry, config, SOLUTION);

    /*--- Compute the root mean square residual ---*/

    SetResidual_RMS(geometry, config);
  }
  SU2_OMP_BARRIER

}

void CTurbSolver::ComputeUnderRelaxationFactor(CSolver **solver, CConfig *config) {

  /* Only apply the turbulent under-relaxation to the SA variants. The
   SA_NEG model is more robust due to allowing for negative nu_tilde,
   so the under-relaxation is not applied to that variant. */

  bool sa_model = ((config->GetKind_Turb_Model() == SA)        ||
                   (config->GetKind_Turb_Model() == SA_E)      ||
                   (config->GetKind_Turb_Model() == SA_COMP)   ||
                   (config->GetKind_Turb_Model() == SA_E_COMP));
  
  bool sst_model = ((config->GetKind_Turb_Model() == SST)      ||
                    (config->GetKind_Turb_Model() == SST_SUST));

  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  const su2double allowableRatio =  0.99;
  const su2double eps = numeric_limits<passivedouble>::epsilon();
  const su2double CFLInc = config->GetCFL_AdaptParam(1);
  const su2double CFLMin = config->GetCFL_AdaptParam(2)*config->GetCFLMaxRedCoeff_Turb();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0; iPoint < nPointDomain; iPoint++) {

    su2double localUnderRelaxation = 1.0;
    for (auto iVar = 0; iVar < nVar; iVar++) {

      /* We impose a limit on the maximum percentage that the
       turbulence variables can change over a nonlinear iteration. */

      const unsigned long index = iPoint * nVar + iVar;
      if (sa_model || sst_model) {
        if (fabs(LinSysSol[index]) > allowableRatio*fabs(nodes->GetSolution(iPoint, iVar))) {
          const su2double invratio = fabs(nodes->GetSolution(iPoint, iVar))/(fabs(LinSysSol[index])+eps);
          localUnderRelaxation = min(allowableRatio*invratio, localUnderRelaxation);
        }
        // if (iVar == 1 && nodes->GetSolution(iPoint, iVar)+LinSysSol[index] < lowerlimit[iVar]) {
        //   const su2double ratio = (nodes->GetSolution(iPoint, iVar)-lowerlimit[iVar])/(fabs(LinSysSol[index])+eps);
        //   localUnderRelaxation = min(ratio, localUnderRelaxation);
        // }
      }

    }

    /* Choose the minimum factor between mean flow and turbulence. */

    localUnderRelaxation = min(localUnderRelaxation, solver[FLOW_SOL]->GetNodes()->GetUnderRelaxation(iPoint));

    /* Threshold the relaxation factor in the event that there is
     a very small value. This helps avoid catastrophic crashes due
     to non-realizable states by canceling the update. */

    if (localUnderRelaxation < 1.0e-10 && nodes->GetLocalCFL(iPoint) > CFLMin*CFLInc) localUnderRelaxation = 0.0;

    /* Store the under-relaxation factor for this point. */

    nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);

  }

}

void CTurbSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver, CConfig *config,
                                       unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  const bool sst_model = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
  const bool implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  const bool first_order = (config->GetTime_Marching() == DT_STEPPING_1ST);
  const bool second_order = (config->GetTime_Marching() == DT_STEPPING_2ND);
  const bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Flow solution, needed to get density. ---*/

  CVariable* flowNodes = solver[FLOW_SOL]->GetNodes();

  /*--- Store the physical time step ---*/

  const su2double TimeStep = config->GetDelta_UnstTimeND();

  /*--- Local variables ---*/

  unsigned short iVar, iMarker, iDim, iNeigh;
  unsigned long iPoint, jPoint, iVertex, iEdge;

  const su2double *U_time_nM1 = nullptr, *U_time_n = nullptr, *U_time_nP1 = nullptr;
  su2double Volume_nM1, Volume_nP1;
  su2double Density_nM1, Density_n, Density_nP1;
  const su2double *Normal = nullptr, *GridVel_i = nullptr, *GridVel_j = nullptr;
  su2double Residual_GCL;

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {

    /*--- Loop over all nodes (excluding halos) ---*/

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

      Volume_nP1 = geometry->node[iPoint]->GetVolume();

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/

      if (sst_model) {

        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        if (incompressible){
          /*--- This is temporary and only valid for constant-density problems:
          density could also be temperature dependent, but as it is not a part
          of the solution vector it's neither stored for previous time steps
          nor updated with the solution at the end of each iteration. */
          Density_nM1 = flowNodes->GetDensity(iPoint);
          Density_n   = flowNodes->GetDensity(iPoint);
          Density_nP1 = flowNodes->GetDensity(iPoint);
        }
        else{
          Density_nM1 = flowNodes->GetSolution_time_n1(iPoint)[0];
          Density_n   = flowNodes->GetSolution_time_n(iPoint,0);
          Density_nP1 = flowNodes->GetSolution(iPoint,0);
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          if (first_order)
            LinSysRes(iPoint,iVar) += ( Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (second_order)
            LinSysRes(iPoint,iVar) += ( 3.0*Density_nP1*U_time_nP1[iVar] - 4.0*Density_n*U_time_n[iVar]
                                       +1.0*Density_nM1*U_time_nM1[iVar] ) * Volume_nP1/(2.0*TimeStep);
        }

      } else {

        for (iVar = 0; iVar < nVar; iVar++) {
          if (first_order)
            LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
          if (second_order)
            LinSysRes(iPoint,iVar) += ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                                       +1.0*U_time_nM1[iVar] ) * Volume_nP1/(2.0*TimeStep);
        }
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
      }

    }

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

      GridVel_i = geometry->node[iPoint]->GetGridVel();
      U_time_n  = nodes->GetSolution_time_n(iPoint);
      Density_n = 1.0;

      if (sst_model) {
        if (incompressible)
          Density_n = flowNodes->GetDensity(iPoint); // Temporary fix
        else
          Density_n = flowNodes->GetSolution_time_n(iPoint,0);
      }

      for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnNeighbor(); iNeigh++) {

        iEdge = geometry->node[iPoint]->GetEdge(iNeigh);
        Normal = geometry->edge[iEdge]->GetNormal();

        jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
        GridVel_j = geometry->node[jPoint]->GetGridVel();

        /*--- Determine whether to consider the normal outward or inward. ---*/
        su2double dir = (geometry->edge[iEdge]->GetNode(0) == iPoint)? 0.5 : -0.5;

        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL += dir*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

        Residual_GCL *= Density_n;

        for (iVar = 0; iVar < nVar; iVar++)
          LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
      }
    }

    /*--- Loop over the boundary edges ---*/

    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

        SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

          /*--- Get the index for node i plus the boundary face normal ---*/

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          /*--- Grid velocities stored at boundary node i ---*/

          GridVel_i = geometry->node[iPoint]->GetGridVel();

          /*--- Compute the GCL term by dotting the grid velocity with the face
           normal. The normal is negated to match the boundary convention. ---*/

          Residual_GCL = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];

          /*--- Compute the GCL component of the source term for node i ---*/

          U_time_n = nodes->GetSolution_time_n(iPoint);

          /*--- Multiply by density at node i for the SST model ---*/

          if (sst_model) {
            if (incompressible)
              Density_n = flowNodes->GetDensity(iPoint); // Temporary fix
            else
              Density_n = flowNodes->GetSolution_time_n(iPoint,0);

            for (iVar = 0; iVar < nVar; iVar++)
              LinSysRes(iPoint,iVar) += Density_n*U_time_n[iVar]*Residual_GCL;
          }
          else {
            for (iVar = 0; iVar < nVar; iVar++)
              LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
          }

        }
      }
    }

    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/

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

      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();

      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/

      if (sst_model) {

        /*--- If this is the SST model, we need to multiply by the density
         in order to get the conservative variables ---*/
        if (incompressible) {
          /*--- This is temporary and only valid for constant-density problems:
          density could also be temperature dependent, but as it is not a part
          of the solution vector it's neither stored for previous time steps
          nor updated with the solution at the end of each iteration. */
          Density_nM1 = flowNodes->GetDensity(iPoint);
          Density_n   = flowNodes->GetDensity(iPoint);
          Density_nP1 = flowNodes->GetDensity(iPoint);
        }
        else {
          Density_nM1 = flowNodes->GetSolution_time_n1(iPoint)[0];
          Density_n   = flowNodes->GetSolution_time_n(iPoint,0);
          Density_nP1 = flowNodes->GetSolution(iPoint,0);
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          if (first_order)
            LinSysRes(iPoint,iVar) += (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (second_order)
            LinSysRes(iPoint,iVar) += (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
                                      + (Density_nM1*U_time_nM1[iVar] - Density_n*U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }

      } else {

        for (iVar = 0; iVar < nVar; iVar++) {
          if (first_order)
            LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
          if (second_order)
            LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
                                       + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
        }
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
      if (implicit) {
        if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
        if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
      }
    }

  } // end dynamic grid

}


void CTurbSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent;
  const su2double *Solution_Fine = nullptr;
  const bool restart_cfl = config->GetRestart_CFL();

  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  SU2_OMP_MASTER
  {

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Skip flow variables ---*/

  unsigned short skipVars = 0;

  if (nDim == 2) skipVars += 6;
  if (nDim == 3) skipVars += 8;

  /*--- Adjust the number of solution variables in the incompressible
   restart. We always carry a space in nVar for the energy equation in the
   mean flow solver, but we only write it to the restart if it is active.
   Therefore, we must reduce skipVars here if energy is inactive so that
   the turbulent variables are read correctly. ---*/

  bool incompressible       = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();

  if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

  /*--- Load data from the restart into correct containers. ---*/

  unsigned long counter = 0, iPoint_Global = 0;
  for (; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {


    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++)
        nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index+iVar]);

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (counter != nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- MPI solution and compute the eddy viscosity ---*/

  solver[MESH_0][TURB_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][TURB_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  } // end SU2_OMP_MASTER, pre and postprocessing are thread-safe.
  SU2_OMP_BARRIER
    
  const bool sst = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
  if (sst) static_cast<CTurbSSTSolver*>(solver[MESH_0][TURB_SOL])->SetPrimitive_Variables(solver[MESH_0]);
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  
  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      su2double Solution_Coarse[MAXNVAR] = {0.0};
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][TURB_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution_Coarse[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][TURB_SOL]->GetNodes()->SetSolution(iPoint,Solution_Coarse);
    }

    SU2_OMP_MASTER
    {
      solver[iMesh][TURB_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver[iMesh][TURB_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    }
    SU2_OMP_BARRIER
    
    if (sst) static_cast<CTurbSSTSolver*>(solver[iMesh][TURB_SOL])->SetPrimitive_Variables(solver[iMesh]);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    solver[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }

  /*--- Store the CFL number if loaded from restart ---*/
  if (restart_cfl)
    for (auto iPoint = 0; iPoint < geometry[MESH_0]->GetnPointDomain(); iPoint++)
      nodes->SetLocalCFL(iPoint, solver[MESH_0][FLOW_SOL]->GetNodes()->GetLocalCFL(iPoint)*config->GetCFLRedCoeff_Turb());

  /*--- Go back to single threaded execution. ---*/
  SU2_OMP_MASTER
  {

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars; Restart_Vars = nullptr;
  delete [] Restart_Data; Restart_Data = nullptr;

  } // end SU2_OMP_MASTER
  SU2_OMP_BARRIER

}
