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

  delete nodes;
}

void CTurbSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver,
                                  CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const unsigned short turbModel = config->GetKind_Turb_Model();
  const bool sst = ((turbModel == SST) || (turbModel == SST_SUST));
  const bool kappa  = config->GetUse_Accurate_Kappa_Jacobians();

  const auto nFlowVarGrad = solver[FLOW_SOL]->GetnPrimVarGrad();

  CVariable* flowNodes = solver[FLOW_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Static arrays of MUSCL-reconstructed flow primitives and turbulence variables (thread safety). ---*/
  su2double turbPrimVar_i[MAXNVAR] = {0.0}, flowPrimVar_i[MAXNVARFLOW] = {0.0};
  su2double turbPrimVar_j[MAXNVAR] = {0.0}, flowPrimVar_j[MAXNVARFLOW] = {0.0};

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring)
  {
  /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
  SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
  for(auto k = 0ul; k < color.size; ++k) {

    const auto iEdge = color.indices[k];

    /*--- Points in edge and normal vectors ---*/

    const auto iPoint = geometry->edge[iEdge]->GetNode(0);
    const auto jPoint = geometry->edge[iEdge]->GetNode(1);

    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Primitive variables w/o reconstruction ---*/

    const auto V_i = flowNodes->GetPrimitive(iPoint);
    const auto V_j = flowNodes->GetPrimitive(jPoint);

    /*--- Turbulent variables w/o reconstruction ---*/

    const auto T_i = sst ? nodes->GetPrimitive(iPoint) : nodes->GetSolution(iPoint);
    const auto T_j = sst ? nodes->GetPrimitive(jPoint) : nodes->GetSolution(jPoint);

    /*--- Grid Movement ---*/

    if (dynamic_grid)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                           geometry->node[jPoint]->GetGridVel());

    bool good_i = true, good_j = true;
    // bool good_i = (!geometry->node[iPoint]->GetPhysicalBoundary());
    // bool good_j = (!geometry->node[jPoint]->GetPhysicalBoundary());
    // bool muscl = (config->GetMUSCL_Turb()) && good_i && good_j;
    bool muscl = (config->GetMUSCL_Turb()) && (good_i || good_j);
    if (muscl) {
      /*--- Reconstruction ---*/

      solver[FLOW_SOL]->ExtrapolateState(solver, geometry, config, iPoint, jPoint, flowPrimVar_i, flowPrimVar_j, 
                                         turbPrimVar_i, turbPrimVar_j, good_i, good_j, nFlowVarGrad, nVar);

      /*--- Check for non-physical solutions after reconstruction. If found, use the
       cell-average value of the solution. This is a locally 1st order approximation,
       which is typically only active during the start-up of a calculation. ---*/

      solver[FLOW_SOL]->CheckExtrapolatedState(config, flowPrimVar_i, flowPrimVar_j, turbPrimVar_i, turbPrimVar_j, nVar, good_i, good_j);

      // muscl = good_i && good_j;

      // /*--- Store the state ---*/

      // numerics->SetPrimitive(muscl ? flowPrimVar_i : V_i, 
      //                        muscl ? flowPrimVar_j : V_j);
      // numerics->SetTurbVar(  muscl ? turbPrimVar_i : T_i, 
                             // muscl ? turbPrimVar_j : T_j);
      muscl = good_i || good_j;

      /*--- Store the state ---*/

      numerics->SetPrimitive(good_i ? flowPrimVar_i : V_i, 
                             good_j ? flowPrimVar_j : V_j);
      numerics->SetTurbVar(  good_i ? turbPrimVar_i : T_i, 
                             good_j ? turbPrimVar_j : T_j);
    }
    else {
      /*--- Nodal values ---*/

      numerics->SetPrimitive(V_i, V_j);
      numerics->SetTurbVar(  T_i, T_j);
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
      if (muscl && kappa) {
        const su2double rho_i = good_i ? flowPrimVar_i[nDim+2] : V_i[nDim+2],
                        rho_j = good_j ? flowPrimVar_j[nDim+2] : V_j[nDim+2];
        SetExtrapolationJacobian(solver, geometry, config,
                                 &rho_i, &rho_j,
                                 residual.jacobian_i, 
                                 residual.jacobian_j,
                                 good_i, good_j,
                                 iPoint, jPoint);
        SetExtrapolationJacobian(solver, geometry, config,
                                 &rho_j, &rho_i,
                                 residual.jacobian_j, 
                                 residual.jacobian_i,
                                 good_j, good_i,
                                 jPoint, iPoint);
      }
      else {
        Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
      }
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

  const auto iPoint = geometry->edge[iEdge]->GetNode(0);
  const auto jPoint = geometry->edge[iEdge]->GetNode(1);

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
    numerics->SetVorticityMag(flowNodes->GetVorticityMag(iPoint),
                              flowNodes->GetVorticityMag(jPoint));
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
      CorrectJacobian(solver, geometry, config, iPoint, jPoint, residual.jacobianWeights_i);
      CorrectJacobian(solver, geometry, config, jPoint, iPoint, residual.jacobianWeights_j);
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

void CTurbSolver::SetExtrapolationJacobian(CSolver             **solver,
                                           const CGeometry     *geometry,
                                           const CConfig       *config,
                                           const su2double     *rho_l, 
                                           const su2double     *rho_r,
                                           const su2double     *const *const dFl_dUl,
                                           const su2double     *const *const dFr_dUr,
                                           const bool          good_i,
                                           const bool          good_j,
                                           const unsigned long iPoint, 
                                           const unsigned long jPoint) {

  const bool wasActive = AD::BeginPassive();

  const bool reconRequired = config->GetReconstructionGradientRequired();
  const unsigned short kindRecon = reconRequired ? config->GetKind_Gradient_Method_Recon()
                                                 : config->GetKind_Gradient_Method();

  const bool gg = (kindRecon == GREEN_GAUSS);
  
  const bool limiter = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (config->GetInnerIter() <= config->GetLimiterIter());

  auto flowNodes = solver[FLOW_SOL]->GetNodes();

  const su2double kappa = config->GetMUSCL_Kappa();
  const su2double sign  = 1.0 - 2.0*(iPoint > jPoint);
  const su2double sign_grad_i = -1.0 + 2.0*(gg);
  const su2double dUl_dVl = *rho_l;
  const su2double dUr_dVr = *rho_r;
  const su2double dVi_dUi = 1.0/flowNodes->GetDensity(iPoint);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Compute the Jacobian terms corresponding to the constant   ---*/
  /*---         term and the difference (0.5*Psi_i*(V_i-V_j)).             ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Store limiters ---*/

  su2double lim_i[MAXNVAR] = {0.0}, lim_j[MAXNVAR] = {0.0};
  for (auto iVar = 0; iVar < nVar; iVar++) {
    lim_i[iVar] = nodes->GetLimiter(iPoint,iVar);
    lim_j[iVar] = nodes->GetLimiter(jPoint,iVar);
  }

  /*--- Store reconstruction weights ---*/

  su2double dVl_dVi[MAXNVAR] = {0.0}, dVr_dVi[MAXNVAR] = {0.0};
  for (auto iVar = 0; iVar < nVar; iVar++) {
    dVl_dVi[iVar] = sign*(1.0 + 0.5*lim_i[iVar]*good_i);
    dVr_dVi[iVar] = sign*(    - 0.5*lim_j[iVar]*good_j);
  }

  for (auto iVar = 0; iVar < nVar; iVar++) {
    for (auto jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] = (dFl_dUl[iVar][jVar]*dUl_dVl*dVl_dVi[jVar]
                              + dFr_dUr[iVar][jVar]*dUr_dVr*dVr_dVi[jVar])*dVi_dUi;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Compute the Jacobian terms corresponding to the nodal      ---*/
  /*---         gradient projection (0.5*Psi_i*gradV_i*dist_ij).           ---*/
  /*--------------------------------------------------------------------------*/

  const auto node_i = geometry->node[iPoint], node_j = geometry->node[jPoint];

  su2double dist_ij[MAXNDIM] = {0.0}, gradWeight[MAXNDIM] = {0.0};
  for (auto iDim = 0; iDim < nDim; iDim++)
    dist_ij[iDim] = node_j->GetCoord(iDim) - node_i->GetCoord(iDim);

  /*--- Green-Gauss surface terms ---*/

  if (gg && node_i->GetPhysicalBoundary()) {

    SetSurfaceGradWeights_GG(gradWeight, geometry, config, iPoint);

    su2double gradWeightDotDist = 0.0;
    for (auto iDim = 0; iDim < nDim; iDim++)
      gradWeightDotDist += gradWeight[iDim]*dist_ij[iDim];

    const su2double factor = sign*gradWeightDotDist*good_i;
    for (auto iVar = 0; iVar < nVar; iVar++)
      dVl_dVi[iVar] = factor*lim_i[iVar];

    for (auto iVar = 0; iVar < nVar; iVar++)
      for (auto jVar = 0; jVar < nVar; jVar++)
        Jacobian_i[iVar][jVar] += dFl_dUl[iVar][jVar]*dUl_dVl*dVl_dVi[jVar]*dVi_dUi;
  }

  /*--- Neighbor node terms ---*/

  for (auto iNeigh = 0; iNeigh < node_i->GetnPoint(); iNeigh++) {
    const auto kPoint = node_i->GetPoint(iNeigh);
    const su2double dVk_dUk = 1.0/flowNodes->GetDensity(kPoint);

    SetGradWeights(gradWeight, solver[TURB_SOL], geometry, config, iPoint, kPoint, reconRequired);

    su2double gradWeightDotDist = 0.0;
    for (auto iDim = 0; iDim < nDim; iDim++)
      gradWeightDotDist += gradWeight[iDim]*dist_ij[iDim];

    const su2double factor = sign*gradWeightDotDist*good_i;
    for (auto iVar = 0; iVar < nVar; iVar++)
      dVl_dVi[iVar] = factor*lim_i[iVar];

    for (auto iVar = 0; iVar < nVar; iVar++) {
      for (auto jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] += dFl_dUl[iVar][jVar]*dUl_dVl*dVl_dVi[jVar]*dVi_dUi*sign_grad_i;
        Jacobian_j[iVar][jVar]  = dFl_dUl[iVar][jVar]*dUl_dVl*dVl_dVi[jVar]*dVk_dUk;
      }
    }

    Jacobian.AddBlock(iPoint, kPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, kPoint, Jacobian_j);
  }

  Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
  Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);

  AD::EndPassive(wasActive);
}

void CTurbSolver::CorrectJacobian(CSolver             **solver,
                                  const CGeometry     *geometry,
                                  const CConfig       *config,
                                  const unsigned long iPoint,
                                  const unsigned long jPoint,
                                  const su2double     *const *const jacobianWeights_i) {

  const bool wasActive = AD::BeginPassive();

  const bool gg  = config->GetKind_Gradient_Method() == GREEN_GAUSS;

  const auto node_i = geometry->node[iPoint];

  CVariable *nodesFlo = solver[FLOW_SOL]->GetNodes();

  const su2double sign = 1.0 - 2.0*(iPoint > jPoint);
  const su2double sign_grad_i = -1.0 + 2.0*(gg);

  for (auto iVar = 0; iVar < nVar; iVar++) {
    for (auto jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] = 0.0;
      Jacobian_j[iVar][jVar] = 0.0;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Compute contributions of surface terms to the Jacobian.    ---*/
  /*---         In Green-Gauss, the weight of the surface node             ---*/
  /*---         contribution is (0.5*Sum(n_v)+n_s)/r. Least squares        ---*/
  /*---         gradients do not have a surface term.                      ---*/
  /*--------------------------------------------------------------------------*/

  
  su2double gradWeight[MAXNDIM] = {0.0};
  if (gg && node_i->GetPhysicalBoundary()) {
    SetSurfaceGradWeights_GG(gradWeight, geometry, config, iPoint);

    for (auto iVar = 0; iVar < nVar; iVar++)
      for (auto iDim = 0; iDim < nDim; iDim++)
        Jacobian_i[iVar][iVar] += sign*jacobianWeights_i[iVar][iDim]*gradWeight[iDim];

  }// physical boundary

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Compute contributions of neighbor nodes to the Jacobian.   ---*/
  /*---         To prevent extra communication overhead, we only consider  ---*/
  /*---         neighbors on the current rank. Note that jacobianWeights_i ---*/
  /*---         is already weighted by 0.5.                                ---*/
  /*--------------------------------------------------------------------------*/

  for (auto iNeigh = 0; iNeigh < node_i->GetnPoint(); iNeigh++) {
    const auto kPoint = node_i->GetPoint(iNeigh);
    const su2double ratio_k = nodesFlo->GetDensity(iPoint)/nodesFlo->GetDensity(kPoint);

    SetGradWeights(gradWeight, solver[TURB_SOL], geometry, config, iPoint, kPoint);
    
    for (auto iVar = 0; iVar < nVar; iVar++) {
      Jacobian_j[iVar][iVar] = 0.0;
      for (auto iDim = 0; iDim < nDim; iDim++) {
        Jacobian_i[iVar][iVar] += sign*jacobianWeights_i[iVar][iDim]*gradWeight[iDim]*sign_grad_i;
        Jacobian_j[iVar][iVar] += sign*jacobianWeights_i[iVar][iDim]*gradWeight[iDim]*ratio_k;
      }
    }

    Jacobian.SubtractBlock(iPoint, kPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, kPoint, Jacobian_j);
  }// iNeigh

  Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
  Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);

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

  const unsigned long InnerIter   = config->GetInnerIter();
  const unsigned long WFStartIter = config->GetWallFunction_Start_Iter();

  const bool adjoint = config->GetContinuous_Adjoint() || (config->GetDiscrete_Adjoint() && config->GetFrozen_Visc_Disc());
  const bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  const bool sst = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);

  const bool disc_adjoint   = config->GetDiscrete_Adjoint();
  const bool restart        = config->GetRestart();
  const bool wall_functions = (config->GetWall_Functions() && ((disc_adjoint) || (InnerIter >= WFStartIter) || (restart)));

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

  /*--- Set residual to 0 if turbulent variables are calculated with wall functions ---*/

  if (wall_functions) {
    for (auto iPoint = 0; iPoint < nPointDomain; iPoint++) {
      const auto node_i = geometry->node[iPoint];
      if (node_i->GetBool_Wall_Neighbor() && nodes->GetTauWall(iPoint) > 0) {

        LinSysRes.SetBlock_Zero(iPoint);
        for (auto iVar = 0; iVar < nVar; iVar++) {
          const unsigned long total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
    }
  }

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

    /*--- Update the turbulent solution. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0; iPoint < nPointDomain; iPoint++)
      for (auto iVar = 0; iVar < nVar; iVar++)
        nodes->AddSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint)*LinSysSol[iPoint*nVar+iVar]);

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

  const bool sa_model = ((config->GetKind_Turb_Model() == SA)        ||
                         (config->GetKind_Turb_Model() == SA_E)      ||
                         (config->GetKind_Turb_Model() == SA_COMP)   ||
                         (config->GetKind_Turb_Model() == SA_E_COMP));
  
  const bool sst_model = ((config->GetKind_Turb_Model() == SST)      ||
                          (config->GetKind_Turb_Model() == SST_SUST));

  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  const su2double allowableDecrease = 0.5;
  const su2double allowableIncrease = 1.0;
  const su2double eps = numeric_limits<passivedouble>::epsilon();
  const su2double CFLInc = config->GetCFL_AdaptParam(1);
  const su2double CFLMin = config->GetCFL_AdaptParam(2)*config->GetCFLMaxRedCoeff_Turb();

  if (sa_model || sst_model) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0; iPoint < nPointDomain; iPoint++) {

      su2double localUnderRelaxation = 1.0;
      for (auto iVar = 0; iVar < nVar; iVar++) {

        /* We impose a limit on the maximum percentage that the
         turbulence variables can change over a nonlinear iteration. */

        const unsigned long index = iPoint * nVar + iVar;
        const su2double allowableRatio = (LinSysSol[index] > 0) ? allowableIncrease : allowableDecrease;
        const su2double allowableChange = allowableRatio*fabs(nodes->GetSolution(iPoint, iVar));
        const su2double change = fabs(LinSysSol[index]);
        if (change > allowableChange) localUnderRelaxation = min(allowableChange/change, localUnderRelaxation);
        
      }

      /* Choose the minimum factor between mean flow and turbulence. */

      // localUnderRelaxation = min(localUnderRelaxation, solver[FLOW_SOL]->GetNodes()->GetUnderRelaxation(iPoint));

      /* Threshold the relaxation factor in the event that there is
       a very small value. This helps avoid catastrophic crashes due
       to non-realizable states by canceling the update. */

      // if (localUnderRelaxation < 0.1) localUnderRelaxation = 0.0;

      /* Store the under-relaxation factor for this point. */

      nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);

    }
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

      /*--- Load local CFL number ---*/

      if (restart_cfl) {
        index += nVar + nDim*(config->GetGrid_Movement()) + 1;
        nodes->SetLocalCFL(iPoint, Restart_Data[index]);
      }

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
