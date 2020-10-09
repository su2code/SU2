/*!
 * \file CNSSolver.cpp
 * \brief Main subrotuines for solving Finite-Volume Navier-Stokes flow problems.
 * \author F. Palacios, T. Economon
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


#include "../../include/solvers/CNSSolver.hpp"
#include "../../include/variables/CNSVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"

CNSSolver::CNSSolver(void) : CEulerSolver() { }

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) :
           CEulerSolver(geometry, config, iMesh, true) {

  /*--- This constructor only allocates/inits what is extra to CEulerSolver. ---*/

  unsigned short iMarker, iDim;
  unsigned long iVertex;

  /*--- Store the values of the temperature and the heat flux density at the boundaries,
   used for coupling with a solid donor cell ---*/
  unsigned short nHeatConjugateVar = 4;

  HeatConjugateVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatConjugateVar[iMarker] = new su2double* [nVertex[iMarker]];
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      HeatConjugateVar[iMarker][iVertex] = new su2double [nHeatConjugateVar]();
      HeatConjugateVar[iMarker][iVertex][0] = config->GetTemperature_FreeStreamND();
    }
  }

  /*--- Allocates a 2D array with variable "outer" sizes and init to 0. ---*/

  auto Alloc2D = [](unsigned long M, const unsigned long* N, su2double**& X) {
    X = new su2double* [M];
    for(unsigned long i = 0; i < M; ++i)
      X[i] = new su2double [N[i]] ();
  };

  /*--- Heat flux in all the markers ---*/

  Alloc2D(nMarker, nVertex, HeatFlux);
  Alloc2D(nMarker, nVertex, HeatFluxTarget);

  /*--- Y plus in all the markers ---*/

  Alloc2D(nMarker, nVertex, YPlus);

  /*--- Skin friction in all the markers ---*/

  CSkinFriction = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSkinFriction[iMarker] = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      CSkinFriction[iMarker][iDim] = new su2double[nVertex[iMarker]] ();
    }
  }

  /*--- Non dimensional aerodynamic coefficients ---*/

  ViscCoeff.allocate(nMarker);
  SurfaceViscCoeff.allocate(config->GetnMarker_Monitoring());

  /*--- Heat flux and buffet coefficients ---*/

  HF_Visc = new su2double[nMarker];
  MaxHF_Visc = new su2double[nMarker];

  Surface_HF_Visc = new su2double[config->GetnMarker_Monitoring()];
  Surface_MaxHF_Visc = new su2double[config->GetnMarker_Monitoring()];

  /*--- Buffet sensor in all the markers and coefficients ---*/

  if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){

    Alloc2D(nMarker, nVertex, Buffet_Sensor);
    Buffet_Metric = new su2double[nMarker];
    Surface_Buffet_Metric = new su2double[config->GetnMarker_Monitoring()];

  }

  /*--- Read farfield conditions from config ---*/

  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  Tke_Inf         = config->GetTke_FreeStreamND();

  /*--- Initialize the seed values for forward mode differentiation. ---*/

  switch(config->GetDirectDiff()) {
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    default:
      /*--- Already done upstream. ---*/
      break;
  }
             
  /*--- Initialize the wall index map and solution matrices on the first iteration. ---*/
  if (config->GetWall_Functions()) {
    unsigned long counter = 0;
    for (auto iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (geometry->node[iPoint]->GetBool_Wall_Neighbor()) {
        nodes->SetWallMap(iPoint,counter);
        counter++;
      }
    }
    if (counter > 0) nodes->InitializeWallSolution(counter);
  }

}

CNSSolver::~CNSSolver(void) {

  unsigned short iMarker, iDim;

  unsigned long iVertex;

  delete [] Buffet_Metric;
  delete [] HF_Visc;
  delete [] MaxHF_Visc;

  delete [] Surface_HF_Visc;
  delete [] Surface_MaxHF_Visc;
  delete [] Surface_Buffet_Metric;

  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        delete [] CSkinFriction[iMarker][iDim];
      }
      delete [] CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }

  if (HeatConjugateVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        delete [] HeatConjugateVar[iMarker][iVertex];
      }
      delete [] HeatConjugateVar[iMarker];
    }
    delete [] HeatConjugateVar;
  }

  if (Buffet_Sensor != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      delete [] Buffet_Sensor[iMarker];
    }
    delete [] Buffet_Sensor;
  }

}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver, CConfig *config, unsigned short iMesh,
                              unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const unsigned long InnerIter   = config->GetInnerIter();
  const unsigned long WFStartIter = config->GetWallFunction_Start_Iter();

  const bool cont_adjoint       = config->GetContinuous_Adjoint();
  const bool disc_adjoint       = config->GetDiscrete_Adjoint();
  const bool limiter_flow       = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  const bool limiter_turb       = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  const bool limiter_adjflow    = (cont_adjoint && (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()));
  const bool edge_limiter_flow  = config->GetEdgeLimiter_Flow();
  const bool edge_limiter_turb  = config->GetEdgeLimiter_Turb();

  const bool restart        = config->GetRestart();
  const bool wall_functions = (config->GetWall_Functions() && ((disc_adjoint) || (InnerIter >= WFStartIter) || (restart)));

  /*--- Common preprocessing steps (implemented by CEulerSolver) ---*/

  CommonPreprocessing(geometry, solver, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  /*--- Compute gradient for MUSCL reconstruction. ---*/

  if (config->GetReconstructionGradientRequired() && (iMesh == MESH_0)) {
    switch (config->GetKind_Gradient_Method_Recon()) {
      case GREEN_GAUSS:
        SetPrimitive_Gradient_GG(geometry, config, true); break;
      case LEAST_SQUARES:
      case WEIGHTED_LEAST_SQUARES:
        SetPrimitive_Gradient_LS(geometry, config, true); break;
      default: break;
    }
  }

  /*--- Compute gradient of the primitive variables ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  else if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  /*--- Compute the limiter in case we need it in the turbulence model or to limit the
   *    viscous terms (check this logic with JST and 2nd order turbulence model) ---*/

  if ((iMesh == MESH_0) && ((limiter_flow && !edge_limiter_flow) || (limiter_turb && !edge_limiter_turb) || limiter_adjflow) && !Output) {
    SetPrimitive_Limiter(geometry, config);
  }

  /*--- Evaluate the vorticity and strain rate magnitude ---*/

  nodes->SetVorticity_StrainMag();

  /*--- Compute the TauWall from the wall functions ---*/

  if (wall_functions) {
    /*--- First reset CFL if needed ---*/
    if ((InnerIter == WFStartIter) && (!restart) && (!disc_adjoint) && !Output) {
      if (rank == MASTER_NODE)
        cout << "---------------------- Switching to wall function -----------------------" << endl;
      ResetCFLAdapt();
      const su2double CFL = config->GetCFL(iMesh);
      for (auto iPoint = 0; iPoint < nPoint; iPoint++)
        nodes->SetLocalCFL(iPoint, CFL);
      Min_CFL_Local = CFL;
      Max_CFL_Local = CFL;
      Avg_CFL_Local = CFL;
    }
    SU2_OMP_MASTER
    ComputeKnoppWallFunction(geometry, solver, config);
    SU2_OMP_BARRIER
  }


}

unsigned long CNSSolver::SetPrimitive_Variables(CSolver **solver, CConfig *config, bool Output) {

  /*--- Number of non-physical points, local to the thread, needs
   *    further reduction if function is called in parallel ---*/
  unsigned long nonPhysicalPoints = 0;

  const unsigned short turb_model = config->GetKind_Turb_Model();
  const bool tkeNeeded = (turb_model == SST) || (turb_model == SST_SUST);

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Retrieve the value of the kinetic energy (if needed). ---*/

    su2double eddy_visc = 0.0, turb_ke = 0.0;

    if (turb_model != NONE && solver[TURB_SOL] != nullptr) {
      eddy_visc = solver[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      if (tkeNeeded) turb_ke = solver[TURB_SOL]->GetNodes()->GetPrimitive(iPoint,0);

      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
        su2double DES_LengthScale = solver[TURB_SOL]->GetNodes()->GetDES_LengthScale(iPoint);
        nodes->SetDES_LengthScale(iPoint, DES_LengthScale);
      }
    }

    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/

    bool physical = static_cast<CNSVariable*>(nodes)->SetPrimVar(iPoint, eddy_visc, turb_ke, GetFluidModel());
    nodes->SetSecondaryVar(iPoint, GetFluidModel());

    /*--- Check for non-realizable states for reporting. ---*/

    nonPhysicalPoints += !physical;

  }

  return nonPhysicalPoints;
}

void CNSSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver,
                                 CNumerics *numerics, CConfig *config) {

  const bool implicit  = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  const bool tkeNeeded = (config->GetKind_Turb_Model() == SST) ||
                         (config->GetKind_Turb_Model() == SST_SUST);

  CVariable* turbNodes = nullptr;
  if (tkeNeeded) turbNodes = solver[TURB_SOL]->GetNodes();

  /*--- Points, coordinates and normal vector in edge ---*/

  const auto edge_i = geometry->edge[iEdge];

  const auto iPoint = edge_i->GetNode(0);
  const auto jPoint = edge_i->GetNode(1);

  numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());

  numerics->SetNormal(edge_i->GetNormal());

  /*--- Primitive and secondary variables. ---*/

  numerics->SetPrimitive(nodes->GetPrimitive(iPoint),
                         nodes->GetPrimitive(jPoint));

  numerics->SetSecondary(nodes->GetSecondary(iPoint),
                         nodes->GetSecondary(jPoint));

  /*--- Gradients. ---*/

  numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                               nodes->GetGradient_Primitive(jPoint));

  /*--- Turbulent kinetic energy. ---*/

  if (tkeNeeded) {
    numerics->SetTurbKineticEnergy(turbNodes->GetPrimitive(iPoint,0),
                                   turbNodes->GetPrimitive(jPoint,0));
    numerics->SetTurbSpecificDissipation(turbNodes->GetPrimitive(iPoint,1),
                                         turbNodes->GetPrimitive(jPoint,1));
    numerics->SetTurbVarGradient(turbNodes->GetGradient(iPoint),
                                 turbNodes->GetGradient(jPoint));
    numerics->SetF1blending(turbNodes->GetF1blending(iPoint),
                            turbNodes->GetF1blending(jPoint));
    numerics->SetF2blending(turbNodes->GetF2blending(iPoint),
                            turbNodes->GetF2blending(jPoint));
    numerics->SetVorticityMag(nodes->GetVorticityMag(iPoint),
                              nodes->GetVorticityMag(jPoint));

  }

  /*--- Gradients. ---*/

  numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                               nodes->GetGradient_Primitive(jPoint));

  /*--- Wall shear stress values (wall functions) ---*/

  numerics->SetTauWall(nodes->GetTauWall(iPoint),
                       nodes->GetTauWall(jPoint));

  /*--- Compute and update residual ---*/

  auto residual = numerics->ComputeResidual(config);

  if (ReducerStrategy) {
    EdgeFluxes.SubtractBlock(iEdge, residual);
    if (implicit)
      Jacobian.UpdateBlocksSub(iEdge, residual.jacobian_i, residual.jacobian_j);
  }
  else {
    LinSysRes.SubtractBlock(iPoint, residual);
    LinSysRes.AddBlock(jPoint, residual);

    if (implicit) {
      Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);

      if (config->GetUse_Accurate_Visc_Jacobians()) {
        CorrectJacobian(solver, geometry, config, iPoint, jPoint, edge_i->GetNormal());
        CorrectJacobian(solver, geometry, config, jPoint, iPoint, edge_i->GetNormal());
      }
    }
  }

}

void CNSSolver::CorrectJacobian(CSolver             **solver,
                                const CGeometry     *geometry,
                                const CConfig       *config,
                                const unsigned long iPoint,
                                const unsigned long jPoint,
                                const su2double     *Normal) {
  
  
  /*--- We're only computing contributions of first neighbors to the Jacobian.
        In Green-Gauss, this contribution is scaled by 0.5*Sum(n_v)/r = 0 for
        volume nodes and (0.5*Sum(n_v)+n_s)/r for surface nodes. So only add to
        the Jacobian if iPoint is on a physical boundary. ---*/

  const bool wasActive = AD::BeginPassive();

  su2double EdgVec[MAXNDIM] = {0.0};
  for (auto iDim = 0; iDim < nDim; iDim++)
    EdgVec[iDim] = geometry->node[jPoint]->GetCoord(iDim)-geometry->node[iPoint]->GetCoord(iDim);

  /*--- Get norm of projection and distance vectors ---*/

  su2double ProjVec = 0.0, Dist2 = 0.0;
  for (auto iDim = 0; iDim < nDim; iDim++){
    ProjVec += Normal[iDim]*EdgVec[iDim];
    Dist2   += EdgVec[iDim]*EdgVec[iDim];
  }

  /*--- Get vector to be multiplied by Jacobian weights ---*/

  su2double Vec[MAXNDIM] = {0.0};
  for (auto iDim = 0; iDim < nDim; iDim++)
    Vec[iDim] = Normal[iDim] - EdgVec[iDim]*ProjVec/Dist2;

  StressTensorJacobian(solver, geometry, config, iPoint, jPoint, Vec);
  HeatFluxJacobian(solver, geometry, config, iPoint, jPoint, Vec);

  AD::EndPassive(wasActive);
}

void CNSSolver::StressTensorJacobian(CSolver             **solver,
                                     const CGeometry     *geometry,
                                     const CConfig       *config,
                                     const unsigned long iPoint,
                                     const unsigned long jPoint,
                                     const su2double     *Vec) {

  const bool gg  = config->GetKind_Gradient_Method() == GREEN_GAUSS;

  const auto node_i = geometry->node[iPoint];

  const su2double sign = 1.0 - 2.0*(iPoint > jPoint);
  const su2double sign_grad_i = -1.0 + 2.0*(gg);
  const su2double delta[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  /*--- Common factors for all Jacobian terms --*/
  const su2double Mean_LaminarVisc = 0.5*(nodes->GetLaminarViscosity(iPoint)+nodes->GetLaminarViscosity(jPoint));
  const su2double Mean_EddyVisc = 0.5*(nodes->GetEddyViscosity(iPoint)+nodes->GetEddyViscosity(jPoint));
  const su2double Mean_Viscosity = Mean_LaminarVisc + Mean_EddyVisc;

  /*--- TODO: Correction with wall function ---*/
  const su2double WF_Factor = nodes->GetTauWallFactor(iPoint);

  const su2double Density_i = nodes->GetDensity(iPoint);
  const su2double Xi_i = WF_Factor*Mean_Viscosity/Density_i;
  

  su2double Mean_Velocity[MAXNDIM] = {0.0};
  for (auto iDim = 0; iDim < nDim; iDim++)
    Mean_Velocity[iDim] = 0.5*(nodes->GetVelocity(iPoint,iDim)+nodes->GetVelocity(jPoint,iDim));

  /*--- Reset first row and last column of Jacobian now so we don't need to later ---*/
  for (auto iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[0][iVar] = Jacobian_i[iVar][nVar-1] = 0.0;
    Jacobian_j[0][iVar] = Jacobian_j[iVar][nVar-1] = 0.0;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Compute contributions of surface terms to the Jacobian.    ---*/
  /*---         In Green-Gauss, the weight of the surface node             ---*/
  /*---         contribution is (0.5*Sum(n_v)+n_s)/r. Least squares        ---*/
  /*---         gradients do not have a surface term.                      ---*/
  /*--------------------------------------------------------------------------*/

  if (gg && node_i->GetPhysicalBoundary()) {

    for (auto iVar = 1; iVar < nVar; iVar++)
      for (auto jVar = 0; jVar < nVar; jVar++)
        Jacobian_i[iVar][jVar] = 0.0;

    const su2double HalfOnVol = 0.5/node_i->GetVolume();
    const su2double factor = -HalfOnVol*sign;
    for (auto iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        const long iVertex = node_i->GetVertex(iMarker);
        if (iVertex != -1) {
          const su2double *gradWeight = geometry->vertex[iMarker][iVertex]->GetNormal();

          /*--- Get projection to be multiplied by divergence terms ---*/
          su2double diagTerm = 0.0;
          for (auto iDim = 0; iDim < nDim; iDim++)
            diagTerm += Vec[iDim]*gradWeight[iDim];

          /*--- Momentum flux Jacobian wrt momentum ---*/
          for (auto iDim = 0; iDim < nDim; iDim++)
            for (auto jDim = 0; jDim < nDim; jDim++)
              Jacobian_i[iDim+1][jDim+1] += factor*Xi_i*(gradWeight[iDim]*Vec[jDim] 
                                          - TWO3*gradWeight[jDim]*Vec[iDim] 
                                          + delta[iDim][jDim]*diagTerm);
        }// iVertex
      }// not send-receive
    }// iMarker

    /*--- Now get density and energy Jacobians for iPoint ---*/
    for (auto iDim = 0; iDim < nDim; iDim++) {
      for (auto jDim = 0; jDim < nDim; jDim++) {
        /*--- Momentum flux Jacobian wrt density ---*/
        Jacobian_i[iDim+1][0] -= Jacobian_i[iDim+1][jDim+1]*nodes->GetVelocity(iPoint,jDim);

        /*--- Energy Jacobian wrt momentum ---*/
        Jacobian_i[nVar-1][iDim+1] += Jacobian_i[jDim+1][iDim+1]*Mean_Velocity[jDim];
      }

      /*--- Energy Jacobian wrt density ---*/
      Jacobian_i[nVar-1][0] -= Jacobian_i[nVar-1][iDim+1]*nodes->GetVelocity(iPoint,iDim);
    }

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
  }// physical boundary

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Compute contributions of neighbor nodes to the Jacobian.   ---*/
  /*---         To prevent extra communication overhead, we only consider  ---*/
  /*---         neighbors on the current rank.                             ---*/
  /*--------------------------------------------------------------------------*/

  su2double gradWeight[MAXNDIM] = {0.0};
  for (auto iNeigh = 0; iNeigh < node_i->GetnPoint(); iNeigh++) {
    const auto kPoint = node_i->GetPoint(iNeigh);
    const su2double Density_k = nodes->GetDensity(kPoint);
    const su2double Xi_k = WF_Factor*Mean_Viscosity/Density_k;
    const su2double factor = 0.5*sign;

    SetGradWeights(gradWeight, solver[FLOW_SOL], geometry, config, iPoint, kPoint);

    /*--- Get projection to be multiplied by divergence terms ---*/
    su2double diagTerm = 0.0;
    for (auto iDim = 0; iDim < nDim; iDim++)
      diagTerm += Vec[iDim]*gradWeight[iDim];

    /*--- Momentum flux Jacobian wrt momentum ---*/
    for (auto iDim = 0; iDim < nDim; iDim++) {
      for (auto jDim = 0; jDim < nDim; jDim++) {
        Jacobian_i[iDim+1][jDim+1] = factor*Xi_i*(gradWeight[iDim]*Vec[jDim] 
                                   - TWO3*gradWeight[jDim]*Vec[iDim] 
                                   + delta[iDim][jDim]*diagTerm)*sign_grad_i;
        Jacobian_j[iDim+1][jDim+1] = factor*Xi_k*(gradWeight[iDim]*Vec[jDim] 
                                   - TWO3*gradWeight[jDim]*Vec[iDim] 
                                   + delta[iDim][jDim]*diagTerm);
      }
    }

    /*--- Now get density and energy Jacobians for kPoint ---*/
    Jacobian_i[nVar-1][0] = Jacobian_j[nVar-1][0] = 0.0;
    for (auto iDim = 0; iDim < nDim; iDim++) {
      Jacobian_i[iDim+1][0] = Jacobian_i[nVar-1][iDim+1] = 0.0;
      Jacobian_j[iDim+1][0] = Jacobian_j[nVar-1][iDim+1] = 0.0;
      for (auto jDim = 0; jDim < nDim; jDim++) {
        /*--- Momentum flux Jacobian wrt density ---*/
        Jacobian_i[iDim+1][0] -= Jacobian_i[iDim+1][jDim+1]*nodes->GetVelocity(iPoint,jDim);
        Jacobian_j[iDim+1][0] -= Jacobian_j[iDim+1][jDim+1]*nodes->GetVelocity(kPoint,jDim);

        /*--- Energy Jacobian wrt momentum ---*/
        Jacobian_i[nVar-1][iDim+1] += Jacobian_i[jDim+1][iDim+1]*Mean_Velocity[jDim];
        Jacobian_j[nVar-1][iDim+1] += Jacobian_j[jDim+1][iDim+1]*Mean_Velocity[jDim];
      }

      /*--- Energy Jacobian wrt density ---*/
      Jacobian_i[nVar-1][0] -= Jacobian_i[nVar-1][iDim+1]*nodes->GetVelocity(iPoint,iDim);
      Jacobian_j[nVar-1][0] -= Jacobian_j[nVar-1][iDim+1]*nodes->GetVelocity(kPoint,iDim);
    }

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, kPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, kPoint, Jacobian_j);
  }// iNeigh

}

void CNSSolver::HeatFluxJacobian(CSolver             **solver,
                                 const CGeometry     *geometry,
                                 const CConfig       *config,
                                 const unsigned long iPoint,
                                 const unsigned long jPoint,
                                 const su2double     *Vec) {

  const bool gg = config->GetKind_Gradient_Method() == GREEN_GAUSS;
  const bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);

  const auto node_i = geometry->node[iPoint];

  CVariable* turbNodes = nullptr;
  if (tkeNeeded) turbNodes = solver[TURB_SOL]->GetNodes();

  const su2double sign = 1.0 - 2.0*(iPoint > jPoint);
  const su2double sign_grad_i = -1.0 + 2.0*(gg);

  /*--- Common factors for all Jacobian terms --*/
  const su2double Mean_LaminarVisc = 0.5*(nodes->GetLaminarViscosity(iPoint)+nodes->GetLaminarViscosity(jPoint));
  const su2double Mean_EddyVisc    = 0.5*(nodes->GetEddyViscosity(iPoint)+nodes->GetEddyViscosity(jPoint));

  const su2double HeatFluxFactor  = Mean_LaminarVisc/Prandtl_Lam + Mean_EddyVisc/Prandtl_Turb;
  const su2double CpOnR           = Gamma/Gamma_Minus_One;
  const su2double ConductivityOnR = CpOnR*HeatFluxFactor;

  const su2double Vel2_i     = nodes->GetVelocity2(iPoint);
  const su2double Density_i  = nodes->GetDensity(iPoint);
  const su2double Pressure_i = nodes->GetPressure(iPoint);
  const su2double Phi_i      = Gamma_Minus_One/Density_i;

  /*--- Additional factor if TKE flux is needed ---*/
  su2double tke_i = 0., tke_visc = 0.;
  if (tkeNeeded) {
    const su2double F1_i = turbNodes->GetF1blending(iPoint);
    const su2double F1_j = turbNodes->GetF1blending(jPoint);
    const su2double sigma_k1 = 0.85, sigma_k2 = 1.0;
    const su2double visc_k_i = 0.5*(F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2)*nodes->GetEddyViscosity(iPoint);
    const su2double visc_k_j = 0.5*(F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2)*nodes->GetEddyViscosity(jPoint);

    tke_visc = (Mean_LaminarVisc + visc_k_i + visc_k_j);
    tke_i = turbNodes->GetPrimitive(iPoint,0);
  }

  /*--- Reset most of Jacobian now so we don't need to later ---*/
  for (auto iVar = 0; iVar < nVar-1; iVar++) {
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

  if (gg && node_i->GetPhysicalBoundary()) {

    for (auto iVar = 0; iVar < nVar; iVar++)
      Jacobian_i[nVar-1][iVar] = 0.0;

    const su2double HalfOnVol = 0.5/node_i->GetVolume();
    for (auto iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        const long iVertex = node_i->GetVertex(iMarker);
        if (iVertex != -1) {
          const su2double *gradWeight = geometry->vertex[iMarker][iVertex]->GetNormal();

          su2double factor = 0.0;
          for (auto iDim = 0; iDim < nDim; iDim++)
            factor += gradWeight[iDim]*Vec[iDim];

          factor *= -HalfOnVol*ConductivityOnR*sign;

          /*--- Density Jacobian ---*/
          Jacobian_i[nVar-1][0] += factor*(-Pressure_i/pow(Density_i,2.0)+0.5*Vel2_i*Phi_i);

          /*--- Momentum Jacobian ---*/
          for (auto jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nVar-1][jDim+1] -= factor*Phi_i*nodes->GetVelocity(iPoint,jDim);

          /*--- Energy Jacobian ---*/
          Jacobian_i[nVar-1][nVar-1] += factor*Phi_i;

          /*--- Tke term ---*/
          if (tkeNeeded) {
            factor *= tke_visc/ConductivityOnR;
            Jacobian_i[nVar-1][0] -= factor*tke_i/Density_i;
          }
        }// iVertex
      }// not send-receive
    }// iMarker

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
  }// physical boundary

  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Compute contributions of neighbor nodes to the Jacobian.   ---*/
  /*---         To prevent extra communication overhead, we only consider  ---*/
  /*---         neighbors on the current rank.                             ---*/
  /*--------------------------------------------------------------------------*/

  su2double gradWeight[MAXNDIM] = {0.0};
  for (auto iNeigh = 0; iNeigh < node_i->GetnPoint(); iNeigh++) {
      
    const auto kPoint = node_i->GetPoint(iNeigh);

    const su2double Vel2_k     = nodes->GetVelocity2(kPoint);
    const su2double Density_k  = nodes->GetDensity(kPoint);
    const su2double Pressure_k = nodes->GetPressure(kPoint);
    const su2double Phi_k      = Gamma_Minus_One/Density_k;

    SetGradWeights(gradWeight, solver[FLOW_SOL], geometry, config, iPoint, kPoint);

    su2double factor = 0.0;
    for (auto iDim = 0; iDim < nDim; iDim++)
      factor += gradWeight[iDim]*Vec[iDim];

    factor *= 0.5*sign*ConductivityOnR;

    /*--- Density Jacobian ---*/
    Jacobian_i[nVar-1][0] = factor*(-Pressure_i/pow(Density_i,2.0)+0.5*Vel2_i*Phi_i)*sign_grad_i;
    Jacobian_j[nVar-1][0] = factor*(-Pressure_k/pow(Density_k,2.0)+0.5*Vel2_k*Phi_k);

    /*--- Momentum Jacobian ---*/
    for (auto jDim = 0; jDim < nDim; jDim++) {
      Jacobian_i[nVar-1][jDim+1] = -factor*Phi_i*nodes->GetVelocity(iPoint,jDim)*sign_grad_i;
      Jacobian_j[nVar-1][jDim+1] = -factor*Phi_k*nodes->GetVelocity(kPoint,jDim);
    }

    /*--- Energy Jacobian ---*/
    Jacobian_i[nVar-1][nVar-1] = factor*Phi_i*sign_grad_i;
    Jacobian_j[nVar-1][nVar-1] = factor*Phi_k;

    /*--- Tke term ---*/
    if (tkeNeeded) {
      const su2double tke_k = turbNodes->GetPrimitive(kPoint,0);
      factor *= tke_visc/ConductivityOnR;
      Jacobian_i[nVar-1][0] -= factor*tke_i/Density_i*sign_grad_i;
      Jacobian_j[nVar-1][0] -= factor*tke_k/Density_k;
    }

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, kPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, kPoint, Jacobian_j);
  }// iNeigh
  
}

void CNSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  su2double Viscosity = 0.0, div_vel, WallDist[3] = {0.0, 0.0, 0.0},
  Area, WallShearStress, TauNormal, factor, RefTemp, RefVel2, RefDensity, GradTemperature, Density = 0.0, WallDistMod, FrictionVel,
  Mach2Vel, Mach_Motion, UnitNormal[3] = {0.0, 0.0, 0.0}, TauElem[3] = {0.0, 0.0, 0.0}, TauTangent[3] = {0.0, 0.0, 0.0},
  Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Cp, thermal_conductivity, MaxNorm = 8.0,
  Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Grad_Temp[3] = {0.0, 0.0, 0.0},
  delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double AxiFactor;
  const su2double *Coord = nullptr, *Coord_Normal = nullptr, *Normal = nullptr;

  string Marker_Tag, Monitoring_Tag;

  su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double RefHeatFlux = config->GetHeat_Flux_Ref();
  su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double *Origin = nullptr;

  if (config->GetnMarker_Monitoring() != 0) { Origin = config->GetRefOriginMoment(0); }

  bool QCR = config->GetQCR();
  bool axisymmetric = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/

  RefTemp = Temperature_Inf;
  RefDensity = Density_Inf;
  if (dynamic_grid) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  } else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*--- Variables initialization ---*/

  AllBoundViscCoeff.setZero();
  SurfaceViscCoeff.setZero();

  AllBound_HF_Visc = 0.0;  AllBound_MaxHF_Visc = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_HF_Visc[iMarker_Monitoring]  = 0.0; Surface_MaxHF_Visc[iMarker_Monitoring]   = 0.0;
  }

  /*--- Loop over the Navier-Stokes markers ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    /*--- Obtain the origin for the moment computation for a particular marker ---*/

    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }

    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL) || (Boundary == HEAT_FLUX) || (Boundary == CHT_WALL_INTERFACE)) {

      /*--- Forces initialization at each Marker ---*/

      ViscCoeff.setZero(iMarker);

      HF_Visc[iMarker] = 0.0;    MaxHF_Visc[iMarker] = 0.0;

      su2double ForceViscous[MAXNDIM] = {0.0}, MomentViscous[MAXNDIM] = {0.0};
      su2double MomentX_Force[MAXNDIM] = {0.0}, MomentY_Force[MAXNDIM] = {0.0}, MomentZ_Force[MAXNDIM] = {0.0};

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        Coord = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPointNormal]->GetCoord();

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(iPoint,iDim+1, jDim);
          }
          Grad_Temp[iDim] = nodes->GetGradient_Primitive(iPoint,0, iDim);
        }

        Viscosity = nodes->GetLaminarViscosity(iPoint);
        Density = nodes->GetDensity(iPoint);

        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);


        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
        }

        /*--- Evaluate Tau ---*/
        
        const su2double wf = nodes->GetTauWallFactor(iPoint);

        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*delta[iDim][jDim];
            Tau[iDim][jDim] *= wf;
          }
        }

        /*--- If necessary evaluate the QCR contribution to Tau ---*/

        if (QCR) {
          su2double den_aux, c_cr1=0.3, O_ik, O_jk;
          unsigned short kDim;

          /*--- Denominator Antisymmetric normalized rotation tensor ---*/

          den_aux = 0.0;
          for (iDim = 0 ; iDim < nDim; iDim++)
            for (jDim = 0 ; jDim < nDim; jDim++)
              den_aux += Grad_Vel[iDim][jDim] * Grad_Vel[iDim][jDim];
          den_aux = sqrt(max(den_aux,1E-10));

          /*--- Adding the QCR contribution ---*/

          for (iDim = 0 ; iDim < nDim; iDim++){
            for (jDim = 0 ; jDim < nDim; jDim++){
              for (kDim = 0 ; kDim < nDim; kDim++){
                O_ik = (Grad_Vel[iDim][kDim] - Grad_Vel[kDim][iDim])/ den_aux;
                O_jk = (Grad_Vel[jDim][kDim] - Grad_Vel[kDim][jDim])/ den_aux;
                Tau[iDim][jDim] -= c_cr1 * (O_ik * Tau[jDim][kDim] + O_jk * Tau[iDim][kDim]);
              }
            }
          }
        }

        /*--- Project Tau in each surface element ---*/

        for (iDim = 0; iDim < nDim; iDim++) {
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++) {
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
          }
        }

        /*--- Compute wall shear stress (using the stress tensor). Compute wall skin friction coefficient, and heat flux on the wall ---*/

        TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];

        WallShearStress = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
          CSkinFriction[iMarker][iDim][iVertex] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);
          WallShearStress += TauTangent[iDim] * TauTangent[iDim];
        }
        WallShearStress = sqrt(WallShearStress);

        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]*UnitNormal[iDim]*UnitNormal[iDim];
        WallDistMod = sqrt(WallDistMod);

        /*--- Compute y+ and non-dimensional velocity ---*/

        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);

        /*--- Compute total and maximum heat flux on the wall ---*/

        GradTemperature = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          GradTemperature -= Grad_Temp[iDim]*UnitNormal[iDim];

        Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
        thermal_conductivity = Cp * Viscosity/Prandtl_Lam;
        HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature*RefHeatFlux;

        /*--- Note that y+, and heat are computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {

          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation ---*/

          su2double Force[MAXNDIM] = {0.0}, MomentDist[MAXNDIM] = {0.0};
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim] * Area * factor * AxiFactor;
            ForceViscous[iDim] += Force[iDim];
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (nDim == 3) {
            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1] += (-Force[1]*Coord[2]);
            MomentX_Force[2] += (Force[2]*Coord[1]);

            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2] += (-Force[2]*Coord[0]);
            MomentY_Force[0] += (Force[0]*Coord[2]);
          }
          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0] += (-Force[0]*Coord[1]);
          MomentZ_Force[1] += (Force[1]*Coord[0]);

        }

        HF_Visc[iMarker]          += HeatFlux[iMarker][iVertex]*Area;
        MaxHF_Visc[iMarker]       += pow(HeatFlux[iMarker][iVertex], MaxNorm);

      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if (Monitoring == YES) {
        if (nDim == 2) {
          ViscCoeff.CD[iMarker]          =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          ViscCoeff.CL[iMarker]          = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          ViscCoeff.CEff[iMarker]        = ViscCoeff.CL[iMarker] / (ViscCoeff.CD[iMarker]+EPS);
          ViscCoeff.CFx[iMarker]         = ForceViscous[0];
          ViscCoeff.CFy[iMarker]         = ForceViscous[1];
          ViscCoeff.CMz[iMarker]         = MomentViscous[2];
          ViscCoeff.CoPx[iMarker]        = MomentZ_Force[1];
          ViscCoeff.CoPy[iMarker]        = -MomentZ_Force[0];
          ViscCoeff.CT[iMarker]          = -ViscCoeff.CFx[iMarker];
          ViscCoeff.CQ[iMarker]          = -ViscCoeff.CMz[iMarker];
          ViscCoeff.CMerit[iMarker]      = ViscCoeff.CT[iMarker] / (ViscCoeff.CQ[iMarker]+EPS);
          MaxHF_Visc[iMarker]            = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          ViscCoeff.CD[iMarker]          =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          ViscCoeff.CL[iMarker]          = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          ViscCoeff.CSF[iMarker]         = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          ViscCoeff.CEff[iMarker]        = ViscCoeff.CL[iMarker]/(ViscCoeff.CD[iMarker] + EPS);
          ViscCoeff.CFx[iMarker]         = ForceViscous[0];
          ViscCoeff.CFy[iMarker]         = ForceViscous[1];
          ViscCoeff.CFz[iMarker]         = ForceViscous[2];
          ViscCoeff.CMx[iMarker]         = MomentViscous[0];
          ViscCoeff.CMy[iMarker]         = MomentViscous[1];
          ViscCoeff.CMz[iMarker]         = MomentViscous[2];
          ViscCoeff.CoPx[iMarker]        = -MomentY_Force[0];
          ViscCoeff.CoPz[iMarker]        = MomentY_Force[2];
          ViscCoeff.CT[iMarker]          = -ViscCoeff.CFz[iMarker];
          ViscCoeff.CQ[iMarker]          = -ViscCoeff.CMz[iMarker];
          ViscCoeff.CMerit[iMarker]      = ViscCoeff.CT[iMarker] / (ViscCoeff.CQ[iMarker] + EPS);
          MaxHF_Visc[iMarker]            = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }

        AllBoundViscCoeff.CD          += ViscCoeff.CD[iMarker];
        AllBoundViscCoeff.CL          += ViscCoeff.CL[iMarker];
        AllBoundViscCoeff.CSF         += ViscCoeff.CSF[iMarker];
        AllBoundViscCoeff.CFx         += ViscCoeff.CFx[iMarker];
        AllBoundViscCoeff.CFy         += ViscCoeff.CFy[iMarker];
        AllBoundViscCoeff.CFz         += ViscCoeff.CFz[iMarker];
        AllBoundViscCoeff.CMx         += ViscCoeff.CMx[iMarker];
        AllBoundViscCoeff.CMy         += ViscCoeff.CMy[iMarker];
        AllBoundViscCoeff.CMz         += ViscCoeff.CMz[iMarker];
        AllBoundViscCoeff.CoPx        += ViscCoeff.CoPx[iMarker];
        AllBoundViscCoeff.CoPy        += ViscCoeff.CoPy[iMarker];
        AllBoundViscCoeff.CoPz        += ViscCoeff.CoPz[iMarker];
        AllBoundViscCoeff.CT          += ViscCoeff.CT[iMarker];
        AllBoundViscCoeff.CQ          += ViscCoeff.CQ[iMarker];
        AllBound_HF_Visc              += HF_Visc[iMarker];
        AllBound_MaxHF_Visc           += pow(MaxHF_Visc[iMarker], MaxNorm);

        /*--- Compute the coefficients per surface ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            SurfaceViscCoeff.CL[iMarker_Monitoring]      += ViscCoeff.CL[iMarker];
            SurfaceViscCoeff.CD[iMarker_Monitoring]      += ViscCoeff.CD[iMarker];
            SurfaceViscCoeff.CSF[iMarker_Monitoring]     += ViscCoeff.CSF[iMarker];
            SurfaceViscCoeff.CEff[iMarker_Monitoring]    += ViscCoeff.CEff[iMarker];
            SurfaceViscCoeff.CFx[iMarker_Monitoring]     += ViscCoeff.CFx[iMarker];
            SurfaceViscCoeff.CFy[iMarker_Monitoring]     += ViscCoeff.CFy[iMarker];
            SurfaceViscCoeff.CFz[iMarker_Monitoring]     += ViscCoeff.CFz[iMarker];
            SurfaceViscCoeff.CMx[iMarker_Monitoring]     += ViscCoeff.CMx[iMarker];
            SurfaceViscCoeff.CMy[iMarker_Monitoring]     += ViscCoeff.CMy[iMarker];
            SurfaceViscCoeff.CMz[iMarker_Monitoring]     += ViscCoeff.CMz[iMarker];
            Surface_HF_Visc[iMarker_Monitoring]          += HF_Visc[iMarker];
            Surface_MaxHF_Visc[iMarker_Monitoring]       += pow(MaxHF_Visc[iMarker],MaxNorm);
          }
        }

      }

    }
  }

  /*--- Update some global coeffients ---*/

  AllBoundViscCoeff.CEff = AllBoundViscCoeff.CL / (AllBoundViscCoeff.CD + EPS);
  AllBoundViscCoeff.CMerit = AllBoundViscCoeff.CT / (AllBoundViscCoeff.CQ + EPS);
  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);


#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    auto Allreduce = [](su2double x) {
      su2double tmp = x; x = 0.0;
      SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      return x;
    };
    AllBoundViscCoeff.CD = Allreduce(AllBoundViscCoeff.CD);
    AllBoundViscCoeff.CL = Allreduce(AllBoundViscCoeff.CL);
    AllBoundViscCoeff.CSF = Allreduce(AllBoundViscCoeff.CSF);
    AllBoundViscCoeff.CEff = AllBoundViscCoeff.CL / (AllBoundViscCoeff.CD + EPS);

    AllBoundViscCoeff.CMx = Allreduce(AllBoundViscCoeff.CMx);
    AllBoundViscCoeff.CMy = Allreduce(AllBoundViscCoeff.CMy);
    AllBoundViscCoeff.CMz = Allreduce(AllBoundViscCoeff.CMz);

    AllBoundViscCoeff.CFx = Allreduce(AllBoundViscCoeff.CFx);
    AllBoundViscCoeff.CFy = Allreduce(AllBoundViscCoeff.CFy);
    AllBoundViscCoeff.CFz = Allreduce(AllBoundViscCoeff.CFz);

    AllBoundViscCoeff.CoPx = Allreduce(AllBoundViscCoeff.CoPx);
    AllBoundViscCoeff.CoPy = Allreduce(AllBoundViscCoeff.CoPy);
    AllBoundViscCoeff.CoPz = Allreduce(AllBoundViscCoeff.CoPz);

    AllBoundViscCoeff.CT = Allreduce(AllBoundViscCoeff.CT);
    AllBoundViscCoeff.CQ = Allreduce(AllBoundViscCoeff.CQ);
    AllBoundViscCoeff.CMerit = AllBoundViscCoeff.CT / (AllBoundViscCoeff.CQ + EPS);

    AllBound_HF_Visc = Allreduce(AllBound_HF_Visc);
    AllBound_MaxHF_Visc = pow(Allreduce(pow(AllBound_MaxHF_Visc, MaxNorm)), 1.0/MaxNorm);

  }

  /*--- Add the forces on the surfaces using all the nodes ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    int nMarkerMon = config->GetnMarker_Monitoring();

    /*--- Use the same buffer for all reductions. We could avoid the copy back into
     *    the original variable by swaping pointers, but it is safer this way... ---*/

    su2double* buffer = new su2double [nMarkerMon];

    auto Allreduce_inplace = [buffer](int size, su2double* x) {
      SU2_MPI::Allreduce(x, buffer, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      for(int i=0; i<size; ++i) x[i] = buffer[i];
    };

    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CL);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CD);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CSF);

    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; iMarker_Monitoring++)
      SurfaceViscCoeff.CEff[iMarker_Monitoring] = SurfaceViscCoeff.CL[iMarker_Monitoring] /
                                                 (SurfaceViscCoeff.CD[iMarker_Monitoring] + EPS);

    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CFx);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CFy);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CFz);

    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CMx);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CMy);
    Allreduce_inplace(nMarkerMon, SurfaceViscCoeff.CMz);

    Allreduce_inplace(nMarkerMon, Surface_HF_Visc);
    Allreduce_inplace(nMarkerMon, Surface_MaxHF_Visc);

    delete [] buffer;

  }

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value)---*/

  TotalCoeff.CD          += AllBoundViscCoeff.CD;
  TotalCoeff.CL          += AllBoundViscCoeff.CL;
  TotalCoeff.CSF         += AllBoundViscCoeff.CSF;
  TotalCoeff.CEff         = TotalCoeff.CL / (TotalCoeff.CD + EPS);
  TotalCoeff.CFx         += AllBoundViscCoeff.CFx;
  TotalCoeff.CFy         += AllBoundViscCoeff.CFy;
  TotalCoeff.CFz         += AllBoundViscCoeff.CFz;
  TotalCoeff.CMx         += AllBoundViscCoeff.CMx;
  TotalCoeff.CMy         += AllBoundViscCoeff.CMy;
  TotalCoeff.CMz         += AllBoundViscCoeff.CMz;
  TotalCoeff.CoPx        += AllBoundViscCoeff.CoPx;
  TotalCoeff.CoPy        += AllBoundViscCoeff.CoPy;
  TotalCoeff.CoPz        += AllBoundViscCoeff.CoPz;
  TotalCoeff.CT          += AllBoundViscCoeff.CT;
  TotalCoeff.CQ          += AllBoundViscCoeff.CQ;
  TotalCoeff.CMerit       = AllBoundViscCoeff.CT / (AllBoundViscCoeff.CQ + EPS);
  Total_Heat         = AllBound_HF_Visc;
  Total_MaxHeat      = AllBound_MaxHF_Visc;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SurfaceCoeff.CL[iMarker_Monitoring]         += SurfaceViscCoeff.CL[iMarker_Monitoring];
    SurfaceCoeff.CD[iMarker_Monitoring]         += SurfaceViscCoeff.CD[iMarker_Monitoring];
    SurfaceCoeff.CSF[iMarker_Monitoring]        += SurfaceViscCoeff.CSF[iMarker_Monitoring];
    SurfaceCoeff.CEff[iMarker_Monitoring]        = SurfaceViscCoeff.CL[iMarker_Monitoring] / (SurfaceCoeff.CD[iMarker_Monitoring] + EPS);
    SurfaceCoeff.CFx[iMarker_Monitoring]        += SurfaceViscCoeff.CFx[iMarker_Monitoring];
    SurfaceCoeff.CFy[iMarker_Monitoring]        += SurfaceViscCoeff.CFy[iMarker_Monitoring];
    SurfaceCoeff.CFz[iMarker_Monitoring]        += SurfaceViscCoeff.CFz[iMarker_Monitoring];
    SurfaceCoeff.CMx[iMarker_Monitoring]        += SurfaceViscCoeff.CMx[iMarker_Monitoring];
    SurfaceCoeff.CMy[iMarker_Monitoring]        += SurfaceViscCoeff.CMy[iMarker_Monitoring];
    SurfaceCoeff.CMz[iMarker_Monitoring]        += SurfaceViscCoeff.CMz[iMarker_Monitoring];
  }

}

void CNSSolver::Buffet_Monitoring(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim;
  su2double *Vel_FS = config->GetVelocity_FreeStream();
  su2double VelMag_FS = 0.0, SkinFrictionMag = 0.0, SkinFrictionDot = 0.0, *Normal, Area, Sref = config->GetRefArea();
  su2double k   = config->GetBuffet_k(),
             lam = config->GetBuffet_lambda();
  string Marker_Tag, Monitoring_Tag;

  for (iDim = 0; iDim < nDim; iDim++){
    VelMag_FS += Vel_FS[iDim]*Vel_FS[iDim];
  }
  VelMag_FS = sqrt(VelMag_FS);

  /*-- Variables initialization ---*/

  Total_Buffet_Metric = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_Buffet_Metric[iMarker_Monitoring] = 0.0;
  }

  /*--- Loop over the Euler and Navier-Stokes markers ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Buffet_Metric[iMarker] = 0.0;

    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL) || (Boundary == HEAT_FLUX) || (Boundary == CHT_WALL_INTERFACE)) {

      /*--- Loop over the vertices to compute the buffet sensor ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Perform dot product of skin friction with freestream velocity ---*/

        SkinFrictionMag = 0.0;
        SkinFrictionDot = 0.0;
        for(iDim = 0; iDim < nDim; iDim++){
          SkinFrictionMag += CSkinFriction[iMarker][iDim][iVertex]*CSkinFriction[iMarker][iDim][iVertex];
          SkinFrictionDot += CSkinFriction[iMarker][iDim][iVertex]*Vel_FS[iDim];
        }
        SkinFrictionMag = sqrt(SkinFrictionMag);

        /*--- Normalize the dot product ---*/

        SkinFrictionDot /= SkinFrictionMag*VelMag_FS;

        /*--- Compute Heaviside function ---*/

        Buffet_Sensor[iMarker][iVertex] = 1./(1. + exp(2.*k*(SkinFrictionDot + lam)));

        /*--- Integrate buffet sensor ---*/

        if(Monitoring == YES){

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0;
          for(iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);

          Buffet_Metric[iMarker] += Buffet_Sensor[iMarker][iVertex]*Area/Sref;

        }

      }

      if(Monitoring == YES){

        Total_Buffet_Metric += Buffet_Metric[iMarker];

        /*--- Per surface buffet metric ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) Surface_Buffet_Metric[iMarker_Monitoring] = Buffet_Metric[iMarker];
        }

      }

    }

  }

#ifdef HAVE_MPI

  /*--- Add buffet metric information using all the nodes ---*/

  su2double MyTotal_Buffet_Metric = Total_Buffet_Metric;
  SU2_MPI::Allreduce(&MyTotal_Buffet_Metric, &Total_Buffet_Metric, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /*--- Add the buffet metric on the surfaces using all the nodes ---*/

  su2double *MySurface_Buffet_Metric = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_Buffet_Metric[iMarker_Monitoring] = Surface_Buffet_Metric[iMarker_Monitoring];
  }

  SU2_MPI::Allreduce(MySurface_Buffet_Metric, Surface_Buffet_Metric, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete [] MySurface_Buffet_Metric;

#endif

}

void CNSSolver::Evaluate_ObjFunc(CConfig *config) {

  unsigned short iMarker_Monitoring, Kind_ObjFunc;
  su2double Weight_ObjFunc;

  /*--- Evaluate objective functions common to Euler and NS solvers ---*/

  CEulerSolver::Evaluate_ObjFunc(config);

  /*--- Evaluate objective functions specific to NS solver ---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {

    Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
    Kind_ObjFunc = config->GetKind_ObjFunc(iMarker_Monitoring);

    switch(Kind_ObjFunc) {
      case BUFFET_SENSOR:
          Total_ComboObj +=Weight_ObjFunc*Surface_Buffet_Metric[iMarker_Monitoring];
          break;
      default:
          break;
    }
  }

}


void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver, CNumerics *conv_numerics,
                                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  su2double *Normal, Area, UnitNormal[3] = {0.0};
  su2double *GridVel;

  const bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary by string name ---*/

  const string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config as well as the
        wall function treatment.---*/

  su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag)/config->GetHeat_Flux_Ref();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (auto iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- If it is a customizable patch, retrieve the specified wall heat flux. ---*/

      if (config->GetMarker_All_PyCustom(val_marker)) Wall_HeatFlux = geometry->GetCustomBoundaryHeatFlux(val_marker, iVertex);

      /*--- Compute dual-grid area and boundary normal ---*/

      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = 0.0;
      for (auto iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (auto iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Initialize the convective & viscous residuals to zero ---*/

      for (auto iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/

      if (dynamic_grid) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (auto iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (auto iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }

      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/

      nodes->SetVelocity_Old(iPoint,Vector);

      for (auto iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);

      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/

      Res_Visc[nDim+1] = Wall_HeatFlux * Area;

      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

      if (dynamic_grid) {

        su2double tau_vel[3] = {0.0}, Grad_Vel[3][3] = {0.0}, tau[3][3] = {0.0},
                  delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

        /*--- Get the grid velocity at the current boundary node ---*/

        su2double ProjGridVel = 0.0;
        for (auto iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;

        /*--- Retrieve other primitive quantities and viscosities ---*/

        const su2double Density  = nodes->GetDensity(iPoint);
        const su2double Pressure = nodes->GetPressure(iPoint);
        const su2double laminar_viscosity = nodes->GetLaminarViscosity(iPoint);
        const su2double eddy_viscosity    = nodes->GetEddyViscosity(iPoint);
        const su2double total_viscosity   = laminar_viscosity + eddy_viscosity;

        for (auto iDim = 0; iDim < nDim; iDim++) {
          for (auto jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(iPoint,iDim+1, jDim);
          }
        }

        /*--- Divergence of the velocity ---*/

        su2double div_vel = 0.0; for (auto iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

        /*--- Compute the viscous stress tensor ---*/

        for (auto iDim = 0; iDim < nDim; iDim++) {
          for (auto jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim]+Grad_Vel[iDim][jDim] ) -
                              TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
        }

        /*--- Dot product of the stress tensor with the grid velocity ---*/

        for (auto iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (auto jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }

        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/

        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (auto iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;

        /*--- Implicit Jacobian contributions due to moving walls ---*/

        if (implicit) {

          /*--- Jacobian contribution related to the pressure term ---*/

          su2double GridVel2 = 0.0;
          for (auto iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (auto iVar = 0; iVar < nVar; iVar++)
            for (auto jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (auto jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;

          /*--- Add the block to the Global Jacobian structure ---*/

          Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

          /*--- Now the Jacobian contribution related to the shear stress ---*/

          for (auto iVar = 0; iVar < nVar; iVar++)
            for (auto jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;

          /*--- Compute closest normal neighbor ---*/

          const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

          /*--- Get coordinates of i & nearest normal and compute distance ---*/

          const auto Coord_i = geometry->node[iPoint]->GetCoord();
          const auto Coord_j = geometry->node[Point_Normal]->GetCoord();

          su2double dist_ij = 0;
          for (auto iDim = 0; iDim < nDim; iDim++)
            dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
          dist_ij = sqrt(dist_ij);

          su2double theta2 = 0.0;
          for (auto iDim = 0; iDim < nDim; iDim++)
            theta2 += UnitNormal[iDim]*UnitNormal[iDim];

          const su2double factor = total_viscosity*Area/(Density*dist_ij);

          if (nDim == 2) {
            const su2double thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            const su2double thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;

            const su2double etaz   = UnitNormal[0]*UnitNormal[1]/3.0;

            const su2double pix = GridVel[0]*thetax + GridVel[1]*etaz;
            const su2double piy = GridVel[0]*etaz   + GridVel[1]*thetay;

            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          } else {
            const su2double thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            const su2double thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            const su2double thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;

            const su2double etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            const su2double etax = UnitNormal[1]*UnitNormal[2]/3.0;
            const su2double etay = UnitNormal[0]*UnitNormal[2]/3.0;

            const su2double pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            const su2double piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            const su2double piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;

            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }

          /*--- Subtract the block from the Global Jacobian structure ---*/

          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

        }
      }

      /*--- Convective contribution to the residual at the wall ---*/

      LinSysRes.AddBlock(iPoint, Res_Conv);

      /*--- Viscous contribution to the residual at the wall ---*/

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

      if (implicit) {
        for (auto iVar = 1; iVar <= nDim; iVar++) {
          const unsigned long total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }

    }
  }
}

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver, CNumerics *conv_numerics,
                                   CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iVertex, iPoint, Point_Normal, total_index;

  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij, theta2;
  su2double Twall, dTdn, dTdrho, thermal_conductivity;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, Pressure = 0.0, Density, Vel2;
  su2double total_viscosity, div_vel, tau_vel[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};
  su2double laminar_viscosity, eddy_viscosity, Grad_Vel[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
  tau[3][3] = {{0.0, 0.0, 0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}, delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature from config
        as well as the wall function treatment.---*/

  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  /*--- Loop over boundary points ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- If it is a customizable patch, retrieve the specified wall temperature. ---*/

      if (config->GetMarker_All_PyCustom(val_marker)) Twall = geometry->GetCustomBoundaryTemperature(val_marker, iVertex);

      /*--- Compute dual-grid area and boundary normal ---*/

      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Calculate useful quantities ---*/

      theta2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        theta2 += UnitNormal[iDim]*UnitNormal[iDim];

      /*--- Compute closest normal neighbor ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Get coordinates of i & nearest normal and compute distance ---*/

      Coord_i = geometry->node[iPoint]->GetCoord();
      Coord_j = geometry->node[Point_Normal]->GetCoord();
      dist_ij = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      dist_ij = sqrt(dist_ij);

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there is grid motion (v = u_wall)---*/

      if (dynamic_grid) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      }
      else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }

      /*--- Initialize the convective & viscous residuals to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }

      /*--- Set the residual, truncation error and velocity value on the boundary ---*/

      nodes->SetVelocity_Old(iPoint,Vector);

      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);

      /*--- Compute the normal gradient in temperature using Twall ---*/

      dTdn = -(nodes->GetTemperature(Point_Normal) - Twall)/dist_ij;

      /*--- Get transport coefficients ---*/

      laminar_viscosity    = nodes->GetLaminarViscosity(iPoint);
      eddy_viscosity       = nodes->GetEddyViscosity(iPoint);
      thermal_conductivity = Cp * ( laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);

      // work in progress on real-gases...
      //thermal_conductivity = nodes->GetThermalConductivity(iPoint);
      //Cp = nodes->GetSpecificHeatCp(iPoint);
      //thermal_conductivity += Cp*eddy_viscosity/Prandtl_Turb;

      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/

      Res_Visc[nDim+1] = thermal_conductivity * dTdn * Area;

      /*--- Calculate Jacobian for implicit time stepping ---*/

      if (implicit) {

        for (iVar = 0; iVar < nVar; iVar ++)
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_i[iVar][jVar] = 0.0;

        /*--- Calculate useful quantities ---*/

        Density = nodes->GetDensity(iPoint);
        Vel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Vel2 += pow(nodes->GetVelocity(iPoint,iDim),2);
        dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );

        /*--- Enforce the no-slip boundary condition in a strong way ---*/

        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }

        /*--- Add contributions to the Jacobian from the weak enforcement of the energy equations ---*/

        Jacobian_i[nDim+1][0]      = -thermal_conductivity*theta2/dist_ij * dTdrho * Area;
        Jacobian_i[nDim+1][nDim+1] = -thermal_conductivity*theta2/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;

        /*--- Subtract the block from the Global Jacobian structure ---*/

        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

      }

      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

      if (dynamic_grid) {

        /*--- Get the grid velocity at the current boundary node ---*/

        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;

        /*--- Retrieve other primitive quantities and viscosities ---*/

        Density  = nodes->GetDensity(iPoint);
        Pressure = nodes->GetPressure(iPoint);
        laminar_viscosity = nodes->GetLaminarViscosity(iPoint);
        eddy_viscosity    = nodes->GetEddyViscosity(iPoint);

        total_viscosity   = laminar_viscosity + eddy_viscosity;

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(iPoint,iDim+1, jDim);
          }
        }

        /*--- Divergence of the velocity ---*/

        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

        /*--- Compute the viscous stress tensor ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim] ) -
                                TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }

        /*--- Dot product of the stress tensor with the grid velocity ---*/

        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }

        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/

        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;

        /*--- Implicit Jacobian contributions due to moving walls ---*/

        if (implicit) {

          /*--- Jacobian contribution related to the pressure term ---*/

          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;

          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;

          /*--- Add the block to the Global Jacobian structure ---*/

          Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

          /*--- Now the Jacobian contribution related to the shear stress ---*/

          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;

          factor = total_viscosity*Area/(Density*dist_ij);

          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;

            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;

            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;

            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          }
          else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;

            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;

            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;

            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }

          /*--- Subtract the block from the Global Jacobian structure ---*/

          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
        }

      }

      /*--- Convective contribution to the residual at the wall ---*/

      LinSysRes.AddBlock(iPoint, Res_Conv);

      /*--- Viscous contribution to the residual at the wall ---*/

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }

    }
  }
}

void CNSSolver::SetRoe_Dissipation(CGeometry *geometry, CConfig *config){

  const unsigned short kind_roe_dissipation = config->GetKind_RoeLowDiss();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0; iPoint < nPoint; iPoint++) {

    if (kind_roe_dissipation == FD || kind_roe_dissipation == FD_DUCROS){

      su2double wall_distance = geometry->node[iPoint]->GetWall_Distance();

      nodes->SetRoe_Dissipation_FD(iPoint, wall_distance);

    } else if (kind_roe_dissipation == NTS || kind_roe_dissipation == NTS_DUCROS) {

      const su2double delta = geometry->node[iPoint]->GetMaxLength();
      assert(delta > 0 && "Delta must be initialized and non-negative");
      nodes->SetRoe_Dissipation_NTS(iPoint, delta, config->GetConst_DES());
    }
  }

}

void CNSSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver, CNumerics *conv_numerics,
                                           CConfig *config, unsigned short val_marker) {

  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iVertex, iPoint, Point_Normal, total_index;

  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij, theta2;
  su2double Twall= 0.0, There, dTdn= 0.0, dTdrho, thermal_conductivity, Tconjugate, HF_FactorHere, HF_FactorConjugate;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, Pressure = 0.0, Density, Vel2;
  su2double total_viscosity, div_vel, tau_vel[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};
  su2double laminar_viscosity, eddy_viscosity, Grad_Vel[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
  tau[3][3] = {{0.0, 0.0, 0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}, delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  su2double Temperature_Ref = config->GetTemperature_Ref();

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Loop over boundary points ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Compute dual-grid area and boundary normal ---*/

      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Calculate useful quantities ---*/

      theta2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        theta2 += UnitNormal[iDim]*UnitNormal[iDim];

      /*--- Compute closest normal neighbor ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Get coordinates of i & nearest normal and compute distance ---*/

      Coord_i = geometry->node[iPoint]->GetCoord();
      Coord_j = geometry->node[Point_Normal]->GetCoord();
      dist_ij = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      dist_ij = sqrt(dist_ij);

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there is grid motion (v = u_wall)---*/

      if (dynamic_grid) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      }
      else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }

      /*--- Initialize the convective & viscous residuals to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }

      /*--- Set the residual, truncation error and velocity value on the boundary ---*/

      nodes->SetVelocity_Old(iPoint,Vector);

      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);

      /*--- Get transport coefficients ---*/

      laminar_viscosity    = nodes->GetLaminarViscosity(iPoint);
      eddy_viscosity       = nodes->GetEddyViscosity(iPoint);
      thermal_conductivity = Cp * ( laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);

      // work in progress on real-gases...
      //thermal_conductivity = nodes->GetThermalConductivity(iPoint);
      //Cp = nodes->GetSpecificHeatCp(iPoint);
      //thermal_conductivity += Cp*eddy_viscosity/Prandtl_Turb;

      /*--- Compute the normal gradient in temperature using Twall ---*/

      There = nodes->GetTemperature(Point_Normal);
      Tconjugate = GetConjugateHeatVariable(val_marker, iVertex, 0)/Temperature_Ref;

      if ((config->GetKind_CHT_Coupling() == AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX) ||
          (config->GetKind_CHT_Coupling() == AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

        /*--- Compute wall temperature from both temperatures ---*/

        HF_FactorHere = thermal_conductivity*config->GetViscosity_Ref()/dist_ij;
        HF_FactorConjugate = GetConjugateHeatVariable(val_marker, iVertex, 2);

        Twall = (There*HF_FactorHere + Tconjugate*HF_FactorConjugate)/(HF_FactorHere + HF_FactorConjugate);
        dTdn = -(There - Twall)/dist_ij;
      }
      else if ((config->GetKind_CHT_Coupling() == DIRECT_TEMPERATURE_NEUMANN_HEATFLUX) ||
              (config->GetKind_CHT_Coupling() == DIRECT_TEMPERATURE_ROBIN_HEATFLUX)) {

        /*--- (Directly) Set wall temperature to conjugate temperature. ---*/

        Twall = Tconjugate;
        dTdn = -(There - Twall)/dist_ij;
      }
      else {
        Twall = dTdn = 0.0;
        SU2_MPI::Error("Unknown CHT coupling method.", CURRENT_FUNCTION);
      }

      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/

      Res_Visc[nDim+1] = thermal_conductivity * dTdn * Area;

      /*--- Calculate Jacobian for implicit time stepping ---*/

      if (implicit) {

        for (iVar = 0; iVar < nVar; iVar ++)
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_i[iVar][jVar] = 0.0;

        /*--- Calculate useful quantities ---*/

        Density = nodes->GetDensity(iPoint);
        Vel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Vel2 += pow(nodes->GetVelocity(iPoint,iDim),2);
        dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );

        /*--- Enforce the no-slip boundary condition in a strong way ---*/

        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }

        /*--- Add contributions to the Jacobian from the weak enforcement of the energy equations ---*/

        Jacobian_i[nDim+1][0]      = -thermal_conductivity*theta2/dist_ij * dTdrho * Area;
        Jacobian_i[nDim+1][nDim+1] = -thermal_conductivity*theta2/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;

        /*--- Subtract the block from the Global Jacobian structure ---*/

        Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

      }

      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

      if (dynamic_grid) {

        /*--- Get the grid velocity at the current boundary node ---*/

        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;

        /*--- Retrieve other primitive quantities and viscosities ---*/

        Density  = nodes->GetDensity(iPoint);
        Pressure = nodes->GetPressure(iPoint);
        laminar_viscosity = nodes->GetLaminarViscosity(iPoint);
        eddy_viscosity    = nodes->GetEddyViscosity(iPoint);

        total_viscosity   = laminar_viscosity + eddy_viscosity;

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(iPoint,iDim+1, jDim);
          }
        }

        /*--- Divergence of the velocity ---*/

        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

        /*--- Compute the viscous stress tensor ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim] )
                              - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }

        /*--- Dot product of the stress tensor with the grid velocity ---*/

        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }

        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/

        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;

        /*--- Implicit Jacobian contributions due to moving walls ---*/

        if (implicit) {

          /*--- Jacobian contribution related to the pressure term ---*/

          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;

          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;

          /*--- Add the block to the Global Jacobian structure ---*/

          Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

          /*--- Now the Jacobian contribution related to the shear stress ---*/

          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;

          factor = total_viscosity*Area/(Density*dist_ij);

          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;

            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;

            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;

            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          }
          else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;

            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;

            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;

            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }

          /*--- Subtract the block from the Global Jacobian structure ---*/

          Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
        }

      }

      /*--- Convective contribution to the residual at the wall ---*/

      LinSysRes.AddBlock(iPoint, Res_Conv);

      /*--- Viscous contribution to the residual at the wall ---*/

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
    }
  }
}

void CNSSolver::ComputeNicholsWallFunction(CGeometry *geometry, CSolver **solver, CConfig *config) {

  unsigned short iDim, jDim, iMarker;
  unsigned long iVertex, iPoint, Point_Normal, counter;

  su2double Area, div_vel, UnitNormal[3], *Normal;
  su2double **grad_primvar, tau[3][3];

  su2double Vel[3] = {0.0, 0.0, 0.0}, VelNormal, VelTang[3], VelTangMod, WallDist[3], WallDistMod;
  su2double T_Normal, P_Normal;
  su2double Density_Wall, T_Wall, P_Wall, Lam_Visc_Wall, Tau_Wall = 0.0;
  su2double *Coord, *Coord_Normal;
  su2double diff, Delta, grad_diff;
  su2double U_Tau, U_Plus, Gam, Beta, Phi, Q, Y_Plus_White, Y_Plus;
  su2double TauElem[3], TauNormal, TauTangent[3], WallShearStress;
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  const unsigned short max_iter = 1000;
  const su2double tol = 1e-12;
  bool converged = true;

  const unsigned short turb_model = config->GetKind_Turb_Model();
  const bool tkeNeeded = (turb_model == SST) || (turb_model == SST_SUST);
  
  /*--- Compute the recovery factor ---*/
  // Double-check: laminar or turbulent Pr for this?
  const su2double Recovery = pow(config->GetPrandtl_Lam(), (1.0/3.0));

  /*--- Typical constants from boundary layer theory ---*/

  const su2double kappa = 0.41;
  const su2double B = 5.0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {

      /*--- Identify the boundary by string name ---*/

      string Marker_Tag = config->GetMarker_All_TagBound(iMarker);

      /*--- Get the specified wall heat flux from config ---*/

      // Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

      /*--- Loop over all of the vertices on this boundary marker ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        /*--- Check if the node belongs to the domain (i.e, not a halo node)
         and the neighbor is not part of the physical boundary ---*/

        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Get coordinates of the current vertex and nearest normal point ---*/

          Coord = geometry->node[iPoint]->GetCoord();
          Coord_Normal = geometry->node[Point_Normal]->GetCoord();

          /*--- Compute dual-grid area and boundary normal ---*/

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt (Area);

          for (iDim = 0; iDim < nDim; iDim++)
            UnitNormal[iDim] = -Normal[iDim]/Area;
          
          if (geometry->vertex[iMarker][iVertex]->GetDonorFound()){
          
            /*--- Get the distance to the exchange location ---*/
            
            const su2double *doubleInfo = config->GetWallFunction_DoubleInfo(Marker_Tag);
            WallDistMod = doubleInfo[0];

            /*--- Get the density, velocity, and energy at the exchange location ---*/
            
            su2double Density_Normal = 0., Energy_Normal = 0., Kine_Normal = 0., VelMod = 0.;
            P_Normal = 0.;
            T_Normal = 0.;
            for (iDim = 0; iDim < nDim; iDim++) Vel[iDim] = 0.;
            
            const unsigned short nDonors = geometry->vertex[iMarker][iVertex]->GetnDonorPoints();
            for (auto iNode = 0; iNode < nDonors; iNode++) {
              const unsigned long donorPoint = geometry->vertex[iMarker][iVertex]->GetInterpDonorPoint(iNode);
              const su2double donorCoeff     = geometry->vertex[iMarker][iVertex]->GetDonorCoeff(iNode);
              
              Density_Normal  += donorCoeff*nodes->GetDensity(donorPoint);
              Energy_Normal   += donorCoeff*nodes->GetEnergy(donorPoint);
              if (tkeNeeded && solver[TURB_SOL] != nullptr)
                Kine_Normal   += donorCoeff*solver[TURB_SOL]->GetNodes()->GetPrimitive(donorPoint,0);
              
              for (iDim = 0; iDim < nDim; iDim++) Vel[iDim] += donorCoeff*nodes->GetVelocity(donorPoint,iDim);
            }

            for (iDim = 0; iDim < nDim; iDim++) VelMod += Vel[iDim]*Vel[iDim];
            Energy_Normal -= Kine_Normal + 0.5*VelMod;
            
            GetFluidModel()->SetTDState_rhoe(Density_Normal, Energy_Normal);
            P_Normal = GetFluidModel()->GetPressure();
            T_Normal = GetFluidModel()->GetTemperature();

            /*--- Compute the wall-parallel velocity at the exchange location ---*/

            VelNormal = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelNormal += Vel[iDim] * UnitNormal[iDim];
            for (iDim = 0; iDim < nDim; iDim++)
              VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

            VelTangMod = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelTangMod += VelTang[iDim]*VelTang[iDim];
            VelTangMod = sqrt(VelTangMod);
          }
          else {
            /*--- Get the velocity, pressure, and temperature at the nearest
             (normal) interior point. ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              Vel[iDim] = nodes->GetVelocity(Point_Normal,iDim);
            P_Normal = nodes->GetPressure(Point_Normal);
            T_Normal = nodes->GetTemperature(Point_Normal);

            /*--- Compute the wall-parallel velocity at first point off the wall ---*/

            VelNormal = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelNormal += Vel[iDim] * UnitNormal[iDim];
            for (iDim = 0; iDim < nDim; iDim++)
              VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

            VelTangMod = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelTangMod += VelTang[iDim]*VelTang[iDim];
            VelTangMod = sqrt(VelTangMod);

            /*--- Compute normal distance of the interior point from the wall ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);

            WallDistMod = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              WallDistMod += WallDist[iDim]*WallDist[iDim];
            WallDistMod = sqrt(WallDistMod);
          }

          /*--- Compute the wall temperature using the Crocco-Buseman equation ---*/

          T_Wall = T_Normal + Recovery*pow(VelTangMod,2.0)/(2.0*Cp);

          /*--- Extrapolate the pressure from the interior & compute the
           wall density using the equation of state ---*/

          P_Wall = P_Normal;
          Density_Wall = P_Wall/(Gas_Constant*T_Wall);

          /*--- Compute the shear stress at the wall in the regular fashion
           by using the stress tensor on the surface ---*/
          
          const su2double StaticEnergy_Wall = T_Wall*Gas_Constant/Gamma_Minus_One;
          GetFluidModel()->SetTDState_rhoe(Density_Wall, StaticEnergy_Wall);
          Lam_Visc_Wall = GetFluidModel()->GetLaminarViscosity();

          grad_primvar  = nodes->GetGradient_Primitive(iPoint);

          div_vel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            div_vel += grad_primvar[iDim+1][iDim];

          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0 ; jDim < nDim; jDim++) {
              Delta = 0.0; if (iDim == jDim) Delta = 1.0;
              tau[iDim][jDim] = Lam_Visc_Wall*(  grad_primvar[jDim+1][iDim]
                                               + grad_primvar[iDim+1][jDim]) -
              TWO3*Lam_Visc_Wall*div_vel*Delta;
            }
            TauElem[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++)
              TauElem[iDim] += tau[iDim][jDim]*UnitNormal[jDim];
          }

          /*--- Compute wall shear stress as the magnitude of the wall-tangential
           component of the shear stress tensor---*/

          TauNormal = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            TauNormal += TauElem[iDim] * UnitNormal[iDim];

          for (iDim = 0; iDim < nDim; iDim++)
            TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];

          WallShearStress = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            WallShearStress += TauTangent[iDim]*TauTangent[iDim];
          WallShearStress = sqrt(WallShearStress);

          /*--- Calculate the quantities from boundary layer theory and
           iteratively solve for a new wall shear stress. Use the current wall
           shear stress as a starting guess for the wall function. ---*/

          counter = 0; diff = 1.0;
          U_Tau = sqrt(WallShearStress/Density_Wall);
          Y_Plus = 0.0; // to avoid warning
          converged = true;

          su2double Y_Plus_Start = Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall;

          /*--- Automatic switch off when y+ < 5 according to Nichols & Nelson (2004) ---*/

          if (Y_Plus_Start < 5.0 || Y_Plus_Start > 1.0e3) {
            nodes->SetTauWall(iPoint,-1.0);
            nodes->SetTauWallFactor(iPoint,1.0);
            continue;
          }

          while (fabs(diff) > tol) {

            /*--- Friction velocity and u+ ---*/

            U_Plus = VelTangMod/U_Tau;

            /*--- Gamma, Beta, Q, and Phi, defined by Nichols & Nelson (2004) ---*/

            Gam  = Recovery*U_Tau*U_Tau/(2.0*Cp*T_Wall);
            Beta = 0.0; // For adiabatic flows only
            Q    = sqrt(Beta*Beta + 4.0*Gam);
            Phi  = asin(-1.0*Beta/Q);

            /*--- Y+ defined by White & Christoph (compressibility and heat transfer) negative value for (2.0*Gam*U_Plus - Beta)/Q ---*/

            Y_Plus_White = exp((kappa/sqrt(Gam))*(asin((2.0*Gam*U_Plus - Beta)/Q) - Phi))*exp(-1.0*kappa*B);

            /*--- Spalding's universal form for the BL velocity with the
             outer velocity form of White & Christoph above. ---*/

            Y_Plus = U_Plus + Y_Plus_White - (exp(-1.0*kappa*B)*
                                              (1.0 + kappa*U_Plus + kappa*kappa*U_Plus*U_Plus/2.0 +
                                               kappa*kappa*kappa*U_Plus*U_Plus*U_Plus/6.0));

            /* --- Define function for Newton method to zero --- */

            diff = (Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall) - Y_Plus;

            /* --- Gradient of function defined above --- */
            grad_diff = Density_Wall * WallDistMod / Lam_Visc_Wall + VelTangMod / (U_Tau * U_Tau) +
                      kappa /(U_Tau * sqrt(Gam)) * asin(U_Plus * sqrt(Gam)) * Y_Plus_White -
                      exp(-1.0 * B * kappa) * (0.5 * pow(VelTangMod * kappa / U_Tau, 3) +
                      pow(VelTangMod * kappa / U_Tau, 2) + VelTangMod * kappa / U_Tau) / U_Tau;

            /* --- Newton Step --- */

            U_Tau = U_Tau - diff / grad_diff;

            counter++;
            if (counter == max_iter) {
              converged = false;
              break;
            }

          }
          /* --- If not converged or Y+ too large, jump to the next point. --- */

          if (!converged || Y_Plus > 1.0e3) {
            nodes->SetTauWall(iPoint,-1.0);
            nodes->SetTauWallFactor(iPoint,1.0);
            continue;
          }

          /*--- Calculate an updated value for the wall shear stress
            using the y+ value, the definition of y+, and the definition of
            the friction velocity. ---*/

          Tau_Wall = (1.0/Density_Wall)*pow(Y_Plus*Lam_Visc_Wall/WallDistMod,2.0);

          /*--- Store this value for the wall shear stress at the node.  ---*/

          nodes->SetTauWall(iPoint,Tau_Wall);
          nodes->SetTauWallFactor(iPoint,Tau_Wall/WallShearStress);
          nodes->SetTemperature(iPoint,T_Wall);
          nodes->SetSolution(iPoint, 0, Density_Wall);
          nodes->SetPrimitive(iPoint, nDim+2, Density_Wall);
          nodes->SetPrimitive(iPoint, nDim+1, P_Wall);
          nodes->SetLaminarViscosity(iPoint, Lam_Visc_Wall);

        }

      }

    }
  }
  
  InitiateComms(geometry, config, PRIMITIVE);
  CompleteComms(geometry, config, PRIMITIVE);
  
  /*--- Communicate values needed for WF ---*/
  InitiateComms(geometry, config, WALL_FUNCTION);
  CompleteComms(geometry, config, WALL_FUNCTION);
  WallFunctionComms(geometry, solver, config);

}

void CNSSolver::ComputeKnoppWallFunction(CGeometry *geometry, CSolver **solver, CConfig *config) {

  unsigned short iDim, jDim, iMarker;
  unsigned long iVertex, iPoint, Point_Normal, counter;

  su2double Area, div_vel, UnitNormal[3], *Normal;
  su2double **grad_primvar, tau[3][3];

  su2double Vel[3] = {0.0, 0.0, 0.0}, VelNormal, VelTang[3], VelTangMod, WallDist[3], WallDistMod;
  su2double T_Normal, P_Normal, Lam_Visc_Normal;
  su2double Density_Wall, T_Wall, P_Wall, Lam_Visc_Wall, Tau_Wall = 0.0;
  su2double *Coord, *Coord_Normal;
  su2double diff, Delta, grad_diff;
  su2double U_Tau, U_Tau_Rei, U_Tau_Spa;
  su2double F_Rei, F_Log, Phi_Rei, Phi_Kno;
  su2double U_Plus, Y_Plus;
  su2double TauElem[3], TauNormal, TauTangent[3], WallShearStress;
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  const unsigned short max_iter = 1000;
  const su2double tol = 1e-12;
  bool converged = true;

  const unsigned short turb_model = config->GetKind_Turb_Model();
  const bool tkeNeeded = (turb_model == SST) || (turb_model == SST_SUST);
  
  /*--- Compute the recovery factor ---*/
  // Double-check: laminar or turbulent Pr for this?
  const su2double Recovery = pow(config->GetPrandtl_Lam(), (1.0/3.0));

  /*--- Typical constants from boundary layer theory ---*/

  const su2double kappa = 0.41;
  const su2double B = 5.2;

  const su2double eps = numeric_limits<passivedouble>::epsilon();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {

      /*--- Identify the boundary by string name ---*/

      string Marker_Tag = config->GetMarker_All_TagBound(iMarker);

      /*--- Get the specified wall heat flux from config ---*/

      // Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

      /*--- Loop over all of the vertices on this boundary marker ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        /*--- Check if the node belongs to the domain (i.e, not a halo node)
         and the neighbor is not part of the physical boundary ---*/

        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Get coordinates of the current vertex and nearest normal point ---*/

          Coord = geometry->node[iPoint]->GetCoord();
          Coord_Normal = geometry->node[Point_Normal]->GetCoord();

          /*--- Compute dual-grid area and boundary normal ---*/

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt (Area);

          for (iDim = 0; iDim < nDim; iDim++)
            UnitNormal[iDim] = -Normal[iDim]/Area;
          
          if (geometry->vertex[iMarker][iVertex]->GetDonorFound()){
          
            /*--- Get the distance to the exchange location ---*/
            
            const su2double *doubleInfo = config->GetWallFunction_DoubleInfo(Marker_Tag);
            WallDistMod = doubleInfo[0];

            /*--- Get the density, velocity, and energy at the exchange location ---*/
            
            su2double Density_Normal = 0., Energy_Normal = 0., Kine_Normal = 0., VelMod = 0.;
            P_Normal = 0.;
            T_Normal = 0.;
            for (iDim = 0; iDim < nDim; iDim++) Vel[iDim] = 0.;
            
            const unsigned short nDonors = geometry->vertex[iMarker][iVertex]->GetnDonorPoints();
            for (auto iNode = 0; iNode < nDonors; iNode++) {
              const unsigned long donorPoint = geometry->vertex[iMarker][iVertex]->GetInterpDonorPoint(iNode);
              const su2double donorCoeff     = geometry->vertex[iMarker][iVertex]->GetDonorCoeff(iNode);
              
              Density_Normal  += donorCoeff*nodes->GetDensity(donorPoint);
              Energy_Normal   += donorCoeff*nodes->GetEnergy(donorPoint);
              if (tkeNeeded && solver[TURB_SOL] != nullptr)
                Kine_Normal   += donorCoeff*solver[TURB_SOL]->GetNodes()->GetPrimitive(donorPoint,0);
              
              for (iDim = 0; iDim < nDim; iDim++) Vel[iDim] += donorCoeff*nodes->GetVelocity(donorPoint,iDim);
            }

            for (iDim = 0; iDim < nDim; iDim++) VelMod += Vel[iDim]*Vel[iDim];
            Energy_Normal -= Kine_Normal + 0.5*VelMod;
            
            GetFluidModel()->SetTDState_rhoe(Density_Normal, Energy_Normal);
            P_Normal = GetFluidModel()->GetPressure();
            T_Normal = GetFluidModel()->GetTemperature();
            Lam_Visc_Normal = GetFluidModel()->GetLaminarViscosity();

            /*--- Compute the wall-parallel velocity at the exchange location ---*/

            VelNormal = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelNormal += Vel[iDim] * UnitNormal[iDim];
            for (iDim = 0; iDim < nDim; iDim++)
              VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

            VelTangMod = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelTangMod += VelTang[iDim]*VelTang[iDim];
            VelTangMod = sqrt(VelTangMod);
          }
          else {
            /*--- Get the velocity, pressure, and temperature at the nearest
             (normal) interior point. ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              Vel[iDim] = nodes->GetVelocity(Point_Normal,iDim);
            P_Normal = nodes->GetPressure(Point_Normal);
            T_Normal = nodes->GetTemperature(Point_Normal);
            Lam_Visc_Normal = nodes->GetLaminarViscosity(Point_Normal);

            /*--- Compute the wall-parallel velocity at first point off the wall ---*/

            VelNormal = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelNormal += Vel[iDim] * UnitNormal[iDim];
            for (iDim = 0; iDim < nDim; iDim++)
              VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

            VelTangMod = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              VelTangMod += VelTang[iDim]*VelTang[iDim];
            VelTangMod = sqrt(VelTangMod);

            /*--- Compute normal distance of the interior point from the wall ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);

            WallDistMod = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              WallDistMod += WallDist[iDim]*WallDist[iDim];
            WallDistMod = sqrt(WallDistMod);
          }

          /*--- Compute the wall temperature using the Crocco-Buseman equation ---*/

          T_Wall = T_Normal + Recovery*pow(VelTangMod,2.0)/(2.0*Cp);

          /*--- Extrapolate the pressure from the interior & compute the
           wall density using the equation of state ---*/

          P_Wall = P_Normal;
          Density_Wall = P_Wall/(Gas_Constant*T_Wall);

          /*--- Compute the shear stress at the wall in the regular fashion
           by using the stress tensor on the surface ---*/
          
          const su2double StaticEnergy_Wall = T_Wall*Gas_Constant/Gamma_Minus_One;
          GetFluidModel()->SetTDState_rhoe(Density_Wall, StaticEnergy_Wall);
          Lam_Visc_Wall = GetFluidModel()->GetLaminarViscosity();

          grad_primvar  = nodes->GetGradient_Primitive(iPoint);

          div_vel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            div_vel += grad_primvar[iDim+1][iDim];

          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0 ; jDim < nDim; jDim++) {
              Delta = 0.0; if (iDim == jDim) Delta = 1.0;
              tau[iDim][jDim] = Lam_Visc_Wall*(  grad_primvar[jDim+1][iDim]
                                               + grad_primvar[iDim+1][jDim]) -
              TWO3*Lam_Visc_Wall*div_vel*Delta;
            }
            TauElem[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++)
              TauElem[iDim] += tau[iDim][jDim]*UnitNormal[jDim];
          }

          /*--- Compute wall shear stress as the magnitude of the wall-tangential
           component of the shear stress tensor---*/

          TauNormal = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            TauNormal += TauElem[iDim] * UnitNormal[iDim];

          for (iDim = 0; iDim < nDim; iDim++)
            TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];

          WallShearStress = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            WallShearStress += TauTangent[iDim]*TauTangent[iDim];
          WallShearStress = sqrt(WallShearStress);

          /*--- Calculate the quantities from boundary layer theory and
           iteratively solve for a new wall shear stress. Use the current wall
           shear stress as a starting guess for the wall function. ---*/

          counter = 0; diff = 1.0;
          U_Tau_Rei = VelTangMod / WallDistMod;
          Y_Plus = Density_Wall * U_Tau_Rei * WallDistMod / Lam_Visc_Wall;
          converged = true;

          /*--- First solve for Riechardt U_Tau ---*/

          if (Y_Plus < 500.) {

            while (fabs(diff) > tol) {

              /*--- y+ ---*/

              Y_Plus = Density_Wall*WallDistMod*U_Tau_Rei/Lam_Visc_Wall;

              /*--- Spalding's universal form for the BL velocity. ---*/

              F_Rei = log(1. + 0.4 * Y_Plus) / kappa + 
                      7.8 * (1. - exp(-Y_Plus / 11.) - 
                      (Y_Plus / 11.) * exp(-Y_Plus / 3.));

              /*--- Blending with classic log law. ---*/

              F_Log = log(Y_Plus) / kappa + 5.1;
              Phi_Rei = tanh(pow(Y_Plus / 27., 4.));
              U_Plus = (1. - Phi_Rei) * F_Rei + Phi_Rei * F_Log;

              /*--- Define function for Newton method to zero ---*/

              diff = U_Plus - (VelTangMod / U_Tau_Rei);

              /*--- Gradient of function defined above ---*/

              grad_diff = VelTangMod /(U_Tau_Rei * U_Tau_Rei) + (1. - Phi_Rei) *
                          0.4 * Density_Wall * WallDistMod / (kappa * Lam_Visc_Wall * (1. + 0.4 * Y_Plus)) +
                          7.8 * (Y_Plus * Y_Plus / (33. * U_Tau_Rei) * exp(-Y_Plus / 3.) - 
                          Y_Plus/(11. * U_Tau_Rei) * exp(-Y_Plus / 11.) - 
                          Y_Plus/(11. * U_Tau_Rei) * exp(-Y_Plus / 3.)) + Phi_Rei *
                          (1. / (kappa * U_Tau_Rei)) + (F_Log - F_Rei) *
                          (4. * pow(Y_Plus / 27., 4.) / U_Tau) * 
                          pow(1. / (cosh(pow(Y_Plus / 27., 4.))), 2.);

              if (grad_diff != grad_diff) {
                // cout << "rank: " << rank << ". node: " << geometry->node[iPoint]->GetGlobalIndex() << ". Y_Plus: " << Y_Plus 
                //      << ". pow(...): " << pow(1. / (cosh(pow(Y_Plus / 27., 4.))), 2.) << endl; 
                converged = false;
                break;
              }

              /* --- Newton Step --- */

              const su2double sign = 1.0 - 2.0*(grad_diff <  0);
              U_Tau_Rei = U_Tau_Rei - diff / (grad_diff + sign*eps);

              counter++;
              if (counter == max_iter) {
                converged = false;
                break;
              }

            }
          }

          /*--- If not converged or y+ too large, jump to the next point. --- */

          if (!converged || Y_Plus > 500.) {
            nodes->SetTauWall(iPoint,-1.0);
            nodes->SetTauWallFactor(iPoint,1.0);
            continue;
          }

          /*--- Now solve for Spalding U_Tau ---*/

          U_Tau_Spa = U_Tau_Rei;

          if (Y_Plus < 300.) {

            counter = 0; diff = 1.0;

            while (fabs(diff) > tol) {

              /*--- u+ ---*/

              U_Plus = VelTangMod/U_Tau_Spa;

              /*--- Spalding's universal form for the BL velocity. ---*/

              Y_Plus = U_Plus + (exp(-kappa * B) *
                       (exp(kappa * U_Plus) - 1. - kappa*U_Plus -
                       kappa*kappa*U_Plus*U_Plus/2.0 -
                       kappa*kappa*kappa*U_Plus*U_Plus*U_Plus/6.0));

              /* --- Define function for Newton method to zero --- */

              diff = Y_Plus - (Density_Wall * U_Tau_Spa * WallDistMod / Lam_Visc_Wall);

              /* --- Gradient of function defined above --- */
              grad_diff = -VelTangMod / (U_Tau_Spa * U_Tau_Spa) + exp(-kappa * B) * 
                          (-kappa * U_Plus / U_Tau_Spa * exp(kappa * U_Plus) + 
                          kappa * U_Plus / U_Tau_Spa +
                          kappa * kappa * U_Plus * U_Plus / U_Tau_Spa +
                          kappa * kappa * kappa * U_Plus * U_Plus * U_Plus / U_Tau_Spa ) - 
                          Density_Wall * WallDistMod / Lam_Visc_Wall;

              if (grad_diff != grad_diff) {
                converged = false;
                break;
              }

              /* --- Newton Step --- */

              const su2double sign = 1.0 - 2.0*(grad_diff <  0);
              U_Tau_Spa = U_Tau_Spa - diff / (grad_diff + sign*eps);

              counter++;
              if (counter == max_iter) {
                converged = false;
                break;
              }

            }

            /*--- If not converged, jump to the next point. --- */

            if (!converged) {
              nodes->SetTauWall(iPoint,-1.0);
              nodes->SetTauWallFactor(iPoint,1.0);
              continue;
            }

            Phi_Kno = tanh(pow(Y_Plus / 50., 2.));
          }

          else Phi_Kno = 1.0;

          U_Tau = (1. - Phi_Kno) * U_Tau_Spa + Phi_Kno * U_Tau_Rei;
          Y_Plus = Density_Wall*WallDistMod*U_Tau/Lam_Visc_Wall;

          /*--- Calculate an updated value for the wall shear stress
            using the y+ value, the definition of y+, and the definition of
            the friction velocity. ---*/

          Tau_Wall = (1.0/Density_Wall)*pow(Y_Plus*Lam_Visc_Wall/WallDistMod,2.0);

          /*--- Store this value for the wall shear stress at the node.  ---*/

          nodes->SetTauWall(iPoint,Tau_Wall);
          nodes->SetTauWallFactor(iPoint,Tau_Wall/WallShearStress);
          nodes->SetTemperature(iPoint,T_Wall);
          nodes->SetSolution(iPoint, 0, Density_Wall);
          nodes->SetPrimitive(iPoint, nDim+2, Density_Wall);
          nodes->SetPrimitive(iPoint, nDim+1, P_Wall);
          nodes->SetLaminarViscosity(iPoint, Lam_Visc_Wall);

        }

      }

    }
  }
  
  InitiateComms(geometry, config, PRIMITIVE);
  CompleteComms(geometry, config, PRIMITIVE);
  
  /*--- Communicate values needed for WF ---*/
  InitiateComms(geometry, config, WALL_FUNCTION);
  CompleteComms(geometry, config, WALL_FUNCTION);
  WallFunctionComms(geometry, solver, config);

  /*--- Correct values in turbulent solver ---*/
  if (solver[TURB_SOL] != nullptr)
    solver[TURB_SOL]->ComputeKnoppWallFunction(geometry, solver, config);

}
