/*!
 * \file CNEMOEulerSolver.cpp
 * \brief Headers of the CNEMOEulerSolver class
 * \author S. R. Copeland, F. Palacios, W. Maier, C. Garbacz
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

#include "../../include/solvers/CNEMOEulerSolver.hpp"
#include "../../include/variables/CNEMONSVariable.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CMutationTCLib.hpp"
#include "../../include/fluid/CSU2TCLib.hpp"
#include "../../include/limiters/CLimiterDetails.hpp"

CNEMOEulerSolver::CNEMOEulerSolver(CGeometry *geometry, CConfig *config,
                           unsigned short iMesh, const bool navier_stokes) :
  CFVMFlowSolverBase<CNEMOEulerVariable, ENUM_REGIME::COMPRESSIBLE>(*geometry, *config) {

  /*--- Based on the navier_stokes boolean, determine if this constructor is
   *    being called by itself, or by its derived class CNEMONSSolver. ---*/
  string description;
  if (navier_stokes) {
    description = "Navier-Stokes";
  }
  else {
    description = "Euler";
  }

  const auto nZone = geometry->GetnZone();
  const bool restart = (config->GetRestart() || config->GetRestart_Flow());
  const auto direct_diff = config->GetDirectDiff();
  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool time_stepping = (config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING);
  const bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();

  int Unst_RestartIter = 0;
  unsigned long iMarker;
  unsigned short nLineLets;
  su2double *Mvec_Inf, Alpha, Beta;

  /*--- A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain ---*/
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the AoA
  before non-dimensionalizing ---*/
  if (!(!restart || (iMesh != MESH_0) || nZone >1 )) {

    /*--- Modify file name for a dual-time unsteady restart ---*/
    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/
    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
    }

    string filename_ = "flow";
    filename_ = config->GetFilename(filename_, ".meta", Unst_RestartIter);

    /*--- Read and store the restart metadata ---*/
    Read_SU2_Restart_Metadata(geometry, config, false, filename_);

  }

  /*--- Define geometric constants in the solver structure ---*/
  nSpecies     = config->GetnSpecies();
  nMarker      = config->GetnMarker_All();
  nDim         = geometry->GetnDim();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Set size of the conserved and primitive vectors ---*/
  //     U: [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  //     V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradV: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // Viscous: append [mu, mu_t]^T
  nVar         = nSpecies + nDim + 2;
  if (navier_stokes) { nPrimVar   = nSpecies + nDim + 10; }
  else {               nPrimVar   = nSpecies +nDim +8;    }
  nPrimVarGrad = nSpecies + nDim + 8;

  /*--- Initialize nVarGrad for deallocation ---*/
  nVarGrad     = nPrimVarGrad;

  /*--- Store the number of vertices on each marker for deallocation ---*/
  nVertex.resize(nMarker);
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  MassFrac_Inf = config->GetGas_Composition();

  /*--- Perform the non-dimensionalization for the flow equations using the
    specified reference values. ---*/
  SetNondimensionalization(config, iMesh);

  /// TODO: This type of variables will be replaced.

  AllocateTerribleLegacyTemporaryVariables();

  /*--- Allocate base class members. ---*/

  Allocate(*config);

  /*--- Allocate Jacobians for implicit time-stepping ---*/
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    /*--- Jacobians and vector  structures for implicit computations ---*/
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
  }
  else {
    if (rank == MASTER_NODE)  cout<< "Explicit Scheme. No Jacobian structure (" << description << "). MG level: " << iMesh <<"."<<endl;
  }

  /*--- Read farfield conditions from the config file ---*/
  Mach_Inf            = config->GetMach();
  Density_Inf         = config->GetDensity_FreeStreamND();
  Pressure_Inf        = config->GetPressure_FreeStreamND();
  Velocity_Inf        = config->GetVelocity_FreeStreamND();
  Temperature_Inf     = config->GetTemperature_FreeStreamND();
  Temperature_ve_Inf  = config->GetTemperature_ve_FreeStreamND();

  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  switch(direct_diff) {
  case NO_DERIVATIVE:
    /*--- Default ---*/
    break;
  case D_DENSITY:
    SU2_TYPE::SetDerivative(Density_Inf, 1.0);
    break;
  case D_PRESSURE:
    SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
    break;
  case D_TEMPERATURE:
    SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
    break;
  case D_MACH: case D_AOA:
  case D_SIDESLIP: case D_REYNOLDS:
  case D_TURB2LAM: case D_DESIGN:
    /*--- Already done in postprocessing of config ---*/
    break;
  default:
    break;
  }

  SetReferenceValues(*config);

  /*--- Vectorize free stream Mach number based on AoA & AoS ---*/
  Mvec_Inf = new su2double[nDim];
  Alpha    = config->GetAoA()*PI_NUMBER/180.0;
  Beta     = config->GetAoS()*PI_NUMBER/180.0;
  if (nDim == 2) {
    Mvec_Inf[0] = cos(Alpha)*Mach_Inf;
    Mvec_Inf[1] = sin(Alpha)*Mach_Inf;
  }
  if (nDim == 3) {
    Mvec_Inf[0] = cos(Alpha)*cos(Beta)*Mach_Inf;
    Mvec_Inf[1] = sin(Beta)*Mach_Inf;
    Mvec_Inf[2] = sin(Alpha)*cos(Beta)*Mach_Inf;
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  if (navier_stokes) {
    nodes      = new CNEMONSVariable    (Pressure_Inf, MassFrac_Inf, Mvec_Inf,
                                         Temperature_Inf, Temperature_ve_Inf,
                                         nPoint, nDim, nVar, nPrimVar, nPrimVarGrad,
                                         config, FluidModel);
    node_infty = new CNEMONSVariable    (Pressure_Inf, MassFrac_Inf, Mvec_Inf,
                                        Temperature_Inf, Temperature_ve_Inf,
                                        1, nDim, nVar, nPrimVar, nPrimVarGrad,
                                        config, FluidModel);
  } else {
    nodes      = new CNEMOEulerVariable(Pressure_Inf, MassFrac_Inf, Mvec_Inf,
                                        Temperature_Inf, Temperature_ve_Inf,
                                        nPoint, nDim, nVar, nPrimVar, nPrimVarGrad,
                                        config, FluidModel);
    node_infty = new CNEMOEulerVariable(Pressure_Inf, MassFrac_Inf, Mvec_Inf,
                                        Temperature_Inf, Temperature_ve_Inf,
                                        1, nDim, nVar, nPrimVar, nPrimVarGrad,
                                        config, FluidModel);
  }
  SetBaseClassPointerToNodes();

  node_infty->SetPrimVar(0, FluidModel);

  /*--- Initial comms. ---*/

  CommunicateInitialState(geometry, config);

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "NEMO.FLOW";

  /*--- Finally, check that the static arrays will be large enough (keep this
   *    check at the bottom to make sure we consider the "final" values). ---*/
  if((nDim > MAXNDIM) || (nPrimVar > MAXNVAR))
    SU2_MPI::Error("Oops! The CNEMOEulerSolver static array sizes are not large enough.",CURRENT_FUNCTION);

   /*--- Deallocate arrays ---*/
  delete [] Mvec_Inf;

}

CNEMOEulerSolver::~CNEMOEulerSolver(void) {

  delete node_infty;
  delete FluidModel;

}

void CNEMOEulerSolver::CommonPreprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                           unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  bool center_jst       = (config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0);
  bool center_jst_ke    = (config->GetKind_Centered_Flow() == JST_KE) && (iMesh == MESH_0);

  /*--- Set the primitive variables ---*/
  ompMasterAssignBarrier(ErrorCounter,0);
  SU2_OMP_ATOMIC
  ErrorCounter += SetPrimitive_Variables(solver_container, config, Output);
  SU2_OMP_BARRIER

  SU2_OMP_MASTER { /*--- Ops that are not OpenMP parallel go in this block. ---*/

    if ((iMesh == MESH_0) && (config->GetComm_Level() == COMM_FULL)) {
      unsigned long tmp = ErrorCounter;
      SU2_MPI::Allreduce(&tmp, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      config->SetNonphysical_Points(ErrorCounter);

      if ((rank == MASTER_NODE) && (ErrorCounter != 0))
        cout << "Warning. The initial solution contains "<< ErrorCounter << " points that are not physical." << endl;
    }
  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER

  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if (center_jst) SetUndivided_Laplacian(geometry, config);
    if (center_jst || center_jst_ke) SetCentered_Dissipation_Sensor(geometry, config);
  }

  /*--- Initialize the Jacobian matrix and residual, not needed for the reducer strategy
   *    as we set blocks (including diagonal ones) and completely overwrite. ---*/

  if(!ReducerStrategy && !Output) {
    LinSysRes.SetValZero();
    if (implicit) Jacobian.SetValZero();
  }
}

void CNEMOEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep,
                                     unsigned short RunTime_EqSystem, bool Output) {
  const unsigned long InnerIter = config->GetInnerIter();
  const bool muscl       = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool limiter     = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  const bool center      = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;
  const bool van_albada  = config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE;

  /*--- Common preprocessing steps ---*/
  CommonPreprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  /*--- Upwind second order reconstruction ---*/
  if (muscl && !center && !Output) {

    /*--- Calculate the gradients ---*/
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config, true);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config, true);
    }

    /*--- Limiter computation ---*/
    if (limiter && !van_albada) {
      SetPrimitive_Limiter(geometry, config);
    }
  }
}

unsigned long CNEMOEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  bool nonphysical = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Incompressible flow, primitive variables ---*/

    nonphysical = nodes->SetPrimVar(iPoint,FluidModel);

    /* Check for non-realizable states for reporting. */

    if (nonphysical) nonPhysicalPoints++;

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return nonPhysicalPoints;
}

void CNEMOEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /*--- Define an object to compute the speed of sound. ---*/
  struct SoundSpeed {
    FORCEINLINE su2double operator() (const CNEMOEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      return 0.5 * (nodes.GetSoundSpeed(iPoint) + nodes.GetSoundSpeed(jPoint));
    }

    FORCEINLINE su2double operator() (const CNEMOEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetSoundSpeed(iPoint);
    }

  } soundSpeed;

  /*--- Define an object to compute the viscous eigenvalue. ---*/
  struct LambdaVisc {
    FORCEINLINE su2double lambda(su2double lamVisc, su2double eddyVisc, su2double rho, su2double k, su2double cv) const {
      /*--- Determine the viscous spectral radius and apply it to the control volume ---*/
      su2double Lambda_1 = (4.0/3.0)*(lamVisc + eddyVisc);
      return (Lambda_1 + k/cv)/rho;
    }

    FORCEINLINE su2double operator() (const CNEMOEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      su2double lamVisc = 0.5*(nodes.GetLaminarViscosity(iPoint) + nodes.GetLaminarViscosity(jPoint));
      su2double eddyVisc = 0.5*(nodes.GetEddyViscosity(iPoint) + nodes.GetEddyViscosity(jPoint));
      su2double thermalCond = 0.5*(nodes.GetThermalConductivity(iPoint) + nodes.GetThermalConductivity(jPoint) +
                                   nodes.GetThermalConductivity_ve(iPoint) + nodes.GetThermalConductivity_ve(jPoint));
      su2double density = 0.5*(nodes.GetDensity(iPoint) + nodes.GetDensity(jPoint));
      su2double cv = 0.5*(nodes.GetRhoCv_tr(iPoint) + nodes.GetRhoCv_ve(iPoint) +
                          nodes.GetRhoCv_tr(jPoint) + nodes.GetRhoCv_ve(jPoint))/ density;
      return lambda(lamVisc, eddyVisc, density, thermalCond, cv);
    }

    FORCEINLINE su2double operator() (const CNEMOEulerVariable& nodes, unsigned long iPoint) const {
      su2double lamVisc = nodes.GetLaminarViscosity(iPoint);
      su2double eddyVisc = nodes.GetEddyViscosity(iPoint);
      su2double thermalCond = nodes.GetThermalConductivity(iPoint) + nodes.GetThermalConductivity_ve(iPoint);
      su2double density = nodes.GetDensity(iPoint);
      su2double cv = (nodes.GetRhoCv_tr(iPoint) + nodes.GetRhoCv_ve(iPoint))/ density;
      return lambda(lamVisc, eddyVisc, density, thermalCond, cv);
    }

  } lambdaVisc;

  /*--- Now instantiate the generic implementation with the two functors above. ---*/

  SetTime_Step_impl(soundSpeed, lambdaVisc, geometry, solver_container, config, iMesh, Iteration);

}

void CNEMOEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {

  /*--- Define an object to compute the speed of sound. ---*/
  struct SoundSpeed {
    FORCEINLINE su2double operator() (const CNEMOEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      return 0.5 * (nodes.GetSoundSpeed(iPoint) + nodes.GetSoundSpeed(jPoint));
    }

    FORCEINLINE su2double operator() (const CNEMOEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetSoundSpeed(iPoint);
    }

  } soundSpeed;

  /*--- Instantiate generic implementation. ---*/

  SetMax_Eigenvalue_impl(soundSpeed, geometry, config);

}

void CNEMOEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;

  CNumerics* numerics = numerics_container[CONV_TERM];

  /*--- Set booleans based on config settings ---*/
  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge, set normal vectors, and number of neighbors ---*/
    iPoint = geometry->edges->GetNode(iEdge, 0);
    jPoint = geometry->edges->GetNode(iEdge, 1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));
    numerics->SetNeighbor(geometry->nodes->GetnNeighbor(iPoint),
                          geometry->nodes->GetnNeighbor(jPoint));

    /*--- Pass conservative & primitive variables w/o reconstruction to CNumerics ---*/
    numerics->SetConservative(nodes->GetSolution(iPoint),  nodes->GetSolution(jPoint));
    numerics->SetPrimitive   (nodes->GetPrimitive(iPoint), nodes->GetPrimitive(jPoint));

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU  ( nodes->GetdPdU(iPoint),   nodes->GetdPdU(jPoint));
    numerics->SetdTdU  ( nodes->GetdTdU(iPoint),   nodes->GetdTdU(jPoint));
    numerics->SetdTvedU( nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));
    numerics->SetEve   ( nodes->GetEve(iPoint),    nodes->GetEve(jPoint));
    numerics->SetCvve  ( nodes->GetCvve(iPoint),   nodes->GetCvve(jPoint));
    numerics->SetGamma ( nodes->GetGamma(iPoint),  nodes->GetGamma(jPoint));

    /*--- Set the largest convective eigenvalue ---*/
    numerics->SetLambda(nodes->GetLambda(iPoint), nodes->GetLambda(jPoint));

    /*--- Compute residuals, and Jacobians ---*/
    auto residual = numerics->ComputeResidual(config);

    /*--- Check residuals/Jacobians --*/
    bool err = CNumerics::CheckResidualNaNs(implicit, nVar, residual);

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.AddBlock(iPoint, residual);
      LinSysRes.SubtractBlock(jPoint, residual);
      if (implicit) {
        Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
      }
    }
  }
}

void CNEMOEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                       CConfig *config, unsigned short iMesh) {

  /*--- Set booleans based on config settings ---*/
  const bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl            = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  const bool limiter          = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE);
  const bool van_albada       = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);

  /*--- Non-physical counter. ---*/
  unsigned long counter_local = 0;
  SU2_OMP_MASTER
  ErrorCounter = 0;
  END_SU2_OMP_MASTER

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[CONV_TERM];

  /*--- Static arrays for MUSCL reconstructed variables ---*/
  su2double     Primitive_i[MAXNVAR] = {0.0},    Primitive_j[MAXNVAR] = {0.0};
  su2double     Conserved_i[MAXNVAR] = {0.0},    Conserved_j[MAXNVAR] = {0.0};
  su2double          dPdU_i[MAXNVAR] = {0.0},         dPdU_j[MAXNVAR] = {0.0};
  su2double          dTdU_i[MAXNVAR] = {0.0},         dTdU_j[MAXNVAR] = {0.0};
  su2double        dTvedU_i[MAXNVAR] = {0.0},       dTvedU_j[MAXNVAR] = {0.0};
  su2double           Eve_i[MAXNVAR] = {0.0},          Eve_j[MAXNVAR] = {0.0};
  su2double          Cvve_i[MAXNVAR] = {0.0},         Cvve_j[MAXNVAR] = {0.0};
  su2double  Project_Grad_i[MAXNVAR] = {0.0}, Project_Grad_j[MAXNVAR] = {0.0};
  su2double Gamma_i = 0.0, Gamma_j = 0.0;

  /*--- Loop over edges and calculate convective fluxes ---*/
  for(unsigned long iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    unsigned short iDim, iVar;

    /*--- Retrieve node numbers and pass edge normal to CNumerics ---*/
    auto iPoint = geometry->edges->GetNode(iEdge, 0);
    auto jPoint = geometry->edges->GetNode(iEdge, 1);

    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    auto Coord_i = geometry->nodes->GetCoord(iPoint);
    auto Coord_j = geometry->nodes->GetCoord(jPoint);

    /*--- Get conserved & primitive variables from CVariable ---*/
    auto U_i = nodes->GetSolution(iPoint);   auto U_j = nodes->GetSolution(jPoint);
    auto V_i = nodes->GetPrimitive(iPoint);  auto V_j = nodes->GetPrimitive(jPoint);

    /*--- Set them with or without high order reconstruction using MUSCL strategy. ---*/
    if (!muscl) {

      numerics->SetPrimitive   (V_i, V_j);
      numerics->SetConservative(U_i, U_j);
      numerics->SetdPdU  (nodes->GetdPdU(iPoint),   nodes->GetdPdU(jPoint));
      numerics->SetdTdU  (nodes->GetdTdU(iPoint),   nodes->GetdTdU(jPoint));
      numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));
      numerics->SetEve   (nodes->GetEve(iPoint),    nodes->GetEve(jPoint));
      numerics->SetCvve  (nodes->GetCvve(iPoint),   nodes->GetCvve(jPoint));
      numerics->SetGamma (nodes->GetGamma(iPoint),  nodes->GetGamma(jPoint));

    } else {

      /*--- High order reconstruction using MUSCL strategy ---*/
      su2double Vector_ij[MAXNDIM] = {0.0};
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_ij[iDim] = 0.5*(Coord_j[iDim] - Coord_i[iDim]);
      }

      /*--- Retrieve gradient information ---*/
      auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
      auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

      /*--- Set and extract limiters ---*/
      su2double *Limiter_i = nullptr, *Limiter_j = nullptr;

      if (limiter && !van_albada){
        Limiter_i = nodes->GetLimiter_Primitive(iPoint);
        Limiter_j = nodes->GetLimiter_Primitive(jPoint);
      }

      su2double lim_i = 2.0;
      su2double lim_j = 2.0;

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Project_Grad_i[iVar] = 0.0; Project_Grad_j[iVar] = 0.0;

        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i[iVar] += Vector_ij[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j[iVar] -= Vector_ij[iDim]*Gradient_j[iVar][iDim];
        }

        if (limiter) {
          if (van_albada) {
            su2double V_ij = V_j[iVar] - V_i[iVar];
            su2double va_lim_i = LimiterHelpers<>::vanAlbadaFunction(Project_Grad_i[iVar], V_ij, EPS);
            su2double va_lim_j = LimiterHelpers<>::vanAlbadaFunction(-Project_Grad_j[iVar], V_ij, EPS);
            lim_i = min(lim_i, va_lim_i);
            lim_j = min(lim_j, va_lim_j);
          } else {
            lim_i = min(lim_i, Limiter_i[iVar]);
            lim_j = min(lim_j, Limiter_j[iVar]);
          }
        } else {
          lim_i = lim_j = 1.0;
        }
      }
      su2double lim_ij = min(lim_i, lim_j);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Primitive_i[iVar] = V_i[iVar] + lim_ij*Project_Grad_i[iVar];
        Primitive_j[iVar] = V_j[iVar] + lim_ij*Project_Grad_j[iVar];
      }

      /*--- Check for non-physical solutions after reconstruction. If found, use the
       cell-average value of the solution. This is a locally 1st order approximation,
       which is typically only active during the start-up of a calculation. ---*/
      bool chk_err_i = CheckNonPhys(Primitive_i);
      bool chk_err_j = CheckNonPhys(Primitive_j);

      nodes->SetNon_Physical(iPoint, chk_err_i);
      nodes->SetNon_Physical(jPoint, chk_err_j);

      /*--- Get updated state, in case the point recovered after the set. ---*/
      chk_err_i = nodes->GetNon_Physical(iPoint);
      chk_err_j = nodes->GetNon_Physical(jPoint);

      counter_local += chk_err_i + chk_err_j;

      /*--- Compute Secondary variables in a thermaodynamically consistent way. ---*/
      if (!chk_err_i) Gamma_i = ComputeConsistentExtrapolation(GetFluidModel(), nSpecies, Primitive_i, dPdU_i, dTdU_i, dTvedU_i, Eve_i, Cvve_i);
      if (!chk_err_j) Gamma_j = ComputeConsistentExtrapolation(GetFluidModel(), nSpecies, Primitive_j, dPdU_j, dTdU_j, dTvedU_j, Eve_j, Cvve_j);

      /*--- Recompute Conserved variables if Roe or MSW scheme ---*/
      if ((config->GetKind_Upwind_Flow() == ROE) || (config->GetKind_Upwind_Flow() == MSW)){
        if (!chk_err_i) RecomputeConservativeVector(Conserved_i, Primitive_i);
        if (!chk_err_j) RecomputeConservativeVector(Conserved_j, Primitive_j);
      }

      /*--- If non-physical, revert to first order ---*/
      numerics->SetConservative(chk_err_i ? U_i : Conserved_i, chk_err_j ? U_j : Conserved_j);
      numerics->SetPrimitive   (chk_err_i ? V_i : Primitive_i, chk_err_j ? V_j : Primitive_j);
      numerics->SetdPdU  (chk_err_i ? nodes->GetdPdU  (iPoint) : dPdU_i,    chk_err_j ? nodes->GetdPdU  (jPoint) : dPdU_j);
      numerics->SetdTdU  (chk_err_i ? nodes->GetdTdU  (iPoint) : dTdU_i,    chk_err_j ? nodes->GetdTdU  (jPoint) : dTdU_j);
      numerics->SetdTvedU(chk_err_i ? nodes->GetdTvedU(iPoint) : dTvedU_i,  chk_err_j ? nodes->GetdTvedU(jPoint) : dTvedU_j);
      numerics->SetEve   (chk_err_i ? nodes->GetEve   (iPoint) : Eve_i,     chk_err_j ? nodes->GetEve   (jPoint) : Eve_j);
      numerics->SetCvve  (chk_err_i ? nodes->GetCvve  (iPoint) : Cvve_i,    chk_err_j ? nodes->GetCvve  (jPoint) : Cvve_j);
      numerics->SetGamma (chk_err_i ? nodes->GetGamma (iPoint) : Gamma_i,   chk_err_j ? nodes->GetGamma (jPoint) : Gamma_j);

    }

    /*--- Compute the residual ---*/
    auto residual = numerics->ComputeResidual(config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err = CNumerics::CheckResidualNaNs(implicit, nVar, residual);

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.AddBlock(iPoint, residual);
      LinSysRes.SubtractBlock(jPoint, residual);
      if (implicit) {
        Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
      }
    }
  }

  /*--- Warning message about non-physical reconstructions. ---*/
  if ((iMesh == MESH_0) && (config->GetComm_Level() == COMM_FULL)) {
    /*--- Add counter results for all threads. ---*/
    SU2_OMP_ATOMIC
    ErrorCounter += counter_local;
    SU2_OMP_BARRIER

    /*--- Add counter results for all ranks. ---*/
    SU2_OMP_MASTER
    {
      counter_local = ErrorCounter;
      SU2_MPI::Reduce(&counter_local, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
      config->SetNonphysical_Reconstr(ErrorCounter);
    }
    END_SU2_OMP_MASTER
    SU2_OMP_BARRIER
  }
}

su2double CNEMOEulerSolver::ComputeConsistentExtrapolation(CNEMOGas *fluidmodel, unsigned short nSpecies, su2double *V,
                                                           su2double* dPdU, su2double* dTdU, su2double* dTvedU,
                                                           su2double* val_eves, su2double *val_Cvves) {

  //NOTE: TODO - this doesnt compute Cvves/ dPdU,etc.yet
  su2double val_gamma;
  vector<su2double> rhos;

  /*--- Rename index information ---*/
  unsigned short T_INDEX   = nSpecies;
  unsigned short TVE_INDEX = nSpecies+1;

  /*--- Rename density vector ---*/
  rhos.resize(nSpecies,0.0);
  for (unsigned short iSpecies=0; iSpecies < nSpecies; iSpecies++ ){
    rhos[iSpecies] = V[iSpecies];
  }

  /*--- Set new fluid state ---*/
  fluidmodel->SetTDStateRhosTTv(rhos, V[T_INDEX], V[TVE_INDEX]);

  /*---Compute the secondary values ---*/
  // This block of code copies a vector to corresponding pointer.
  auto it = val_eves;
  auto& ref = fluidmodel->ComputeSpeciesEve(V[TVE_INDEX]);
  for (auto v : ref) {
    *it = v;  ++it;
  }

  val_gamma = fluidmodel->ComputeGamma();

  return val_gamma;
}

void CNEMOEulerSolver::RecomputeConservativeVector(su2double *U, const su2double *V) const {

  /*---Useful variables ---*/
  vector<su2double> rhos;
  rhos.resize(nSpecies,0.0);

  /*--- Set Indices ---*/
  //Make these in a general location
  unsigned short RHO_INDEX = nodes->GetRhoIndex();
  unsigned short T_INDEX   = nodes->GetTIndex();
  unsigned short TVE_INDEX = nodes->GetTveIndex();
  unsigned short VEL_INDEX = nodes->GetVelIndex();

  /*--- Set densities and mass fraction ---*/
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    U[iSpecies]    = V[iSpecies];
    rhos[iSpecies] = V[iSpecies];
  }

  /*--- Set momentum and compute v^2 ---*/
  //TODO: geometry toolbox
  su2double sqvel = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    U[nSpecies+iDim] = V[RHO_INDEX]*V[VEL_INDEX+iDim];
    sqvel           += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];
  }

  /*--- Set the fluidmodel and recompute energies ---*/
  FluidModel->SetTDStateRhosTTv( rhos, V[T_INDEX], V[TVE_INDEX]);
  const auto& Energies = FluidModel->ComputeMixtureEnergies();

  /*--- Set conservative energies ---*/
  U[nSpecies+nDim]   = V[RHO_INDEX]*(Energies[0]+0.5*sqvel);
  U[nSpecies+nDim+1] = V[RHO_INDEX]*(Energies[1]);

}

bool CNEMOEulerSolver::CheckNonPhys(const su2double *V) const {

  su2double Tmin, Tmax, Tvemin, Tvemax;

  /*--- Set booleans ---*/
  bool nonPhys = false;

  /*--- Set Indices ---*/
  //Make these in a general location
  unsigned short RHOS_INDEX = nodes->GetRhosIndex();
  unsigned short T_INDEX    = nodes->GetTIndex();
  unsigned short TVE_INDEX  = nodes->GetTveIndex();
  unsigned short P_INDEX    = nodes->GetPIndex();
  unsigned short A_INDEX    = nodes->GetAIndex();

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0; Tmax   = 8E4;
  Tvemin = 50.0; Tvemax = 8E4;

  /*--- Check whether state makes sense ---*/
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    if (V[RHOS_INDEX+iSpecies] < 0.0) nonPhys = true;

  if (V[P_INDEX] < 0.0) nonPhys = true;

  if (V[T_INDEX] < Tmin || V[T_INDEX] > Tmax) nonPhys = true;

  if (V[TVE_INDEX] < Tvemin || V[TVE_INDEX] > Tvemax) nonPhys = true;

  if (V[A_INDEX] < 0.0 ) nonPhys = true;

  return nonPhys;

}

void CNEMOEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  unsigned short iVar;
  unsigned long iPoint;

  /*--- Assign booleans ---*/
  bool err        = false;
  bool implicit   = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool frozen     = config->GetFrozen();
  bool monoatomic = config->GetMonoatomic();
  bool axisymm    = config->GetAxisymmetric();
  bool viscous    = config->GetViscous();
  bool rans       = (config->GetKind_Turb_Model() != TURB_MODEL::NONE);

  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];

  /*--- Initialize the error counter ---*/
  unsigned long eAxi_local = 0;
  unsigned long eChm_local = 0;
  unsigned long eVib_local = 0;

  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

  /*--- Preprocess viscous axisymm variables (if necessary) ---*/
  if (axisymm && viscous) {
    ComputeAxisymmetricAuxGradients(geometry,config);
  }

  AD::StartNoSharedReading();

  /*--- loop over interior points ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Set conserved & primitive variables  ---*/
    numerics->SetConservative(nodes->GetSolution(iPoint),  nullptr);
    numerics->SetPrimitive   (nodes->GetPrimitive(iPoint), nullptr);

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(nodes->GetdPdU(iPoint),     nullptr);
    numerics->SetdTdU(nodes->GetdTdU(iPoint),     nullptr);
    numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nullptr);
    numerics->SetEve(nodes->GetEve(iPoint),       nullptr);
    numerics->SetCvve(nodes->GetCvve(iPoint),     nullptr);

    /*--- Set volume of the dual grid cell ---*/
    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint), nullptr);

    /*--- Compute finite rate chemistry ---*/

    if(!monoatomic){
      if(!frozen){
        /*--- Compute the non-equilibrium chemistry ---*/
        auto residual = numerics->ComputeChemistry(config);

        /*--- Check for errors before applying source to the linear system ---*/
        err = CNumerics::CheckResidualNaNs(implicit, nVar, residual);

        /*--- Apply the chemical sources to the linear system ---*/
        if (!err) {
          LinSysRes.SubtractBlock(iPoint, residual);
          if (implicit)
            Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
        } else
          eChm_local++;
      }
    }

    /*--- Compute vibrational energy relaxation ---*/
    /// NOTE: Jacobians don't account for relaxation time derivatives

    if (!monoatomic){
      auto residual = numerics->ComputeVibRelaxation(config);

      /*--- Check for errors before applying source to the linear system ---*/
      err = CNumerics::CheckResidualNaNs(implicit, nVar, residual);

      /*--- Apply the vibrational relaxation terms to the linear system ---*/
      if (!err) {
        LinSysRes.SubtractBlock(iPoint, residual);
        if (implicit)
          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
      } else
        eVib_local++;
    }

    /*--- Compute axisymmetric source terms (if needed) ---*/
    if (axisymm) {

      /*--- If necessary, set variables needed for viscous computation ---*/
      if (viscous) {

        /*--- Set gradient of primitive variables ---*/
        numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nullptr);

        /*--- Set gradient of auxillary variables ---*/
        numerics->SetAuxVarGrad(nodes->GetAuxVarGradient(iPoint), nullptr);

        /*--- Set diffusion coefficient ---*/
        numerics->SetDiffusionCoeff(nodes->GetDiffusionCoeff(iPoint), nullptr);

        /*--- Laminar viscosity ---*/
        numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint), 0.0);

        /*--- Eddy viscosity ---*/
        numerics->SetEddyViscosity(nodes->GetEddyViscosity(iPoint), 0.0);

        /*--- Thermal conductivity ---*/
        numerics->SetThermalConductivity(nodes->GetThermalConductivity(iPoint), 0.0);

        /*--- Vib-el. thermal conductivity ---*/
        numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint), 0.0);

        /*--- Set turbulence kinetic energy ---*/
        if (rans){
          CVariable* turbNodes = solver_container[TURB_SOL]->GetNodes();
          numerics->SetTurbKineticEnergy(turbNodes->GetSolution(iPoint,0), 0.0);
        }
      }

      auto residual = numerics->ComputeAxisymmetric(config);

      /*--- Check for errors before applying source to the linear system ---*/
      err = CNumerics::CheckResidualNaNs(implicit, nVar, residual);

      /*--- Apply the update to the linear system ---*/
      if (!err) {
        LinSysRes.AddBlock(iPoint, residual);
        if (implicit)
          Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
      } else
        eAxi_local++;
    }
  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  /*--- Checking for NaN ---*/
  unsigned long eAxi_global = eAxi_local;
  unsigned long eChm_global = eChm_local;
  unsigned long eVib_global = eVib_local;

  //THIS IS NO FUN
  if ((eAxi_global != 0) ||
      (eChm_global != 0) ||
      (eVib_global != 0)) {
    cout << "Warning!! Instances of NaN in the following source terms: " << endl;
    cout << "Axisymmetry: " << eAxi_global << endl;
    cout << "Chemical:    " << eChm_global << endl;
    cout << "Vib. Relax:  " << eVib_global << endl;
  }
}

void CNEMOEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                            CConfig *config, unsigned short iRKStep) {

  Explicit_Iteration<RUNGE_KUTTA_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CNEMOEulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                              CConfig *config, unsigned short iRKStep) {

  Explicit_Iteration<CLASSICAL_RK4_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CNEMOEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  Explicit_Iteration<EULER_EXPLICIT>(geometry, solver_container, config, 0);
}

void CNEMOEulerSolver::PrepareImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {

  struct DummyPrec {
    const bool active = false;
    FORCEINLINE su2double** operator() (const CConfig*, unsigned long, su2double) const { return nullptr; }
  } precond;

  PrepareImplicitIteration_impl(precond, geometry, config);
}

void CNEMOEulerSolver::CompleteImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {

  CompleteImplicitIteration_impl<true>(geometry, config);
}

void CNEMOEulerSolver::ComputeUnderRelaxationFactor(const CConfig *config) {

  /* Loop over the solution update given by relaxing the linear
   system for this nonlinear iteration. */

  const su2double allowableRatio = 0.2;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    su2double localUnderRelaxation = 1.0;

    su2double num = 0.0;
    su2double denom = 0.0;

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      /* We impose a limit on the maximum percentage that the
       density (sum of all species) and energy can change over a nonlinear iteration. */

      const unsigned long index = iPoint * nVar + iVar;
      if (iVar < config->GetnSpecies()) {
        num   += fabs(LinSysSol[index]);
        denom += fabs(nodes->GetSolution(iPoint, iVar));

        /*--- If final density/species, compute Under-relaxation ---*/
        if (iVar == (config ->GetnSpecies()-1)){
          su2double ratio = (num/(denom+EPS));
          if (ratio > allowableRatio) {
            localUnderRelaxation = min(allowableRatio / ratio, localUnderRelaxation);
          }
        }

        /*--- Energy ---*/
        if (iVar == (nVar-2)){
          su2double ratio = fabs(LinSysSol[index]) / (fabs(nodes->GetSolution(iPoint, iVar)) + EPS);
          if (ratio > allowableRatio) {
            localUnderRelaxation = min(allowableRatio / ratio, localUnderRelaxation);
          }
        }
      }
    }

    /* Threshold the relaxation factor in the event that there is
     a very small value. This helps avoid catastrophic crashes due
     to non-realizable states by canceling the update. */

    if (localUnderRelaxation < 1e-10) localUnderRelaxation = 0.0;

    /* Store the under-relaxation factor for this point. */

    nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);
  }
  END_SU2_OMP_FOR
}

void CNEMOEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {

  su2double
  Temperature_FreeStream = 0.0, Temperature_ve_FreeStream = 0.0, Mach2Vel_FreeStream         = 0.0,
  ModVel_FreeStream      = 0.0, Energy_FreeStream         = 0.0, ModVel_FreeStreamND         = 0.0,
  Velocity_Reynolds      = 0.0, Omega_FreeStream          = 0.0, Omega_FreeStreamND          = 0.0,
  Viscosity_FreeStream   = 0.0, Density_FreeStream        = 0.0, Pressure_FreeStream         = 0.0,
  Tke_FreeStream         = 0.0, Length_Ref                = 0.0, Density_Ref                 = 0.0,
  Pressure_Ref           = 0.0, Velocity_Ref              = 0.0, Temperature_Ref             = 0.0,
  Temperature_ve_Ref     = 0.0, Time_Ref                  = 0.0, Omega_Ref                   = 0.0,
  Force_Ref              = 0.0, Gas_Constant_Ref          = 0.0, Viscosity_Ref               = 0.0,
  Conductivity_Ref       = 0.0, Energy_Ref                = 0.0, Pressure_FreeStreamND       = 0.0,
  Energy_FreeStreamND    = 0.0, Temperature_FreeStreamND  = 0.0, Temperature_ve_FreeStreamND = 0.0,
  Gas_ConstantND         = 0.0, Viscosity_FreeStreamND    = 0.0, sqvel                       = 0.0,
  Tke_FreeStreamND       = 0.0, Total_UnstTimeND          = 0.0, Delta_UnstTimeND            = 0.0,
  soundspeed             = 0.0, GasConstant_Inf           = 0.0, Froude                      = 0.0,
  Density_FreeStreamND   = 0.0, Heat_Flux_Ref             = 0.0;

  su2double Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0};

  unsigned short iDim;

  /*--- Local variables ---*/
  su2double Alpha         = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta          = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach          = config->GetMach();
  su2double Reynolds      = config->GetReynolds();

  bool unsteady           = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool viscous            = config->GetViscous();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = false;
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == TURB_MODEL::SST));
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);

  /*--- Instatiate the fluid model ---*/
  switch (config->GetKind_FluidModel()) {
  case MUTATIONPP:
   #if defined(HAVE_MPP) && !defined(CODI_REVERSE_TYPE) && !defined(CODI_FORWARD_TYPE)
     FluidModel = new CMutationTCLib(config, nDim);
   #else
     SU2_MPI::Error(string("Either 1) Mutation++ has not been configured/compiled (add '-Denable-mpp=true' to your meson string) or 2) CODI must be deactivated since it is not compatible with Mutation++."),
     CURRENT_FUNCTION);
   #endif
   break;
  case SU2_NONEQ:
   FluidModel = new CSU2TCLib(config, nDim, viscous);
   break;
  }

  /*--- Compute the Free Stream Pressure, Temperatrue, and Density ---*/
  Pressure_FreeStream        = config->GetPressure_FreeStream();
  Temperature_FreeStream     = config->GetTemperature_FreeStream();
  Temperature_ve_FreeStream  = config->GetTemperature_ve_FreeStream();

  /*---                                     ---*/
  /*--- Compressible non dimensionalization ---*/
  /*---                                     ---*/

  /*--- Set mixture state based on pressure, mass fractions and temperatures ---*/
  FluidModel->SetTDStatePTTv(Pressure_FreeStream, MassFrac_Inf,
                             Temperature_FreeStream, Temperature_ve_FreeStream);

  /*--- Compute Gas Constant ---*/
  GasConstant_Inf = FluidModel->ComputeGasConstant();
  config->SetGas_Constant(GasConstant_Inf);

  /*--- Compute the freestream density, soundspeed ---*/
  Density_FreeStream = FluidModel->GetDensity();
  soundspeed         = FluidModel->ComputeSoundSpeed();
  Gamma              = FluidModel->ComputeGamma();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*soundspeed;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*soundspeed;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*soundspeed;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*soundspeed;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*soundspeed;
  }

  /*--- Compute the modulus of the free stream velocity ---*/
  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++){
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  }
  sqvel = ModVel_FreeStream;
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Calculate energies ---*/
  const auto& energies = FluidModel->ComputeMixtureEnergies();

  /*--- Viscous initialization ---*/
  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
          To accomplish this, simply set the non-dimensional coefficients to the
          dimensional ones. This will be overruled later.---*/
    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());
    config->SetMu_ConstantND(config->GetMu_Constant());

    /*--- First, check if there is mesh motion. If yes, use the Mach
         number relative to the body to initialize the flow. ---*/

    if (dynamic_grid) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
    else Velocity_Reynolds = ModVel_FreeStream;

    if (!reynolds_init) {

      /*--- Thermodynamics quantities based initialization ---*/
      Viscosity_FreeStream = FluidModel->GetViscosity();
      Energy_FreeStream    = energies[0] + 0.5*sqvel;

    } else {

      /*--- Reynolds based initialization not present in NEMO ---*/
      SU2_MPI::Error("Only thermodynamics quantities based initialization: set pressure, temperatures and flag INIT_OPTION= TD_CONDITIONS." , CURRENT_FUNCTION);
    }

    config->SetViscosity_FreeStream(Viscosity_FreeStream);

    /*--- Compute Reynolds number ---*/
    Reynolds = (Density_FreeStream*Velocity_Reynolds*config->GetLength_Reynolds())/Viscosity_FreeStream;
    config->SetReynolds(Reynolds);

    /*--- Turbulence kinetic energy ---*/
    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  } else {

    /*--- For inviscid flow, energy is calculated from the specified
       FreeStream quantities using the proper gas law. ---*/
    Energy_FreeStream    = energies[0] + 0.5*sqvel;

  }

  config->SetDensity_FreeStream(Density_FreeStream);

  /*-- Compute the freestream energy. ---*/
  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
     Lref is one because we have converted the grid to meters. ---*/

  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref       = 1.0;
    Density_Ref        = 1.0;
    Temperature_Ref    = 1.0;
    Temperature_ve_Ref = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref       = Pressure_FreeStream;       // Pressure_FreeStream = 1.0
    Density_Ref        = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref    = Temperature_FreeStream;    // Temperature_FreeStream = 1.0
    Temperature_ve_Ref = Temperature_ve_FreeStream; // Temperature_ve_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref       = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref        = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref    = Temperature_FreeStream;    // Temp_FreeStream = 1.0
    Temperature_ve_Ref = Temperature_ve_FreeStream; // Temp_ve_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref       = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref        = Density_FreeStream;                  // Density_FreeStream = 1.0
    Temperature_Ref    = Temperature_FreeStream;              // Temp_FreeStream = 1.0
    Temperature_ve_Ref = Temperature_ve_FreeStream;           // Temp_ve_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);
  config->SetTemperature_ve_Ref(Temperature_ve_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = config->GetDensity_Ref()*Velocity_Ref*Velocity_Ref*Length_Ref*Length_Ref; config->SetForce_Ref(Force_Ref);
  Heat_Flux_Ref     = Density_Ref*Velocity_Ref*Velocity_Ref*Velocity_Ref;           config->SetHeat_Flux_Ref(Heat_Flux_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDARD_GRAVITY*Length_Ref);         config->SetFroude(Froude);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/

  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND    = Temperature_FreeStream/config->GetTemperature_Ref();       config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  Temperature_ve_FreeStreamND = Temperature_ve_FreeStream/config->GetTemperature_ve_Ref(); config->SetTemperature_ve_FreeStreamND(Temperature_ve_FreeStreamND);
  Gas_ConstantND              = config->GetGas_Constant()/Gas_Constant_Ref;                config->SetGas_ConstantND(Gas_ConstantND);

  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);

  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);

  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStream(Omega_FreeStream);

  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/(Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);

  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

  Energy_FreeStreamND = energies[0] + 0.5*ModVel_FreeStreamND *ModVel_FreeStreamND;

  if (viscous) {

    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

    /*--- Sutherland's model ---*/
    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());

    /* constant thermal conductivity model */
    config->SetThermal_Conductivity_ConstantND(config->GetThermal_Conductivity_Constant()/Conductivity_Ref);

  }

  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is the master node and first domain ---*/

  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {

    cout.precision(6);

    if (viscous) {
      if (reynolds_init){
        cout << "Viscous flow: Computing pressure using the equation of state for multi-species and multi-temperatures" << endl;
        cout << "based on the free-stream temperatures and a density computed" << endl;
        cout << "from the Reynolds number." << endl;
      } else {
        cout << "Viscous flow: Computing density using the equation of state for multi-species and multi-temperatures" << endl;
        cout << "based on the free-stream temperatures and pressure." << endl;
      }
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "and pressure using the the equation of state for multi-species and multi-temperatures." << endl;
    }

    if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;

    stringstream NonDimTableOut, ModelTableOut;
    stringstream Unit;

    cout << endl;
    PrintingToolbox::CTablePrinter ModelTable(&ModelTableOut);
    ModelTableOut <<"-- Models:"<< endl;

    ModelTable.AddColumn("Mixture", 25);
    ModelTable.AddColumn("Fluid Model", 25);
    ModelTable.AddColumn("Transport Model", 25);
    ModelTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
    ModelTable.PrintHeader();

    PrintingToolbox::CTablePrinter NonDimTable(&NonDimTableOut);
    NonDimTable.AddColumn("Name", 22);
    NonDimTable.AddColumn("Dim. value", 14);
    NonDimTable.AddColumn("Ref. value", 14);
    NonDimTable.AddColumn("Unit", 10);
    NonDimTable.AddColumn("Non-dim. value", 14);
    NonDimTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);

    NonDimTableOut <<"-- Fluid properties:"<< endl;

    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << Gamma;
    Unit.str("");

    ModelTable << config->GetGasModel();

    switch(config->GetKind_FluidModel()){
    case SU2_NONEQ:
      ModelTable << "SU2 NonEq";
      break;
    case MUTATIONPP:
      ModelTable << "Mutation++ NonEq";
      break;
    }

    if (viscous) {

      switch(config->GetKind_TransCoeffModel()){
      case TRANSCOEFFMODEL::WILKE:
        ModelTable << "Wilke-Blottner-Eucken";
        NonDimTable.PrintFooter();
        break;

      case TRANSCOEFFMODEL::GUPTAYOS:
        ModelTable << "Gupta-Yos";
        NonDimTable.PrintFooter();
        break;

      case TRANSCOEFFMODEL::CHAPMANN_ENSKOG:
        ModelTable << "CHAPMANN-ENSKOG_LDLT";
        NonDimTable.PrintFooter();
        break;

      default:
        break;
      }
    } else {
      ModelTable << "-" ;
    }

    NonDimTable.PrintFooter();
    NonDimTableOut <<"-- Initial and free-stream conditions:"<< endl;
    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Static Pressure" << config->GetPressure_FreeStream() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "kg/m^3";
    else if (config->GetSystemMeasurements() == US) Unit << "slug/ft^3";
    NonDimTable << "Density" << config->GetDensity_FreeStream() << config->GetDensity_Ref() << Unit.str() << config->GetDensity_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "K";
    else if (config->GetSystemMeasurements() == US) Unit << "R";
    NonDimTable << " T-R Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "K";
    else if (config->GetSystemMeasurements() == US) Unit << "R";
    NonDimTable << " V-E Temperature" << config->GetTemperature_ve_FreeStream() << config->GetTemperature_ve_Ref() << Unit.str() << config->GetTemperature_ve_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
    else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
    NonDimTable << "Total Energy" << config->GetEnergy_FreeStream() << config->GetEnergy_Ref() << Unit.str() << config->GetEnergy_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m/s";
    else if (config->GetSystemMeasurements() == US) Unit << "ft/s";
    NonDimTable << "Velocity-X" << config->GetVelocity_FreeStream()[0] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[0];
    NonDimTable << "Velocity-Y" << config->GetVelocity_FreeStream()[1] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[1];
    if (nDim == 3){
      NonDimTable << "Velocity-Z" << config->GetVelocity_FreeStream()[2] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[2];
    }
    NonDimTable << "Velocity Magnitude" << config->GetModVel_FreeStream() << config->GetVelocity_Ref() << Unit.str() << config->GetModVel_FreeStreamND();
    Unit.str("");

    if (viscous){
      NonDimTable.PrintFooter();
      if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
      NonDimTable << "Viscosity" << config->GetViscosity_FreeStream() << config->GetViscosity_Ref() << Unit.str() << config->GetViscosity_FreeStreamND();
      Unit.str("");
      if (turbulent){
        if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
        else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
        NonDimTable << "Turb. Kin. Energy" << config->GetTke_FreeStream() << config->GetTke_FreeStream()/config->GetTke_FreeStreamND() << Unit.str() << config->GetTke_FreeStreamND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "1/s";
        else if (config->GetSystemMeasurements() == US) Unit << "1/s";
        NonDimTable << "Spec. Dissipation" << config->GetOmega_FreeStream() << config->GetOmega_FreeStream()/config->GetOmega_FreeStreamND() << Unit.str() << config->GetOmega_FreeStreamND();
        Unit.str("");
      }
    }

    NonDimTable.PrintFooter();
    NonDimTable << "Mach Number" << "-" << "-" << "-" << config->GetMach();
    if (viscous){
      NonDimTable << "Reynolds Number" << "-" << "-" << "-" << config->GetReynolds();
    }
    if (gravity) {
      NonDimTable << "Froude Number" << "-" << "-" << "-" << Froude;
      NonDimTable << "Wave Length"   << "-" << "-" << "-" << 2.0*PI_NUMBER*Froude*Froude;
    }
    NonDimTable.PrintFooter();
    ModelTable.PrintFooter();

    if (unsteady){
      NonDimTableOut << "-- Unsteady conditions" << endl;
      NonDimTable.PrintHeader();
      NonDimTable << "Total Time" << config->GetMax_Time() << config->GetTime_Ref() << "s" << config->GetMax_Time()/config->GetTime_Ref();
      Unit.str("");
      NonDimTable << "Time Step" << config->GetTime_Step() << config->GetTime_Ref() << "s" << config->GetDelta_UnstTimeND();
      Unit.str("");
      NonDimTable.PrintFooter();
    }

    cout << ModelTableOut.str();
    cout << NonDimTableOut.str();

  }
}

void CNEMOEulerSolver::SetReferenceValues(const CConfig& config) {

  DynamicPressureRef = 0.5 * Density_Inf * GeometryToolbox::SquaredNorm(nDim, Velocity_Inf);
  AeroCoeffForceRef =  DynamicPressureRef * config.GetRefArea();

}

void CNEMOEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim, jDim, iSpecies, iVar, jVar;
  unsigned long iPoint, iVertex;

  /*--- Allocate the necessary vector structures ---*/
  su2double Normal[MAXNDIM] = {0.0}, UnitNormal[MAXNDIM] = {0.0};

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negative for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      /*--- Calculate parameters from the geometry ---*/
      su2double Area = GeometryToolbox::Norm(nDim, Normal);

      for (iDim = 0; iDim < nDim; iDim++){
        UnitNormal[iDim] = -Normal[iDim]/Area;
      }

      /*--- Retrieve the pressure on the vertex ---*/
      su2double P = nodes->GetPressure(iPoint);

      /*--- Apply the flow-tangency b.c. to the convective flux ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        Residual[iSpecies] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        Residual[nSpecies+iDim] = P * UnitNormal[iDim] * Area;
      }
      Residual[nSpecies+nDim]   = 0.0;
      Residual[nSpecies+nDim+1] = 0.0;

      /*--- Add value to the residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- If using implicit time-stepping, calculate b.c. contribution to Jacobian ---*/
      if (implicit) {

        /*--- Allocate arrays needed for implicit solver ---*/
        su2double Velocity[MAXNDIM] = {0.0};

        /*--- Get species molar mass ---*/
        auto& Ms = FluidModel->GetSpeciesMolarMass();

        /*--- Initialize Jacobian ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;

        /*--- Calculate state i ---*/
        su2double rho     = nodes->GetDensity(iPoint);
        su2double rhoE    = nodes->GetSolution(iPoint,nSpecies+nDim);
        su2double rhoEve  = nodes->GetSolution(iPoint,nSpecies+nDim+1);
        auto dPdU    = nodes->GetdPdU(iPoint);
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);

        su2double conc = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          su2double cs = nodes->GetMassFraction(iPoint,iSpecies);
          conc += cs * rho/Ms[iSpecies];

          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_i[nSpecies+iDim][iSpecies] = dPdU[iSpecies] * UnitNormal[iDim];
            Jacobian_i[iSpecies][nSpecies+iDim] = cs * UnitNormal[iDim];
          }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            Jacobian_i[nSpecies+iDim][nSpecies+jDim] = Velocity[iDim]*UnitNormal[jDim]
                + dPdU[nSpecies+jDim]*UnitNormal[iDim];
          }
          Jacobian_i[nSpecies+iDim][nSpecies+nDim]   = dPdU[nSpecies+nDim]  *UnitNormal[iDim];
          Jacobian_i[nSpecies+iDim][nSpecies+nDim+1] = dPdU[nSpecies+nDim+1]*UnitNormal[iDim];

          Jacobian_i[nSpecies+nDim][nSpecies+iDim]   = (rhoE+P)/rho * UnitNormal[iDim];
          Jacobian_i[nSpecies+nDim+1][nSpecies+iDim] = rhoEve/rho   * UnitNormal[iDim];
        }

        /*--- Integrate over the dual-grid area ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = Jacobian_i[iVar][jVar] * Area;

        /*--- Apply the contribution to the system ---*/
        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

      }
    }
  }
}

void CNEMOEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double *V_infty, *V_domain, *U_infty, *U_domain;

  /*--- Set booleans from configuration parameters ---*/
  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool viscous  = config->GetViscous();

  /*--- Allocate arrays ---*/
  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Allocate the value at the infinity ---*/
    V_infty = GetCharacPrimVar(val_marker, iVertex);

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Retrieve index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor(); //only used for implicit

      /*--- Pass boundary node normal to CNumerics ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at the boundary node & free-stream ---*/
      U_domain = nodes->GetSolution(iPoint);
      V_domain = nodes->GetPrimitive(iPoint);
      U_infty  = node_infty->GetSolution(0);
      V_infty  = node_infty->GetPrimitive(0);

      /*--- Pass conserved & primitive variables to CNumerics ---*/
      conv_numerics->SetConservative(U_domain, U_infty);
      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Pass supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU  (nodes->GetdPdU(iPoint),   node_infty->GetdPdU(0));
      conv_numerics->SetdTdU  (nodes->GetdTdU(iPoint),   node_infty->GetdTdU(0));
      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));
      conv_numerics->SetEve   (nodes->GetEve(iPoint),    node_infty->GetEve(0));
      conv_numerics->SetCvve  (nodes->GetCvve(iPoint),   node_infty->GetCvve(0));
      conv_numerics->SetGamma (nodes->GetGamma(iPoint),  node_infty->GetGamma(0));

      /*--- Compute the convective residual (and Jacobian) ---*/
      // Note: This uses the specified boundary num. method specified in driver_structure.cpp
      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Apply contribution to the linear system ---*/
      LinSysRes.AddBlock(iPoint, residual);

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      /*--- Viscous contribution ---*/
      if (viscous) {
        su2double Coord_Reflected[MAXNDIM];
        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected );
        visc_numerics->SetNormal(Normal);

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetConservative(nodes->GetSolution(iPoint),
                                       node_infty->GetSolution(0) );
        visc_numerics->SetPrimitive(nodes->GetPrimitive(iPoint),
                                    node_infty->GetPrimitive(0) );
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          node_infty->GetGradient_Primitive(0) );

        /*--- Pass supplementary information to CNumerics ---*/
        visc_numerics->SetdPdU  (nodes->GetdPdU(iPoint),   node_infty->GetdPdU(0));
        visc_numerics->SetdTdU  (nodes->GetdTdU(iPoint),   node_infty->GetdTdU(0));
        visc_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));
        visc_numerics->SetEve   (nodes->GetEve(iPoint),    node_infty->GetEve(0));
        visc_numerics->SetCvve  (nodes->GetCvve(iPoint),   node_infty->GetCvve(0));

        /*--- Species diffusion coefficients ---*/
        visc_numerics->SetDiffusionCoeff(nodes->GetDiffusionCoeff(iPoint),
                                         nodes->GetDiffusionCoeff(iPoint));

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint),
                                           nodes->GetLaminarViscosity(iPoint));

        /*--- Eddy viscosity ---*/
        visc_numerics->SetEddyViscosity(nodes->GetEddyViscosity(iPoint),
                                        nodes->GetEddyViscosity(iPoint));

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(nodes->GetThermalConductivity(iPoint),
                                              nodes->GetThermalConductivity(iPoint));

        /*--- Vib-el. thermal conductivity ---*/
        visc_numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint),
                                                 nodes->GetThermalConductivity_ve(iPoint));

        /*--- Compute and update residual ---*/
        auto residual = visc_numerics->ComputeResidual(config);

        LinSysRes.SubtractBlock(iPoint, residual);
        if (implicit) {
          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;
}

void CNEMOEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  SU2_MPI::Error("BC_INLET: Not operational in NEMO.", CURRENT_FUNCTION);

  unsigned short iVar, iDim, iSpecies, RHO_INDEX, nSpecies;

  unsigned long iVertex, iPoint;
  su2double  T_Total, P_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3] = {0.0};

  const su2double *Flow_Dir;
  bool implicit             = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  bool dynamic_grid         = config->GetGrid_Movement();
  su2double Two_Gamma_M1    = 2.0/Gamma_Minus_One;
  su2double Gas_Constant    = config->GetGas_ConstantND();
  INLET_TYPE Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  su2double *U_domain = new su2double[nVar];      su2double *U_inlet = new su2double[nVar];
  su2double *V_domain = new su2double[nPrimVar];  su2double *V_inlet = new su2double[nPrimVar];
  su2double *Normal   = new su2double[nDim];

  nSpecies = config->GetnSpecies();
  su2double *Spec_Density = new su2double[nSpecies];
  for(iSpecies=0; iSpecies<nSpecies; iSpecies++)
    Spec_Density[iSpecies] = 0.0;               /*--- To avoid a compiler warning. ---*/

  RHO_INDEX = nodes->GetRhoIndex();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Retrieve solution at this boundary node ---*/
      for (iVar = 0; iVar < nVar; iVar++)     U_domain[iVar] = nodes->GetSolution(iPoint, iVar);
      for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Build the fictitious intlet state based on characteristics ---*/

      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
       therefore we can specify all but one state variable at the inlet.
       The outgoing Riemann invariant provides the final piece of info.
       Adapted from an original implementation in the Stanford University
       multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
       written by Edwin van der Weide, last modified 04-20-2009. ---*/

      switch (Kind_Inlet) {

      /*--- Total properties have been specified at the inlet. ---*/
      case INLET_TYPE::TOTAL_CONDITIONS:

        /*--- Retrieve the specified total conditions for this inlet. ---*/
        P_Total  = config->GetInlet_Ptotal(Marker_Tag);
        T_Total  = config->GetInlet_Ttotal(Marker_Tag);
        Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_Total /= config->GetPressure_Ref();
        T_Total /= config->GetTemperature_Ref();

        /*--- Store primitives and set some variables for clarity. ---*/
        Density = V_domain[RHO_INDEX];
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = U_domain[nSpecies+iDim]/Density;

        Velocity2   = GeometryToolbox::SquaredNorm(nDim, Velocity);
        Energy      = U_domain[nVar-2]/Density;
        Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
        H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
        SoundSpeed2 = Gamma*Pressure/Density;

        /*--- Compute the acoustic Riemann invariant that is extrapolated
           from the domain interior. ---*/
        Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
        for (iDim = 0; iDim < nDim; iDim++)
          Riemann += Velocity[iDim]*UnitNormal[iDim];

        /*--- Total speed of sound ---*/
        SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

        /*--- Dot product of normal and flow direction. This should
           be negative due to outward facing boundary normal convention. ---*/
        alpha = GeometryToolbox::DotProduct(nDim, UnitNormal, Flow_Dir);

        /*--- Coefficients in the quadratic equation for the velocity ---*/
        aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
        bb = -1.0*Gamma_Minus_One*alpha*Riemann;
        cc =  0.5*Gamma_Minus_One*Riemann*Riemann
            -2.0*SoundSpeed_Total2/Gamma_Minus_One;

        /*--- Solve quadratic equation for velocity magnitude. Value must
           be positive, so the choice of root is clear. ---*/
        dd = bb*bb - 4.0*aa*cc;
        dd = sqrt(max(0.0,dd));
        Vel_Mag   = (-bb + dd)/(2.0*aa);
        Vel_Mag   = max(0.0,Vel_Mag);
        Velocity2 = Vel_Mag*Vel_Mag;

        /*--- Compute speed of sound from total speed of sound eqn. ---*/
        SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

        /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
        Mach2 = Velocity2/SoundSpeed2;
        Mach2 = min(1.0,Mach2);
        Velocity2   = Mach2*SoundSpeed2;
        Vel_Mag     = sqrt(Velocity2);
        SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

        /*--- Compute new velocity vector at the inlet ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

        /*--- Static temperature from the speed of sound relation ---*/
        Temperature = SoundSpeed2/(Gamma*Gas_Constant);
        //NEED TVE AS WELL

        /*--- Static pressure using isentropic relation at a point ---*/
        Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);

        /*--- Density at the inlet from the gas law ---*/
        Density = Pressure/(Gas_Constant*Temperature);
        //NEED SPECIES DENSITIES

        /*--- Using pressure, density, & velocity, compute the energy ---*/
        Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
        //NEED EVE AS WELL

        /*--- Conservative variables, using the derived quantities ---*/
        for (iSpecies=0; iSpecies<nSpecies; iSpecies++)
          U_inlet[iSpecies] = Spec_Density[iSpecies];
        for (iDim = 0; iDim < nDim; iDim++)
          U_inlet[nSpecies+iDim] = Velocity[iDim]*Density;
        U_inlet[nVar-2] = Energy*Density;
        //U_inlet[nVar-1]=Eve

        /*--- Primitive variables, using the derived quantities ---*/
        for (iSpecies=0; iSpecies<nSpecies; iSpecies++)
          V_inlet[iSpecies] = Spec_Density[iSpecies];
        V_inlet[nSpecies] = Temperature;
        //V_inlet[nSpecies+1] = Tve
        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[nSpecies+2] = Velocity[iDim];
        V_inlet[nSpecies+nDim+2] = Pressure;
        V_inlet[RHO_INDEX] = Density;
        //V_inlet[H_INDEX] = H;
        //V_inlet[A_INDEX] = A;
        //V_inlet[RHO_CVTR_INDEX] = rcvtr;
        //V_inlet[RHO_CVVE_INDEX] = rcvve;

        break;

        /*--- Mass flow has been specified at the inlet. ---*/
      case INLET_TYPE::MASS_FLOW:

        /*--- Retrieve the specified mass flow for the inlet. ---*/
        Density  = config->GetInlet_Ttotal(Marker_Tag);
        Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
        Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        Density /= config->GetDensity_Ref();
        Vel_Mag /= config->GetVelocity_Ref();

        /*--- Get primitives from current inlet state. ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = nodes->GetVelocity(iPoint, iDim);
        Pressure    = nodes->GetPressure(iPoint);
        SoundSpeed2 = Gamma*Pressure/U_domain[0];

        /*--- Compute the acoustic Riemann invariant that is extrapolated
           from the domain interior. ---*/
        Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
        for (iDim = 0; iDim < nDim; iDim++)
          Riemann += Velocity[iDim]*UnitNormal[iDim];

        /*--- Speed of sound squared for fictitious inlet state ---*/
        SoundSpeed2 = Riemann;
        for (iDim = 0; iDim < nDim; iDim++)
          SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitNormal[iDim];

        SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
        SoundSpeed2 = SoundSpeed2*SoundSpeed2;

        /*--- Pressure for the fictitious inlet state ---*/
        Pressure = SoundSpeed2*Density/Gamma;

        /*--- Energy for the fictitious inlet state ---*/
        Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Vel_Mag*Vel_Mag;

        /*--- Conservative variables, using the derived quantities ---*/
        U_inlet[0] = Density;
        for (iDim = 0; iDim < nDim; iDim++)
          U_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim]*Density;
        U_inlet[nDim+1] = Energy*Density;

        /*--- Primitive variables, using the derived quantities ---*/
        V_inlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
        V_inlet[nDim+1] = Pressure;
        V_inlet[nDim+2] = Density;

        break;

      default:
        SU2_MPI::Error("Unsupported INLET_TYPE.", CURRENT_FUNCTION);
        break;
      }

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetConservative(U_domain, U_inlet);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/
      auto residual = conv_numerics->ComputeResidual(config);

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

//      /*--- Viscous contribution ---*/
//      if (viscous) {

//        /*--- Set the normal vector and the coordinates ---*/
//        visc_numerics->SetNormal(Normal);
//        su2double Coord_Reflected[MAXNDIM];
//        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//        visc_numerics->SetCoord(geometry->nodes->GetCoord(), Coord_Reflected);

//        /*--- Primitive variables, and gradient ---*/
//        visc_numerics->SetPrimitive(V_domain, V_inlet);
//        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(), nodes->GetGradient_Primitive());

//        /*--- Laminar viscosity ---*/
//        visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(), nodes->GetLaminarViscosity());

//        /*--- Compute and update residual ---*/
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);

//        /*--- Jacobian contribution for implicit integration ---*/
//        if (implicit)
//          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] U_domain;
  delete [] U_inlet;
  delete [] V_domain;
  delete [] V_inlet;
  delete [] Normal;
  delete [] Spec_Density;
}

void CNEMOEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iVar, iDim, iSpecies;
  unsigned long iVertex, iPoint;
  su2double Pressure, P_Exit, Velocity[3], Temperature, Tve, Velocity2, Entropy, Density,
  Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit, Area, UnitNormal[3];
  vector<su2double> rhos;

  rhos.resize(nSpecies,0.0);

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool dynamic_grid = config->GetGrid_Movement();
  bool gravity      = config->GetGravityForce();
  bool implicit     = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  su2double *U_domain;
  su2double *U_outlet = new su2double[nVar];
  su2double *V_domain, *V_outlet;
  su2double *Normal   = new su2double[nDim];
  su2double *Ys       = new su2double[nSpecies];

  unsigned short T_INDEX       = nodes->GetTIndex();
  unsigned short TVE_INDEX     = nodes->GetTveIndex();
  unsigned short VEL_INDEX     = nodes->GetVelIndex();
  unsigned short P_INDEX       = nodes->GetPIndex();
  unsigned short RHO_INDEX     = nodes->GetRhoIndex();
  unsigned short H_INDEX       = nodes->GetHIndex();
  unsigned short A_INDEX       = nodes->GetAIndex();
  unsigned short RHOCVTR_INDEX = nodes->GetRhoCvtrIndex();
  unsigned short RHOCVVE_INDEX = nodes->GetRhoCvveIndex();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(val_marker, iVertex);

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Current solution at this boundary node ---*/
      U_domain = nodes->GetSolution(iPoint);
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Initialize solution at outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++)     U_outlet[iVar] = 0.0;
      for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = 0.0;

      /*--- Build the fictitious intlet state based on characteristics ---*/

      /*--- Compute Gamma using domain state ---*/
      Gamma = nodes->GetGamma(iPoint);
      Gamma_Minus_One = Gamma - 1.0;

      /*--- Retrieve the specified back pressure for this outlet. ---*/
      if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->nodes->GetCoord(iPoint, nDim-1)*STANDARD_GRAVITY;
      else         P_Exit = config->GetOutlet_Pressure(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
       of boundary update depends on this. ---*/
      Density = V_domain[RHO_INDEX];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[VEL_INDEX+iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Temperature = V_domain[T_INDEX];
      Tve         = V_domain[TVE_INDEX];
      Pressure    = V_domain[P_INDEX];
      SoundSpeed  = V_domain[A_INDEX];
      Mach_Exit   = sqrt(Velocity2)/SoundSpeed;

      /*--- Compute Species Concentrations ---*/
      for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
        Ys[iSpecies] = V_domain[iSpecies]/Density;
      }

      /*--- Recompute boundary state depending Mach number ---*/
      if (Mach_Exit >= 1.0) {

        /*--- Supersonic exit flow: there are no incoming characteristics,
         so no boundary condition is necessary. Set outlet state to current
         state so that upwinding handles the direction of propagation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)     U_outlet[iVar] = U_domain[iVar];
        for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];

      } else {

        /*--- Subsonic exit flow: there is one incoming characteristic,
         therefore one variable can be specified (back pressure) and is used
         to update the conservative variables. Compute the entropy and the
         acoustic Riemann variable. These invariants, as well as the
         tangential velocity components, are extrapolated. The Temperatures
         (T and Tve) and species concentraition are also assumed to be extrapolated.
         ---*/

        Entropy = Pressure*pow(1.0/Density,Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        //     U: [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
        //     V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        Pressure   = P_Exit;
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }

        /*--- Primitive variables, using the derived quantities ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
          V_outlet[iSpecies] = Ys[iSpecies]*Density;
          rhos[iSpecies]     = V_outlet[iSpecies];
        }

        V_outlet[T_INDEX]     = V_domain[T_INDEX];
        V_outlet[TVE_INDEX]   = V_domain[TVE_INDEX];

        for (iDim = 0; iDim < nDim; iDim++){
          V_outlet[VEL_INDEX+iDim] = Velocity[iDim];
        }

        V_outlet[P_INDEX]     = Pressure;
        V_outlet[RHO_INDEX]   = Density;
        V_outlet[A_INDEX]     = SoundSpeed;

        /*--- Set mixture state and compute quantities ---*/
        FluidModel->SetTDStateRhosTTv(rhos, Temperature, Tve);
        V_outlet[RHOCVTR_INDEX] = FluidModel->ComputerhoCvtr();
        V_outlet[RHOCVVE_INDEX] = FluidModel->ComputerhoCvve();

        const auto& energies = FluidModel->ComputeMixtureEnergies();

        /*--- Conservative variables, using the derived quantities ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
          U_outlet[iSpecies] = V_outlet[iSpecies];
        }

        for (iDim = 0; iDim < nDim; iDim++)
          U_outlet[nSpecies+iDim] = Velocity[iDim]*Density;

        U_outlet[nVar-2] = (energies[0] + 0.5*Velocity2) * Density;
        U_outlet[nVar-1] = energies[1] * Density;

      }

      /*--- Setting Last remaining variables ---*/
      V_outlet[H_INDEX]= (U_outlet[nVar-2]+Pressure)/Density;

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetConservative(U_domain, U_outlet);
      conv_numerics->SetPrimitive(V_domain,V_outlet);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Passing supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU  (nodes->GetdPdU(iPoint),   node_infty->GetdPdU(0));
      conv_numerics->SetdTdU  (nodes->GetdTdU(iPoint),   node_infty->GetdTdU(0));
      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));
      conv_numerics->SetEve   (nodes->GetEve(iPoint),    node_infty->GetEve(0));
      conv_numerics->SetCvve  (nodes->GetCvve(iPoint),   node_infty->GetCvve(0));
      conv_numerics->SetGamma (nodes->GetGamma(iPoint),  node_infty->GetGamma(0));

      /*--- Compute the residual using an upwind scheme ---*/
      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      /*--- Viscous contribution ---*/
//      if (viscous) {

//        /*--- Set the normal vector and the coordinates ---*/
//        visc_numerics->SetNormal(Normal);
//        su2double Coord_Reflected[MAXNDIM];
//        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                                 geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//        visc_numerics->SetCoord(geometry->nodes->GetCoord(), Coord_Reflected);

//        /*--- Primitive variables, and gradient ---*/
//        visc_numerics->SetPrimitive(V_domain, V_outlet);
//        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(), nodes->GetGradient_Primitive());

//        /*--- Conservative variables, and gradient ---*/
//        visc_numerics->SetConservative(U_domain, U_outlet);

//        /*--- Pass supplementary information to CNumerics ---*/
//        visc_numerics->SetdPdU(nodes->GetdPdU(), node_infty->GetdPdU());
//        visc_numerics->SetdTdU(nodes->GetdTdU(), node_infty->GetdTdU());
//        visc_numerics->SetdTvedU(nodes->GetdTvedU(), node_infty->GetdTvedU());

//        /*--- Species diffusion coefficients ---*/
//        visc_numerics->SetDiffusionCoeff(nodes->GetDiffusionCoeff(),
//                                         node_infty->GetDiffusionCoeff() );

//        /*--- Laminar viscosity ---*/
//        visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(),
//                                           node_infty->GetLaminarViscosity() );

//        /*--- Thermal conductivity ---*/
//        visc_numerics->SetThermalConductivity(nodes->GetThermalConductivity(),
//                                              node_infty->GetThermalConductivity());

//        /*--- Vib-el. thermal conductivity ---*/
//        visc_numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(),
//                                                 node_infty->GetThermalConductivity_ve() );

//        /*--- Laminar viscosity ---*/
//        visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(), nodes->GetLaminarViscosity());

//        /*--- Compute and update residual ---*/
//        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//        LinSysRes.SubtractBlock(iPoint, Residual);

//        /*--- Jacobian contribution for implicit integration ---*/
//        if (implicit)
//          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
//      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] U_outlet;
  delete [] Normal;
  delete [] Ys;

}

void CNEMOEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                           CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

SU2_MPI::Error("BC_SUPERSONIC_INLET: Not operational in NEMO.", CURRENT_FUNCTION);

//  unsigned short iDim, iVar;
//  unsigned long iVertex, iPoint, Point_Normal;
//  su2double Density, Pressure, Temperature, Temperature_ve, Energy, *Velocity, Velocity2, soundspeed;
//  su2double Gas_Constant = config->GetGas_ConstantND();
//
//  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
//  bool dynamic_grid  = config->GetGrid_Movement();
//  bool viscous              = config->GetViscous();
//  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
//
//  su2double RuSI  = UNIVERSAL_GAS_CONSTANT;
//  su2double Ru = 1000.0*RuSI;
//
//  su2double *U_inlet = new su2double[nVar];     su2double *U_domain = new su2double[nVar];
//  su2double *V_inlet = new su2double[nPrimVar]; su2double *V_domain = new su2double[nPrimVar];
//  su2double *Normal = new su2double[nDim];

/* ----------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------- */
/* The block of code commented below needs to be updated to use Fluidmodel class */
/* ----------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------- */

//  /*--- Supersonic inlet flow: there are no outgoing characteristics,
//   so all flow variables can be imposed at the inlet.
//   First, retrieve the specified values for the primitive variables. ---*/
//  //ASSUME TVE = T for the time being
//  auto Mass_Frac      = config->GetInlet_MassFrac(Marker_Tag);
//  Temperature    = config->GetInlet_Temperature(Marker_Tag);
//  Pressure       = config->GetInlet_Pressure(Marker_Tag);
//  Velocity       = config->GetInlet_Velocity(Marker_Tag);
//  Temperature_ve = Temperature;
//
//  /*--- Compute Density and Species Densities ---*/
//  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
//    denom += Mass_Frac[iSpecies] * (Ru/Ms[iSpecies]) * Temperature;
//  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
//    denom += Mass_Frac[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Temperature_ve;
//  Density = Pressure / denom;
//
//  /*--- Compute Soundspeed and Velocity squared ---*/
//  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
//    conc += Mass_Frac[iSpecies]*Density/Ms[iSpecies];
//    rhoCvtr += Density*Mass_Frac[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
//  }
//  soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * Pressure/Density);
//
//  Velocity2 = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    Velocity2 += Velocity[iDim]*Velocity[iDim];
//
//  /*--- Non-dim. the inputs if necessary. ---*/
//  // Need to update this portion
//  //Temperature = Temperature/config->GetTemperature_Ref();
//  //Pressure    = Pressure/config->GetPressure_Ref();
//  //Density     = Density/config->GetDensity_Ref();
//  //for (iDim = 0; iDim < nDim; iDim++)
//  //  Velocity[iDim] = Velocity[iDim]/config->GetVelocity_Ref();
//
//  /*--- Compute energy (RRHO) from supplied primitive quanitites ---*/
//  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
//
//    // Species density
//    rhos = Mass_Frac[iSpecies]*Density;
//
//    // Species formation energy
//    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
//
//    // Species vibrational energy
//    if (thetav[iSpecies] != 0.0)
//      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Temperature_ve)-1.0);
//    else
//      Ev = 0.0;
//
//    // Species electronic energy
//    num = 0.0;
//    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Temperature_ve);
//    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
//      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Temperature_ve);
//      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Temperature_ve);
//    }
//    Ee = Ru/Ms[iSpecies] * (num/denom);
//
//    // Mixture total energy
//    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (Temperature-Tref[iSpecies])
//                    + Ev + Ee + Ef + 0.5*Velocity2);
//
//    // Mixture vibrational-electronic energy
//    rhoEve += rhos * (Ev + Ee);
//  }
//
//  /*--- Setting Conservative Variables ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    U_inlet[iSpecies] = Mass_Frac[iSpecies]*Density;
//  for (iDim = 0; iDim < nDim; iDim++)
//    U_inlet[nSpecies+iDim] = Density*Velocity[iDim];
//  U_inlet[nVar-2] = rhoE;
//  U_inlet[nVar-1] = rhoEve;
//
//  /*--- Setting Primitive Vaariables ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    V_inlet[iSpecies] = Mass_Frac[iSpecies]*Density;
//  V_inlet[nSpecies] = Temperature;
//  V_inlet[nSpecies+1] = Temperature_ve;
//  for (iDim = 0; iDim < nDim; iDim++)
//    V_inlet[nSpecies+2+iDim] = Velocity[iDim];
//  V_inlet[nSpecies+2+nDim] = Pressure;
//  V_inlet[nSpecies+3+nDim] = Density;
//  V_inlet[nSpecies+4+nDim] = rhoE+Pressure/Density;
//  V_inlet[nSpecies+5+nDim] = soundspeed;
//  V_inlet[nSpecies+6+nDim] = rhoCvtr;

/* ----------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------- */
/* The block of code that needs to be updated to use Fluidmodel class end here   */
/* ----------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------- */

//  //This requires Newtown Raphson.....So this is not currently operational (See Deathstar)
//  //V_inlet[nSpecies+7+nDim] = rhoCvve;
//
//  /*--- Loop over all the vertices on this boundary marker ---*/
//  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//    if (geometry->nodes->GetDomain(iPoint)) {
//
//      /*--- Index of the closest interior node ---*/
//      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
//
//      /*--- Current solution at this boundary node ---*/
//      for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = nodes->GetSolution(iPoint,iVar);
//      for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = nodes->GetPrimitive(iPoint,iVar);
//
//      /*--- Normal vector for this vertex (negate for outward convention) ---*/
//      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
//      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
//
//      su2double Area = 0.0;
//      for (iDim = 0; iDim < nDim; iDim++)
//        Area += Normal[iDim]*Normal[iDim];
//      Area = sqrt (Area);
//
//      /*--- Set various quantities in the solver class ---*/
//      conv_numerics->SetNormal(Normal);
//      conv_numerics->SetConservative(U_domain, U_inlet);
//      conv_numerics->SetPrimitive(V_domain, V_inlet);
//
//      /*--- Pass supplementary info to CNumerics ---*/
//      conv_numerics->SetdPdU(nodes->GetdPdU(iPoint), node_infty->GetdPdU(0));
//      conv_numerics->SetdTdU(nodes->GetdTdU(iPoint), node_infty->GetdTdU(0));
//      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));
//      conv_numerics->SetEve(nodes->GetEve(iPoint), node_infty->GetEve(0));
//      conv_numerics->SetCvve(nodes->GetCvve(iPoint), node_infty->GetCvve(0));
//
//      if (dynamic_grid)
//        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
//                                  geometry->nodes->GetGridVel(iPoint));
//
//      /*--- Compute the residual using an upwind scheme ---*/
//      auto residual = conv_numerics->ComputeResidual(config);
//      LinSysRes.AddBlock(iPoint, residual);
//
//      /*--- Jacobian contribution for implicit integration ---*/
//      //if (implicit)
//      //  Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
//
//      /*--- Viscous contribution ---*/
//      if (viscous) {
//
//        /*--- Set the normal vector and the coordinates ---*/
//        visc_numerics->SetNormal(Normal);
//        su2double Coord_Reflected[MAXNDIM];
//        GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
//                                           geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);
//
//        /*--- Primitive variables, and gradient ---*/
//        visc_numerics->SetPrimitive(V_domain, V_inlet);
//        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
//
//        /*--- Laminar viscosity ---*/
//        visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint), nodes->GetLaminarViscosity(iPoint));
//
//        /*--- Compute and update residual ---*/
//        auto residual = visc_numerics->ComputeResidual(config);
//        LinSysRes.SubtractBlock(iPoint, residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//        //if (implicit)
//        //  Jacobian.SubtractBlock2Diag(iPoint, residual.Jacobian_i);
//      }
//
//    }
//  }
//
//  /*--- Free locally allocated memory ---*/
//  delete [] U_domain;
//  delete [] U_inlet;
//  delete [] V_domain;
//  delete [] V_inlet;
//  delete [] Normal;

}

void CNEMOEulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_outlet, *V_domain;
  su2double *U_outlet, *U_domain;

  bool implicit     = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double *Normal = new su2double[nDim];

  /*--- Supersonic outlet flow: there are no ingoing characteristics,
   so all flow variables can should be interpolated from the domain. ---*/

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Current solution at this boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);
      U_domain = nodes->GetSolution(iPoint);

      /*--- Allocate the value at the outlet ---*/
      V_outlet = GetCharacPrimVar(val_marker, iVertex);
      V_outlet = V_domain;
      U_outlet = U_domain;

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      conv_numerics->SetConservative(U_domain, U_outlet);

      /*--- Pass supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU  (nodes->GetdPdU(iPoint),   nodes->GetdPdU(iPoint));
      conv_numerics->SetdTdU  (nodes->GetdTdU(iPoint),   nodes->GetdTdU(iPoint));
      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(iPoint));
      conv_numerics->SetEve   (nodes->GetEve(iPoint),    nodes->GetEve(iPoint));
      conv_numerics->SetCvve  (nodes->GetCvve(iPoint),   nodes->GetCvve(iPoint));
      conv_numerics->SetGamma (nodes->GetGamma(iPoint),  nodes->GetGamma(iPoint));

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/
      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}
