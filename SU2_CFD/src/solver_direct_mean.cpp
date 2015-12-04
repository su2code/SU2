/*!
 * \file solution_direct_mean.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 4.0.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/solver_structure.hpp"

CEulerSolver::CEulerSolver(void) : CSolver() {
  
  /*--- Basic array initialization ---*/
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  
  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  
  /*--- Surface based array initialization ---*/
  
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL; Surface_CSideForce_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  
  Surface_CLift = NULL; Surface_CDrag = NULL; Surface_CSideForce = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;
  
  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;
  
  /*--- Supersonic simulation array initialization ---*/
  
  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;
  
  /*--- Engine simulation array initialization ---*/
  
  Inflow_MassFlow = NULL;   Inflow_Pressure = NULL;
  Inflow_Mach = NULL;       Inflow_Area = NULL;
  Bleed_MassFlow = NULL;    Bleed_Pressure = NULL;
  Bleed_Temperature = NULL; Inflow_Area = NULL;
  Exhaust_Pressure = NULL;  Exhaust_Temperature = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  
  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  
  /*--- Fixed CL mode initialization (cauchy criteria) ---*/
  
  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  
}

CEulerSolver::CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double StaticEnergy, Density, Velocity2, Pressure, Temperature, dull_val;
  int Unst_RestartIter;
  ifstream restart_file;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool adjoint = config->GetAdjoint();
  string filename = config->GetSolution_FlowFileName();
  
  unsigned short direct_diff = config->GetDirectDiff();
  unsigned short nMarkerTurboPerf = config->Get_nMarkerTurboPerf();
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Array initialization ---*/
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL; CEff_Inv = NULL;
  CMx_Inv = NULL;   CMy_Inv = NULL;   CMz_Inv = NULL;
  CFx_Inv = NULL;   CFy_Inv = NULL;   CFz_Inv = NULL;
  
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL; Surface_CSideForce_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL;   Surface_CFy_Inv = NULL;   Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL;   Surface_CMy_Inv = NULL;   Surface_CMz_Inv = NULL;
  
  Surface_CLift = NULL; Surface_CDrag = NULL; Surface_CSideForce = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL;   Surface_CFy = NULL;   Surface_CFz = NULL;
  Surface_CMx = NULL;   Surface_CMy = NULL;   Surface_CMz = NULL;
  
  ForceInviscid = NULL;  MomentInviscid = NULL;
  CPressure = NULL;      CPressureTarget = NULL; HeatFlux = NULL;
  HeatFluxTarget = NULL; YPlus = NULL;
  
  CMerit_Inv = NULL; CT_Inv = NULL; CQ_Inv = NULL;
  
  CEquivArea_Inv = NULL; CNearFieldOF_Inv = NULL;
  
  Inflow_MassFlow = NULL; Exhaust_MassFlow = NULL; Exhaust_Area = NULL;      Exhaust_Pressure = NULL;
  Inflow_Pressure = NULL; Inflow_Mach = NULL;      Inflow_Area = NULL;       Exhaust_Temperature = NULL;
  Bleed_MassFlow = NULL;  Bleed_Pressure = NULL;   Bleed_Temperature = NULL; Bleed_Area = NULL;
  
  iPoint_UndLapl = NULL;  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  Secondary = NULL; Secondary_i = NULL; Secondary_j = NULL;
  CharacPrimVar = NULL;
  Cauchy_Serie = NULL;
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp)
   Incompressible flow, primitive variables (P, vx, vy, vz, rho, beta, lamMu, EddyMu).
   FreeSurface Incompressible flow, primitive variables (P, vx, vy, vz, rho, beta, lamMu, EddyMu, LevelSet, Dist).
   ---*/
  
  nDim = geometry->GetnDim();
  
  if (incompressible) { nVar = nDim+1; nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  if (freesurface)    { nVar = nDim+2; nPrimVar = nDim+7; nPrimVarGrad = nDim+6; }
  if (compressible)   { nVar = nDim+2;
    nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
    nSecondaryVar = 2; nSecondaryVarGrad = 2;
  }
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);
  
  /*--- Allocate the node variables ---*/
  
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  
  Primitive   = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the Secondary solution ---*/
  
  Secondary   = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  
  if (roe_turkel) {
    LowMach_Precontioner = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new su2double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  }
  
  else {
    if (rank == MASTER_NODE) cout << "Explicit scheme. No Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      cvector[iVar] = new su2double [nDim];
    
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  CharacPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Force definition and coefficient arrays for all of the markers ---*/
  
  CPressure = new su2double* [nMarker];
  CPressureTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];
    CPressureTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Non-dimensional coefficients ---*/
  
  ForceInviscid     = new su2double[nDim];
  MomentInviscid    = new su2double[3];
  CDrag_Inv         = new su2double[nMarker];
  CLift_Inv         = new su2double[nMarker];
  CSideForce_Inv    = new su2double[nMarker];
  CMx_Inv           = new su2double[nMarker];
  CMy_Inv           = new su2double[nMarker];
  CMz_Inv           = new su2double[nMarker];
  CEff_Inv          = new su2double[nMarker];
  CFx_Inv           = new su2double[nMarker];
  CFy_Inv           = new su2double[nMarker];
  CFz_Inv           = new su2double[nMarker];
  
  Surface_CLift_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce_Inv = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CLift          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Rotorcraft coefficients ---*/
  
  CT_Inv           = new su2double[nMarker];
  CQ_Inv           = new su2double[nMarker];
  CMerit_Inv       = new su2double[nMarker];
  
  /*--- Supersonic coefficients ---*/
  
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];
  
  /*--- Engine simulation ---*/
  
  Inflow_MassFlow     = new su2double[nMarker];
  Inflow_Pressure     = new su2double[nMarker];
  Inflow_Mach         = new su2double[nMarker];
  Inflow_Area         = new su2double[nMarker];
  
  Exhaust_MassFlow    = new su2double[nMarker];
  Exhaust_Pressure    = new su2double[nMarker];
  Exhaust_Temperature = new su2double[nMarker];
  Exhaust_Area        = new su2double[nMarker];
  
  Bleed_MassFlow      = new su2double[nMarker];
  Bleed_Pressure      = new su2double[nMarker];
  Bleed_Temperature   = new su2double[nMarker];
  Bleed_Area          = new su2double[nMarker];
  
  /*--- Init total coefficients ---*/
  
  Total_CDrag   = 0.0;	Total_CLift        = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy          = 0.0;  Total_CMz          = 0.0;
  Total_CEff    = 0.0;	Total_CEquivArea   = 0.0;  Total_CNearFieldOF = 0.0;
  Total_CFx     = 0.0;	Total_CFy          = 0.0;  Total_CFz          = 0.0;
  Total_CT      = 0.0;	Total_CQ           = 0.0;  Total_CMerit       = 0.0;
  Total_MaxHeat = 0.0;  Total_Heat         = 0.0;
  Total_CpDiff  = 0.0;  Total_HeatFluxDiff = 0.0;
  
  /*--- Read farfield conditions ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();
  
  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  
  switch(direct_diff){
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
  
  
  /*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;
    
    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;
    
    Bleed_MassFlow[iMarker]      = 0.0;
    Bleed_Temperature[iMarker]   = Temperature_Inf;
    Bleed_Pressure[iMarker]      = Pressure_Inf;
    Bleed_Area[iMarker]          = 0.0;
  }
  
  /*--- Initializate quantities for the mixing process ---*/
  
  AveragedVelocity = new su2double* [nMarker];
  AveragedNormal = new su2double* [nMarker];
  AveragedGridVel = new su2double* [nMarker];
  AveragedFlux = new su2double* [nMarker];
  TotalFlux = new su2double* [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedVelocity[iMarker] = new su2double [nDim];
    AveragedNormal[iMarker] = new su2double [nDim];
    AveragedGridVel[iMarker] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      AveragedVelocity[iMarker][iDim] = 0.0;
      AveragedNormal[iMarker][iDim] = 0.0;
      AveragedGridVel [iMarker][iDim] = 0.0;
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedFlux[iMarker] = new su2double [nVar];
    TotalFlux[iMarker] = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      AveragedFlux[iMarker][iVar] = 0.0;
      TotalFlux[iMarker][iVar] = 0.0;
    }
  }
  
  AveragedNormalVelocity = new su2double[nMarker];
  AveragedTangVelocity = new su2double[nMarker];
  ExtAveragedNormalVelocity = new su2double[nMarker];
  ExtAveragedTangVelocity = new su2double[nMarker];
  MassFlow= new su2double[nMarker];
  FlowAngle= new su2double[nMarker];
  AveragedEnthalpy  = new su2double[nMarker];
  AveragedPressure  = new su2double[nMarker];
  AveragedTotPressure  = new su2double[nMarker];
  AveragedTotTemperature  = new su2double[nMarker];
  ExtAveragedTotPressure  = new su2double[nMarker];
  ExtAveragedTotTemperature  = new su2double[nMarker];
  AveragedDensity   = new su2double[nMarker];
  ExtAveragedPressure  = new su2double[nMarker];
  ExtAveragedDensity   = new su2double[nMarker];
  AveragedSoundSpeed= new su2double[nMarker];
  AveragedEntropy   = new su2double[nMarker];
  AveragedTangGridVelocity = new su2double[nMarker];
  AveragedMach = new su2double[nMarker];
  AveragedNormalMach = new su2double[nMarker];
  AveragedTangMach = new su2double[nMarker];
  
  
  /*--- Initializate quantities for turboperformace ---*/
  
  TotalStaticEfficiency = new su2double[nMarkerTurboPerf];
  TotalTotalEfficiency = new su2double[nMarkerTurboPerf];
  KineticEnergyLoss= new su2double[nMarkerTurboPerf];
  TotalPressureLoss= new su2double[nMarkerTurboPerf];
  MassFlowIn= new su2double[nMarkerTurboPerf];
  MassFlowOut= new su2double[nMarkerTurboPerf];
  FlowAngleIn= new su2double[nMarkerTurboPerf];
  FlowAngleOut= new su2double[nMarkerTurboPerf];
  EulerianWork= new su2double[nMarkerTurboPerf];
  TotalEnthalpyIn= new su2double[nMarkerTurboPerf];
  PressureRatio= new su2double[nMarkerTurboPerf];
  PressureOut= new su2double[nMarkerTurboPerf];
  EnthalpyOut= new su2double[nMarkerTurboPerf];
  MachIn= new su2double[nMarkerTurboPerf];
  MachOut= new su2double[nMarkerTurboPerf];
  NormalMachIn= new su2double[nMarkerTurboPerf];
  NormalMachOut= new su2double[nMarkerTurboPerf];
  VelocityOutIs= new su2double[nMarkerTurboPerf];
  
  for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
    TotalStaticEfficiency[iMarker]= 0.0;
    TotalTotalEfficiency[iMarker]= 0.0;
    KineticEnergyLoss[iMarker]= 0.0;
    TotalPressureLoss[iMarker]= 0.0;
    MassFlowIn[iMarker]= 0.0;
    MassFlowOut[iMarker]= 0.0;
    FlowAngleIn[iMarker]= 0.0;
    FlowAngleOut[iMarker]= 0.0;
    EulerianWork[iMarker]= 0.0;
    TotalEnthalpyIn[iMarker]= 0.0;
    PressureRatio[iMarker]= 0.0;
    PressureOut[iMarker]= 0.0;
    EnthalpyOut[iMarker]= 0.0;
    MachIn[iMarker]= 0.0;
    MachOut[iMarker]= 0.0;
    NormalMachIn[iMarker]= 0.0;
    NormalMachOut[iMarker]= 0.0;
    VelocityOutIs[iMarker]= 0.0;
  }
  
  
  /*--- Initialize the cauchy critera array for fixed CL mode ---*/
  
  if (config->GetFixed_CL_Mode())
    
    Cauchy_Serie = new su2double [config->GetCauchy_Elems()+1];
  
  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  
  if (!restart || (iMesh != MESH_0)) {
    
    /*--- Restart the solution from the free-stream state ---*/
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
    
  } else {
    
    /*--- Modify file name for an unsteady restart ---*/
    
    if (dual_time) {
      
      if (adjoint) { Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter()) - 1; }
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }
    if (nZone >1)
      filename= config->GetRestart_FlowFileName(filename, iZone);
    
    /*--- Open the restart file, throw an error if this fails. ---*/
    
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    
    long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    
    /*--- First, set all indices to a negative value by default ---*/
    
    for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
      Global2Local[iPoint] = -1;
    
    /*--- Now fill array with the transform values only for local points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    
    /*--- Read all lines in the restart file ---*/
    
    long iPoint_Local;
    unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
    
    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    while (getline (restart_file, text_line)) {
      istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
      
      iPoint_Local = Global2Local[iPoint_Global];
      
      /*--- Load the solution for this node. Note that the first entry
       on the restart file line is the global index, followed by the
       node coordinates, and then the conservative variables. ---*/
      
      if (iPoint_Local >= 0) {
        if (compressible) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        }
        if (incompressible) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        }
        if (freesurface) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        }
        node[iPoint_Local] = new CEulerVariable(Solution, nDim, nVar, config);
        iPoint_Global_Local++;
      }
      iPoint_Global++;
    }
    
    /*--- Detect a wrong solution file ---*/
    
    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
    
#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
        cout << "It could be empty lines at the end of the file." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CEulerVariable(Solution, nDim, nVar, config);
    
    /*--- Close the restart file ---*/
    
    restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    
    delete [] Global2Local;
    
  }
  
  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  
  if (compressible) {
    
    counter_local = 0;
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      Density = node[iPoint]->GetSolution(0);
      
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);
      
      StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;
      
      FluidModel->SetTDState_rhoe(Density, StaticEnergy);
      Pressure= FluidModel->GetPressure();
      Temperature= FluidModel->GetTemperature();
      
      /*--- Use the values at the infinity ---*/
      
      if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
        Solution[0] = Density_Inf;
        for (iDim = 0; iDim < nDim; iDim++)
          Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
        Solution[nDim+1] = Energy_Inf*Density_Inf;
        node[iPoint]->SetSolution(Solution);
        node[iPoint]->SetSolution_Old(Solution);
        counter_local++;
      }
      
    }
    
    /*--- Warning message about non-physical points ---*/
    
    if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
      SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
      counter_global = counter_local;
#endif
      if ((rank == MASTER_NODE) && (counter_global != 0))
        cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
    }
    
  }
  
  /*--- Define solver parameters needed for execution of destructor ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;
  
  /*--- Perform the MPI communication of the solution ---*/
//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);
  
}

CEulerSolver::~CEulerSolver(void) {
  unsigned short iVar, iMarker;
  
  /*--- Array deallocation ---*/
  if (CDrag_Inv != NULL)         delete [] CDrag_Inv;
  if (CLift_Inv != NULL)         delete [] CLift_Inv;
  if (CSideForce_Inv != NULL)    delete [] CSideForce_Inv;
  if (CMx_Inv != NULL)           delete [] CMx_Inv;
  if (CMy_Inv != NULL)           delete [] CMy_Inv;
  if (CMz_Inv != NULL)           delete [] CMz_Inv;
  if (CFx_Inv != NULL)           delete [] CFx_Inv;
  if (CFy_Inv != NULL)           delete [] CFy_Inv;
  if (CFz_Inv != NULL)           delete [] CFz_Inv;
  if (Surface_CLift_Inv != NULL) delete[] Surface_CLift_Inv;
  if (Surface_CDrag_Inv != NULL) delete[] Surface_CDrag_Inv;
  if (Surface_CSideForce_Inv != NULL) delete[] Surface_CSideForce_Inv;
  if (Surface_CEff_Inv != NULL) delete[] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;
  if (Surface_CLift != NULL)    delete [] Surface_CLift;
  if (Surface_CDrag != NULL)    delete [] Surface_CDrag;
  if (Surface_CSideForce != NULL) delete [] Surface_CSideForce;
  if (Surface_CEff != NULL) delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)          delete [] CEff_Inv;
  if (CMerit_Inv != NULL)        delete [] CMerit_Inv;
  if (CT_Inv != NULL)            delete [] CT_Inv;
  if (CQ_Inv != NULL)            delete [] CQ_Inv;
  if (CEquivArea_Inv != NULL)    delete [] CEquivArea_Inv;
  if (CNearFieldOF_Inv != NULL)  delete [] CNearFieldOF_Inv;
  if (ForceInviscid != NULL)     delete [] ForceInviscid;
  if (MomentInviscid != NULL)    delete [] MomentInviscid;
  if (Inflow_MassFlow != NULL)  delete [] Inflow_MassFlow;
  if (Exhaust_MassFlow != NULL)  delete [] Exhaust_MassFlow;
  if (Exhaust_Area != NULL)      delete [] Exhaust_Area;
  if (Inflow_Pressure != NULL)  delete [] Inflow_Pressure;
  if (Inflow_Mach != NULL)      delete [] Inflow_Mach;
  if (Inflow_Area != NULL)      delete [] Inflow_Area;
  if (Bleed_Pressure != NULL)  delete [] Bleed_Pressure;
  if (Bleed_Temperature != NULL)      delete [] Bleed_Temperature;
  if (Exhaust_Pressure != NULL)  delete [] Exhaust_Pressure;
  if (Exhaust_Temperature != NULL)      delete [] Exhaust_Temperature;
  if (Bleed_Area != NULL)      delete [] Bleed_Area;
  if (iPoint_UndLapl != NULL)       delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL)       delete [] jPoint_UndLapl;
  if (Primitive != NULL)        delete [] Primitive;
  if (Primitive_i != NULL)      delete [] Primitive_i;
  if (Primitive_j != NULL)      delete [] Primitive_j;
  //  if (Secondary != NULL)        delete [] Secondary;
  if (Secondary_i != NULL)      delete [] Secondary_i;
  if (Secondary_j != NULL)      delete [] Secondary_j;
  
  if (LowMach_Precontioner != NULL) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete LowMach_Precontioner[iVar];
    delete [] LowMach_Precontioner;
  }
  
  if (CPressure != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete CPressure[iMarker];
    delete [] CPressure;
  }
  
  if (CPressureTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete CPressureTarget[iMarker];
    delete [] CPressureTarget;
  }
  
  //  if (CharacPrimVar != NULL) {
  //    for (iMarker = 0; iMarker < nMarker; iMarker++) {
  //      for (iVertex = 0; iVertex < nVertex; iVertex++) {
  //        delete CharacPrimVar[iMarker][iVertex];
  //      }
  //    }
  //    delete [] CharacPrimVar;
  //  }
  
  if (HeatFlux != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete HeatFlux[iMarker];
    }
    delete [] HeatFlux;
  }
  
  if (HeatFluxTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete HeatFluxTarget[iMarker];
    }
    delete [] HeatFluxTarget;
  }
  
  if (YPlus != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete YPlus[iMarker];
    }
    delete [] YPlus;
  }
  
  if (Cauchy_Serie != NULL)
    delete [] Cauchy_Serie;
  
}

void CEulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CEulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new su2double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;
      
    }
    
  }
  
}

void CEulerSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR, *Buffer_Receive_Neighbor = NULL, *Buffer_Send_Neighbor = NULL;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      Buffer_Receive_Neighbor = new unsigned short [nBufferR_Vector];
      Buffer_Send_Neighbor = new unsigned short[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetLambda();
        Buffer_Send_Neighbor[iVertex] = geometry->node[iPoint]->GetnPoint();
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_Neighbor, nBufferS_Vector, MPI_UNSIGNED_SHORT, send_to, 1,
                        Buffer_Receive_Neighbor, nBufferR_Vector, MPI_UNSIGNED_SHORT, receive_from, 1, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
        Buffer_Receive_Neighbor[iVertex] = Buffer_Send_Neighbor[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      delete [] Buffer_Send_Neighbor;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetLambda(Buffer_Receive_Lambda[iVertex]);
        geometry->node[iPoint]->SetnNeighbor(Buffer_Receive_Neighbor[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      delete [] Buffer_Receive_Neighbor;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetSensor(Buffer_Receive_Lambda[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      
    }
    
  }
}

void CEulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Limiter[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}

void CEulerSolver::Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient_Primitive(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
  }
  
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CEulerSolver::Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nPrimVarGrad];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad;        nBufferR_Vector = nVertexR*nPrimVarGrad;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Primitive(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Limiter[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Limiter[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          node[iPoint]->SetLimiter_Primitive(iVar, Limiter[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
  
  delete [] Limiter;
  
}

//void CEulerSolver::Set_MPI_Secondary_Gradient(CGeometry *geometry, CConfig *config) {
//  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
//  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
//  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
//  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
//  int send_to, receive_from;
//
//  su2double **Gradient = new su2double* [nSecondaryVarGrad];
//  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//    Gradient[iVar] = new su2double[nDim];
//
//#ifdef HAVE_MPI
//  MPI_Status status;
//#endif
//
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//
//    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
//        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
//
//      MarkerS = iMarker;  MarkerR = iMarker+1;
//
//      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
//      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
//
//      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
//      nBufferS_Vector = nVertexS*nSecondaryVarGrad*nDim;        nBufferR_Vector = nVertexR*nSecondaryVarGrad*nDim;
//
//      /*--- Allocate Receive and send buffers  ---*/
//      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
//      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
//
//      /*--- Copy the solution old that should be sended ---*/
//      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
//        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Buffer_Send_Gradient[iDim*nSecondaryVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Secondary(iVar, iDim);
//      }
//
//#ifdef HAVE_MPI
//
//      /*--- Send/Receive information using Sendrecv ---*/
//      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
//                   Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
//
//#else
//
//      /*--- Receive information without MPI ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Buffer_Receive_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//      }
//
//#endif
//
//      /*--- Deallocate send buffer ---*/
//      delete [] Buffer_Send_Gradient;
//
//      /*--- Do the coordinate transformation ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//
//        /*--- Find point and its type of transformation ---*/
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
//
//        /*--- Retrieve the supplied periodic information. ---*/
//        angles = config->GetPeriodicRotation(iPeriodic_Index);
//
//        /*--- Store angles separately for clarity. ---*/
//        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
//        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
//        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
//
//        /*--- Compute the rotation matrix. Note that the implicit
//         ordering is rotation about the x-axis, y-axis,
//         then z-axis. Note that this is the transpose of the matrix
//         used during the preprocessing stage. ---*/
//        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
//        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
//        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
//
//        /*--- Copy conserved variables before performing transformation. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//
//        /*--- Need to rotate the gradients for all conserved variables. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//          if (nDim == 2) {
//            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//          }
//          else {
//            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
//          }
//        }
//
//        /*--- Store the received information ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            node[iPoint]->SetGradient_Secondary(iVar, iDim, Gradient[iVar][iDim]);
//
//      }
//
//      /*--- Deallocate receive buffer ---*/
//      delete [] Buffer_Receive_Gradient;
//
//    }
//
//  }
//
//  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//    delete [] Gradient[iVar];
//  delete [] Gradient;
//
//}

//void CEulerSolver::Set_MPI_Secondary_Limiter(CGeometry *geometry, CConfig *config) {
//  unsigned short iVar, iMarker, MarkerS, MarkerR;
//  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
//  su2double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
//  int send_to, receive_from;
//
//  su2double *Limiter = new su2double [nSecondaryVarGrad];
//
//#ifdef HAVE_MPI
//  MPI_Status status;
//#endif
//
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//
//    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
//        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
//
//      MarkerS = iMarker;  MarkerR = iMarker+1;
//
//      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
//      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
//
//      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
//      nBufferS_Vector = nVertexS*nSecondaryVarGrad;        nBufferR_Vector = nVertexR*nSecondaryVarGrad;
//
//      /*--- Allocate Receive and send buffers  ---*/
//      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
//      Buffer_Send_Limit = new su2double[nBufferS_Vector];
//
//      /*--- Copy the solution old that should be sended ---*/
//      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
//        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Secondary(iVar);
//      }
//
//#ifdef HAVE_MPI
//
//      /*--- Send/Receive information using Sendrecv ---*/
//      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
//                   Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
//
//#else
//
//      /*--- Receive information without MPI ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
//      }
//
//#endif
//
//      /*--- Deallocate send buffer ---*/
//      delete [] Buffer_Send_Limit;
//
//      /*--- Do the coordinate transformation ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//
//        /*--- Find point and its type of transformation ---*/
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//
//        /*--- Copy conserved variables before performing transformation. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
//
//        /*--- Copy transformed conserved variables back into buffer. ---*/
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          node[iPoint]->SetLimiter_Secondary(iVar, Limiter[iVar]);
//
//      }
//
//      /*--- Deallocate receive buffer ---*/
//      delete [] Buffer_Receive_Limit;
//
//    }
//
//  }
//
//  delete [] Limiter;
//
//}

void CEulerSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {
  
  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
  Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0, Velocity_Reynolds = 0.0,
  Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
  Temperature_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref= 0.0,
  Froude = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0, TgammaR = 0.0;
  
  unsigned short iDim;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Local variables ---*/
  
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool compressible       = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface        = (config->GetKind_Regime() == FREESURFACE);
  bool unsteady           = (config->GetUnsteady_Simulation() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool standard_air       = (config->GetKind_FluidModel() == STANDARD_AIR);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);
  bool aeroelastic        = config->GetAeroelastic_Simulation();
  
  /*--- Set temperature via the flutter speed index ---*/
  if (aeroelastic) {
    su2double vf             = config->GetAeroelastic_Flutter_Speed_Index();
    su2double w_alpha        = config->GetAeroelastic_Frequency_Pitch();
    su2double b              = config->GetLength_Reynolds()/2.0; // airfoil semichord, Reynolds length is by defaul 1.0
    su2double mu             = config->GetAeroelastic_Airfoil_Mass_Ratio();
    // The temperature times gamma times the gas constant. Depending on the FluidModel temp is calculated below.
    TgammaR = ((vf*vf)*(b*b)*(w_alpha*w_alpha)*mu) / (Mach*Mach);
  }
  
  /*--- Compressible non dimensionalization ---*/
  
  if (compressible) {
    
    /*--- Compute the Free Stream velocity, using the Mach number ---*/
    
    Pressure_FreeStream = config->GetPressure_FreeStream();
    Density_FreeStream  = config->GetDensity_FreeStream();
    Temperature_FreeStream  = config->GetTemperature_FreeStream();
    
    switch (config->GetKind_FluidModel()) {
        
      case STANDARD_AIR:
        
        if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
        else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);
        
        FluidModel = new CIdealGas(1.4, config->GetGas_Constant());
        if (free_stream_temp) {
          if (aeroelastic) {
            Temperature_FreeStream = TgammaR / (config->GetGas_Constant()*1.4);
            config->SetTemperature_FreeStream(Temperature_FreeStream);
          }
          FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
          Density_FreeStream = FluidModel->GetDensity();
          config->SetDensity_FreeStream(Density_FreeStream);
        }
        else {
          FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
          Temperature_FreeStream = FluidModel->GetTemperature();
          config->SetTemperature_FreeStream(Temperature_FreeStream);
        }
        break;
        
      case IDEAL_GAS:
        
        FluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
        if (free_stream_temp) {
          FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
          Density_FreeStream = FluidModel->GetDensity();
          config->SetDensity_FreeStream(Density_FreeStream);
        }
        else {
          FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
          Temperature_FreeStream = FluidModel->GetTemperature();
          config->SetTemperature_FreeStream(Temperature_FreeStream);
        }
        break;
        
      case VW_GAS:
        
        FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                         config->GetPressure_Critical(), config->GetTemperature_Critical());
        if (free_stream_temp) {
          FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
          Density_FreeStream = FluidModel->GetDensity();
          config->SetDensity_FreeStream(Density_FreeStream);
        }
        else {
          FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
          Temperature_FreeStream = FluidModel->GetTemperature();
          config->SetTemperature_FreeStream(Temperature_FreeStream);
        }
        break;
        
      case PR_GAS:
        
        FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                       config->GetTemperature_Critical(), config->GetAcentric_Factor());
        if (free_stream_temp) {
          FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
          Density_FreeStream = FluidModel->GetDensity();
          config->SetDensity_FreeStream(Density_FreeStream);
        }
        else {
          FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
          Temperature_FreeStream = FluidModel->GetTemperature();
          config->SetTemperature_FreeStream(Temperature_FreeStream);
        }
        break;
        
    }
    
    Mach2Vel_FreeStream = FluidModel->GetSoundSpeed();
    
    /*--- Compute the Free Stream velocity, using the Mach number ---*/
    
    if (nDim == 2) {
      config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
      config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
    }
    if (nDim == 3) {
      config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
      config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
      config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    }
    
    /*--- Compute the modulus of the free stream velocity ---*/
    
    ModVel_FreeStream = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
    ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);
    
    /*--- Viscous initialization ---*/
    
    if (viscous) {
      
      /*--- Reynolds based initialization ---*/
      
      if (reynolds_init) {
        
        /*--- First, check if there is mesh motion. If yes, use the Mach
         number relative to the body to initialize the flow. ---*/
        
        if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
        else Velocity_Reynolds = ModVel_FreeStream;
        
        /*--- Change of measurement system, hard coded value working only with STANDAR AIR model ---*/
        
        if (standard_air) {
          if (config->GetSystemMeasurements() == SI) {
            config->SetMu_RefND(1.716E-5);
            config->SetMu_SND(110.4);
            config->SetMu_Temperature_RefND(273.15);
          }
          if (config->GetSystemMeasurements() == US) {
            config->SetMu_RefND(3.62E-7);
            config->SetMu_SND(198.72);
            config->SetMu_Temperature_RefND(518.7);
          }
        }
        
        /*--- For viscous flows, pressure will be computed from a density
         that is found from the Reynolds number. The viscosity is computed
         from the dimensional version of Sutherland's law ---*/
        
        FluidModel->SetLaminarViscosityModel(config);
        
        Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
        config->SetViscosity_FreeStream(Viscosity_FreeStream);
        
        Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
        config->SetDensity_FreeStream(Density_FreeStream);
        FluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
        Pressure_FreeStream = FluidModel->GetPressure();
        config->SetPressure_FreeStream(Pressure_FreeStream);
        Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
        
      }
      
      /*--- Thermodynamics quantities based initialization ---*/
      
      else {
        
        FluidModel->SetLaminarViscosityModel(config);
        Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
        config->SetViscosity_FreeStream(Viscosity_FreeStream);
        Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
        
      }
      
      /*--- Turbulence kinetic energy ---*/
      
      Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
      
    }
    else {
      
      /*--- For inviscid flow, energy is calculated from the specified
       FreeStream quantities using the proper gas law. ---*/
      
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
      
    }
    
    /*-- Compute the freestream energy. ---*/
    
    if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);
    
    /*--- Compute non dimensional quantities. By definition,
     Lref is one because we have converted the grid to meters. ---*/
    
    if (config->GetRef_NonDim() == DIMENSIONAL) {
      Pressure_Ref      = 1.0;
      Density_Ref       = 1.0;
      Temperature_Ref   = 1.0;
    }
    else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
      Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
      Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
      Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
    }
    else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
      Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
      Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
      Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
    }
    else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
      Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
      Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
      Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
    }
    config->SetPressure_Ref(Pressure_Ref);
    config->SetDensity_Ref(Density_Ref);
    config->SetTemperature_Ref(Temperature_Ref);
    
    Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
    Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
    Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
    Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
    Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;                        config->SetForce_Ref(Force_Ref);
    Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
    Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
    Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
    Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);         config->SetFroude(Froude);
    
  }
  
  /*--- Incompressible non dimensionalization ---*/
  
  else {
    
    /*--- Reference length = 1 (by default)
     Reference density   = liquid density or freestream
     Reference viscosity = liquid viscosity or freestream
     Reference velocity  = liquid velocity or freestream
     Reference pressure  = Reference density * Reference velocity * Reference velocity
     Reynolds number based on the liquid or reference viscosity ---*/
    
    Pressure_FreeStream = 0.0; config->SetPressure_FreeStream(Pressure_FreeStream);
    Density_FreeStream  = config->GetDensity_FreeStream();
    ModVel_FreeStream   = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
    ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);
    
    /*--- Additional reference values defined by Pref, Tref, Rho_ref. By definition,
     Lref is one because we have converted the grid to meters.---*/
    
    Length_Ref   = config->GetLength_Reynolds();              config->SetLength_Ref(Length_Ref);
    Density_Ref  = Density_FreeStream;                        config->SetDensity_Ref(Density_Ref);
    Velocity_Ref = ModVel_FreeStream;                         config->SetVelocity_Ref(Velocity_Ref);
    Pressure_Ref = Density_Ref*(Velocity_Ref*Velocity_Ref);   config->SetPressure_Ref(Pressure_Ref);
    
    if (viscous) {
      Viscosity_FreeStream = config->GetViscosity_FreeStream();
      Reynolds = Density_Ref*Velocity_Ref*Length_Ref / Viscosity_FreeStream; config->SetReynolds(Reynolds);
      Viscosity_Ref = Viscosity_FreeStream * Reynolds;                       config->SetViscosity_Ref(Viscosity_Ref);
    }
    
    /*--- Compute Mach number ---*/
    
    Mach = ModVel_FreeStream / sqrt(config->GetBulk_Modulus()/Density_FreeStream);   config->SetMach(Mach);
    
    /*--- Compute Alpha angle ---*/
    
    if (nDim == 2) Alpha = atan(config->GetVelocity_FreeStream()[1]/config->GetVelocity_FreeStream()[0])*180.0/PI_NUMBER;
    else Alpha = atan(config->GetVelocity_FreeStream()[2]/config->GetVelocity_FreeStream()[0])*180.0/PI_NUMBER;
    config->SetAoA(Alpha);
    
    /*--- Compute Beta angle ---*/
    
    if (nDim == 2) Beta = 0.0;
    else Beta = asin(config->GetVelocity_FreeStream()[1]/ModVel_FreeStream)*180.0/PI_NUMBER;
    config->SetAoS(Beta);
    
    /*--- Compute Froude ---*/
    
    Froude = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);   config->SetFroude(Froude);
    Time_Ref = Length_Ref/Velocity_Ref;                             config->SetTime_Ref(Time_Ref);
    
  }
  
  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  
  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }
  
  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  
  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);
  
  
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
  
  switch (config->GetKind_FluidModel()) {
      
    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND);
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case VW_GAS:
      FluidModel = new CVanDerWaalsGas(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                       config->GetTemperature_Critical()/config->GetTemperature_Ref());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
    case PR_GAS:
      FluidModel = new CPengRobinson(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                     config->GetTemperature_Critical()/config->GetTemperature_Ref(), config->GetAcentric_Factor());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;
      
  }
  
  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
  
  if (viscous) {
    
    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_ConstantND()/Viscosity_Ref);
    
    /*--- Sutherland's model ---*/
    
    config->SetMu_RefND(config->GetMu_RefND()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_SND()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_RefND()/config->GetTemperature_Ref());
    
    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_ConstantND()/Conductivity_Ref);
    
    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);
    
  }
  
  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);
  
  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);
  
  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);
  
  /*--- Write output to the console if this is the master node and first domain ---*/
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {
    
    cout.precision(6);
    
    if (compressible) {
      if (viscous) {
        cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
        cout << "based on the free-stream temperature and a density computed" << endl;
        cout << "from the Reynolds number." << endl;
      } else {
        cout << "Inviscid flow: Computing density based on free-stream" << endl;
        cout << "temperature and pressure using the ideal gas law." << endl;
      }
    }
    
    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;
    
    if (incompressible || freesurface) {
      cout << "Viscous and Inviscid flow: rho_ref, and vel_ref" << endl;
      cout << "are based on the free-stream values, p_ref = rho_ref*vel_ref^2." << endl;
      cout << "The free-stream value of the pressure is 0." << endl;
      cout << "Mach number: "<< config->GetMach() << ", computed using the Bulk modulus." << endl;
      cout << "Angle of attack (deg): "<< config->GetAoA() << ", computed using the the free-stream velocity." << endl;
      cout << "Side slip angle (deg): "<< config->GetAoS() << ", computed using the the free-stream velocity." << endl;
      if (viscous) cout << "Reynolds number: " << config->GetReynolds() << ", computed using free-stream values."<< endl;
      cout << "Only dimensional computation, the grid should be dimensional." << endl;
    }
    
    cout <<"-- Input conditions:"<< endl;
    
    if (compressible) {
      switch (config->GetKind_FluidModel()) {
          
        case STANDARD_AIR:
          cout << "Fluid Model: STANDARD_AIR "<< endl;
          cout << "Specific gas constant: " << config->GetGas_Constant();
          if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
          cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
          cout << "Specific Heat Ratio: "<< Gamma << endl;
          break;
          
        case IDEAL_GAS:
          cout << "Fluid Model: IDEAL_GAS "<< endl;
          cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
          cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
          cout << "Specific Heat Ratio: "<< Gamma << endl;
          break;
          
        case VW_GAS:
          cout << "Fluid Model: Van der Waals "<< endl;
          cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
          cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
          cout << "Specific Heat Ratio: "<< Gamma << endl;
          cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
          cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
          cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
          cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
          break;
          
        case PR_GAS:
          cout << "Fluid Model: Peng-Robinson "<< endl;
          cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
          cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< endl;
          cout << "Specific Heat Ratio: "<< Gamma << endl;
          cout << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << endl;
          cout << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << endl;
          cout << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << endl;
          cout << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << endl;
          break;
          
      }
      if (viscous) {
        switch (config->GetKind_ViscosityModel()) {
            
          case CONSTANT_VISCOSITY:
            cout << "Viscosity Model: CONSTANT_VISCOSITY  "<< endl;
            cout << "Laminar Viscosity: " << config->GetMu_ConstantND()*Viscosity_Ref;
            if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
            else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
            cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
            break;
            
          case SUTHERLAND:
            cout << "Viscosity Model: SUTHERLAND "<< endl;
            cout << "Ref. Laminar Viscosity: " << config->GetMu_RefND()*Viscosity_Ref;
            if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
            else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
            cout << "Ref. Temperature: " << config->GetMu_Temperature_RefND()*config->GetTemperature_Ref();
            if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
            else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
            cout << "Sutherland Constant: "<< config->GetMu_SND()*config->GetTemperature_Ref();
            if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
            else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
            cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
            cout << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< endl;
            cout << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< endl;
            break;
            
        }
        switch (config->GetKind_ConductivityModel()) {
            
          case CONSTANT_PRANDTL:
            cout << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
            cout << "Prandtl: " << config->GetPrandtl_Lam()<< endl;
            break;
            
          case CONSTANT_CONDUCTIVITY:
            cout << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
            cout << "Molecular Conductivity: " << config->GetKt_ConstantND()*Conductivity_Ref<< " W/m^2.K." << endl;
            cout << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< endl;
            break;
            
        }
      }
    }
    
    if (incompressible || freesurface) {
      cout << "Bulk modulus: " << config->GetBulk_Modulus();
      if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
      cout << "Artificial compressibility factor: " << config->GetArtComp_Factor();
      if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    }
    
    cout << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+Mach*Mach*0.5*(Gamma-1.0), Gamma/(Gamma-1.0) );
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    if (compressible) {
      cout << "Free-stream temperature: " << config->GetTemperature_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    }
    
    cout << "Free-stream density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      cout << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) cout << " m/s. ";
    else if (config->GetSystemMeasurements() == US) cout << " ft/s. ";
    
    cout << "Magnitude: "	<< config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;
    
    if (compressible) {
      cout << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    }
    
    if (viscous) {
      cout << "Free-stream viscosity: " << config->GetViscosity_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
        cout << "Free-stream specific dissipation: " << config->GetOmega_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " 1/s." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " 1/s." << endl;
      }
    }
    
    if (unsteady) { cout << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << endl; }
    
    /*--- Print out reference values. ---*/
    
    cout <<"-- Reference values:"<< endl;
    
    if (compressible) {
      cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
    }
    
    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    if (compressible) {
      cout << "Reference temperature: " << config->GetTemperature_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    }
    
    cout << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    cout << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;
    
    if (compressible) {
      cout << "Reference energy per unit mass: " << config->GetEnergy_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
    }
    
    if (incompressible || freesurface) {
      cout << "Reference length: " << config->GetLength_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " in." << endl;
    }
    
    if (viscous) {
      cout << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      if (compressible){
        cout << "Reference conductivity: " << config->GetConductivity_Ref();
        if (config->GetSystemMeasurements() == SI) cout << " W/m^2.K." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " lbf/ft.s.R." << endl;
      }
    }
    
    
    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;
    
    /*--- Print out resulting non-dim values here. ---*/
    
    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;
    if (viscous) {
      cout << "Reynolds number (non-dim): " << config->GetReynolds() <<". Re length: " << config->GetLength_Reynolds();
      if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " ft." << endl;
    }
    if (gravity) {
      cout << "Froude number (non-dim): " << Froude << endl;
      cout << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*Froude*Froude << endl;
    }
    
    if (compressible) {
      cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
      cout << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;
    }
    
    cout << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;
    
    cout << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << endl;
    
    if (nDim == 2) {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "	 << config->GetModVel_FreeStreamND() << endl;
    
    if (compressible)
      cout << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << endl;
    
    if (viscous) {
      cout << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << endl;
      if (turbulent) {
        cout << "Free-stream turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << endl;
        cout << "Free-stream specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << endl;
      }
    }
    
    if (unsteady) {
      cout << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << endl;
      cout << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << endl;
    }
    
    cout << endl;
    
  }
  
}

void CEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar, iDim;
  su2double Density, Pressure, yFreeSurface, PressFreeSurface, Froude, yCoord, Velx, Vely, Velz, RhoVelx, RhoVely, RhoVelz, YCoord = 0.0,
  ZCoord = 0.0, DensityInc, ViscosityInc, Heaviside, LevelSet, lambda, Area_Children, Area_Parent, LevelSet_Fine, epsilon,
  *Solution_Fine, *Solution, PressRef, yCoordRef;
  
  unsigned short nDim = geometry[MESH_0]->GetnDim();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool rans = ((config->GetKind_Solver() == RANS) ||
               (config->GetKind_Solver() == ADJ_RANS) ||
               (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool gravity = (config->GetGravityForce() == YES);
  bool engine_intake = config->GetEngine_Intake();
  
  
  /*--- Set the location and value of the free-surface ---*/
  
  if (freesurface) {
    
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        
        /*--- Set initial boundary condition at iter 0 ---*/
        
        if ((ExtIter == 0) && (!restart)) {
          
          /*--- Compute the level set value in all the MG levels (basic case, distance to
           the Y/Z plane, and interpolate the solution to the coarse levels ---*/
          
          if (iMesh == MESH_0) {
            YCoord = geometry[iMesh]->node[iPoint]->GetCoord(1);
            if (nDim == 2) LevelSet = YCoord - config->GetFreeSurface_Zero();
            else {
              ZCoord = geometry[iMesh]->node[iPoint]->GetCoord(2);
              LevelSet = ZCoord - config->GetFreeSurface_Zero();
            }
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(nDim+1, LevelSet);
          }
          else {
            Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
            LevelSet = 0.0;
            for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
              Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
              Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
              LevelSet_Fine = solver_container[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution(nDim+1);
              LevelSet += LevelSet_Fine*Area_Children/Area_Parent;
            }
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(nDim+1, LevelSet);
          }
          
          /*--- Compute the flow solution using the level set value. ---*/
          
          epsilon = config->GetFreeSurface_Thickness();
          Heaviside = 0.0;
          if (LevelSet < -epsilon) Heaviside = 1.0;
          if (fabs(LevelSet) <= epsilon) Heaviside = 1.0 - (0.5*(1.0+(LevelSet/epsilon)+(1.0/PI_NUMBER)*sin(PI_NUMBER*LevelSet/epsilon)));
          if (LevelSet > epsilon) Heaviside = 0.0;
          
          /*--- Set the value of the incompressible density for free surface flows (density ratio g/l) ---*/
          
          lambda = config->GetRatioDensity();
          DensityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetDensity_FreeStreamND();
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetDensityInc(DensityInc);
          
          /*--- Set the value of the incompressible viscosity for free surface flows (viscosity ratio g/l) ---*/
          
          lambda = config->GetRatioViscosity();
          ViscosityInc = (lambda + (1.0 - lambda)*Heaviside)*config->GetViscosity_FreeStreamND();
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetLaminarViscosityInc(ViscosityInc);
          
          /*--- Update solution with the new pressure ---*/
          
          yFreeSurface = config->GetFreeSurface_Zero();
          PressFreeSurface = solver_container[iMesh][FLOW_SOL]->GetPressure_Inf();
          Density = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetDensityInc();
          Froude = config->GetFroude();
          yCoord = geometry[iMesh]->node[iPoint]->GetCoord(nDim-1);
          Pressure = PressFreeSurface + Density*((yFreeSurface-yCoord)/(Froude*Froude));
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Pressure);
          
          /*--- Update solution with the new velocity ---*/
          
          Velx = solver_container[iMesh][FLOW_SOL]->GetVelocity_Inf(0);
          Vely = solver_container[iMesh][FLOW_SOL]->GetVelocity_Inf(1);
          RhoVelx = Velx * Density; RhoVely = Vely * Density;
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(1, RhoVelx);
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(2, RhoVely);
          if (nDim == 3) {
            Velz = solver_container[iMesh][FLOW_SOL]->GetVelocity_Inf(2);
            RhoVelz = Velz * Density;
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(3, RhoVelz);
          }
          
        }
        
      }
      
      /*--- Set the MPI communication ---*/
      
      solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
      solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution_Old(geometry[iMesh], config);
      
    }
    
  }
  
  /*--- Set the pressure value in simulations with gravity ---*/
  
  
  if (incompressible && gravity ) {
    
    
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        
        /*--- Set initial boundary condition at iter 0 ---*/
        
        if ((ExtIter == 0) && (!restart)) {
          
          /*--- Update solution with the new pressure ---*/
          
          PressRef = solver_container[iMesh][FLOW_SOL]->GetPressure_Inf();
          Density = solver_container[iMesh][FLOW_SOL]->GetDensity_Inf();
          yCoordRef = 0.0;
          yCoord = geometry[iMesh]->node[iPoint]->GetCoord(nDim-1);
          Pressure = PressRef + Density*((yCoordRef-yCoord)/(config->GetFroude()*config->GetFroude()));
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Pressure);
          
        }
        
      }
      
      /*--- Set the MPI communication ---*/
      
      solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
      solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution_Old(geometry[iMesh], config);
      
    }
    
  }
  
  
  /*--- Set subsonic initial condition for engine intakes ---*/
  
  if (engine_intake) {
    
    /*--- Set initial boundary condition at iteration 0 ---*/
    
    if ((ExtIter == 0) && (!restart)) {
      
      su2double Velocity_Box[3] = {0.0, 0.0, 0.0}, Velocity_BoxND[3] = {0.0, 0.0, 0.0}, Viscosity_Box,
      Density_Box, Density_BoxND, Pressure_Box, Pressure_BoxND, ModVel_Box, ModVel_BoxND, Energy_BoxND,
      T_ref = 0.0, S = 0.0, Mu_ref = 0.0, *Coord, MinCoordValues[3],
      MaxCoordValues[3], *Subsonic_Engine_Box;
      
      su2double Mach = 0.40;
      su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
      su2double Beta  = config->GetAoS()*PI_NUMBER/180.0;
      
      su2double Gamma_Minus_One = Gamma - 1.0;
      su2double Gas_Constant = config->GetGas_Constant();
      
      su2double Temperature_Box = config->GetTemperature_FreeStream();
      su2double Mach2Vel_Box = sqrt(Gamma*Gas_Constant*Temperature_Box);
      
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          
          Velocity_Box[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_Box;
          Velocity_Box[1] = sin(Beta)*Mach*Mach2Vel_Box;
          Velocity_Box[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_Box;
          
          ModVel_Box = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            ModVel_Box += Velocity_Box[iDim]*Velocity_Box[iDim];
          }
          ModVel_Box = sqrt(ModVel_Box);
          
          if (config->GetViscous()) {
            if (config->GetSystemMeasurements() == SI) { T_ref = 273.15; S = 110.4; Mu_ref = 1.716E-5; }
            if (config->GetSystemMeasurements() == US) { T_ref = 518.7; S = 198.72; Mu_ref = 3.62E-7; }
            Viscosity_Box = Mu_ref*(pow(Temperature_Box/T_ref, 1.5) * (T_ref+S)/(Temperature_Box+S));
            Density_Box   = config->GetReynolds()*Viscosity_Box/(ModVel_Box*config->GetLength_Reynolds());
            Pressure_Box  = Density_Box*Gas_Constant*Temperature_Box;
          }
          else {
            Pressure_Box = config->GetPressure_FreeStream();
            Density_Box = Pressure_Box/(Gas_Constant*Temperature_Box);
          }
          
          Density_BoxND  = Density_Box/config->GetDensity_Ref();
          Pressure_BoxND = Pressure_Box/config->GetPressure_Ref();
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_BoxND[iDim] = Velocity_Box[iDim]/config->GetVelocity_Ref();
          }
          
          ModVel_BoxND = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            ModVel_BoxND += Velocity_BoxND[iDim]*Velocity_BoxND[iDim];
          }
          ModVel_BoxND = sqrt(ModVel_BoxND);
          
          Energy_BoxND = Pressure_BoxND/(Density_BoxND*Gamma_Minus_One)+0.5*ModVel_BoxND*ModVel_BoxND;
          
          Coord = geometry[iMesh]->node[iPoint]->GetCoord();
          
          Subsonic_Engine_Box = config->GetSubsonic_Engine_Box();
          
          MinCoordValues[0] = Subsonic_Engine_Box[0]; MinCoordValues[1] = Subsonic_Engine_Box[1]; MinCoordValues[2] = Subsonic_Engine_Box[2];
          MaxCoordValues[0] = Subsonic_Engine_Box[3]; MaxCoordValues[1] = Subsonic_Engine_Box[4]; MaxCoordValues[2] = Subsonic_Engine_Box[5];
          
          if (((Coord[0] >= MinCoordValues[0]) && (Coord[0] <= MaxCoordValues[0])) &&
              ((Coord[1] >= MinCoordValues[1]) && (Coord[1] <= MaxCoordValues[1])) &&
              ((Coord[2] >= MinCoordValues[2]) && (Coord[2] <= MaxCoordValues[2]))) {
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Density_BoxND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(iDim+1, Density_BoxND*Velocity_BoxND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(nVar-1, Density_BoxND*Energy_BoxND);
            
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(0, Density_BoxND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(iDim+1, Density_BoxND*Velocity_BoxND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(nVar-1, Density_BoxND*Energy_BoxND);
            
          }
          
        }
        
        /*--- Set the MPI communication ---*/
        
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
        solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution_Old(geometry[iMesh], config);
        
      }
      
    }
    
  }
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  
  if (restart && (ExtIter == 0)) {
    
    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
      }
      solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
    delete [] Solution;
    
    /*--- Interpolate the turblence variable also, if needed ---*/
    
    if (rans) {
      
      unsigned short nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
      Solution = new su2double[nVar_Turb];
      for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
          for (iVar = 0; iVar < nVar_Turb; iVar++) Solution[iVar] = 0.0;
          for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
            Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
            Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
            Solution_Fine = solver_container[iMesh-1][TURB_SOL]->node[Point_Fine]->GetSolution();
            for (iVar = 0; iVar < nVar_Turb; iVar++) {
              Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
            }
          }
          solver_container[iMesh][TURB_SOL]->node[iPoint]->SetSolution(Solution);
        }
        solver_container[iMesh][TURB_SOL]->Set_MPI_Solution(geometry[iMesh], config);
        solver_container[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver_container[iMesh], config, iMesh);
      }
      delete [] Solution;
    }
    
  }
  
  /*--- The value of the solution for the first iteration of the dual time ---*/
  
  if (dual_time && (ExtIter == 0 || (restart && (long)ExtIter == config->GetUnst_RestartIter()))) {
    
    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/
    
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
        if (rans) {
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
        }
      }
    }
    
    if ((restart && (long)ExtIter == config->GetUnst_RestartIter()) &&
        (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
      /*--- Load an additional restart file for a 2nd-order restart ---*/
      
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1));
      
      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1));
      
      /*--- Push back this new solution to time level N. ---*/
      
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          if (rans) {
            solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }
  }
}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long ErrorCounter = 0;
  
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned long ExtIter = config->GetExtIter();
  bool adjoint          = config->GetAdjoint();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool low_fidelity     = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool second_order     = ((config->GetSpatialOrder_Flow() == SECOND_ORDER) || (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
  bool limiter          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (!low_fidelity) && (ExtIter <= config->GetLimiterIter()));
  bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
  bool freesurface      = (config->GetKind_Regime() == FREESURFACE);
  bool engine           = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineBleed() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk    = ((config->GetnMarker_ActDisk_Inlet() != 0) || (config->GetnMarker_ActDisk_Outlet() != 0));
  bool fixed_cl         = config->GetFixed_CL_Mode();
  
  /*--- Compute the engine properties ---*/
  
  if (engine) { GetEngine_Properties(geometry, config, iMesh, Output); }
  
  /*--- Compute the actuator disk properties ---*/
  
  if (actuator_disk) { GetActuatorDisk_Properties(geometry, config, iMesh, Output); }
  
  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/
  
  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }
  
  /*--- Compute distance function to zero level set (Set LevelSet and Distance primitive variables)---*/
  
  if (freesurface) { SetFreeSurface_Distance(geometry, config); }
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Upwind second order reconstruction ---*/
  
  if ((second_order && !center) && ((iMesh == MESH_0) || low_fidelity) && !Output) {
    
    /*--- Gradient computation ---*/
    
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
      //    	if (compressible && !ideal_gas) SetSecondary_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
      //    	if (compressible && !ideal_gas) SetSecondary_Gradient_LS(geometry, config);
    }
    
    
    /*--- Limiter computation ---*/
    
    if ((limiter) && (iMesh == MESH_0) && !Output) {
      SetPrimitive_Limiter(geometry, config);
      //    	if (compressible && !ideal_gas) SetSecondary_Limiter(geometry, config);
    }
    
  }
  
  /*--- Artificial dissipation ---*/
  
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && ((iMesh == MESH_0) || low_fidelity)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !config->GetDiscrete_Adjoint()) Jacobian.SetValZero();
  
  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) { }

unsigned long CEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;
  
  bool compressible         = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible       = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface          = (config->GetKind_Regime() == FREESURFACE);
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Incompressible flow, primitive variables nDim+3, (P, vx, vy, vz, rho, beta),
     FreeSurface Incompressible flow, primitive variables nDim+4, (P, vx, vy, vz, rho, beta, dist),
     Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
    if (compressible) {
      RightSol = node[iPoint]->SetPrimVar_Compressible(FluidModel);
      node[iPoint]->SetSecondaryVar_Compressible(FluidModel);
    }
    
    if (incompressible) {
      RightSol = node[iPoint]->SetPrimVar_Incompressible(Density_Inf, config);
    }
    
    
    if (freesurface){
      RightSol = node[iPoint]->SetPrimVar_FreeSurface(config);
    }
    
    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
}
void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {
  
  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Mean_BetaInc2, Lambda, Local_Delta_Time, Mean_DensityInc, Mean_LevelSet,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j, Delta = 0.0, a, b, c, e, f;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  su2double epsilon = config->GetFreeSurface_Thickness();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetMax_Lambda_Inv(0.0);
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    if (compressible) {
      Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
      Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    }
    if (incompressible) {
      Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
      Mean_BetaInc2 = 0.5 * (node[iPoint]->GetBetaInc2() + node[jPoint]->GetBetaInc2());
      Mean_DensityInc = 0.5 * (node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
      Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
    }
    if (freesurface) {
      Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
      Mean_BetaInc2 = 0.5 * (node[iPoint]->GetBetaInc2() + node[jPoint]->GetBetaInc2());
      Mean_DensityInc = 0.5 * (node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
      Mean_LevelSet = 0.5 * (node[iPoint]->GetLevelSet() + node[jPoint]->GetLevelSet());
      
      if (Mean_LevelSet < -epsilon) Delta = 0.0;
      if (fabs(Mean_LevelSet) <= epsilon) Delta = 0.5*(1.0+cos(PI_NUMBER*Mean_LevelSet/epsilon))/epsilon;
      if (Mean_LevelSet > epsilon) Delta = 0.0;
      
      a = Mean_BetaInc2/Mean_DensityInc, b = Mean_LevelSet/Mean_DensityInc;
      c = (1.0 - config->GetRatioDensity())*Delta*config->GetDensity_FreeStreamND();
      e = (2.0*fabs(Mean_ProjVel) + b*c*fabs(Mean_ProjVel)), f = sqrt(4.0*a*Area*Area + e*e);
      Mean_SoundSpeed = 0.5*f;
      Mean_ProjVel = 0.5*e;
      
    }
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      if (compressible) {
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      }
      if (incompressible) {
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_BetaInc2 = node[iPoint]->GetBetaInc2();
        Mean_DensityInc = node[iPoint]->GetDensityInc();
        Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
      }
      if (freesurface) {
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_BetaInc2 = node[iPoint]->GetBetaInc2();
        Mean_DensityInc = node[iPoint]->GetDensityInc();
        Mean_LevelSet = node[iPoint]->GetLevelSet();
        
        if (Mean_LevelSet < -epsilon) Delta = 0.0;
        if (fabs(Mean_LevelSet) <= epsilon) Delta = 0.5*(1.0+cos(PI_NUMBER*Mean_LevelSet/epsilon))/epsilon;
        if (Mean_LevelSet > epsilon) Delta = 0.0;
        
        a = Mean_BetaInc2/Mean_DensityInc; b = Mean_LevelSet/Mean_DensityInc;
        c = (1.0 - config->GetRatioDensity())*Delta*config->GetDensity_FreeStreamND();
        e = (2.0*fabs(Mean_ProjVel) + b*c*fabs(Mean_ProjVel)); f = sqrt(4.0*a*Area*Area + e*e);
        Mean_SoundSpeed = 0.5*f;
        Mean_ProjVel = 0.5*e;
      }
      
      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;
    
    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool second_order = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
  bool low_fidelity = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, set normal vectors, and number of neighbors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
    /*--- Set primitive variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    
    /*--- Set the largest convective eigenvalue ---*/
    
    numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());
    
    /*--- Set undivided laplacian an pressure based sensor ---*/
    
    if ((second_order || low_fidelity)) {
      numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
      numerics->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor());
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- Compute residuals, and Jacobians ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    
    /*--- Update convective and artificial dissipation residuals ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit computation ---*/
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
  }
  
}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, RoeVelocity[3] = {0.0,0.0,0.0}, R, sq_vel, RoeEnthalpy,
  *V_i, *V_j, *S_i, *S_j, *Limiter_i = NULL, *Limiter_j = NULL, YDistance, GradHidrosPress, sqvel, Non_Physical = 1.0;
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;
  bool neg_density_i = false, neg_density_j = false, neg_pressure_i = false, neg_pressure_j = false, neg_sound_speed = false;
  
  
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool low_fidelity     = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool second_order     = (((config->GetSpatialOrder_Flow() == SECOND_ORDER) || (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)) && ((iMesh == MESH_0) || low_fidelity));
  bool limiter          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && !low_fidelity);
  bool freesurface      = (config->GetKind_Regime() == FREESURFACE);
  bool compressible     = (config->GetKind_Regime() == COMPRESSIBLE);
  bool grid_movement    = config->GetGrid_Movement();
  bool roe_turkel       = (config->GetKind_Upwind_Flow() == TURKEL);
  bool ideal_gas        = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );
  
  /*--- Loop over all the edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Roe Turkel preconditioning ---*/
    
    if (roe_turkel) {
      sqvel = 0.0;
      for (iDim = 0; iDim < nDim; iDim ++)
        sqvel += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
      numerics->SetVelocity2_Inf(sqvel);
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    /*--- Get primitive variables ---*/
    
    V_i = node[iPoint]->GetPrimitive(); V_j = node[jPoint]->GetPrimitive();
    S_i = node[iPoint]->GetSecondary(); S_j = node[jPoint]->GetSecondary();
    
    /*--- The zero order reconstruction includes the gradient
     of the hydrostatic pressure constribution ---*/
    
    if (freesurface) {
      
      YDistance = 0.5*(geometry->node[jPoint]->GetCoord(nDim-1)-geometry->node[iPoint]->GetCoord(nDim-1));
      GradHidrosPress = node[iPoint]->GetDensityInc()/(config->GetFroude()*config->GetFroude());
      Primitive_i[0] = V_i[0] - GradHidrosPress*YDistance;
      GradHidrosPress = node[jPoint]->GetDensityInc()/(config->GetFroude()*config->GetFroude());
      Primitive_j[0] = V_j[0] + GradHidrosPress*YDistance;
      
      for (iVar = 1; iVar < nPrimVar; iVar++) {
        Primitive_i[iVar] = V_i[iVar]+EPS;
        Primitive_j[iVar] = V_j[iVar]+EPS;
      }
      
    }
    
    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (second_order) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter_Primitive();
        Limiter_j = node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        Non_Physical = node[iPoint]->GetNon_Physical()*node[jPoint]->GetNon_Physical();
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim]*Non_Physical;
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim]*Non_Physical;
        }
        if (limiter) {
          Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      /*--- Recompute the extrapolated quantities in a
       thermodynamic consistent way  ---*/
      
      if (!ideal_gas) { ComputeConsExtrapolation(config); }
      
      /*--- Check for non-physical solutions after reconstruction. If found,
       use the cell-average value of the solution. This results in a locally
       first-order approximation, but this is typically only active
       during the start-up of a calculation. If non-physical, use the
       cell-averaged state. ---*/
      
      if (compressible) {
        
        neg_pressure_i = (Primitive_i[nDim+1] < 0.0); neg_pressure_j = (Primitive_j[nDim+1] < 0.0);
        neg_density_i  = (Primitive_i[nDim+2] < 0.0); neg_density_j  = (Primitive_j[nDim+2] < 0.0);
        
        R = sqrt(fabs(Primitive_j[nDim+2]/Primitive_i[nDim+2]));
        sq_vel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          RoeVelocity[iDim] = (R*Primitive_j[iDim+1]+Primitive_i[iDim+1])/(R+1);
          sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
        }
        RoeEnthalpy = (R*Primitive_j[nDim+3]+Primitive_i[nDim+3])/(R+1);
        neg_sound_speed = ((Gamma-1)*(RoeEnthalpy-0.5*sq_vel) < 0.0);
        
      }
      
      if (neg_sound_speed) {
        for (iVar = 0; iVar < nPrimVar; iVar++) {
          Primitive_i[iVar] = V_i[iVar];
          Primitive_j[iVar] = V_j[iVar]; }
        if (compressible) {
          Secondary_i[0] = S_i[0]; Secondary_i[1] = S_i[1];
          Secondary_j[0] = S_i[0]; Secondary_j[1] = S_i[1]; }
        counter_local++;
      }
      
      if (neg_density_i || neg_pressure_i) {
        for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = V_i[iVar];
        if (compressible) { Secondary_i[0] = S_i[0]; Secondary_i[1] = S_i[1]; }
        counter_local++;
      }
      
      if (neg_density_j || neg_pressure_j) {
        for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = V_j[iVar];
        if (compressible) { Secondary_j[0] = S_j[0]; Secondary_j[1] = S_j[1]; }
        counter_local++;
      }
      
      numerics->SetPrimitive(Primitive_i, Primitive_j);
      numerics->SetSecondary(Secondary_i, Secondary_j);
      
    }
    else {
      
      /*--- Set conservative variables without reconstruction ---*/
      
      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);
      
      if (freesurface) {
        numerics->SetPrimitive(Primitive_i, Primitive_j);
      }
      
    }
    
    /*--- Compute the residual ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    
    /*--- Update residual value ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit Jacobians ---*/
    
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
    
    /*--- Roe Turkel preconditioning, set the value of beta ---*/
    
    if (roe_turkel) {
      node[iPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
      node[jPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
    }
    
  }
  
  /*--- Warning message about non-physical reconstructions ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);
  }
  
}

void CEulerSolver::ComputeConsExtrapolation(CConfig *config) {
  
  unsigned short iDim;
  
  su2double density_i = Primitive_i[nDim+2];
  su2double pressure_i = Primitive_i[nDim+1];
  su2double velocity2_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    velocity2_i += Primitive_i[iDim+1]*Primitive_i[iDim+1];
  }
  
  FluidModel->SetTDState_Prho(pressure_i, density_i);
  
  Primitive_i[0]= FluidModel->GetTemperature();
  Primitive_i[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_i[nDim+1]/Primitive_i[nDim+2] + 0.5*velocity2_i;
  Primitive_i[nDim+4]= FluidModel->GetSoundSpeed();
  Secondary_i[0]=FluidModel->GetdPdrho_e();
  Secondary_i[1]=FluidModel->GetdPde_rho();
  
  
  su2double density_j = Primitive_j[nDim+2];
  su2double pressure_j = Primitive_j[nDim+1];
  su2double velocity2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    velocity2_j += Primitive_j[iDim+1]*Primitive_j[iDim+1];
  }
  
  FluidModel->SetTDState_Prho(pressure_j, density_j);
  
  Primitive_j[0]= FluidModel->GetTemperature();
  Primitive_j[nDim+3]= FluidModel->GetStaticEnergy() + Primitive_j[nDim+1]/Primitive_j[nDim+2] + 0.5*velocity2_j;
  Primitive_j[nDim+4]=FluidModel->GetSoundSpeed();
  Secondary_j[0]=FluidModel->GetdPdrho_e();
  Secondary_j[1]=FluidModel->GetdPde_rho();
  
}

void CEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar, jVar;
  unsigned long iPoint;
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  bool axisymmetric   = config->GetAxisymmetric();
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool gravity        = (config->GetGravityForce() == YES);
  bool time_spectral  = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool windgust       = config->GetWind_Gust();
  
  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  
  if (rotating_frame) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the conservative variables ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
      if (incompressible || freesurface) {
        /*--- Set incompressible density  ---*/
        numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());
      }
      
      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }
  
  if (gravity) {
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution  ---*/
      numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
      /*--- Set incompressible density  ---*/
      if (incompressible || freesurface) {
        numerics->SetDensityInc(node[iPoint]->GetDensityInc(), node[iPoint]->GetDensityInc());
      }
      
      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
    
  }
  
  if (freesurface) {
    
    unsigned long iPoint;
    su2double Vol, x_o, x_od, x, z, levelset, DampingFactor;
    su2double factor = config->GetFreeSurface_Damping_Length();
    
    x_o = config->GetFreeSurface_Outlet();
    x_od = x_o - factor*2.0*PI_NUMBER*config->GetFroude()*config->GetFroude();
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      Vol = geometry->node[iPoint]->GetVolume();
      x = geometry->node[iPoint]->GetCoord()[0];
      z = geometry->node[iPoint]->GetCoord()[nDim-1]-config->GetFreeSurface_Zero();
      levelset = node[iPoint]->GetSolution(nDim+1);
      
      DampingFactor = 0.0;
      if (x >= x_od)
        DampingFactor = config->GetFreeSurface_Damping_Coeff()*pow((x-x_od)/(x_o-x_od), 2.0);
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
      }
      
      Residual[nDim+1] = Vol*(levelset-z)*DampingFactor;
      Jacobian_i[nDim+1][nDim+1] = Vol*DampingFactor;
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
    
  }
  
  if (time_spectral) {
    
    su2double Volume, Source;
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Get control volume ---*/
      Volume = geometry->node[iPoint]->GetVolume();
      
      /*--- Get stored time spectral source term ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Source = node[iPoint]->GetTimeSpectral_Source(iVar);
        Residual[iVar] = Source*Volume;
      }
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }
  
  if (windgust) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the wind gust ---*/
      numerics->SetWindGust(node[iPoint]->GetWindGust(), node[iPoint]->GetWindGust());
      
      /*--- Load the wind gust derivatives ---*/
      numerics->SetWindGustDer(node[iPoint]->GetWindGustDer(), node[iPoint]->GetWindGustDer());
      
      /*--- Load the primitive variables ---*/
      numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[iPoint]->GetPrimitive());
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}

void CEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */
  
  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = node[iPoint]->GetSolution() */
  /* Get the volume (or area in 2D) associated with iPoint = node[iPoint]->GetVolume() */
  /* Get the vector of geometric coordinates of point iPoint = node[iPoint]->GetCoord() */
  
}

void CEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
  
  su2double *Normal, Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Mean_BetaInc2, Lambda, Mean_DensityInc,
  ProjVel, ProjVel_i, ProjVel_j, *GridVel, *GridVel_i, *GridVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement = config->GetGrid_Movement();
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetLambda(0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    if (compressible) {
      Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
      Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    }
    if (incompressible || freesurface) {
      Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
      Mean_BetaInc2 = 0.5 * (node[iPoint]->GetBetaInc2() + node[jPoint]->GetBetaInc2());
      Mean_DensityInc = 0.5 * (node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
      Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
    }
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      if (compressible) {
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      }
      if (incompressible || freesurface) {
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_BetaInc2 = node[iPoint]->GetBetaInc2();
        Mean_DensityInc = node[iPoint]->GetDensityInc();
        Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
      }
      
      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddLambda(Lambda);
      }
      
    }
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_MaxEigenvalue(geometry, config);
  
}

void CEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, jPoint, iEdge;
  su2double Pressure_i = 0, Pressure_j = 0, *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  
  Diff = new su2double[nVar];
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Solution differences ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
    /*--- Correction for compressible flows which use the enthalpy ---*/
    
    if (compressible) {
      Pressure_i = node[iPoint]->GetPressure();
      Pressure_j = node[jPoint]->GetPressure();
      Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar-1) + Pressure_j);
    }
    
    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Undivided_Laplacian(geometry, config);
  
  delete [] Diff;
  
}

void CEulerSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  su2double Pressure_i = 0.0, Pressure_j = 0.0;
  bool boundary_i, boundary_j;
  
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  
  /*--- Reset variables to store the undivided pressure ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    iPoint_UndLapl[iPoint] = 0.0;
    jPoint_UndLapl[iPoint] = 0.0;
  }
  
  /*--- Evaluate the pressure sensor ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Get the pressure, or density for incompressible solvers ---*/
    
    if (compressible) {
      Pressure_i = node[iPoint]->GetPressure();
      Pressure_j = node[jPoint]->GetPressure();
    }
    if (incompressible || freesurface) {
      Pressure_i = node[iPoint]->GetDensityInc();
      Pressure_j = node[jPoint]->GetDensityInc();
    }
    
    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both on the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    
  }
  
  /*--- Set pressure switch for each point ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetSensor(fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Dissipation_Switch(geometry, config);
  
}

void CEulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double Pressure = 0.0, *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord, Area,
  factor, NFPressOF, RefVel2, RefTemp, RefDensity, RefPressure, Mach2Vel, Mach_Motion,
  Force[3] = {0.0,0.0,0.0};
  string Marker_Tag, Monitoring_Tag;
  
#ifdef HAVE_MPI
  su2double MyAllBound_CDrag_Inv, MyAllBound_CLift_Inv, MyAllBound_CSideForce_Inv, MyAllBound_CMx_Inv, MyAllBound_CMy_Inv, MyAllBound_CMz_Inv, MyAllBound_CFx_Inv, MyAllBound_CFy_Inv, MyAllBound_CFz_Inv, MyAllBound_CT_Inv, MyAllBound_CQ_Inv, MyAllBound_CNearFieldOF_Inv, *MySurface_CLift_Inv = NULL, *MySurface_CDrag_Inv = NULL, *MySurface_CSideForce_Inv = NULL, *MySurface_CEff_Inv = NULL, *MySurface_CFx_Inv = NULL, *MySurface_CFy_Inv = NULL, *MySurface_CFz_Inv = NULL, *MySurface_CMx_Inv = NULL, *MySurface_CMy_Inv = NULL, *MySurface_CMz_Inv = NULL;
#endif
  
  su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefAreaCoeff    = config->GetRefAreaCoeff();
  su2double RefLengthMoment = config->GetRefLengthMoment();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double *Origin         = config->GetRefOriginMoment(0);
  bool grid_movement        = config->GetGrid_Movement();
  bool compressible         = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible       = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface          = (config->GetKind_Regime() == FREESURFACE);
  
  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream
   values, which is the standard convention. ---*/
  
  RefTemp     = Temperature_Inf;
  RefDensity  = Density_Inf;
  RefPressure = Pressure_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  /*-- Variables initialization ---*/
  
  Total_CDrag = 0.0;        Total_CLift = 0.0; Total_CSideForce = 0.0; Total_CEff = 0.0;
  Total_CMx = 0.0;          Total_CMy = 0.0;   Total_CMz = 0.0;
  Total_CFx = 0.0;          Total_CFy = 0.0;   Total_CFz = 0.0;
  Total_CT = 0.0;           Total_CQ = 0.0;    Total_CMerit = 0.0;
  Total_CNearFieldOF = 0.0; Total_Heat = 0.0;  Total_MaxHeat = 0.0;
  
  AllBound_CDrag_Inv = 0.0;        AllBound_CLift_Inv = 0.0; AllBound_CSideForce_Inv = 0.0;
  AllBound_CMx_Inv = 0.0;          AllBound_CMy_Inv = 0.0;   AllBound_CMz_Inv = 0.0;
  AllBound_CFx_Inv = 0.0;          AllBound_CFy_Inv = 0.0;   AllBound_CFz_Inv = 0.0;
  AllBound_CT_Inv = 0.0;           AllBound_CQ_Inv = 0.0;    AllBound_CMerit_Inv = 0.0;
  AllBound_CNearFieldOF_Inv = 0.0; AllBound_CEff_Inv = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift_Inv[iMarker_Monitoring]      = 0.0; Surface_CDrag_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CSideForce_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0; Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0; Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0; Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CLift[iMarker_Monitoring]          = 0.0; Surface_CDrag[iMarker_Monitoring]          = 0.0;
    Surface_CSideForce[iMarker_Monitoring]     = 0.0; Surface_CEff[iMarker_Monitoring]           = 0.0;
    Surface_CFx[iMarker_Monitoring]            = 0.0; Surface_CFy[iMarker_Monitoring]            = 0.0;
    Surface_CFz[iMarker_Monitoring]            = 0.0; Surface_CMx[iMarker_Monitoring]            = 0.0;
    Surface_CMy[iMarker_Monitoring]            = 0.0; Surface_CMz[iMarker_Monitoring]            = 0.0;
  }
  
  /*--- Loop over the Euler and Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
        (Boundary == ISOTHERMAL) || (Boundary == NEARFIELD_BOUNDARY)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CDrag_Inv[iMarker] = 0.0;        CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;
      CMx_Inv[iMarker] = 0.0;          CMy_Inv[iMarker] = 0.0;   CMz_Inv[iMarker] = 0.0;
      CFx_Inv[iMarker] = 0.0;          CFy_Inv[iMarker] = 0.0;   CFz_Inv[iMarker] = 0.0;
      CT_Inv[iMarker] = 0.0;           CQ_Inv[iMarker] = 0.0;    CMerit_Inv[iMarker] = 0.0;
      CNearFieldOF_Inv[iMarker] = 0.0; CEff_Inv[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
      MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
      NFPressOF = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (compressible)   Pressure = node[iPoint]->GetPressure();
        if (incompressible || freesurface) Pressure = node[iPoint]->GetPressureInc();
        
        CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefAreaCoeff;
        
        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          
          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/
          
          NFPressOF += 0.5*(Pressure - Pressure_Inf)*(Pressure - Pressure_Inf)*Normal[nDim-1];
          
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) {
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }
          
          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf)*Normal[iDim]*factor;
            ForceInviscid[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (nDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLengthMoment;
            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLengthMoment;
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        
        if (Boundary != NEARFIELD_BOUNDARY) {
          if (nDim == 2) {
            CDrag_Inv[iMarker]  =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            CLift_Inv[iMarker]  = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
            CEff_Inv[iMarker]   = CLift_Inv[iMarker] / (CDrag_Inv[iMarker]+EPS);
            CMz_Inv[iMarker]    = MomentInviscid[2];
            CFx_Inv[iMarker]    = ForceInviscid[0];
            CFy_Inv[iMarker]    = ForceInviscid[1];
            CT_Inv[iMarker]     = -CFx_Inv[iMarker];
            CQ_Inv[iMarker]     = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker] = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }
          if (nDim == 3) {
            CDrag_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
            CLift_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
            CSideForce_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            CEff_Inv[iMarker]       = CLift_Inv[iMarker] / (CDrag_Inv[iMarker] + EPS);
            CMx_Inv[iMarker]        = MomentInviscid[0];
            CMy_Inv[iMarker]        = MomentInviscid[1];
            CMz_Inv[iMarker]        = MomentInviscid[2];
            CFx_Inv[iMarker]        = ForceInviscid[0];
            CFy_Inv[iMarker]        = ForceInviscid[1];
            CFz_Inv[iMarker]        = ForceInviscid[2];
            CT_Inv[iMarker]         = -CFz_Inv[iMarker];
            CQ_Inv[iMarker]         = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker]     = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }
          
          AllBound_CDrag_Inv        += CDrag_Inv[iMarker];
          AllBound_CLift_Inv        += CLift_Inv[iMarker];
          AllBound_CSideForce_Inv   += CSideForce_Inv[iMarker];
          AllBound_CEff_Inv          = AllBound_CLift_Inv / (AllBound_CDrag_Inv + EPS);
          AllBound_CMx_Inv          += CMx_Inv[iMarker];
          AllBound_CMy_Inv          += CMy_Inv[iMarker];
          AllBound_CMz_Inv          += CMz_Inv[iMarker];
          AllBound_CFx_Inv          += CFx_Inv[iMarker];
          AllBound_CFy_Inv          += CFy_Inv[iMarker];
          AllBound_CFz_Inv          += CFz_Inv[iMarker];
          AllBound_CT_Inv           += CT_Inv[iMarker];
          AllBound_CQ_Inv           += CQ_Inv[iMarker];
          AllBound_CMerit_Inv        = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
          
          /*--- Compute the coefficients per surface ---*/
          
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
            Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            if (Marker_Tag == Monitoring_Tag) {
              Surface_CLift_Inv[iMarker_Monitoring]      += CLift_Inv[iMarker];
              Surface_CDrag_Inv[iMarker_Monitoring]      += CDrag_Inv[iMarker];
              Surface_CSideForce_Inv[iMarker_Monitoring] += CSideForce_Inv[iMarker];
              Surface_CEff_Inv[iMarker_Monitoring]        = CLift_Inv[iMarker] / (CDrag_Inv[iMarker] + EPS);
              Surface_CFx_Inv[iMarker_Monitoring]        += CFx_Inv[iMarker];
              Surface_CFy_Inv[iMarker_Monitoring]        += CFy_Inv[iMarker];
              Surface_CFz_Inv[iMarker_Monitoring]        += CFz_Inv[iMarker];
              Surface_CMx_Inv[iMarker_Monitoring]        += CMx_Inv[iMarker];
              Surface_CMy_Inv[iMarker_Monitoring]        += CMy_Inv[iMarker];
              Surface_CMz_Inv[iMarker_Monitoring]        += CMz_Inv[iMarker];
            }
          }
          
        }
        
        /*--- At the Nearfield SU2 only cares about the pressure coeffient ---*/
        
        else {
          CNearFieldOF_Inv[iMarker] = NFPressOF;
          AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
        }
        
      }
      
      
    }
  }
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CDrag_Inv        = AllBound_CDrag_Inv;        AllBound_CDrag_Inv = 0.0;
  MyAllBound_CLift_Inv        = AllBound_CLift_Inv;        AllBound_CLift_Inv = 0.0;
  MyAllBound_CSideForce_Inv   = AllBound_CSideForce_Inv;   AllBound_CSideForce_Inv = 0.0;
  AllBound_CEff_Inv = 0.0;
  MyAllBound_CMx_Inv          = AllBound_CMx_Inv;          AllBound_CMx_Inv = 0.0;
  MyAllBound_CMy_Inv          = AllBound_CMy_Inv;          AllBound_CMy_Inv = 0.0;
  MyAllBound_CMz_Inv          = AllBound_CMz_Inv;          AllBound_CMz_Inv = 0.0;
  MyAllBound_CFx_Inv          = AllBound_CFx_Inv;          AllBound_CFx_Inv = 0.0;
  MyAllBound_CFy_Inv          = AllBound_CFy_Inv;          AllBound_CFy_Inv = 0.0;
  MyAllBound_CFz_Inv          = AllBound_CFz_Inv;          AllBound_CFz_Inv = 0.0;
  MyAllBound_CT_Inv           = AllBound_CT_Inv;           AllBound_CT_Inv = 0.0;
  MyAllBound_CQ_Inv           = AllBound_CQ_Inv;           AllBound_CQ_Inv = 0.0;
  AllBound_CMerit_Inv = 0.0;
  MyAllBound_CNearFieldOF_Inv = AllBound_CNearFieldOF_Inv; AllBound_CNearFieldOF_Inv = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CDrag_Inv, &AllBound_CDrag_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CLift_Inv, &AllBound_CLift_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSideForce_Inv, &AllBound_CSideForce_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Inv = AllBound_CLift_Inv / (AllBound_CDrag_Inv + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Inv, &AllBound_CMx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Inv, &AllBound_CMy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Inv, &AllBound_CMz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Inv, &AllBound_CFx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Inv, &AllBound_CFy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Inv, &AllBound_CFz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Inv, &AllBound_CT_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Inv, &AllBound_CQ_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Inv = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CNearFieldOF_Inv, &AllBound_CNearFieldOF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CLift_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CDrag_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSideForce_Inv = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CLift_Inv[iMarker_Monitoring]      = Surface_CLift_Inv[iMarker_Monitoring];
    MySurface_CDrag_Inv[iMarker_Monitoring]      = Surface_CDrag_Inv[iMarker_Monitoring];
    MySurface_CSideForce_Inv[iMarker_Monitoring] = Surface_CSideForce_Inv[iMarker_Monitoring];
    MySurface_CEff_Inv[iMarker_Monitoring]       = Surface_CEff_Inv[iMarker_Monitoring];
    MySurface_CFx_Inv[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    MySurface_CFy_Inv[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    MySurface_CFz_Inv[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    MySurface_CMx_Inv[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    MySurface_CMy_Inv[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    MySurface_CMz_Inv[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
    
    Surface_CLift_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CDrag_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CSideForce_Inv[iMarker_Monitoring] = 0.0;
    Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CLift_Inv, Surface_CLift_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CDrag_Inv, Surface_CDrag_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSideForce_Inv, Surface_CSideForce_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Inv[iMarker_Monitoring] = Surface_CLift_Inv[iMarker_Monitoring] / (Surface_CDrag_Inv[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Inv, Surface_CFx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Inv, Surface_CFy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Inv, Surface_CFz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Inv, Surface_CMx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Inv, Surface_CMy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Inv, Surface_CMz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  delete [] MySurface_CLift_Inv; delete [] MySurface_CDrag_Inv; delete [] MySurface_CSideForce_Inv;
  delete [] MySurface_CEff_Inv;  delete [] MySurface_CFx_Inv;   delete [] MySurface_CFy_Inv;
  delete [] MySurface_CFz_Inv;   delete [] MySurface_CMx_Inv;   delete [] MySurface_CMy_Inv;
  delete [] MySurface_CMz_Inv;
  
#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  
  Total_CDrag         = AllBound_CDrag_Inv;
  Total_CLift         = AllBound_CLift_Inv;
  Total_CSideForce    = AllBound_CSideForce_Inv;
  Total_CEff          = Total_CLift / (Total_CDrag + EPS);
  Total_CMx           = AllBound_CMx_Inv;
  Total_CMy           = AllBound_CMy_Inv;
  Total_CMz           = AllBound_CMz_Inv;
  Total_CFx           = AllBound_CFx_Inv;
  Total_CFy           = AllBound_CFy_Inv;
  Total_CFz           = AllBound_CFz_Inv;
  Total_CT            = AllBound_CT_Inv;
  Total_CQ            = AllBound_CQ_Inv;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);
  Total_CNearFieldOF  = AllBound_CNearFieldOF_Inv;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift[iMarker_Monitoring]      = Surface_CLift_Inv[iMarker_Monitoring];
    Surface_CDrag[iMarker_Monitoring]      = Surface_CDrag_Inv[iMarker_Monitoring];
    Surface_CSideForce[iMarker_Monitoring] = Surface_CSideForce_Inv[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CLift_Inv[iMarker_Monitoring] / (Surface_CDrag_Inv[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
  }
  
}

void CEulerSolver::TurboPerformance(CSolver *solver, CConfig *config, unsigned short inMarker,  unsigned short outMarker, unsigned short Kind_TurboPerf, unsigned short inMarkerTP ){
  
  su2double  avgVel2In, avgVel2Out,avgVelRel2In, avgVelRel2Out, avgGridVel2In, avgGridVel2Out, avgTotalEnthalpyIn= 0.0,avgTotalRothalpyIn,
  avgTotalEnthalpyOut, avgTotalRothalpyOut, avgTotalEnthalpyOutIs, avgEnthalpyOut, avgEnthalpyOutIs,
  avgPressureOut, avgTotalRelPressureIn, avgTotalRelPressureOut, avgEntropyIn, avgEntropyOut;
  unsigned short iDim;
  
  
  /*--- compute or retrieve inlet information ---*/
  avgVelRel2In= 0.0;
  avgGridVel2In= 0.0;
  avgVel2In= 0.0;
  for (iDim = 0; iDim < nDim; iDim++){
    avgVelRel2In +=( AveragedVelocity[inMarker][iDim] - AveragedGridVel[inMarker][iDim])*( AveragedVelocity[inMarker][iDim] - AveragedGridVel[inMarker][iDim]);
    avgGridVel2In += AveragedGridVel[inMarker][iDim]*AveragedGridVel[inMarker][iDim];
    avgVel2In += AveragedVelocity[inMarker][iDim]*AveragedVelocity[inMarker][iDim];
  }
  
  avgTotalRothalpyIn = AveragedEnthalpy[inMarker] + 0.5*avgVelRel2In - 0.5*avgGridVel2In;
  avgTotalEnthalpyIn = AveragedEnthalpy[inMarker] + 0.5*avgVel2In;
  avgEntropyIn = AveragedEntropy[inMarker];
  FluidModel->SetTDState_hs(avgTotalRothalpyIn, avgEntropyIn);
  avgTotalRelPressureIn  = FluidModel->GetPressure();

  /*--- compute or retrieve outlet information ---*/
  avgVelRel2Out = 0.0;
  avgGridVel2Out = 0.0;
  avgVel2Out = 0.0;
  for (iDim = 0; iDim < nDim; iDim++){
    avgVelRel2Out += (solver->GetAveragedVelocity(outMarker)[iDim]- solver->GetAveragedGridVelocity(outMarker)[iDim])*(solver->GetAveragedVelocity(outMarker)[iDim]- solver->GetAveragedGridVelocity(outMarker)[iDim]);
    avgGridVel2Out += solver->GetAveragedGridVelocity(outMarker)[iDim]*solver->GetAveragedGridVelocity(outMarker)[iDim];
    avgVel2Out += solver->GetAveragedVelocity(outMarker)[iDim]*solver->GetAveragedVelocity(outMarker)[iDim];
  }
  avgTotalRothalpyOut = solver->GetAveragedEnthalpy(outMarker) + 0.5*avgVelRel2Out - 0.5*avgGridVel2Out;
  avgTotalEnthalpyOut = solver->GetAveragedEnthalpy(outMarker) + 0.5*avgVel2Out;
  avgEntropyOut = solver->GetAveragedEntropy(outMarker);
  avgEnthalpyOut = solver->GetAveragedEnthalpy(outMarker);
  FluidModel->SetTDState_hs(avgTotalRothalpyOut, avgEntropyOut);
  avgTotalRelPressureOut  =  FluidModel->GetPressure();
  avgPressureOut= solver->GetAveragedPressure(outMarker);
  
  /*--- compute outlet isoentropic conditions ---*/
  FluidModel->SetTDState_Ps(avgPressureOut, avgEntropyIn);
  avgEnthalpyOutIs = FluidModel->GetStaticEnergy() + avgPressureOut/FluidModel->GetDensity();
  avgTotalEnthalpyOutIs = avgEnthalpyOutIs + 0.5*avgVel2Out;
  
  /*--- store turboperformance informations ---*/
  PressureOut[inMarkerTP] = avgPressureOut;
  PressureRatio[inMarkerTP] = avgTotalRelPressureIn/avgPressureOut;
  
  switch(Kind_TurboPerf){
    case BLADE:
      
      TotalPressureLoss[inMarkerTP] = (avgTotalRelPressureIn - avgTotalRelPressureOut)/(avgTotalRelPressureOut - avgPressureOut) ;
      KineticEnergyLoss[inMarkerTP] = (avgEnthalpyOut - avgEnthalpyOutIs)/(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out);
      EulerianWork[inMarkerTP] = avgTotalEnthalpyIn - avgTotalEnthalpyOut;
      TotalEnthalpyIn[inMarkerTP] = avgTotalRothalpyIn;
      FlowAngleIn[inMarkerTP]= FlowAngle[inMarker];
      FlowAngleOut[inMarkerTP]= solver->GetFlowAngle(outMarker);
      MassFlowIn[inMarkerTP]= MassFlow[inMarker];
      MassFlowOut[inMarkerTP]= solver->GetMassFlow(outMarker);
      MachIn[inMarkerTP]= AveragedMach[inMarker];
      MachOut[inMarkerTP]= solver->GetAveragedMach(outMarker);
      NormalMachIn[inMarkerTP]= AveragedNormalMach[inMarker];
      NormalMachOut[inMarkerTP]= solver->GetAveragedNormalMach(outMarker);
      EnthalpyOut[inMarkerTP]= avgEnthalpyOut;
      VelocityOutIs[inMarkerTP]=sqrt(2.0*(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out));
      break;
      
    case STAGE: case TURBINE:
      
      TotalTotalEfficiency[inMarkerTP] = (avgTotalEnthalpyIn - avgTotalEnthalpyOut)/(avgTotalEnthalpyIn - avgTotalEnthalpyOutIs);
      TotalStaticEfficiency[inMarkerTP] = (avgTotalEnthalpyIn - avgTotalEnthalpyOut)/(avgTotalEnthalpyIn - avgEnthalpyOutIs);
      TotalEnthalpyIn[inMarkerTP]= avgTotalEnthalpyIn;
      EnthalpyOut[inMarkerTP] = avgTotalEnthalpyOut;
      break;
      
    default:
      cout << "Warning! Invalid TurboPerformance option!" << endl;
      exit(EXIT_FAILURE);
      break;
  }
  
  
  
  
  
}

void CEulerSolver::MPITurboPerformance(CConfig *config){

	unsigned short iMarker, iMarkerTP;
	su2double  avgVel2In, avgVel2Out,avgVelRel2In, avgVelRel2Out, avgGridVel2In, avgGridVel2Out, avgTotalEnthalpyIn= 0.0,avgTotalRothalpyIn,
	avgTotalEnthalpyOut, avgTotalRothalpyOut, avgTotalEnthalpyOutIs, avgEnthalpyOut, avgEnthalpyOutIs,
	avgPressureOut, avgTotalRelPressureIn, avgTotalRelPressureOut, avgEntropyIn, avgEntropyOut, flowAngleIn, massFlowIn, machIn, normalMachIn, 	flowAngleOut,
	massFlowOut, machOut, normalMachOut;

	unsigned short iDim, i, n1, n2, n1t,n2t;

	int rank = MASTER_NODE, rankIn = MASTER_NODE, rankOut= MASTER_NODE;
	int size = SINGLE_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  su2double *TurbPerfIn= NULL,*TurbPerfOut= NULL;
  su2double *TotTurbPerfIn = NULL,*TotTurbPerfOut = NULL;
  n1  = 8;
  n2  = 10;
  n1t = n1*size;
  n2t = n2*size;
  TurbPerfIn = new su2double[n1];
  TurbPerfOut = new su2double[n2];

  for (i=0;i<n1;i++)
  	TurbPerfIn[i]= -1.0;
  for (i=0;i<n2;i++)
		TurbPerfOut[i]= -1.0;
#endif

  avgTotalRothalpyIn     = -1.0;
  avgTotalEnthalpyIn		 = -1.0;
  avgEntropyIn           = -1.0;
  avgTotalRelPressureIn  = -1.0;
	flowAngleIn						 = -1.0;
	massFlowIn						 = -1.0;
	machIn								 = -1.0;
	normalMachIn					 = -1.0;
  avgTotalRothalpyOut    = -1.0;
	avgTotalEnthalpyOut    = -1.0;
	avgTotalRelPressureOut = -1.0;
	avgPressureOut				 = -1.0;
	avgEnthalpyOut				 = -1.0;
	avgGridVel2Out				 = -1.0;
	flowAngleOut					 = -1.0;
	massFlowOut						 = -1.0;
	machOut								 = -1.0;
	normalMachOut					 = -1.0;

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iMarkerTP=1; iMarkerTP < config->Get_nMarkerTurboPerf()+1; iMarkerTP++)
			if (config->GetMarker_All_TurboPerformance(iMarker) == iMarkerTP){

				/*--- compute or retrieve inlet information ---*/
				if (config->GetMarker_All_TurboPerformanceFlag(iMarker) == INFLOW){
					avgVelRel2In= 0.0;
					avgGridVel2In= 0.0;
					avgVel2In= 0.0;
					for (iDim = 0; iDim < nDim; iDim++){
						avgVelRel2In +=( AveragedVelocity[iMarker][iDim] - AveragedGridVel[iMarker][iDim])*( AveragedVelocity[iMarker][iDim] - AveragedGridVel[iMarker][iDim]);
						avgGridVel2In += AveragedGridVel[iMarker][iDim]*AveragedGridVel[iMarker][iDim];
						avgVel2In += AveragedVelocity[iMarker][iDim]*AveragedVelocity[iMarker][iDim];
					}

					avgTotalRothalpyIn 			= AveragedEnthalpy[iMarker] + 0.5*avgVelRel2In - 0.5*avgGridVel2In;
					avgTotalEnthalpyIn 			= AveragedEnthalpy[iMarker] + 0.5*avgVel2In;
					avgEntropyIn 						= AveragedEntropy[iMarker];
					FluidModel->SetTDState_hs(avgTotalRothalpyIn, avgEntropyIn);
					avgTotalRelPressureIn   = FluidModel->GetPressure();
					flowAngleIn							= FlowAngle[iMarker];
					massFlowIn							= MassFlow[iMarker];
					machIn									= AveragedMach[iMarker];
					normalMachIn						= AveragedNormalMach[iMarker];
#ifdef HAVE_MPI
					TurbPerfIn[0] = avgTotalRothalpyIn;
					TurbPerfIn[1] = avgTotalEnthalpyIn;
					TurbPerfIn[2] = avgEntropyIn;
					TurbPerfIn[3] = avgTotalRelPressureIn;
					TurbPerfIn[4] = flowAngleIn;
					TurbPerfIn[5] = massFlowIn;
					TurbPerfIn[6] = machIn;
					TurbPerfIn[7] =	normalMachIn;
#endif
				}

				/*--- compute or retrieve outlet information ---*/
				if (config->GetMarker_All_TurboPerformanceFlag(iMarker) == OUTFLOW){
					avgVelRel2Out = 0.0;
					avgGridVel2Out = 0.0;
					avgVel2Out = 0.0;
					for (iDim = 0; iDim < nDim; iDim++){
						avgVelRel2Out += ( AveragedVelocity[iMarker][iDim] - AveragedGridVel[iMarker][iDim])*( AveragedVelocity[iMarker][iDim] - AveragedGridVel[iMarker][iDim]);
						avgGridVel2Out += AveragedGridVel[iMarker][iDim]*AveragedGridVel[iMarker][iDim];
						avgVel2Out += AveragedVelocity[iMarker][iDim]*AveragedVelocity[iMarker][iDim];
					}
					avgTotalRothalpyOut       = AveragedEnthalpy[iMarker] + 0.5*avgVelRel2Out - 0.5*avgGridVel2Out;
					avgTotalEnthalpyOut       = AveragedEnthalpy[iMarker] + 0.5*avgVel2Out;
					avgEntropyOut             = AveragedEntropy[iMarker];
					avgEnthalpyOut            = AveragedEnthalpy[iMarker];
					FluidModel->SetTDState_hs(avgTotalRothalpyOut, avgEntropyOut);
					avgTotalRelPressureOut    =  FluidModel->GetPressure();
					avgPressureOut= AveragedPressure[iMarker];
					flowAngleOut							= FlowAngle[iMarker];
					massFlowOut							  = MassFlow[iMarker];
					machOut									  = AveragedMach[iMarker];
					normalMachOut						  = AveragedNormalMach[iMarker];

#ifdef HAVE_MPI
					TurbPerfOut[0] = avgTotalRothalpyOut;
					TurbPerfOut[1] = avgTotalEnthalpyOut;
					TurbPerfOut[2] = avgTotalRelPressureOut;
					TurbPerfOut[3] = avgPressureOut;
					TurbPerfOut[4] = avgEnthalpyOut;
					TurbPerfOut[5] = avgGridVel2Out;
					TurbPerfOut[6] = flowAngleOut;
					TurbPerfOut[7] = massFlowOut;
					TurbPerfOut[8] = machOut;
					TurbPerfOut[9] = normalMachOut;
#endif

				}
			}


#ifdef HAVE_MPI
if (rank == MASTER_NODE){
  TotTurbPerfIn = new su2double[n1t];
  TotTurbPerfOut = new su2double[n2t];
  for (i=0;i<n1t;i++)
  	TotTurbPerfIn[i]= -1.0;
  for (i=0;i<n2t;i++)
		TotTurbPerfOut[i]= -1.0;
	}
	SU2_MPI::Gather(TurbPerfIn, n1, MPI_DOUBLE, TotTurbPerfIn, n1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
	SU2_MPI::Gather(TurbPerfOut, n2, MPI_DOUBLE,TotTurbPerfOut, n2, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

	delete [] TurbPerfIn, delete [] TurbPerfOut;

if (rank == MASTER_NODE){
		for (i=0;i<size;i++){

			if(TotTurbPerfIn[n1*i] > 0.0){
				avgTotalRothalpyIn 		 = 0.0;
				avgTotalRothalpyIn 		 = TotTurbPerfIn[n1*i];
				avgTotalEnthalpyIn 		 = 0.0;
				avgTotalEnthalpyIn 		 = TotTurbPerfIn[n1*i+1];
				avgEntropyIn 					 = 0.0;
				avgEntropyIn 			     = TotTurbPerfIn[n1*i+2];
				avgTotalRelPressureIn  = 0.0;
				avgTotalRelPressureIn  = TotTurbPerfIn[n1*i+3];
				flowAngleIn						 = 0.0;
				flowAngleIn						 = TotTurbPerfIn[n1*i+4];
				massFlowIn						 = 0.0;
				massFlowIn						 = TotTurbPerfIn[n1*i+5];
				machIn								 = 0.0;
				machIn								 = TotTurbPerfIn[n1*i+6];
				normalMachIn					 = 0.0;
				normalMachIn					 = TotTurbPerfIn[n1*i+7];
			}

			if(TotTurbPerfOut[n2*i] > 0.0){
				avgTotalRothalpyOut    = 0.0;
				avgTotalRothalpyOut    = TotTurbPerfOut[n2*i];
				avgTotalEnthalpyOut    = 0.0;
				avgTotalEnthalpyOut    = TotTurbPerfOut[n2*i+1];
				avgTotalRelPressureOut = 0.0;
				avgTotalRelPressureOut = TotTurbPerfOut[n2*i+2];
				avgPressureOut				 = 0.0;
				avgPressureOut				 = TotTurbPerfOut[n2*i+3];
				avgEnthalpyOut				 = 0.0;
				avgEnthalpyOut				 = TotTurbPerfOut[n2*i+4];
				avgGridVel2Out				 = 0.0;
				avgGridVel2Out				 = TotTurbPerfOut[n2*i+5];
				flowAngleOut					 = 0.0;
				flowAngleOut					 = TotTurbPerfOut[n2*i+6];
				massFlowOut						 = 0.0;
				massFlowOut 					 = TotTurbPerfOut[n2*i+7];
				machOut								 = 0.0;
				machOut								 = TotTurbPerfOut[n2*i+8];
				normalMachOut					 = 0.0;
				normalMachOut					 = TotTurbPerfOut[n2*i+9];
			}
		}

		delete [] TotTurbPerfIn, delete [] TotTurbPerfOut;
}

#endif
	if (rank == MASTER_NODE){

				/*--- compute outlet isoentropic conditions ---*/
				FluidModel->SetTDState_Ps(avgPressureOut, avgEntropyIn);
				avgEnthalpyOutIs = FluidModel->GetStaticEnergy() + avgPressureOut/FluidModel->GetDensity();
				avgTotalEnthalpyOutIs = avgEnthalpyOutIs + 0.5*avgVel2Out;

				/*--- store turboperformance informations ---*/
				PressureOut[0] = avgPressureOut;
				PressureRatio[0] = avgTotalRelPressureIn/avgPressureOut;

				switch(config->GetKind_TurboPerf(0)){
					case BLADE:

						TotalPressureLoss[0] = (avgTotalRelPressureIn - avgTotalRelPressureOut)/(avgTotalRelPressureOut - avgPressureOut) ;
						KineticEnergyLoss[0] = (avgEnthalpyOut - avgEnthalpyOutIs)/(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out);
						EulerianWork[0]      = avgTotalEnthalpyIn - avgTotalEnthalpyOut;
						TotalEnthalpyIn[0]   = avgTotalRothalpyIn;
						FlowAngleIn[0]       = flowAngleIn;
						FlowAngleOut[0]      = flowAngleOut;
						MassFlowIn[0]        = massFlowIn;
						MassFlowOut[0]       = massFlowOut;
						MachIn[0]            = machIn;
						MachOut[0]           = machOut;
						NormalMachIn[0]      = normalMachIn;
						NormalMachOut[0]     = normalMachOut;
						EnthalpyOut[0]       = avgEnthalpyOut;
						VelocityOutIs[0]     =sqrt(2.0*(avgTotalRothalpyIn - avgEnthalpyOut + 0.5*avgGridVel2Out));
						break;

					case STAGE: case TURBINE:

						TotalTotalEfficiency[0] = (avgTotalEnthalpyIn - avgTotalEnthalpyOut)/(avgTotalEnthalpyIn - avgTotalEnthalpyOutIs);
						TotalStaticEfficiency[0] = (avgTotalEnthalpyIn - avgTotalEnthalpyOut)/(avgTotalEnthalpyIn - avgEnthalpyOutIs);
						TotalEnthalpyIn[0]= avgTotalEnthalpyIn;
						EnthalpyOut[0] = avgTotalEnthalpyOut;
						break;

					default:
						cout << "Warning! Invalid TurboPerformance option!" << endl;
						exit(EXIT_FAILURE);
						break;
				}


	}
}


void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetAdjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);
    
    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
  
}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  bool adjoint = config->GetAdjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);
    
    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
    
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar, jVar;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Delta, *local_Res_TruncError, Vol;
  
  bool adjoint = config->GetAdjoint();
  bool roe_turkel = config->GetKind_Upwind_Flow() == TURKEL;
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    
    if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
      if (roe_turkel) {
        SetPreconditioner(config, iPoint);
        for (iVar = 0; iVar < nVar; iVar ++ )
          for (jVar = 0; jVar < nVar; jVar ++ )
            LowMach_Precontioner[iVar][jVar] = Delta*LowMach_Precontioner[iVar][jVar];
        Jacobian.AddBlock(iPoint, iPoint, LowMach_Precontioner);
      }
      else {
        Jacobian.AddVal2Diag(iPoint, Delta);
      }
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }
  
  /*--- Solve or smooth the linear system ---*/
  
  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint]->AddSolution(iVar, config->GetRelaxation_Factor_Flow()*LinSysSol[iPoint*nVar+iVar]);
      }
    }
  }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CEulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal;
  
  /*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
  PrimVar_Vertex = new su2double [nPrimVarGrad];
  PrimVar_i = new su2double [nPrimVarGrad];
  PrimVar_j = new su2double [nPrimVarGrad];
  
  /*--- Set Gradient_Primitive to zero ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
      PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
    }
    
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
      }
    }
  }
  
  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          PrimVar_Vertex[iVar] = node[iPoint]->GetPrimitive(iVar);
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
            node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
          }
      }
    }
  }
  
  /*--- Update gradient value ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar, iDim) / (geometry->node[iPoint]->GetVolume());
        node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
      }
    }
  }
  
  delete [] PrimVar_Vertex;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
  Set_MPI_Primitive_Gradient(geometry, config);
  
}

void CEulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get primitives from CVariable ---*/
    
    PrimVar_i = node[iPoint]->GetPrimitive();
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        cvector[iVar][iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      PrimVar_j = node[jPoint]->GetPrimitive();
      
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (weight != 0.0) {
        
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        
        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
        
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    
    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }
    
    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
        }
        
        node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
      }
    }
    
  }
  
  Set_MPI_Primitive_Gradient(geometry, config);
  
}

void CEulerSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Primitive_i, *Primitive_j,
  dave, LimK, eps2, eps1, dm, dp, du, y, limiter;
  
  /*--- Initialize solution max and solution min and the limiter in the entire domain --*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      node[iPoint]->SetSolution_Max(iVar, -EPS);
      node[iPoint]->SetSolution_Min(iVar, EPS);
      node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
    }
  }
  
  /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Get the primitive variables ---*/
    
    Primitive_i = node[iPoint]->GetPrimitive();
    Primitive_j = node[jPoint]->GetPrimitive();
    
    /*--- Compute the maximum, and minimum values for nodes i & j ---*/
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      du = (Primitive_j[iVar] - Primitive_i[iVar]);
      node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
      node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
      node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
      node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
    }
    
  }
  
  
  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
          else dp = node[iPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar))
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
          else dp = node[jPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar))
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
        
      }
      
    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
      }
    }
    
  }
  
  /*--- Venkatakrishnan limiter ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    
    dave = config->GetRefElemLength();
    LimK = config->GetLimiterCoeff();
    eps1 = LimK*dave;
    eps2 = eps1*eps1*eps1;
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar))
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar))
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
      }
      
    }
    
  }
  
  /*--- Limiter MPI ---*/
  
  Set_MPI_Primitive_Limiter(geometry, config);
  
}

//void CEulerSolver::SetSecondary_Gradient_GG(CGeometry *geometry, CConfig *config) {
//  unsigned long iPoint, jPoint, iEdge, iVertex;
//  unsigned short iDim, iVar, iMarker;
//  su2double *SecondaryVar_Vertex, *SecondaryVar_i, *SecondaryVar_j, SecondaryVar_Average,
//  Partial_Gradient, Partial_Res, *Normal;
//
//  /*--- Gradient Secondary variables compressible (temp, vx, vy, vz, P, rho)
//   Gradient Secondary variables incompressible (rho, vx, vy, vz, beta) ---*/
//  SecondaryVar_Vertex = new su2double [nSecondaryVarGrad];
//  SecondaryVar_i = new su2double [nSecondaryVarGrad];
//  SecondaryVar_j = new su2double [nSecondaryVarGrad];
//
//  /*--- Set Gradient_Secondary to zero ---*/
//  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
//    node[iPoint]->SetGradient_SecondaryZero(nSecondaryVarGrad);
//
//  /*--- Loop interior edges ---*/
//  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//    iPoint = geometry->edge[iEdge]->GetNode(0);
//    jPoint = geometry->edge[iEdge]->GetNode(1);
//
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      SecondaryVar_i[iVar] = node[iPoint]->GetSecondary(iVar);
//      SecondaryVar_j[iVar] = node[jPoint]->GetSecondary(iVar);
//    }
//
//    Normal = geometry->edge[iEdge]->GetNormal();
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      SecondaryVar_Average =  0.5 * ( SecondaryVar_i[iVar] + SecondaryVar_j[iVar] );
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Partial_Res = SecondaryVar_Average*Normal[iDim];
//        if (geometry->node[iPoint]->GetDomain())
//          node[iPoint]->AddGradient_Secondary(iVar, iDim, Partial_Res);
//        if (geometry->node[jPoint]->GetDomain())
//          node[jPoint]->SubtractGradient_Secondary(iVar, iDim, Partial_Res);
//      }
//    }
//  }
//
//  /*--- Loop boundary edges ---*/
//  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
//    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
//      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
//      if (geometry->node[iPoint]->GetDomain()) {
//
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          SecondaryVar_Vertex[iVar] = node[iPoint]->GetSecondary(iVar);
//
//        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++) {
//            Partial_Res = SecondaryVar_Vertex[iVar]*Normal[iDim];
//            node[iPoint]->SubtractGradient_Secondary(iVar, iDim, Partial_Res);
//          }
//      }
//    }
//  }
//
//  /*--- Update gradient value ---*/
//  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Partial_Gradient = node[iPoint]->GetGradient_Secondary(iVar, iDim) / (geometry->node[iPoint]->GetVolume());
//        node[iPoint]->SetGradient_Secondary(iVar, iDim, Partial_Gradient);
//      }
//    }
//  }
//
//  delete [] SecondaryVar_Vertex;
//  delete [] SecondaryVar_i;
//  delete [] SecondaryVar_j;
//
//  Set_MPI_Secondary_Gradient(geometry, config);
//
//}

//void CEulerSolver::SetSecondary_Gradient_LS(CGeometry *geometry, CConfig *config) {
//
//  unsigned short iVar, iDim, jDim, iNeigh;
//  unsigned long iPoint, jPoint;
//  su2double *SecondaryVar_i, *SecondaryVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
//  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
//  bool singular;
//
//  /*--- Loop over points of the grid ---*/
//
//  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
//
//    /*--- Set the value of the singular ---*/
//    singular = false;
//
//    /*--- Get coordinates ---*/
//
//    Coord_i = geometry->node[iPoint]->GetCoord();
//
//    /*--- Get Secondarys from CVariable ---*/
//
//    SecondaryVar_i = node[iPoint]->GetSecondary();
//
//    /*--- Inizialization of variables ---*/
//
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//      for (iDim = 0; iDim < nDim; iDim++)
//        cvector[iVar][iDim] = 0.0;
//
//    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
//    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0; detR2 = 0.0;
//
//    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
//      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
//      Coord_j = geometry->node[jPoint]->GetCoord();
//
//      SecondaryVar_j = node[jPoint]->GetSecondary();
//
//      weight = 0.0;
//      for (iDim = 0; iDim < nDim; iDim++)
//        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
//
//      /*--- Sumations for entries of upper triangular matrix R ---*/
//
//      if (weight != 0.0) {
//
//        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
//        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
//        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
//
//        if (nDim == 3) {
//          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
//          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
//          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
//          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
//        }
//
//        /*--- Entries of c:= transpose(A)*b ---*/
//
//        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(SecondaryVar_j[iVar]-SecondaryVar_i[iVar])/weight;
//
//      }
//
//    }
//
//    /*--- Entries of upper triangular matrix R ---*/
//
//    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
//    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
//    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
//
//    if (nDim == 3) {
//      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
//      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
//      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
//    }
//
//    /*--- Compute determinant ---*/
//
//    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
//    else detR2 = (r11*r22*r33)*(r11*r22*r33);
//
//    /*--- Detect singular matrices ---*/
//
//    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
//
//    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
//
//    if (singular) {
//      for (iDim = 0; iDim < nDim; iDim++)
//        for (jDim = 0; jDim < nDim; jDim++)
//          Smatrix[iDim][jDim] = 0.0;
//    }
//    else {
//      if (nDim == 2) {
//        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
//        Smatrix[0][1] = -r11*r12/detR2;
//        Smatrix[1][0] = Smatrix[0][1];
//        Smatrix[1][1] = r11*r11/detR2;
//      }
//      else {
//        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
//        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
//        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
//        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
//        Smatrix[0][2] = (z13*z33)/detR2;
//        Smatrix[1][0] = Smatrix[0][1];
//        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
//        Smatrix[1][2] = (z23*z33)/detR2;
//        Smatrix[2][0] = Smatrix[0][2];
//        Smatrix[2][1] = Smatrix[1][2];
//        Smatrix[2][2] = (z33*z33)/detR2;
//      }
//    }
//
//    /*--- Computation of the gradient: S*c ---*/
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        product = 0.0;
//        for (jDim = 0; jDim < nDim; jDim++) {
//          product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
//        }
//
//        node[iPoint]->SetGradient_Secondary(iVar, iDim, product);
//      }
//    }
//
//  }
//
//  Set_MPI_Secondary_Gradient(geometry, config);
//
//}

//void CEulerSolver::SetSecondary_Limiter(CGeometry *geometry, CConfig *config) {
//
//  unsigned long iEdge, iPoint, jPoint;
//  unsigned short iVar, iDim;
//  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Secondary_i, *Secondary_j,
//  dave, LimK, eps2, dm, dp, du, limiter;
//
//  /*--- Initialize solution max and solution min in the entire domain --*/
//  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      node[iPoint]->SetSolution_Max(iVar, -EPS);
//      node[iPoint]->SetSolution_Min(iVar, EPS);
//    }
//  }
//
//  /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
//  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//
//    /*--- Point identification, Normal vector and area ---*/
//    iPoint = geometry->edge[iEdge]->GetNode(0);
//    jPoint = geometry->edge[iEdge]->GetNode(1);
//
//    /*--- Get the conserved variables ---*/
//    Secondary_i = node[iPoint]->GetSecondary();
//    Secondary_j = node[jPoint]->GetSecondary();
//
//    /*--- Compute the maximum, and minimum values for nodes i & j ---*/
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      du = (Secondary_j[iVar] - Secondary_i[iVar]);
//      node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
//      node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
//      node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
//      node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
//    }
//  }
//
//  /*--- Initialize the limiter --*/
//  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
//    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//      node[iPoint]->SetLimiter_Secondary(iVar, 2.0);
//    }
//  }
//
//  /*--- Venkatakrishnan limiter ---*/
//
//  if (config->GetKind_SlopeLimit() == VENKATAKRISHNAN) {
//
//    /*-- Get limiter parameters from the configuration file ---*/
//    dave = config->GetRefElemLength();
//    LimK = config->GetLimiterCoeff();
//    eps2 = pow((LimK*dave), 3.0);
//
//    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//
//      iPoint     = geometry->edge[iEdge]->GetNode(0);
//      jPoint     = geometry->edge[iEdge]->GetNode(1);
//      Gradient_i = node[iPoint]->GetGradient_Secondary();
//      Gradient_j = node[jPoint]->GetGradient_Secondary();
//      Coord_i    = geometry->node[iPoint]->GetCoord();
//      Coord_j    = geometry->node[jPoint]->GetCoord();
//
//      for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
//
//        /*--- Calculate the interface left gradient, delta- (dm) ---*/
//        dm = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
//
//        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
//        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
//        else dp = node[iPoint]->GetSolution_Min(iVar);
//
//        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
//
//        if (limiter < node[iPoint]->GetLimiter_Secondary(iVar))
//          if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiter_Secondary(iVar, limiter);
//
//        /*-- Repeat for point j on the edge ---*/
//        dm = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
//
//        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
//        else dp = node[jPoint]->GetSolution_Min(iVar);
//
//        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
//
//        if (limiter < node[jPoint]->GetLimiter_Secondary(iVar))
//          if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiter_Secondary(iVar, limiter);
//      }
//    }
//  }
//
//  /*--- Limiter MPI ---*/
//  Set_MPI_Secondary_Limiter(geometry, config);
//
//}

void CEulerSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) {
  unsigned short iDim, jDim, iVar, jVar;
  su2double Beta, local_Mach, Beta2, rho, enthalpy, soundspeed, sq_vel;
  su2double *U_i = NULL;
  su2double Beta_min = config->GetminTurkelBeta();
  su2double Beta_max = config->GetmaxTurkelBeta();
  
  
  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();
  Beta 		    = max(Beta_min, min(local_Mach, Beta_max));
  Beta2 		    = Beta*Beta;
  
  U_i = node[iPoint]->GetSolution();
  
  rho = U_i[0];
  enthalpy = node[iPoint]->GetEnthalpy();
  soundspeed = node[iPoint]->GetSoundSpeed();
  sq_vel = node[iPoint]->GetVelocity2();
  
  /*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
  LowMach_Precontioner[0][0] = 0.5*sq_vel;
  LowMach_Precontioner[0][nVar-1] = 1.0;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[0][1+iDim] = -1.0*U_i[iDim+1]/rho;
  
  for (iDim = 0; iDim < nDim; iDim ++) {
    LowMach_Precontioner[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
    LowMach_Precontioner[iDim+1][nVar-1] = U_i[iDim+1]/rho;
    for (jDim = 0; jDim < nDim; jDim ++) {
      LowMach_Precontioner[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
    }
  }
  
  LowMach_Precontioner[nVar-1][0] = 0.5*sq_vel*enthalpy;
  LowMach_Precontioner[nVar-1][nVar-1] = enthalpy;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;
  
  
  for (iVar = 0; iVar < nVar; iVar ++ ) {
    for (jVar = 0; jVar < nVar; jVar ++ ) {
      LowMach_Precontioner[iVar][jVar] = (1.0/(Beta2+EPS) - 1.0) * (Gamma-1.0)/(soundspeed*soundspeed)*LowMach_Precontioner[iVar][jVar];
      if (iVar == jVar)
        LowMach_Precontioner[iVar][iVar] += 1.0;
    }
  }
  
}

void CEulerSolver::GetEngine_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iDim, iMarker, iMarker_EngineInflow, iMarker_EngineBleed, iMarker_EngineExhaust, iVar;
  unsigned long iVertex, iPoint;
  su2double Pressure, Temperature, Velocity[3], Velocity2, MassFlow, Density, Energy, Area,
  Mach, SoundSpeed, Flow_Dir[3], alpha;
  
  su2double Gas_Constant                  = config->GetGas_ConstantND();
  unsigned short nMarker_EngineInflow  = config->GetnMarker_EngineInflow();
  unsigned short nMarker_EngineBleed   = config->GetnMarker_EngineBleed();
  unsigned short nMarker_EngineExhaust = config->GetnMarker_EngineExhaust();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    Inflow_MassFlow[iMarker] = 0.0;
    Inflow_Mach[iMarker] = 0.0;
    Inflow_Pressure[iMarker] = 0.0;
    Inflow_Area[iMarker] = 0.0;
    
    Bleed_MassFlow[iMarker] = 0.0;
    Bleed_Pressure[iMarker] = 0.0;
    Bleed_Temperature[iMarker] = 0.0;
    Bleed_Area[iMarker] = 0.0;
    
    Exhaust_MassFlow[iMarker] = 0.0;
    Exhaust_Pressure[iMarker] = 0.0;
    Exhaust_Temperature[iMarker] = 0.0;
    Exhaust_Area[iMarker] = 0.0;
    
    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          Density = node[iPoint]->GetSolution(0);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Area += Vector[iDim]*Vector[iDim];
            Velocity[iDim] = node[iPoint]->GetSolution(iDim+1)/Density;
            Velocity2 += Velocity[iDim]*Velocity[iDim];
            MassFlow -= Vector[iDim]*node[iPoint]->GetSolution(iDim+1);
          }
          
          Area       = sqrt (Area);
          Energy     = node[iPoint]->GetSolution(nVar-1)/Density;
          Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
          SoundSpeed = sqrt(Gamma*Pressure/Density);
          Mach       = sqrt(Velocity2)/SoundSpeed;
          
          /*--- Compute the Inflow_MassFlow, Inflow_Pressure, Inflow_Mach, and Inflow_Area ---*/
          
          Inflow_MassFlow[iMarker] += MassFlow;
          Inflow_Pressure[iMarker] += Pressure*Area;
          Inflow_Mach[iMarker] += Mach*Area;
          Inflow_Area[iMarker] += Area;
          
        }
      }
      
    }
    
    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_BLEED) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          Density = node[iPoint]->GetSolution(0);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Area += Vector[iDim]*Vector[iDim];
            Velocity[iDim] = node[iPoint]->GetSolution(iDim+1)/Density;
            Velocity2 += Velocity[iDim]*Velocity[iDim];
            MassFlow += Vector[iDim]*node[iPoint]->GetSolution(iDim+1);
          }
          
          Area       = sqrt (Area);
          Energy     = node[iPoint]->GetSolution(nVar-1)/Density;
          Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
          Temperature = Pressure / (Gas_Constant * Density);
          
          /*--- Compute the Bleed_MassFlow, Bleed_Pressure, Bleed_Temperature, and Bleed_Area ---*/
          
          Bleed_MassFlow[iMarker] += MassFlow;
          Bleed_Pressure[iMarker] += Pressure*Area;
          Bleed_Temperature[iMarker] += Temperature*Area;
          Bleed_Area[iMarker] += Area;
          
        }
      }
      
    }
    
    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          Density = node[iPoint]->GetSolution(0);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Area += Vector[iDim]*Vector[iDim];
            Velocity[iDim] = node[iPoint]->GetSolution(iDim+1)/Density;
            Velocity2 += Velocity[iDim]*Velocity[iDim];
            MassFlow += Vector[iDim]*node[iPoint]->GetSolution(iDim+1);
          }
          
          Area       = sqrt (Area);
          Energy     = node[iPoint]->GetSolution(nVar-1)/Density;
          Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
          Temperature = Pressure / (Gas_Constant * Density);
          
          /*--- Compute the mass Exhaust_MassFlow ---*/
          
          Exhaust_MassFlow[iMarker] += MassFlow;
          Exhaust_Pressure[iMarker] += Pressure*Area;
          Exhaust_Temperature[iMarker] += Temperature*Area;
          Exhaust_Area[iMarker] += Area;
          
        }
      }
      
    }
    
  }
  
  /*--- Copy to the appropriate structure ---*/
  
  su2double *Inflow_MassFlow_Local = new su2double [nMarker_EngineInflow];
  su2double *Inflow_Mach_Local = new su2double [nMarker_EngineInflow];
  su2double *Inflow_Pressure_Local = new su2double [nMarker_EngineInflow];
  su2double *Inflow_Area_Local = new su2double [nMarker_EngineInflow];
  
  su2double *Inflow_MassFlow_Total = new su2double [nMarker_EngineInflow];
  su2double *Inflow_Mach_Total = new su2double [nMarker_EngineInflow];
  su2double *Inflow_Pressure_Total = new su2double [nMarker_EngineInflow];
  su2double *Inflow_Area_Total = new su2double [nMarker_EngineInflow];
  
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
    Inflow_MassFlow_Local[iMarker_EngineInflow] = 0.0;
    Inflow_Mach_Local[iMarker_EngineInflow] = 0.0;
    Inflow_Pressure_Local[iMarker_EngineInflow] = 0.0;
    Inflow_Area_Local[iMarker_EngineInflow] = 0.0;
    
    Inflow_MassFlow_Total[iMarker_EngineInflow] = 0.0;
    Inflow_Mach_Total[iMarker_EngineInflow] = 0.0;
    Inflow_Pressure_Total[iMarker_EngineInflow] = 0.0;
    Inflow_Area_Total[iMarker_EngineInflow] = 0.0;
  }
  
  su2double *Bleed_MassFlow_Local = new su2double [nMarker_EngineBleed];
  su2double *Bleed_Temperature_Local = new su2double [nMarker_EngineBleed];
  su2double *Bleed_Pressure_Local = new su2double [nMarker_EngineBleed];
  su2double *Bleed_Area_Local = new su2double [nMarker_EngineBleed];
  
  su2double *Bleed_MassFlow_Total = new su2double [nMarker_EngineBleed];
  su2double *Bleed_Temperature_Total = new su2double [nMarker_EngineBleed];
  su2double *Bleed_Pressure_Total = new su2double [nMarker_EngineBleed];
  su2double *Bleed_Area_Total = new su2double [nMarker_EngineBleed];
  
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++) {
    Bleed_MassFlow_Local[iMarker_EngineBleed] = 0.0;
    Bleed_Temperature_Local[iMarker_EngineBleed] = 0.0;
    Bleed_Pressure_Local[iMarker_EngineBleed] = 0.0;
    Bleed_Area_Local[iMarker_EngineBleed] = 0.0;
    
    Bleed_MassFlow_Total[iMarker_EngineBleed] = 0.0;
    Bleed_Temperature_Total[iMarker_EngineBleed] = 0.0;
    Bleed_Pressure_Total[iMarker_EngineBleed] = 0.0;
    Bleed_Area_Total[iMarker_EngineBleed] = 0.0;
  }
  
  su2double *Exhaust_MassFlow_Local = new su2double [nMarker_EngineExhaust];
  su2double *Exhaust_Temperature_Local = new su2double [nMarker_EngineExhaust];
  su2double *Exhaust_Pressure_Local = new su2double [nMarker_EngineExhaust];
  su2double *Exhaust_Area_Local = new su2double [nMarker_EngineExhaust];
  
  su2double *Exhaust_MassFlow_Total = new su2double [nMarker_EngineExhaust];
  su2double *Exhaust_Temperature_Total = new su2double [nMarker_EngineExhaust];
  su2double *Exhaust_Pressure_Total = new su2double [nMarker_EngineExhaust];
  su2double *Exhaust_Area_Total = new su2double [nMarker_EngineExhaust];
  
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
    Exhaust_MassFlow_Local[iMarker_EngineExhaust] = 0.0;
    Exhaust_Temperature_Local[iMarker_EngineExhaust] = 0.0;
    Exhaust_Pressure_Local[iMarker_EngineExhaust] = 0.0;
    Exhaust_Area_Local[iMarker_EngineExhaust] = 0.0;
    
    Exhaust_MassFlow_Total[iMarker_EngineExhaust] = 0.0;
    Exhaust_Temperature_Total[iMarker_EngineExhaust] = 0.0;
    Exhaust_Pressure_Total[iMarker_EngineExhaust] = 0.0;
    Exhaust_Area_Total[iMarker_EngineExhaust] = 0.0;
  }
  
  /*--- Compute the numerical fan face Mach number, mach number, temperature and the total area ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW) {
      
      for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
        
        /*--- Add the Inflow_MassFlow, Inflow_Mach, Inflow_Pressure and Inflow_Area to the particular boundary ---*/
        
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_EngineInflow(iMarker_EngineInflow)) {
          Inflow_MassFlow_Local[iMarker_EngineInflow] += Inflow_MassFlow[iMarker];
          Inflow_Mach_Local[iMarker_EngineInflow] += Inflow_Mach[iMarker];
          Inflow_Pressure_Local[iMarker_EngineInflow] += Inflow_Pressure[iMarker];
          Inflow_Area_Local[iMarker_EngineInflow] += Inflow_Area[iMarker];
        }
        
      }
      
    }
    
    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_BLEED) {
      
      for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++) {
        
        /*--- Add the Bleed_MassFlow, Bleed_Temperature, Bleed_Pressure and Bleed_Area to the particular boundary ---*/
        
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_EngineBleed(iMarker_EngineBleed)) {
          Bleed_MassFlow_Local[iMarker_EngineBleed] += Bleed_MassFlow[iMarker];
          Bleed_Temperature_Local[iMarker_EngineBleed] += Bleed_Temperature[iMarker];
          Bleed_Pressure_Local[iMarker_EngineBleed] += Bleed_Pressure[iMarker];
          Bleed_Area_Local[iMarker_EngineBleed] += Bleed_Area[iMarker];
        }
        
      }
      
    }
    
    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST) {
      
      for (iMarker_EngineExhaust= 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
        
        /*--- Add the Exhaust_MassFlow, and Exhaust_Area to the particular boundary ---*/
        
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_EngineExhaust(iMarker_EngineExhaust)) {
          Exhaust_MassFlow_Local[iMarker_EngineExhaust] += Exhaust_MassFlow[iMarker];
          Exhaust_Temperature_Local[iMarker_EngineExhaust] += Exhaust_Temperature[iMarker];
          Exhaust_Pressure_Local[iMarker_EngineExhaust] += Exhaust_Pressure[iMarker];
          Exhaust_Area_Local[iMarker_EngineExhaust] += Exhaust_Area[iMarker];
        }
        
      }
      
    }
    
  }
  
#ifdef HAVE_MPI
  
  SU2_MPI::Allreduce(Inflow_MassFlow_Local, Inflow_MassFlow_Total, nMarker_EngineInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Inflow_Mach_Local, Inflow_Mach_Total, nMarker_EngineInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Inflow_Pressure_Local, Inflow_Pressure_Total, nMarker_EngineInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Inflow_Area_Local, Inflow_Area_Total, nMarker_EngineInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  SU2_MPI::Allreduce(Bleed_MassFlow_Local, Bleed_MassFlow_Total, nMarker_EngineBleed, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Bleed_Temperature_Local, Bleed_Temperature_Total, nMarker_EngineBleed, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Bleed_Pressure_Local, Bleed_Pressure_Total, nMarker_EngineBleed, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Bleed_Area_Local, Bleed_Area_Total, nMarker_EngineBleed, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  SU2_MPI::Allreduce(Exhaust_MassFlow_Local, Exhaust_MassFlow_Total, nMarker_EngineExhaust, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Exhaust_Temperature_Local, Exhaust_Temperature_Total, nMarker_EngineExhaust, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Exhaust_Pressure_Local, Exhaust_Pressure_Total, nMarker_EngineExhaust, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Exhaust_Area_Local, Exhaust_Area_Total, nMarker_EngineExhaust, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#else
  
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
    Inflow_MassFlow_Total[iMarker_EngineInflow]   = Inflow_MassFlow_Local[iMarker_EngineInflow];
    Inflow_Mach_Total[iMarker_EngineInflow]       = Inflow_Mach_Local[iMarker_EngineInflow];
    Inflow_Pressure_Total[iMarker_EngineInflow]   = Inflow_Pressure_Local[iMarker_EngineInflow];
    Inflow_Area_Total[iMarker_EngineInflow]       = Inflow_Area_Local[iMarker_EngineInflow];
  }
  
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++) {
    Bleed_MassFlow_Total[iMarker_EngineBleed]    = Bleed_MassFlow_Local[iMarker_EngineBleed];
    Bleed_Temperature_Total[iMarker_EngineBleed] = Bleed_Temperature_Local[iMarker_EngineBleed];
    Bleed_Pressure_Total[iMarker_EngineBleed]    = Bleed_Pressure_Local[iMarker_EngineBleed];
    Bleed_Area_Total[iMarker_EngineBleed]        = Bleed_Area_Local[iMarker_EngineBleed];
  }
  
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
    Exhaust_MassFlow_Total[iMarker_EngineExhaust]  = Exhaust_MassFlow_Local[iMarker_EngineExhaust];
    Exhaust_Temperature_Total[iMarker_EngineExhaust] = Exhaust_Temperature_Local[iMarker_EngineExhaust];
    Exhaust_Pressure_Total[iMarker_EngineExhaust]   = Exhaust_Pressure_Local[iMarker_EngineExhaust];
    Exhaust_Area_Total[iMarker_EngineExhaust]      = Exhaust_Area_Local[iMarker_EngineExhaust];
  }
  
#endif
  
  /*--- Compute the value of Inflow_Area_Total, and Inflow_Pressure_Total, and
   set the value in the config structure for future use ---*/
  
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
    if (Inflow_Area_Total[iMarker_EngineInflow] != 0.0) Inflow_Mach_Total[iMarker_EngineInflow] /= Inflow_Area_Total[iMarker_EngineInflow];
    else Inflow_Mach_Total[iMarker_EngineInflow] = 0.0;
    if (Inflow_Area_Total[iMarker_EngineInflow] != 0.0) Inflow_Pressure_Total[iMarker_EngineInflow] /= Inflow_Area_Total[iMarker_EngineInflow];
    else Inflow_Pressure_Total[iMarker_EngineInflow] = 0.0;
    
    if (iMesh == MESH_0) {
      config->SetInflow_Mach(iMarker_EngineInflow, Inflow_Mach_Total[iMarker_EngineInflow]);
      config->SetInflow_Pressure(iMarker_EngineInflow, Inflow_Pressure_Total[iMarker_EngineInflow]);
    }
    
  }
  
  /*--- Compute the value of Bleed_Area_Total, and Bleed_Pressure_Total, and
   set the value in the config structure for future use ---*/
  
  for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++) {
    if (Bleed_Area_Total[iMarker_EngineBleed] != 0.0) Bleed_Temperature_Total[iMarker_EngineBleed] /= Bleed_Area_Total[iMarker_EngineBleed];
    else Bleed_Temperature_Total[iMarker_EngineBleed] = 0.0;
    if (Bleed_Area_Total[iMarker_EngineBleed] != 0.0) Bleed_Pressure_Total[iMarker_EngineBleed] /= Bleed_Area_Total[iMarker_EngineBleed];
    else Bleed_Pressure_Total[iMarker_EngineBleed] = 0.0;
    
    if (iMesh == MESH_0) {
      config->SetBleed_Temperature(iMarker_EngineBleed, Bleed_Temperature_Total[iMarker_EngineBleed]);
      config->SetBleed_MassFlow(iMarker_EngineBleed, Bleed_MassFlow_Total[iMarker_EngineBleed]);
      config->SetBleed_Pressure(iMarker_EngineBleed, Bleed_Pressure_Total[iMarker_EngineBleed]);
    }
    
  }
  
  /*--- Compute the value of Exhaust_Area_Total, and Exhaust_Pressure_Total, and
   set the value in the config structure for future use ---*/
  
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
    if (Exhaust_Area_Total[iMarker_EngineExhaust] != 0.0) Exhaust_Temperature_Total[iMarker_EngineExhaust] /= Exhaust_Area_Total[iMarker_EngineExhaust];
    else Exhaust_Temperature_Total[iMarker_EngineExhaust] = 0.0;
    if (Exhaust_Area_Total[iMarker_EngineExhaust] != 0.0) Exhaust_Pressure_Total[iMarker_EngineExhaust] /= Exhaust_Area_Total[iMarker_EngineExhaust];
    else Exhaust_Pressure_Total[iMarker_EngineExhaust] = 0.0;
    
    if (iMesh == MESH_0) {
      config->SetExhaust_Temperature(iMarker_EngineExhaust, Exhaust_Temperature_Total[iMarker_EngineExhaust]);
      config->SetExhaust_Pressure(iMarker_EngineExhaust, Exhaust_Pressure_Total[iMarker_EngineExhaust]);
    }
    
  }
  
  bool write_heads = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output) {
    
    cout.precision(4);
    cout.setf(ios::fixed, ios::floatfield);
    
    cout << endl << "---------------------------- Engine properties --------------------------" << endl;
    for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
      cout << "Engine inflow ("<< config->GetMarker_EngineInflow(iMarker_EngineInflow);
      if (config->GetSystemMeasurements() == SI) cout << "): Mass flow (kg/s): ";
      else if (config->GetSystemMeasurements() == US) cout << "): Mass flow (slug/s): ";
      cout << Inflow_MassFlow_Total[iMarker_EngineInflow] * config->GetDensity_Ref() * config->GetVelocity_Ref()
      << ", Mach: " << Inflow_Mach_Total[iMarker_EngineInflow];
      if (config->GetSystemMeasurements() == SI) cout << ", Area (m^2): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Area (ft^2): ";
      cout << Inflow_Area_Total[iMarker_EngineInflow] <<"."<< endl;
    }
    
    for (iMarker_EngineBleed = 0; iMarker_EngineBleed < nMarker_EngineBleed; iMarker_EngineBleed++) {
      cout << "Engine bleed ("<< config->GetMarker_EngineBleed(iMarker_EngineBleed);
      if (config->GetSystemMeasurements() == SI) cout << "): Mass flow (kg/s): ";
      else if (config->GetSystemMeasurements() == US) cout << "): Mass flow (slug/s): ";
      cout << Bleed_MassFlow_Total[iMarker_EngineBleed] * config->GetDensity_Ref() * config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Temp (K): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Temp (R): ";
      cout << Bleed_Temperature_Total[iMarker_EngineBleed] * config->GetTemperature_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Area (m^2): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Area (ft^2): ";
      cout << Bleed_Area_Total[iMarker_EngineBleed] <<"."<< endl;
    }
    
    for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
      cout << "Engine exhaust ("<< config->GetMarker_EngineExhaust(iMarker_EngineExhaust);
      if (config->GetSystemMeasurements() == SI) cout << "): Mass flow (kg/s): ";
      else if (config->GetSystemMeasurements() == US) cout << "): Mass flow (slug/s): ";
      cout << Exhaust_MassFlow_Total[iMarker_EngineExhaust] * config->GetDensity_Ref() * config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Temp (K): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Temp (R): ";
      cout << Exhaust_Temperature_Total[iMarker_EngineExhaust] * config->GetTemperature_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Pressure (Pa): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Pressure (psf): ";
      cout << Exhaust_Pressure_Total[iMarker_EngineExhaust] * config->GetPressure_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Area (m^2): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Area (ft^2): ";
      cout << Exhaust_Area_Total[iMarker_EngineExhaust] <<"."<< endl;
    }
    cout << "-------------------------------------------------------------------------" << endl;
    
  }
  
  /*--- Check the flow orientation in the engine ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW) ||
        (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST) ||
        (config->GetMarker_All_KindBC(iMarker) == ENGINE_BLEED)) {
      
      /*--- Loop over all the vertices on this boundary marker ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Normal vector for this vertex (negate for outward convention) ---*/
        
        geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
        
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
        
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Area += Vector[iDim]*Vector[iDim];
        Area = sqrt (Area);
        
        /*--- Compute unitary vector ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Vector[iDim] /= Area;
        
        /*--- The flow direction is defined by the local velocity on the surface ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Flow_Dir[iDim] = node[iPoint]->GetSolution(iDim+1) / node[iPoint]->GetSolution(0);
        
        /*--- Dot product of normal and flow direction. ---*/
        
        alpha = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          alpha += Vector[iDim]*Flow_Dir[iDim];
        
        /*--- Flow in the wrong direction. ---*/
        
        if (((config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST) ||
             (config->GetMarker_All_KindBC(iMarker) == ENGINE_BLEED)) && (alpha > 0.0)) {
          
          /*--- Copy the old solution ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_Old(iVar));
          
        }
        
        if ((config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW) && (alpha < 0.0)) {
          
          /*--- Copy the old solution ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_Old(iVar));
          
        }
        
      }
      
    }
    
  }
  
  
  delete [] Inflow_MassFlow_Local;
  delete [] Inflow_Mach_Local;
  delete [] Inflow_Pressure_Local;
  delete [] Inflow_Area_Local;
  
  delete [] Inflow_MassFlow_Total;
  delete [] Inflow_Mach_Total;
  delete [] Inflow_Pressure_Total;
  delete [] Inflow_Area_Total;
  
  delete [] Exhaust_MassFlow_Local;
  delete [] Exhaust_Temperature_Local;
  delete [] Exhaust_Pressure_Local;
  delete [] Exhaust_Area_Local;
  
  delete [] Exhaust_MassFlow_Total;
  delete [] Exhaust_Temperature_Total;
  delete [] Exhaust_Pressure_Total;
  delete [] Exhaust_Area_Total;
  
  delete [] Bleed_MassFlow_Local;
  delete [] Bleed_Temperature_Local;
  delete [] Bleed_Pressure_Local;
  delete [] Bleed_Area_Local;
  
  delete [] Bleed_MassFlow_Total;
  delete [] Bleed_Temperature_Total;
  delete [] Bleed_Pressure_Total;
  delete [] Bleed_Area_Total;
  
}

void CEulerSolver::GetActuatorDisk_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iDim, iMarker;
  unsigned long iVertex, iPoint;
  su2double Pressure, Temperature, Velocity[3], Velocity2, MassFlow, Density, Energy, Area;
  unsigned short iMarker_ActDiskInlet, iMarker_ActDiskOutlet;
  
  su2double Gas_Constant = config->GetGas_ConstantND();
  
  unsigned short nMarker_ActDiskInlet = config->GetnMarker_ActDisk_Inlet();
  unsigned short nMarker_ActDiskOutlet = config->GetnMarker_ActDisk_Outlet();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  su2double *Inlet_MassFlow = new su2double [config->GetnMarker_All()];
  su2double *Inlet_Pressure = new su2double [config->GetnMarker_All()];
  su2double *Inlet_Temperature = new su2double [config->GetnMarker_All()];
  su2double *Inlet_Area = new su2double [config->GetnMarker_All()];
  
  su2double *Outlet_MassFlow = new su2double [config->GetnMarker_All()];
  su2double *Outlet_Pressure = new su2double [config->GetnMarker_All()];
  su2double *Outlet_Temperature = new su2double [config->GetnMarker_All()];
  su2double *Outlet_Area = new su2double [config->GetnMarker_All()];
  
  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    Inlet_MassFlow[iMarker] = 0.0;
    Inlet_Pressure[iMarker] = 0.0;
    Inlet_Temperature[iMarker] = 0.0;
    Inlet_Area[iMarker] = 0.0;
    
    Outlet_MassFlow[iMarker] = 0.0;
    Outlet_Pressure[iMarker] = 0.0;
    Outlet_Temperature[iMarker] = 0.0;
    Outlet_Area[iMarker] = 0.0;
    
    if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          Density = node[iPoint]->GetSolution(0);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Area += Vector[iDim]*Vector[iDim];
            Velocity[iDim] = node[iPoint]->GetSolution(iDim+1)/Density;
            Velocity2 += Velocity[iDim]*Velocity[iDim];
            MassFlow += Vector[iDim]*node[iPoint]->GetSolution(iDim+1);
          }
          
          Area       = sqrt (Area);
          Energy     = node[iPoint]->GetSolution(nVar-1)/Density;
          Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
          Temperature = Pressure / (Gas_Constant * Density);
          
          /*--- Compute the Inlet_MassFlow, Inlet_Pressure, Inlet_Temperature, and Inlet_Area ---*/
          
          Inlet_MassFlow[iMarker] += MassFlow;
          Inlet_Pressure[iMarker] += Pressure*Area;
          Inlet_Temperature[iMarker] += Temperature*Area;
          Inlet_Area[iMarker] += Area;
          
        }
      }
      
    }
    
    if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          Density = node[iPoint]->GetSolution(0);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Area += Vector[iDim]*Vector[iDim];
            Velocity[iDim] = node[iPoint]->GetSolution(iDim+1)/Density;
            Velocity2 += Velocity[iDim]*Velocity[iDim];
            MassFlow += Vector[iDim]*node[iPoint]->GetSolution(iDim+1);
          }
          
          Area       = sqrt (Area);
          Energy     = node[iPoint]->GetSolution(nVar-1)/Density;
          Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
          Temperature = Pressure / (Gas_Constant * Density);
          
          /*--- Compute the mass Outlet_MassFlow ---*/
          
          Outlet_MassFlow[iMarker] += MassFlow;
          Outlet_Pressure[iMarker] += Pressure*Area;
          Outlet_Temperature[iMarker] += Temperature*Area;
          Outlet_Area[iMarker] += Area;
          
        }
      }
      
    }
    
  }
  
  /*--- Copy to the appropriate structure ---*/
  
  su2double *Inlet_MassFlow_Local = new su2double [nMarker_ActDiskInlet];
  su2double *Inlet_Temperature_Local = new su2double [nMarker_ActDiskInlet];
  su2double *Inlet_Pressure_Local = new su2double [nMarker_ActDiskInlet];
  su2double *Inlet_Area_Local = new su2double [nMarker_ActDiskInlet];
  
  su2double *Inlet_MassFlow_Total = new su2double [nMarker_ActDiskInlet];
  su2double *Inlet_Temperature_Total = new su2double [nMarker_ActDiskInlet];
  su2double *Inlet_Pressure_Total = new su2double [nMarker_ActDiskInlet];
  su2double *Inlet_Area_Total = new su2double [nMarker_ActDiskInlet];
  
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
    Inlet_MassFlow_Local[iMarker_ActDiskInlet] = 0.0;
    Inlet_Temperature_Local[iMarker_ActDiskInlet] = 0.0;
    Inlet_Pressure_Local[iMarker_ActDiskInlet] = 0.0;
    Inlet_Area_Local[iMarker_ActDiskInlet] = 0.0;
    
    Inlet_MassFlow_Total[iMarker_ActDiskInlet] = 0.0;
    Inlet_Temperature_Total[iMarker_ActDiskInlet] = 0.0;
    Inlet_Pressure_Total[iMarker_ActDiskInlet] = 0.0;
    Inlet_Area_Total[iMarker_ActDiskInlet] = 0.0;
  }
  
  su2double *Outlet_MassFlow_Local = new su2double [nMarker_ActDiskOutlet];
  su2double *Outlet_Temperature_Local = new su2double [nMarker_ActDiskOutlet];
  su2double *Outlet_Pressure_Local = new su2double [nMarker_ActDiskOutlet];
  su2double *Outlet_Area_Local = new su2double [nMarker_ActDiskOutlet];
  
  su2double *Outlet_MassFlow_Total = new su2double [nMarker_ActDiskOutlet];
  su2double *Outlet_Temperature_Total = new su2double [nMarker_ActDiskOutlet];
  su2double *Outlet_Pressure_Total = new su2double [nMarker_ActDiskOutlet];
  su2double *Outlet_Area_Total = new su2double [nMarker_ActDiskOutlet];
  
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
    Outlet_MassFlow_Local[iMarker_ActDiskOutlet] = 0.0;
    Outlet_Temperature_Local[iMarker_ActDiskOutlet] = 0.0;
    Outlet_Pressure_Local[iMarker_ActDiskOutlet] = 0.0;
    Outlet_Area_Local[iMarker_ActDiskOutlet] = 0.0;
    
    Outlet_MassFlow_Total[iMarker_ActDiskOutlet] = 0.0;
    Outlet_Temperature_Total[iMarker_ActDiskOutlet] = 0.0;
    Outlet_Pressure_Total[iMarker_ActDiskOutlet] = 0.0;
    Outlet_Area_Total[iMarker_ActDiskOutlet] = 0.0;
  }
  
  /*--- Compute the numerical fan face Mach number, mach number, temperature and the total area ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) {
      
      for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
        
        /*--- Add the Inlet_MassFlow, Inlet_Temperature, Inlet_Pressure and Inlet_Area to the particular boundary ---*/
        
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_ActDisk_Inlet(iMarker_ActDiskInlet)) {
          Inlet_MassFlow_Local[iMarker_ActDiskInlet] += Inlet_MassFlow[iMarker];
          Inlet_Temperature_Local[iMarker_ActDiskInlet] += Inlet_Temperature[iMarker];
          Inlet_Pressure_Local[iMarker_ActDiskInlet] += Inlet_Pressure[iMarker];
          Inlet_Area_Local[iMarker_ActDiskInlet] += Inlet_Area[iMarker];
        }
        
      }
      
    }
    
    if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) {
      
      for (iMarker_ActDiskOutlet= 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
        
        /*--- Add the Outlet_MassFlow, and Outlet_Area to the particular boundary ---*/
        
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_ActDisk_Outlet(iMarker_ActDiskOutlet)) {
          Outlet_MassFlow_Local[iMarker_ActDiskOutlet] += Outlet_MassFlow[iMarker];
          Outlet_Temperature_Local[iMarker_ActDiskOutlet] += Outlet_Temperature[iMarker];
          Outlet_Pressure_Local[iMarker_ActDiskOutlet] += Outlet_Pressure[iMarker];
          Outlet_Area_Local[iMarker_ActDiskOutlet] += Outlet_Area[iMarker];
        }
        
      }
      
    }
    
  }
  
#ifdef HAVE_MPI
  
  SU2_MPI::Allreduce(Inlet_MassFlow_Local, Inlet_MassFlow_Total, nMarker_ActDiskInlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Inlet_Temperature_Local, Inlet_Temperature_Total, nMarker_ActDiskInlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Inlet_Pressure_Local, Inlet_Pressure_Total, nMarker_ActDiskInlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Inlet_Area_Local, Inlet_Area_Total, nMarker_ActDiskInlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  SU2_MPI::Allreduce(Outlet_MassFlow_Local, Outlet_MassFlow_Total, nMarker_ActDiskOutlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Outlet_Temperature_Local, Outlet_Temperature_Total, nMarker_ActDiskOutlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Outlet_Pressure_Local, Outlet_Pressure_Total, nMarker_ActDiskOutlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Outlet_Area_Local, Outlet_Area_Total, nMarker_ActDiskOutlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#else
  
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
    Inlet_MassFlow_Total[iMarker_ActDiskInlet]    = Inlet_MassFlow_Local[iMarker_ActDiskInlet];
    Inlet_Temperature_Total[iMarker_ActDiskInlet] = Inlet_Temperature_Local[iMarker_ActDiskInlet];
    Inlet_Pressure_Total[iMarker_ActDiskInlet]    = Inlet_Pressure_Local[iMarker_ActDiskInlet];
    Inlet_Area_Total[iMarker_ActDiskInlet]        = Inlet_Area_Local[iMarker_ActDiskInlet];
  }
  
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
    Outlet_MassFlow_Total[iMarker_ActDiskOutlet]  = Outlet_MassFlow_Local[iMarker_ActDiskOutlet];
    Outlet_Temperature_Total[iMarker_ActDiskOutlet] = Outlet_Temperature_Local[iMarker_ActDiskOutlet];
    Outlet_Pressure_Total[iMarker_ActDiskOutlet]   = Outlet_Pressure_Local[iMarker_ActDiskOutlet];
    Outlet_Area_Total[iMarker_ActDiskOutlet]      = Outlet_Area_Local[iMarker_ActDiskOutlet];
  }
  
#endif
  
  /*--- Compute the value of Inlet_Area_Total, and Inlet_Pressure_Total, and
   set the value in the config structure for future use ---*/
  
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
    if (Inlet_Area_Total[iMarker_ActDiskInlet] != 0.0) Inlet_Temperature_Total[iMarker_ActDiskInlet] /= Inlet_Area_Total[iMarker_ActDiskInlet];
    else Inlet_Temperature_Total[iMarker_ActDiskInlet] = 0.0;
    if (Inlet_Area_Total[iMarker_ActDiskInlet] != 0.0) Inlet_Pressure_Total[iMarker_ActDiskInlet] /= Inlet_Area_Total[iMarker_ActDiskInlet];
    else Inlet_Pressure_Total[iMarker_ActDiskInlet] = 0.0;
  }
  
  /*--- Compute the value of Outlet_Area_Total, and Outlet_Pressure_Total, and
   set the value in the config structure for future use ---*/
  
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
    if (Outlet_Area_Total[iMarker_ActDiskOutlet] != 0.0) Outlet_Temperature_Total[iMarker_ActDiskOutlet] /= Outlet_Area_Total[iMarker_ActDiskOutlet];
    else Outlet_Temperature_Total[iMarker_ActDiskOutlet] = 0.0;
    if (Outlet_Area_Total[iMarker_ActDiskOutlet] != 0.0) Outlet_Pressure_Total[iMarker_ActDiskOutlet] /= Outlet_Area_Total[iMarker_ActDiskOutlet];
    else Outlet_Pressure_Total[iMarker_ActDiskOutlet] = 0.0;
  }
  
  bool write_heads = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output) {
    
    cout.precision(4);
    cout.setf(ios::fixed, ios::floatfield);
    
    cout << endl << "------------------------ Actuator Disk properties -----------------------" << endl;
    
    for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
      cout << "Actuator Disk Inlet ("<< config->GetMarker_ActDisk_Inlet(iMarker_ActDiskInlet);
      if (config->GetSystemMeasurements() == SI) cout << "): Mass flow (kg/s): ";
      else if (config->GetSystemMeasurements() == US) cout << "): Mass flow (slug/s): ";
      cout << -Inlet_MassFlow_Total[iMarker_ActDiskInlet] * config->GetDensity_Ref() * config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Temp (K): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Temp (R): ";
      cout << Inlet_Temperature_Total[iMarker_ActDiskInlet] * config->GetTemperature_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Pressure (Pa): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Pressure (psf): ";
      cout << Inlet_Pressure_Total[iMarker_ActDiskInlet] * config->GetPressure_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Area (m^2): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Area (ft^2): ";
      cout << Inlet_Area_Total[iMarker_ActDiskInlet] <<"."<< endl;
    }
    
    for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
      cout << "Actuator Disk Outlet ("<< config->GetMarker_ActDisk_Outlet(iMarker_ActDiskOutlet);
      if (config->GetSystemMeasurements() == SI) cout << "): Mass flow (kg/s): ";
      else if (config->GetSystemMeasurements() == US) cout << "): Mass flow (slug/s): ";
      cout << Outlet_MassFlow_Total[iMarker_ActDiskOutlet] * config->GetDensity_Ref() * config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Temp (K): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Temp (R): ";
      cout << Outlet_Temperature_Total[iMarker_ActDiskOutlet] * config->GetTemperature_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Pressure (Pa): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Pressure (psf): ";
      cout << Outlet_Pressure_Total[iMarker_ActDiskOutlet] * config->GetPressure_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ", Area (m^2): ";
      else if (config->GetSystemMeasurements() == US) cout << ", Area (ft^2): ";
      cout << Outlet_Area_Total[iMarker_ActDiskOutlet] <<"."<< endl;
    }
    cout << "-------------------------------------------------------------------------" << endl;
    
  }
  
  delete [] Outlet_MassFlow_Local;
  delete [] Outlet_Temperature_Local;
  delete [] Outlet_Pressure_Local;
  delete [] Outlet_Area_Local;
  
  delete [] Outlet_MassFlow_Total;
  delete [] Outlet_Temperature_Total;
  delete [] Outlet_Pressure_Total;
  delete [] Outlet_Area_Total;
  
  delete [] Inlet_MassFlow_Local;
  delete [] Inlet_Temperature_Local;
  delete [] Inlet_Pressure_Local;
  delete [] Inlet_Area_Local;
  
  delete [] Inlet_MassFlow_Total;
  delete [] Inlet_Temperature_Total;
  delete [] Inlet_Pressure_Total;
  delete [] Inlet_Area_Total;
  
  
  delete [] Inlet_MassFlow;
  delete [] Inlet_Pressure;
  delete [] Inlet_Temperature;
  delete [] Inlet_Area;
  
  delete [] Outlet_MassFlow;
  delete [] Outlet_Pressure;
  delete [] Outlet_Temperature;
  delete [] Outlet_Area;
  
}

void CEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iDim, iCounter;
  bool Update_AoA = false;
  su2double Target_CL, AoA_inc, AoA;
  su2double DampingFactor = config->GetDamp_Fixed_CL();
  unsigned long Iter_Fixed_CL = config->GetIter_Fixed_CL();
  unsigned long ExtIter = config->GetExtIter();
  su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
  su2double Vel_Infty[3], Vel_Infty_Mag;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Only the fine mesh level should check the convergence criteria ---*/
  
  if (iMesh == MESH_0) {
    
    /*--- Initialize the update flag to false ---*/
    
    Update_AoA = false;
    
    /*--- Reevaluate Angle of Attack at a fix number of iterations ---*/
    
    if (ExtIter % Iter_Fixed_CL == 0) { Update_AoA = true; };
    
    /*--- Store the update boolean for use on other mesh levels in the MG ---*/
    
    config->SetUpdate_AoA(Update_AoA);
    
  }
  
  else {
    Update_AoA = config->GetUpdate_AoA();
  }
  
  /*--- If we are within two digits of convergence in the CL coefficient,
   compute an updated value for the AoA at the farfield. We are iterating
   on the AoA in order to match the specified fixed lift coefficient. ---*/
  
  if (Update_AoA) {
    
    /*--- Retrieve the specified target CL value. ---*/
    
    Target_CL = config->GetTarget_CL();
    
    /*--- Retrieve the old AoA. ---*/
    
    AoA_old = config->GetAoA()*PI_NUMBER/180.0;
    
    /*--- Estimate the increment in AoA based on a 2*pi lift curve slope ---*/
    
    AoA_inc = (1.0/(2.0*PI_NUMBER))*(Target_CL - Total_CLift);
    
    /*--- Compute a new value for AoA on the fine mesh only ---*/
    
    if (iMesh == MESH_0)
      AoA = (1.0 - DampingFactor)*AoA_old + DampingFactor * (AoA_old + AoA_inc);
    else
      AoA = config->GetAoA()*PI_NUMBER/180.0;
    
    /*--- Update the freestream velocity vector at the farfield ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty[iDim] = GetVelocity_Inf(iDim);
    
    /*--- Compute the magnitude of the free stream velocity ---*/
    
    Vel_Infty_Mag = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      Vel_Infty_Mag += Vel_Infty[iDim]*Vel_Infty[iDim];
    Vel_Infty_Mag = sqrt(Vel_Infty_Mag);
    
    /*--- Compute the new freestream velocity with the updated AoA ---*/
    
    if (nDim == 2) {
      Vel_Infty[0] = cos(AoA)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(AoA)*Vel_Infty_Mag;
    }
    if (nDim == 3) {
      Vel_Infty[0] = cos(AoA)*cos(Beta)*Vel_Infty_Mag;
      Vel_Infty[1] = sin(Beta)*Vel_Infty_Mag;
      Vel_Infty[2] = sin(AoA)*cos(Beta)*Vel_Infty_Mag;
    }
    
    /*--- Store the new freestream velocity vector for the next iteration ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_Inf[iDim] = Vel_Infty[iDim];
    }
    
    /*--- Only the fine mesh stores the updated values for AoA in config ---*/
    if (iMesh == MESH_0) {
      for (iDim = 0; iDim < nDim; iDim++)
        config->SetVelocity_FreeStreamND(Vel_Infty[iDim], iDim);
      config->SetAoA(AoA*180.0/PI_NUMBER);
    }
    
    /*--- Reset the local cauchy criteria ---*/
    Cauchy_Value = 0.0;
    Cauchy_Counter = 0;
    for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
      Cauchy_Serie[iCounter] = 0.0;
  }
  
  /*--- Output some information to the console with the headers ---*/
  
  bool write_heads = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output) {
    cout.precision(7);
    cout.setf(ios::fixed, ios::floatfield);
    cout << endl << "----------------------------- Fixed CL Mode -----------------------------" << endl;
    cout << "Target CL: " << config->GetTarget_CL();
    cout << ", Current CL: " << Total_CLift;
    cout << ", Current AoA: " << config->GetAoA() << " deg." << endl;
    cout << "-------------------------------------------------------------------------" << endl;
  }
  
  
}

void CEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar, kVar, jDim;
  unsigned long iPoint, iVertex;
  su2double Density = 0.0, Pressure = 0.0, *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3], *NormalArea,
  ProjGridVel = 0.0, turb_ke;
  su2double Density_b, StaticEnergy_b, Enthalpy_b, *Velocity_b, Kappa_b, Chi_b, Energy_b, VelMagnitude2_b, Pressure_b;
  su2double Density_i, *Velocity_i, ProjVelocity_i = 0.0, Energy_i, VelMagnitude2_i;
  su2double **Jacobian_b, **DubDu;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  
  Normal = new su2double[nDim];
  NormalArea = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_i = new su2double[nDim];
  Jacobian_b = new su2double*[nVar];
  DubDu = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_b[iVar] = new su2double[nVar];
    DubDu[iVar] = new su2double[nVar];
  }
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negative for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++) {
        NormalArea[iDim] = -Normal[iDim];
        UnitNormal[iDim] = -Normal[iDim]/Area;
      }
      
      /*--- Compressible solver ---*/
      
      if (compressible) {
        
        /*--- Get the state i ---*/
        
        VelMagnitude2_i = 0.0; ProjVelocity_i = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
          ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
          VelMagnitude2_i += Velocity_i[iDim]*Velocity_i[iDim];
        }
        Density_i = node[iPoint]->GetDensity();
        Energy_i = node[iPoint]->GetEnergy();
        
        /*--- Compute the boundary state b ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity_b[iDim] = Velocity_i[iDim] - ProjVelocity_i * UnitNormal[iDim]; //Force the velocity to be tangential to the surface.
        
        if (grid_movement) {
          GridVel = geometry->node[iPoint]->GetGridVel();
          ProjGridVel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
          for (iDim = 0; iDim < nDim; iDim++) Velocity_b[iDim] += GridVel[iDim] - ProjGridVel * UnitNormal[iDim];
        }
        
        VelMagnitude2_b = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          VelMagnitude2_b += Velocity_b[iDim] * Velocity_b[iDim];
        
        /*--- Compute the residual ---*/
        
        turb_ke = 0.0;
        if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        
        Density_b = Density_i;
        StaticEnergy_b = Energy_i - 0.5 * VelMagnitude2_i - turb_ke;
        Energy_b = StaticEnergy_b + 0.5 * VelMagnitude2_b + turb_ke;
        
        FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
        Kappa_b = FluidModel->GetdPde_rho() / Density_b;
        Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
        Pressure_b = FluidModel->GetPressure();
        Enthalpy_b = Energy_b + Pressure_b/Density_b;
        
        numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, NormalArea, Residual);
        
        /*--- Grid velocity correction to the energy term ---*/
        if (grid_movement) {
          GridVel = geometry->node[iPoint]->GetGridVel();
          ProjGridVel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
          Residual[nVar-1] += Pressure_b*ProjGridVel*Area;
        }
        
        /*--- Add the Reynolds stress tensor contribution ---*/
        
        if (tkeNeeded) {
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[iDim+1] += (2.0/3.0)*Density_b*turb_ke*NormalArea[iDim];
        }
        
      }
      
      /*--- Incompressible solver ---*/
      
      if (incompressible || freesurface) {
        
        /*--- Compute the residual ---*/
        
        Pressure = node[iPoint]->GetPressureInc();
        Density = node[iPoint]->GetPressureInc();
        
        Residual[0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] = Pressure*NormalArea[iDim];
        if (compressible || freesurface) {
          Residual[nVar-1] = 0.0;
        }
        
        /*--- Add the Reynolds stress tensor contribution ---*/
        
        if (tkeNeeded) {
          turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
          for (iDim = 0; iDim < nDim; iDim++)
            Residual[iDim+1] += (2.0/3.0)*Density*turb_ke*NormalArea[iDim];
        }
        
        /*--- Adjustment to energy equation due to grid motion ---*/
        
        if (grid_movement) {
          ProjGridVel = 0.0;
          GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel += GridVel[iDim]*UnitNormal[iDim];
          Residual[nVar-1] = Pressure*ProjGridVel*Area;
        }
        
      }
      
      /*--- Add value to the residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Form Jacobians for implicit computations ---*/
      
      if (implicit) {
        
        /*--- Initialize Jacobian ---*/
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
        
        /*--- Compressible solver ---*/
        
        if (compressible) {
          
          /*--- Compute DubDu ---*/
          
          for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++)
              DubDu[iVar][jVar]= 0.0;
            DubDu[iVar][iVar]= 1.0;
          }
          
          for (iDim = 0; iDim < nDim; iDim++)
            for (jDim = 0; jDim<nDim; jDim++)
              DubDu[iDim+1][jDim+1] -= UnitNormal[iDim]*UnitNormal[jDim];
          DubDu[nVar-1][0] += 0.5*ProjVelocity_i*ProjVelocity_i;
          for (iDim = 0; iDim < nDim; iDim++) {
            DubDu[nVar-1][iDim+1] -= ProjVelocity_i*UnitNormal[iDim];
          }
          
          /*--- Compute flux Jacobian in state b ---*/
          
          numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, NormalArea, 1, Jacobian_b);
          
          // Check for grid movement, should be already considered since Jacobian b is computed from u_b
          // if (grid_movement) {
          // Jacobian_b[nVar-1][0] += 0.5*ProjGridVel*ProjGridVel;
          // for (iDim = 0; iDim < nDim; iDim++)
          // Jacobian_b[nVar-1][iDim+1] -= ProjGridVel * UnitNormal[iDim];
          // }
          
          /*--- Compute numerical flux Jacobian at node i ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              for (kVar = 0; kVar < nVar; kVar++)
                Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
          
          /*--- Add the Jacobian to the sparse matrix ---*/
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
        }
        
        /*--- Incompressible solver ---*/
        
        if (incompressible || freesurface) {
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_i[iDim+1][0] = -Normal[iDim];
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        }
        
      }
    }
  }
  
  delete [] Normal;
  delete [] NormalArea;
  delete [] Velocity_b;
  delete [] Velocity_i;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_b[iVar];
    delete [] DubDu[iVar];
  }
  delete [] Jacobian_b;
  delete [] DubDu;
  
}

void CEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  
  su2double *GridVel;
  su2double Area, UnitNormal[3] = {0.0,0.0,0.0};
  su2double Density, Pressure, Energy,  Velocity[3] = {0.0,0.0,0.0};
  su2double Density_Bound, Pressure_Bound, Vel_Bound[3] = {0.0,0.0,0.0};
  su2double Density_Infty, Pressure_Infty, Vel_Infty[3] = {0.0,0.0,0.0};
  su2double SoundSpeed, Entropy, Velocity2, Vn;
  su2double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound;
  su2double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  su2double RiemannPlus, RiemannMinus;
  su2double *V_infty, *V_domain;
  
  su2double Gas_Constant     = config->GetGas_ConstantND();
  
  bool implicit       = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool grid_movement  = config->GetGrid_Movement();
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool viscous        = config->GetViscous();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS ) ||
                     (config->GetKind_Solver() == DISC_ADJ_RANS))
                    && (config->GetKind_Turb_Model() == SST));
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Allocate the value at the infinity ---*/
    V_infty = GetCharacPrimVar(val_marker, iVertex);
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Construct solution state at infinity (far-field) ---*/
      
      if (compressible) {
        
        /*--- Construct solution state at infinity for compressible flow by
         using Riemann invariants, and then impose a weak boundary condition
         by computing the flux using this new state for U. See CFD texts by
         Hirsch or Blazek for more detail. Adapted from an original
         implementation in the Stanford University multi-block (SUmb) solver
         in the routine bcFarfield.f90 written by Edwin van der Weide,
         last modified 06-12-2005. First, compute the unit normal at the
         boundary nodes. ---*/
        
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);
        
        for (iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = Normal[iDim]/Area;
        
        /*--- Store primitive variables (density, velocities, velocity squared,
         energy, pressure, and sound speed) at the boundary node, and set some
         other quantities for clarity. Project the current flow velocity vector
         at this boundary node into the local normal direction, i.e. compute
         v_bound.n.  ---*/
        
        Density_Bound = V_domain[nDim+2];
        Vel2_Bound = 0.0; Vn_Bound = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Vel_Bound[iDim] = V_domain[iDim+1];
          Vel2_Bound     += Vel_Bound[iDim]*Vel_Bound[iDim];
          Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
        }
        Pressure_Bound   = node[iPoint]->GetPressure();
        SoundSpeed_Bound = sqrt(Gamma*Pressure_Bound/Density_Bound);
        Entropy_Bound    = pow(Density_Bound, Gamma)/Pressure_Bound;
        
        /*--- Store the primitive variable state for the freestream. Project
         the freestream velocity vector into the local normal direction,
         i.e. compute v_infty.n. ---*/
        
        Density_Infty = GetDensity_Inf();
        Vel2_Infty = 0.0; Vn_Infty = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Vel_Infty[iDim] = GetVelocity_Inf(iDim);
          Vel2_Infty     += Vel_Infty[iDim]*Vel_Infty[iDim];
          Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
        }
        Pressure_Infty   = GetPressure_Inf();
        SoundSpeed_Infty = sqrt(Gamma*Pressure_Infty/Density_Infty);
        Entropy_Infty    = pow(Density_Infty, Gamma)/Pressure_Infty;
        
        /*--- Adjust the normal freestream velocity for grid movement ---*/
        
        Qn_Infty = Vn_Infty;
        if (grid_movement) {
          GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            Qn_Infty -= GridVel[iDim]*UnitNormal[iDim];
        }
        
        /*--- Compute acoustic Riemann invariants: R = u.n +/- 2c/(gamma-1).
         These correspond with the eigenvalues (u+c) and (u-c), respectively,
         which represent the acoustic waves. Positive characteristics are
         incoming, and a physical boundary condition is imposed (freestream
         state). This occurs when either (u.n+c) > 0 or (u.n-c) > 0. Negative
         characteristics are leaving the domain, and numerical boundary
         conditions are required by extrapolating from the interior state
         using the Riemann invariants. This occurs when (u.n+c) < 0 or
         (u.n-c) < 0. Note that grid movement is taken into account when
         checking the sign of the eigenvalue. ---*/
        
        /*--- Check whether (u.n+c) is greater or less than zero ---*/
        
        if (Qn_Infty > -SoundSpeed_Infty) {
          /*--- Subsonic inflow or outflow ---*/
          RiemannPlus = Vn_Bound + 2.0*SoundSpeed_Bound/Gamma_Minus_One;
        } else {
          /*--- Supersonic inflow ---*/
          RiemannPlus = Vn_Infty + 2.0*SoundSpeed_Infty/Gamma_Minus_One;
        }
        
        /*--- Check whether (u.n-c) is greater or less than zero ---*/
        
        if (Qn_Infty > SoundSpeed_Infty) {
          /*--- Supersonic outflow ---*/
          RiemannMinus = Vn_Bound - 2.0*SoundSpeed_Bound/Gamma_Minus_One;
        } else {
          /*--- Subsonic outflow ---*/
          RiemannMinus = Vn_Infty - 2.0*SoundSpeed_Infty/Gamma_Minus_One;
        }
        
        /*--- Compute a new value for the local normal velocity and speed of
         sound from the Riemann invariants. ---*/
        
        Vn = 0.5 * (RiemannPlus + RiemannMinus);
        SoundSpeed = 0.25 * (RiemannPlus - RiemannMinus)*Gamma_Minus_One;
        
        /*--- Construct the primitive variable state at the boundary for
         computing the flux for the weak boundary condition. The values
         that we choose to construct the solution (boundary or freestream)
         depend on whether we are at an inflow or outflow. At an outflow, we
         choose boundary information (at most one characteristic is incoming),
         while at an inflow, we choose infinity values (at most one
         characteristic is outgoing). ---*/
        
        if (Qn_Infty > 0.0)   {
          /*--- Outflow conditions ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Bound[iDim] + (Vn-Vn_Bound)*UnitNormal[iDim];
          Entropy = Entropy_Bound;
        } else  {
          /*--- Inflow conditions ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Infty[iDim] + (Vn-Vn_Infty)*UnitNormal[iDim];
          Entropy = Entropy_Infty;
        }
        
        /*--- Recompute the primitive variables. ---*/
        
        Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
        Velocity2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
        Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Store new primitive state for computing the flux. ---*/
        
        V_infty[0] = Pressure/(Gas_Constant*Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_infty[iDim+1] = Velocity[iDim];
        V_infty[nDim+1] = Pressure;
        V_infty[nDim+2] = Density;
        V_infty[nDim+3] = Energy + Pressure/Density;
        
      }
      if (incompressible) {
        
        /*--- All the values computed from the infinity ---*/
        
        V_infty[0] = GetPressure_Inf();
        for (iDim = 0; iDim < nDim; iDim++)
          V_infty[iDim+1] = GetVelocity_Inf(iDim);
        V_infty[nDim+1] = GetDensity_Inf();
        V_infty[nDim+2] = config->GetArtComp_Factor();
        
      }
      if (freesurface) {
        
        /*--- All the values computed from the infinity ---*/
        
        V_infty[0] = GetPressure_Inf();
        for (iDim = 0; iDim < nDim; iDim++)
          V_infty[iDim+1] = GetVelocity_Inf(iDim);
        V_infty[nDim+1] = GetDensity_Inf();
        V_infty[nDim+2] = config->GetArtComp_Factor();
        
      }
      
      /*--- Set various quantities in the numerics class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      if (grid_movement) {
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      }
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Convective Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous residual contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        if (compressible) {
          V_infty[nDim+5] = node[iPoint]->GetLaminarViscosity();
          V_infty[nDim+6] = node[iPoint]->GetEddyViscosity();
        }
        if (incompressible) {
          V_infty[nDim+3] = node[iPoint]->GetLaminarViscosityInc();
          V_infty[nDim+4] = node[iPoint]->GetEddyViscosityInc();
        }
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                          node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                              solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update viscous residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Viscous Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container,
                              CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, kVar;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double P_Total, T_Total, P_static, T_static, Rho_static, *Mach, *Flow_Dir, Area, UnitNormal[3];
  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
  su2double *Velocity_e, Velocity2_e, VelMag_e, Enthalpy_e, Entropy_e, Energy_e = 0.0, StaticEnthalpy_e, StaticEnergy_e, Density_e = 0.0, Pressure_e;
  su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  su2double ProjVelocity_i;
  su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_e, *u_i, *u_b;
  su2double *gridVel;
  su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  
  su2double *Normal, *FlowDirMix, TangVelocity, NormalVelocity;
  Normal = new su2double[nDim];
  su2double ext_flow_angle;
  
  Velocity_i = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_e = new su2double[nDim];
  FlowDirMix = new su2double[nDim];
  Lambda_i = new su2double[nVar];
  u_i = new su2double[nVar];
  u_e = new su2double[nVar];
  u_b = new su2double[nVar];
  dw = new su2double[nVar];
  
  S_boundary = new su2double[8];
  
  P_Tensor = new su2double*[nVar];
  invP_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new su2double[nVar];
    invP_Tensor[iVar] = new su2double[nVar];
  }
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    V_boundary= GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /* --- Compute the internal state u_i --- */
      Velocity2_i = 0;
      for (iDim=0; iDim < nDim; iDim++)
      {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      
      
      Density_i = node[iPoint]->GetDensity();
      
      Energy_i = node[iPoint]->GetEnergy();
      StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
      
      FluidModel->SetTDState_rhoe(Density_i, StaticEnergy_i);
      
      Pressure_i = FluidModel->GetPressure();
      Enthalpy_i = Energy_i + Pressure_i/Density_i;
      
      SoundSpeed_i = FluidModel->GetSoundSpeed();
      
      Kappa_i = FluidModel->GetdPde_rho() / Density_i;
      Chi_i = FluidModel->GetdPdrho_e() - Kappa_i * StaticEnergy_i;
      
      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      
      /*--- Build the external state u_e from boundary data and internal node ---*/
      
      switch(config->GetKind_Data_Riemann(Marker_Tag))
      {
          //TODO(turbo), generilize for 3D case
          //TODO(turbo), generilize for Inlet and Outlet in for backflow treatment
          //TODO(turbo), implement not uniform inlet and radial equilibrium for the outlet
        case TOTAL_CONDITIONS_PT:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          if (gravity) P_Total = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
          else P_Total  = config->GetRiemann_Var1(Marker_Tag);
          T_Total  = config->GetRiemann_Var2(Marker_Tag);
          Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();
          
          /* --- Computes the total state --- */
          FluidModel->SetTDState_PT(P_Total, T_Total);
          Enthalpy_e = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
          Entropy_e = FluidModel->GetEntropy();
          
          /* --- Compute the boundary state u_e --- */
          Velocity2_e = Velocity2_i;
          if (nDim == 2){
            NormalVelocity= -sqrt(Velocity2_e)*Flow_Dir[0];
            TangVelocity= -sqrt(Velocity2_e)*Flow_Dir[1];
            Velocity_e[0]= UnitNormal[0]*NormalVelocity - UnitNormal[1]*TangVelocity;
            Velocity_e[1]= UnitNormal[1]*NormalVelocity + UnitNormal[0]*TangVelocity;
          }else{
            for (iDim = 0; iDim < nDim; iDim++)
            	Velocity_e[iDim] = sqrt(Velocity2_e)*Flow_Dir[iDim];

//						NormalVelocity= -sqrt(Velocity2_e)*Flow_Dir[0];
//						TangVelocity= -sqrt(Velocity2_e)*Flow_Dir[1];
//						Velocity_e[0]= UnitNormal[0]*NormalVelocity - UnitNormal[1]*TangVelocity;
//						Velocity_e[1]= UnitNormal[1]*NormalVelocity + UnitNormal[0]*TangVelocity;
//						Velocity_e[2] = sqrt(Velocity2_e)*Flow_Dir[2];
          }
          StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
          FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case STATIC_SUPERSONIC_INFLOW_PT:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
          else P_static  = config->GetRiemann_Var1(Marker_Tag);
          T_static  = config->GetRiemann_Var2(Marker_Tag);
          Mach = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          P_static /= config->GetPressure_Ref();
          T_static /= config->GetTemperature_Ref();
          
          /* --- Computes the total state --- */
          FluidModel->SetTDState_PT(P_static, T_static);
          
          /* --- Compute the boundary state u_e --- */
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Mach[iDim]*FluidModel->GetSoundSpeed();
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case STATIC_SUPERSONIC_INFLOW_PD:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          
          if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;/// check in which case is true (only freesurface?)
          else P_static  = config->GetRiemann_Var1(Marker_Tag);
          Rho_static  = config->GetRiemann_Var2(Marker_Tag);
          Mach = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          P_static /= config->GetPressure_Ref();
          Rho_static /= config->GetDensity_Ref();
          
          /* --- Computes the total state --- */
          FluidModel->SetTDState_Prho(P_static, Rho_static);
          
          /* --- Compute the boundary state u_e --- */
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Mach[iDim]*FluidModel->GetSoundSpeed();
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case MIXING_IN:
          
          /*--- Retrieve the specified total conditions for this boundary. ---*/
          P_Total = ExtAveragedTotPressure[val_marker];
          T_Total = ExtAveragedTotTemperature[val_marker];
          ext_flow_angle = atan(ExtAveragedTangVelocity[val_marker]/ExtAveragedNormalVelocity[val_marker]);
          FlowDirMix[0] = cos(ext_flow_angle);
          FlowDirMix[1] = sin(ext_flow_angle);
          
          /* --- Computes the total state --- */
          FluidModel->SetTDState_PT(P_Total, T_Total);
          Enthalpy_e = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
          Entropy_e = FluidModel->GetEntropy();
          
          /* --- Compute the boundary state u_e --- */
          Velocity2_e = Velocity2_i;
          if (nDim == 2){
            NormalVelocity= -sqrt(Velocity2_e)*FlowDirMix[0];
            TangVelocity= -sqrt(Velocity2_e)*FlowDirMix[1];
            Velocity_e[0]= UnitNormal[0]*NormalVelocity - UnitNormal[1]*TangVelocity;
            Velocity_e[1]= UnitNormal[1]*NormalVelocity + UnitNormal[0]*TangVelocity;
          }else{
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity_e[iDim] = sqrt(Velocity2_e)*FlowDirMix[iDim];
          }
          StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
          FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
          Density_e = FluidModel->GetDensity();
          StaticEnergy_e = FluidModel->GetStaticEnergy();
          Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
          if (tkeNeeded) Energy_e += GetTke_Inf();
          break;
          
        case DENSITY_VELOCITY:
          
          /*--- Retrieve the specified density and velocity magnitude ---*/
          Density_e  = config->GetRiemann_Var1(Marker_Tag);
          VelMag_e   = config->GetRiemann_Var2(Marker_Tag);
          Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          Density_e /= config->GetDensity_Ref();
          VelMag_e /= config->GetVelocity_Ref();
          
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity_e[iDim] = VelMag_e*Flow_Dir[iDim];
          Energy_e = Energy_i;
          break;
          
        case MIXING_OUT:
          
          /*--- Retrieve the staic pressure for this boundary. ---*/
          Pressure_e = ExtAveragedPressure[val_marker];
          Density_e = Density_i;
          
          /* --- Compute the boundary state u_e --- */
          FluidModel->SetTDState_Prho(Pressure_e, Density_e);
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Velocity_i[iDim];
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
          break;
          
        case STATIC_PRESSURE:
          
          /*--- Retrieve the staic pressure for this boundary. ---*/
          Pressure_e = config->GetRiemann_Var1(Marker_Tag);
          Pressure_e /= config->GetPressure_Ref();
          Density_e = Density_i;
          
          /* --- Compute the boundary state u_e --- */
          FluidModel->SetTDState_Prho(Pressure_e, Density_e);
          Velocity2_e = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity_e[iDim] = Velocity_i[iDim];
            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
          }
          Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
          break;
          
        default:
          cout << "Warning! Invalid Riemann input!" << endl;
          exit(EXIT_FAILURE);
          break;
          
      }
      
      /*--- Compute P (matrix of right eigenvectors) ---*/
      conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);
      
      /*--- Compute inverse P (matrix of left eigenvectors)---*/
      conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);
      
      /*--- eigenvalues contribution due to grid motion ---*/
      if (grid_movement){
        gridVel = geometry->node[iPoint]->GetGridVel();
        
        su2double ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel   += gridVel[iDim]*UnitNormal[iDim];
        ProjVelocity_i -= ProjGridVel;
      }
      
      /*--- Flow eigenvalues ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Lambda_i[iDim] = ProjVelocity_i;
      Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
      Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;
      
      /* --- Compute the boundary state u_e --- */
      u_e[0] = Density_e;
      for (iDim = 0; iDim < nDim; iDim++)
        u_e[iDim+1] = Velocity_e[iDim]*Density_e;
      u_e[nVar-1] = Energy_e*Density_e;
      
      /* --- Compute the boundary state u_i --- */
      u_i[0] = Density_i;
      for (iDim = 0; iDim < nDim; iDim++)
        u_i[iDim+1] = Velocity_i[iDim]*Density_i;
      u_i[nVar-1] = Energy_i*Density_i;
      
      /*--- Compute the characteristic jumps ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dw[iVar] = 0;
        for (jVar = 0; jVar < nVar; jVar++)
          dw[iVar] += invP_Tensor[iVar][jVar] * (u_e[jVar] - u_i[jVar]);
        
      }
      
      /*--- Compute the boundary state u_b using characteristics ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        u_b[iVar] = u_i[iVar];
        
        for (jVar = 0; jVar < nVar; jVar++)
        {
          if (Lambda_i[jVar] < 0)
          {
            u_b[iVar] += P_Tensor[iVar][jVar]*dw[jVar];
            
          }
        }
      }
      
      
      /*--- Compute the thermodynamic state in u_b ---*/
      Density_b = u_b[0];
      Velocity2_b = 0;
      for (iDim = 0; iDim < nDim; iDim++)
      {
        Velocity_b[iDim] = u_b[iDim+1]/Density_b;
        Velocity2_b += Velocity_b[iDim]*Velocity_b[iDim];
      }
      Energy_b = u_b[nVar-1]/Density_b;
      StaticEnergy_b = Energy_b - 0.5*Velocity2_b;
      FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);
      Pressure_b = FluidModel->GetPressure();
      Temperature_b = FluidModel->GetTemperature();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      
      /*--- Compute the residuals ---*/
      conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);
      
      /*--- Residual contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        su2double projVelocity = 0.0;
        
        for (iDim = 0; iDim < nDim; iDim++)
          projVelocity +=  gridVel[iDim]*Normal[iDim];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] -= projVelocity *(u_b[iVar]);
      }
      
      if (implicit) {
        
        Jacobian_b = new su2double*[nVar];
        DubDu = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++)
        {
          Jacobian_b[iVar] = new su2double[nVar];
          DubDu[iVar] = new su2double[nVar];
        }
        
        /*--- Initialize DubDu to unit matrix---*/
        
        for (iVar = 0; iVar < nVar; iVar++)
        {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0;
          
          DubDu[iVar][iVar]= 1;
        }
        
        /*--- Compute DubDu -= RNL---*/
        for (iVar=0; iVar<nVar; iVar++)
        {
          for (jVar=0; jVar<nVar; jVar++)
          {
            for (kVar=0; kVar<nVar; kVar++)
            {
              if (Lambda_i[kVar]<0)
                DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
            }
          }
        }
        
        /*--- Compute flux Jacobian in state b ---*/
        conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);
        
        /*--- Jacobian contribution due to grid motion ---*/
        if (grid_movement)
        {
          gridVel = geometry->node[iPoint]->GetGridVel();
          su2double projVelocity = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            projVelocity +=  gridVel[iDim]*Normal[iDim];
          for (iVar = 0; iVar < nVar; iVar++){
            Residual[iVar] -= projVelocity *(u_b[iVar]);
            Jacobian_b[iVar][iVar] -= projVelocity;
          }
          
        }
        
        /*--- initiate Jacobian_i to zero matrix ---*/
        for (iVar=0; iVar<nVar; iVar++)
          for (jVar=0; jVar<nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        
        /*--- Compute numerical flux Jacobian at node i ---*/
        for (iVar=0; iVar<nVar; iVar++) {
          for (jVar=0; jVar<nVar; jVar++) {
            for (kVar=0; kVar<nVar; kVar++) {
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
            }
          }
        }
        
        for (iVar = 0; iVar < nVar; iVar++) {
          delete [] Jacobian_b[iVar];
          delete [] DubDu[iVar];
        }
        delete [] Jacobian_b;
        delete [] DubDu;
      }
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Primitive variables, using the derived quantities ---*/
        V_boundary[0] = Temperature_b;
        for (iDim = 0; iDim < nDim; iDim++)
          V_boundary[iDim+1] = Velocity_b[iDim];
        V_boundary[nDim+1] = Pressure_b;
        V_boundary[nDim+2] = Density_b;
        V_boundary[nDim+3] = Enthalpy_b;
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_boundary[nDim+5] = FluidModel->GetLaminarViscosity();
        V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();
        V_boundary[nDim+7] = FluidModel->GetThermalConductivity();
        V_boundary[nDim+8] = FluidModel->GetCp();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_boundary);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Secondary variables ---*/
        S_domain = node[iPoint]->GetSecondary();
        
        /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
        
        S_boundary[0]= FluidModel->GetdPdrho_e();
        S_boundary[1]= FluidModel->GetdPde_rho();
        
        S_boundary[2]= FluidModel->GetdTdrho_e();
        S_boundary[3]= FluidModel->GetdTde_rho();
        
        /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
        
        S_boundary[4]= FluidModel->Getdmudrho_T();
        S_boundary[5]= FluidModel->GetdmudT_rho();
        
        S_boundary[6]= FluidModel->Getdktdrho_T();
        S_boundary[7]= FluidModel->GetdktdT_rho();
        
        visc_numerics->SetSecondary(S_domain, S_boundary);
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] Velocity_e;
  delete [] Velocity_b;
  delete [] Velocity_i;
  delete [] FlowDirMix;
  
  delete [] S_boundary;
  delete [] Lambda_i;
  delete [] u_i;
  delete [] u_e;
  delete [] u_b;
  delete [] dw;
  
  
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}


void CEulerSolver::Mixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker) {
  
  unsigned long iVertex, iPoint, nVert;
  unsigned short iDim, iVar;
  unsigned short mixing_process = config->GetKind_MixingProcess();
  su2double Pressure = 0.0, Density = 0.0, Enthalpy = 0.0,  *Velocity = NULL, *Normal, *gridVel,
  Area, TotalArea, TotalAreaPressure, TotalAreaDensity, *TotalAreaVelocity, UnitNormal[3];
  string Marker_Tag, Monitoring_Tag;
  su2double val_init_pressure;
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool grid_movement        = config->GetGrid_Movement();
  su2double TotalDensity, TotalPressure, *TotalVelocity, TotalNormal, avgVel2, avgTotalEnthaply;
  int rank = MASTER_NODE;
  /*-- Variables declaration and allocation ---*/
  Velocity = new su2double[nDim];
  Normal = new su2double[nDim];
  TotalVelocity = new su2double[nDim];
  TotalAreaVelocity = new su2double[nDim];
  
  for (iDim=0; iDim<nDim; iDim++) {
    TotalVelocity[iDim]=0;
    TotalAreaVelocity[iDim]=0;
  }
  
  TotalDensity = 0.0;
  TotalPressure = 0.0;
  TotalAreaPressure=0.0;
  TotalAreaDensity=0.0;
  TotalArea = 0.0;
  TotalNormal=0.0;
  
  /*--- Forces initialization for Marker vector ---*/
  AveragedPressure[val_Marker] = 0.0;
  AveragedEnthalpy[val_Marker] = 0.0;
  AveragedDensity[val_Marker] = 0.0;
  AveragedSoundSpeed[val_Marker] = 0.0;
  
  for (iDim=0;iDim < nDim;iDim++){
    AveragedVelocity[val_Marker][iDim] = 0.0;
    AveragedNormal[val_Marker][iDim] = 0.0;
    AveragedGridVel[val_Marker][iDim] = 0.0;
  }
  
  for (iVar=0;iVar<nVar;iVar++)
    TotalFlux[val_Marker][iVar]= 0.0;
  

  /*--- Loop over the vertices to compute the averaged quantities ---*/
  nVert = 0;
  for (iVertex = 0; iVertex < geometry->GetnVertex(val_Marker); iVertex++) {
    
    iPoint = geometry->vertex[val_Marker][iVertex]->GetNode();
    
    /*--- Compute the integral fluxes for the boundaries ---*/
    if (compressible) {
      Pressure = node[iPoint]->GetPressure();
      Density = node[iPoint]->GetDensity();
      Enthalpy = node[iPoint]->GetEnthalpy();
    }
    else {
      cout << "!!! Mixing process for incompressible and freesurface does not available yet !!! " << endl;
      cout << "Press any key to exit..." << endl;
      cin.get();
      exit(1);
    }
    
    /*--- Note that the fluxes from halo cells are discarded ---*/
    if ( (geometry->node[iPoint]->GetDomain())  ) {
      nVert++;
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_Marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      su2double VelNormal = 0.0, VelSq = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        UnitNormal[iDim] = Normal[iDim]/Area;
        Velocity[iDim] = node[iPoint]->GetPrimitive(iDim+1);
        VelNormal += UnitNormal[iDim]*Velocity[iDim];
        VelSq += Velocity[iDim]*Velocity[iDim];
      }
      
      
      /*--- Compute the integral fluxes for the boundary of interest ---*/
      
      if ((mixing_process == AREA_AVERAGE) || (mixing_process == MIXEDOUT_AVERAGE)){
        
        TotalFlux[val_Marker][0] += Area*(Density*VelNormal );
        for (iDim = 1; iDim < nDim+1; iDim++)
          TotalFlux[val_Marker][iDim] += Area*(Density*VelNormal*Velocity[iDim -1] + Pressure*UnitNormal[iDim -1] );
        TotalFlux[val_Marker][nDim+1] += Area*(Density*VelNormal*Enthalpy );
        
        TotalArea += Area;
        TotalAreaPressure += Area*Pressure;
        TotalAreaDensity  += Area*Density;
        for (iDim = 0; iDim < nDim; iDim++)
          TotalAreaVelocity[iDim] += Area*Velocity[iDim];
        
      }else{
        
        TotalDensity += Density;
        TotalPressure += Pressure;
        for (iDim = 0; iDim < nDim; iDim++)
          TotalVelocity[iDim] += Velocity[iDim];
        
        
      }
      for (iDim = 0; iDim < nDim; iDim++) AveragedNormal[val_Marker][iDim] +=Normal[iDim];
      if (grid_movement){
        gridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          AveragedGridVel[val_Marker][iDim] +=gridVel[iDim];
      }
    }
  }
  

  /*--- Compute the averaged value for the boundary of interest ---*/
  for (iDim = 0; iDim < nDim; iDim++){
    AveragedNormal[val_Marker][iDim] /=nVert;
    TotalNormal+= AveragedNormal[val_Marker][iDim]*AveragedNormal[val_Marker][iDim];
  }
  for (iDim = 0; iDim < nDim; iDim++){ AveragedNormal[val_Marker][iDim] /=sqrt(TotalNormal);
//  cout <<" normal vector "<< AveragedNormal[val_Marker][iDim] << " comp "<< iDim <<endl;
  }
  if (grid_movement){
    for (iDim = 0; iDim < nDim; iDim++)
      AveragedGridVel[val_Marker][iDim] /=nVert;
  }

  switch(mixing_process){
      
      
    case ALGEBRAIC_AVERAGE:
      AveragedDensity[val_Marker] = TotalDensity / nVert;
      AveragedPressure[val_Marker] = TotalPressure / nVert;
      for (iDim = 0; iDim < nDim; iDim++)
        AveragedVelocity[val_Marker][iDim] = TotalVelocity[iDim] / nVert;
      break;
      
    case AREA_AVERAGE:
      AveragedDensity[val_Marker] = TotalAreaDensity / TotalArea;
      AveragedPressure[val_Marker] = TotalAreaPressure / TotalArea;
      for (iDim = 0; iDim < nDim; iDim++)
        AveragedVelocity[val_Marker][iDim] = TotalAreaVelocity[iDim] / TotalArea;
      break;
      
    case MIXEDOUT_AVERAGE:
      for (iVar = 0; iVar<nVar; iVar++){
        AveragedFlux[val_Marker][iVar] = TotalFlux[val_Marker][iVar]/TotalArea;
      }
      val_init_pressure = TotalAreaPressure/TotalArea;
      
      if (abs(AveragedFlux[val_Marker][0])<(10.0e-9)*TotalAreaDensity) {
        cout << "Mass flux is 0.0 so a Area Averaged algorithm is used for the Mixing Procees" << endl;
        AveragedDensity[val_Marker] = TotalAreaDensity / TotalArea;
        AveragedPressure[val_Marker] = TotalAreaPressure / TotalArea;
        for (iDim = 0; iDim < nDim; iDim++)
          AveragedVelocity[val_Marker][iDim] = TotalAreaVelocity[iDim] / TotalArea;
        
      }else {
        MixedOut_Average (val_init_pressure, AveragedFlux[val_Marker], AveragedNormal[val_Marker], &AveragedPressure[val_Marker], &AveragedDensity[val_Marker]);
        for (iDim = 1; iDim < nDim +1;iDim++)
          AveragedVelocity[val_Marker][iDim-1]= ( AveragedFlux[val_Marker][iDim] - AveragedPressure[val_Marker]*AveragedNormal[val_Marker][iDim-1] ) / AveragedFlux[val_Marker][0];
      }
      break;
      
      
    default:
      cout << "Warning! Invalid MIXING_PROCESS input!" << endl;
      exit(EXIT_FAILURE);
      break;
  }
  
  /* --- compute static averaged quantities ---*/
  FluidModel->SetTDState_Prho(AveragedPressure[val_Marker], AveragedDensity[val_Marker]);
  AveragedEnthalpy[val_Marker] = FluidModel->GetStaticEnergy() + AveragedPressure[val_Marker]/AveragedDensity[val_Marker];
  AveragedSoundSpeed[val_Marker] = FluidModel->GetSoundSpeed();
  AveragedEntropy[val_Marker] = FluidModel->GetEntropy();
  AveragedNormalVelocity[val_Marker]= AveragedNormal[val_Marker][0]*AveragedVelocity[val_Marker][0] + AveragedNormal[val_Marker][1]*AveragedVelocity[val_Marker][1];
  AveragedTangVelocity[val_Marker]= AveragedNormal[val_Marker][0]*AveragedVelocity[val_Marker][1] - AveragedNormal[val_Marker][1]*AveragedVelocity[val_Marker][0];
  MassFlow[val_Marker]= AveragedDensity[val_Marker]*AveragedNormalVelocity[val_Marker]*TotalArea;
  FlowAngle[val_Marker]= atan(AveragedTangVelocity[val_Marker]/AveragedNormalVelocity[val_Marker]);
  
  /* --- compute total averaged quantities ---*/
  avgVel2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) avgVel2 += AveragedVelocity[val_Marker][iDim]*AveragedVelocity[val_Marker][iDim];
  
  avgTotalEnthaply = AveragedEnthalpy[val_Marker] + 0.5*avgVel2;
  FluidModel->SetTDState_hs(avgTotalEnthaply,AveragedEntropy[val_Marker]);
  AveragedTotTemperature[val_Marker] = FluidModel->GetTemperature();
  AveragedTotPressure[val_Marker] = FluidModel->GetPressure();
  
  if(grid_movement){
    AveragedTangGridVelocity[val_Marker] = AveragedNormal[val_Marker][0]*AveragedGridVel[val_Marker][1]-AveragedNormal[val_Marker][1]*AveragedGridVel[val_Marker][0];
    AveragedMach[val_Marker] = sqrt(AveragedNormalVelocity[val_Marker]*AveragedNormalVelocity[val_Marker] + (AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker])*(AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker]));
    AveragedMach[val_Marker] /= AveragedSoundSpeed[val_Marker];
    AveragedTangMach[val_Marker] = (AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker])/AveragedSoundSpeed[val_Marker];
    FlowAngle[val_Marker]= atan((AveragedTangVelocity[val_Marker] - AveragedTangGridVelocity[val_Marker])/AveragedNormalVelocity[val_Marker]);
    
  }else{
    AveragedMach[val_Marker] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      AveragedMach[val_Marker] += AveragedVelocity[val_Marker][iDim]*AveragedVelocity[val_Marker][iDim];
    }
    AveragedMach[val_Marker] = sqrt(AveragedMach[val_Marker])/AveragedSoundSpeed[val_Marker];
    AveragedTangMach[val_Marker] = AveragedTangVelocity[val_Marker]/AveragedSoundSpeed[val_Marker];
    
  }
  
  AveragedNormalMach[val_Marker] = AveragedNormalVelocity[val_Marker]/AveragedSoundSpeed[val_Marker];
  
  
  if ((AveragedDensity[val_Marker]!= AveragedDensity[val_Marker]) or (AveragedEnthalpy[val_Marker]!=AveragedEnthalpy[val_Marker]))
    cout<<"nan in mixing process in boundary "<<config->GetMarker_All_TagBound(val_Marker)<< endl;
  
  /*--- Free locally allocated memory ---*/
  delete [] Velocity;
  delete [] Normal;
  delete [] TotalVelocity;
  delete [] TotalAreaVelocity;
}


void CEulerSolver::MPIMixing_Process(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short marker_flag) {

  unsigned long iVertex, iPoint, nVert;
  unsigned short iDim, iVar, iMarker, iMarkerTP;
  unsigned short mixing_process = config->GetKind_MixingProcess();
  su2double Pressure = 0.0, Density = 0.0, Enthalpy = 0.0,  *Velocity = NULL, *Normal, *gridVel,
  Area, TotalArea, TotalAreaPressure, TotalAreaDensity, *TotalAreaVelocity, UnitNormal[3];
  string Marker_Tag, Monitoring_Tag;
  su2double val_init_pressure;
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool grid_movement        = config->GetGrid_Movement();
  su2double TotalDensity, TotalPressure, *TotalVelocity, *TotalNormal, avgVel2, avgTotalEnthaply, *TotalFluxes, *TotalGridVel;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  /*-- Variables declaration and allocation ---*/
  Velocity = new su2double[nDim];
  Normal = new su2double[nDim];
  TotalVelocity = new su2double[nDim];
  TotalAreaVelocity = new su2double[nDim];
  TotalNormal = new su2double[nDim];
  TotalGridVel = new su2double[nDim];
  TotalFluxes = new su2double[nVar];





#ifdef HAVE_MPI
  su2double MyTotalDensity, MyTotalPressure, MyTotalAreaDensity, MyTotalAreaPressure, *MyTotalFluxes = NULL;
  su2double MyTotalArea, *MyTotalNormal= NULL,*MyTotalGridVel= NULL, *MyTotalVelocity = NULL, *MyTotalAreaVelocity = NULL;
  unsigned long My_nVert;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  /*--- Forces initialization for contenitors ---*/
  for (iVar=0;iVar<nVar;iVar++)
    TotalFluxes[iVar]= 0.0;
  for (iDim=0; iDim<nDim; iDim++) {
      TotalVelocity[iDim]=0.0;
      TotalAreaVelocity[iDim]=0.0;
      TotalNormal[iDim]=0.0;
      TotalGridVel[iDim]=0.0;
  }

  TotalDensity = 0.0;
  TotalPressure = 0.0;
  TotalAreaPressure=0.0;
  TotalAreaDensity=0.0;
  TotalArea = 0.0;
  nVert = 0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    /*--- Loop over the vertices to compute the averaged quantities ---*/
    for (iMarkerTP=1; iMarkerTP < config->Get_nMarkerTurboPerf()+1; iMarkerTP++)
      if (config->GetMarker_All_TurboPerformance(iMarker) == iMarkerTP)
        if (config->GetMarker_All_TurboPerformanceFlag(iMarker) == marker_flag){

      	for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      		iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					/*--- Compute the integral fluxes for the boundaries ---*/
					if (compressible) {
						Pressure = node[iPoint]->GetPressure();
						Density = node[iPoint]->GetDensity();
						Enthalpy = node[iPoint]->GetEnthalpy();
					}
					else {
						cout << "!!! Mixing process for incompressible and freesurface does not available yet !!! " << endl;
						cout << "Press any key to exit..." << endl;
						cin.get();
						exit(1);
					}

					/*--- Note that the fluxes from halo cells are discarded ---*/
					if ( (geometry->node[iPoint]->GetDomain())  ) {
						nVert++;

						/*--- Normal vector for this vertex (negate for outward convention) ---*/
						geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
						for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
						Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
						su2double VelNormal = 0.0, VelSq = 0.0;
						for (iDim = 0; iDim < nDim; iDim++) {
							UnitNormal[iDim] = Normal[iDim]/Area;
							Velocity[iDim] = node[iPoint]->GetPrimitive(iDim+1);
							VelNormal += UnitNormal[iDim]*Velocity[iDim];
							VelSq += Velocity[iDim]*Velocity[iDim];
						}


						/*--- Compute the integral fluxes for the boundary of interest ---*/

						if ((mixing_process == AREA_AVERAGE) || (mixing_process == MIXEDOUT_AVERAGE)){

							TotalFluxes[0] += Area*(Density*VelNormal );
							for (iDim = 1; iDim < nDim+1; iDim++)
								TotalFluxes[iDim] += Area*(Density*VelNormal*Velocity[iDim -1] + Pressure*UnitNormal[iDim -1] );
							TotalFluxes[nDim+1] += Area*(Density*VelNormal*Enthalpy );

							TotalArea += Area;
							TotalAreaPressure += Area*Pressure;
							TotalAreaDensity  += Area*Density;
							for (iDim = 0; iDim < nDim; iDim++)
								TotalAreaVelocity[iDim] += Area*Velocity[iDim];

						}else{

							TotalDensity += Density;
							TotalPressure += Pressure;
							for (iDim = 0; iDim < nDim; iDim++)
								TotalVelocity[iDim] += Velocity[iDim];


						}
						for (iDim = 0; iDim < nDim; iDim++) TotalNormal[iDim] +=Normal[iDim];
						if (grid_movement){
							gridVel = geometry->node[iPoint]->GetGridVel();
							for (iDim = 0; iDim < nDim; iDim++)
								TotalGridVel[iDim] +=gridVel[iDim];
						}
					}
				}

			}

  }

#ifdef HAVE_MPI

	/*--- Add information using all the nodes ---*/

	MyTotalDensity       = TotalDensity; 							TotalDensity         = 0;
	MyTotalPressure      = TotalPressure;  					  TotalPressure        = 0;
	MyTotalAreaDensity   = TotalAreaDensity; 					TotalAreaDensity     = 0;
	MyTotalAreaPressure  = TotalAreaPressure;         TotalAreaPressure    = 0;
	MyTotalArea          = TotalArea;                 TotalArea            = 0;
	My_nVert						 = nVert;											nVert								 = 0;

	SU2_MPI::Allreduce(&MyTotalDensity, &TotalDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(&MyTotalPressure, &TotalPressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(&MyTotalAreaDensity, &TotalAreaDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(&MyTotalAreaPressure, &TotalAreaPressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(&MyTotalArea, &TotalArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(&My_nVert, &nVert, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);



	MyTotalFluxes					 = new su2double[nVar];
	MyTotalVelocity        = new su2double[nDim];
	MyTotalAreaVelocity    = new su2double[nDim];
	MyTotalNormal          = new su2double[nDim];
	MyTotalGridVel				 = new su2double[nDim];

	for (iVar = 0; iVar < nVar; iVar++) {
		MyTotalFluxes[iVar]  = TotalFluxes[iVar];
		TotalFluxes[iVar]    = 0.0;
	}

	for (iDim = 0; iDim < nDim; iDim++) {
		MyTotalVelocity[iDim]      			  = TotalVelocity[iDim];
		TotalVelocity[iDim]        			  = 0.0;
		MyTotalAreaVelocity[iDim]  			  = TotalAreaVelocity[iDim];
		TotalAreaVelocity[iDim]    				= 0.0;
		MyTotalNormal[iDim]				 				= TotalNormal[iDim];
		TotalNormal[iDim] 								= 0.0;
		MyTotalGridVel[iDim] 							= TotalGridVel[iDim];
		TotalGridVel[iDim]								= 0.0;
		}

	SU2_MPI::Allreduce(MyTotalFluxes, TotalFluxes, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(MyTotalVelocity, TotalVelocity, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(MyTotalAreaVelocity, TotalAreaVelocity, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(MyTotalNormal, TotalNormal, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	SU2_MPI::Allreduce(MyTotalGridVel, TotalGridVel, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	delete [] MyTotalFluxes; delete [] MyTotalVelocity; delete [] MyTotalAreaVelocity;
	delete [] MyTotalNormal; delete [] MyTotalGridVel;

#endif


	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iMarkerTP=1; iMarkerTP < config->Get_nMarkerTurboPerf()+1; iMarkerTP++)
			if (config->GetMarker_All_TurboPerformance(iMarker) == iMarkerTP)
				if (config->GetMarker_All_TurboPerformanceFlag(iMarker) == marker_flag){


					/*--- Compute the averaged value for the boundary of interest ---*/
					su2double Normal2 = 0.0;
					for (iDim = 0; iDim < nDim; iDim++){
						TotalNormal[iDim] /=nVert;
						Normal2 += TotalNormal[iDim]*TotalNormal[iDim];
					}
					for (iDim = 0; iDim < nDim; iDim++)AveragedNormal[iMarker][iDim] = TotalNormal[iDim]/sqrt(Normal2);

					if (grid_movement){
						for (iDim = 0; iDim < nDim; iDim++)AveragedGridVel[iMarker][iDim] =TotalGridVel[iDim]/nVert;
					}

					switch(mixing_process){
						case ALGEBRAIC_AVERAGE:
							AveragedDensity[iMarker] = TotalDensity / nVert;
							AveragedPressure[iMarker] = TotalPressure / nVert;
							for (iDim = 0; iDim < nDim; iDim++)
								AveragedVelocity[iMarker][iDim] = TotalVelocity[iDim] / nVert;
							break;

						case AREA_AVERAGE:
							AveragedDensity[iMarker] = TotalAreaDensity / TotalArea;
							AveragedPressure[iMarker] = TotalAreaPressure / TotalArea;
							for (iDim = 0; iDim < nDim; iDim++)
								AveragedVelocity[iMarker][iDim] = TotalAreaVelocity[iDim] / TotalArea;
							break;

						case MIXEDOUT_AVERAGE:
							for (iVar = 0; iVar<nVar; iVar++){
								AveragedFlux[iMarker][iVar] = TotalFluxes[iVar]/TotalArea;
							}
							val_init_pressure = TotalAreaPressure/TotalArea;

							if (abs(AveragedFlux[iMarker][0])<(10.0e-9)*TotalAreaDensity) {
								cout << "Mass flux is 0.0 so a Area Averaged algorithm is used for the Mixing Procees" << endl;
								AveragedDensity[iMarker] = TotalAreaDensity / TotalArea;
								AveragedPressure[iMarker] = TotalAreaPressure / TotalArea;
								for (iDim = 0; iDim < nDim; iDim++)
									AveragedVelocity[iMarker][iDim] = TotalAreaVelocity[iDim] / TotalArea;

							}else {
								MixedOut_Average (val_init_pressure, AveragedFlux[iMarker], AveragedNormal[iMarker], &AveragedPressure[iMarker], &AveragedDensity[iMarker]);
								for (iDim = 1; iDim < nDim +1;iDim++)
									AveragedVelocity[iMarker][iDim-1]= ( AveragedFlux[iMarker][iDim] - AveragedPressure[iMarker]*AveragedNormal[iMarker][iDim-1] ) / AveragedFlux[iMarker][0];
							}
							break;


						default:
							cout << "Warning! Invalid MIXING_PROCESS input!" << endl;
							exit(EXIT_FAILURE);
							break;
					}

					/* --- compute static averaged quantities ---*/
					FluidModel->SetTDState_Prho(AveragedPressure[iMarker], AveragedDensity[iMarker]);
					AveragedEnthalpy[iMarker] = FluidModel->GetStaticEnergy() + AveragedPressure[iMarker]/AveragedDensity[iMarker];
					AveragedSoundSpeed[iMarker] = FluidModel->GetSoundSpeed();
					AveragedEntropy[iMarker] = FluidModel->GetEntropy();
					AveragedNormalVelocity[iMarker]= AveragedNormal[iMarker][0]*AveragedVelocity[iMarker][0] + AveragedNormal[iMarker][1]*AveragedVelocity[iMarker][1];
					AveragedTangVelocity[iMarker]= AveragedNormal[iMarker][0]*AveragedVelocity[iMarker][1] - AveragedNormal[iMarker][1]*AveragedVelocity[iMarker][0];
					MassFlow[iMarker]= AveragedDensity[iMarker]*AveragedNormalVelocity[iMarker]*TotalArea;
					FlowAngle[iMarker]= atan(AveragedTangVelocity[iMarker]/AveragedNormalVelocity[iMarker]);

					/* --- compute total averaged quantities ---*/
					avgVel2 = 0.0;
					for (iDim = 0; iDim < nDim; iDim++) avgVel2 += AveragedVelocity[iMarker][iDim]*AveragedVelocity[iMarker][iDim];

					avgTotalEnthaply = AveragedEnthalpy[iMarker] + 0.5*avgVel2;
					FluidModel->SetTDState_hs(avgTotalEnthaply,AveragedEntropy[iMarker]);
					AveragedTotTemperature[iMarker] = FluidModel->GetTemperature();
					AveragedTotPressure[iMarker] = FluidModel->GetPressure();

					if(grid_movement){
						AveragedTangGridVelocity[iMarker] = AveragedNormal[iMarker][0]*AveragedGridVel[iMarker][1]-AveragedNormal[iMarker][1]*AveragedGridVel[iMarker][0];
						AveragedMach[iMarker] = sqrt(AveragedNormalVelocity[iMarker]*AveragedNormalVelocity[iMarker] + (AveragedTangVelocity[iMarker] - AveragedTangGridVelocity[iMarker])*(AveragedTangVelocity[iMarker] - AveragedTangGridVelocity[iMarker]));
						AveragedMach[iMarker] /= AveragedSoundSpeed[iMarker];
						AveragedTangMach[iMarker] = (AveragedTangVelocity[iMarker] - AveragedTangGridVelocity[iMarker])/AveragedSoundSpeed[iMarker];
						FlowAngle[iMarker]= atan((AveragedTangVelocity[iMarker] - AveragedTangGridVelocity[iMarker])/AveragedNormalVelocity[iMarker]);

					}else{
						AveragedMach[iMarker] = 0.0;
						for (iDim = 0; iDim < nDim; iDim++) {
							AveragedMach[iMarker] += AveragedVelocity[iMarker][iDim]*AveragedVelocity[iMarker][iDim];
						}
						AveragedMach[iMarker] = sqrt(AveragedMach[iMarker])/AveragedSoundSpeed[iMarker];
						AveragedTangMach[iMarker] = AveragedTangVelocity[iMarker]/AveragedSoundSpeed[iMarker];

					}

					AveragedNormalMach[iMarker] = AveragedNormalVelocity[iMarker]/AveragedSoundSpeed[iMarker];


#ifdef HAVE_MPI
					if ((AveragedDensity[iMarker]!= AveragedDensity[iMarker]) || (AveragedEnthalpy[iMarker]!=AveragedEnthalpy[iMarker])){
						if(size > 1 && rank == MASTER_NODE) cout<<"nan in mixing process in boundary "<<config->GetMarker_All_TagBound(iMarker)<< endl;
						else cout<<"nan in mixing process in boundary "<<config->GetMarker_All_TagBound(iMarker)<< endl;
					}
#else
					if ((AveragedDensity[iMarker]!= AveragedDensity[iMarker])||(AveragedEnthalpy[iMarker]!=AveragedEnthalpy[iMarker]))
						cout<<"nan in mixing process in boundary "<<config->GetMarker_All_TagBound(iMarker)<< endl;
#endif
				}
  /*--- Free locally allocated memory ---*/
//  MPI_Comm_free(&marker_comm);
  delete [] Velocity;
  delete [] Normal;
  delete [] TotalVelocity;
  delete [] TotalAreaVelocity;
  delete [] TotalFluxes;
  delete [] TotalNormal;
  delete [] TotalGridVel;
}


void CEulerSolver::MixedOut_Average (su2double val_init_pressure, su2double *val_Averaged_Flux, su2double *val_normal,
                                     su2double *pressure_mix, su2double *density_mix) {
  
  unsigned short maxiter = 10;
  unsigned short iter = 0;
  su2double toll = 1.0e-07;
  su2double resdl = 0.0;
  
  su2double *val_func = new su2double, *val_right_func = new su2double, *val_left_func = new su2double;
  su2double deltaP, *p_mix = new su2double, *p_mix_right = new su2double, *p_mix_left = new su2double;
  su2double epsilon = 1.0e-04;
  su2double relax_factor = 1;
  
  *pressure_mix = val_init_pressure;
  
  /*--- Newton-Raphson's method with central difference formula ---*/
  
  while ( iter <= maxiter ) {
    deltaP = 2*epsilon*(*pressure_mix);
    *p_mix_right = *pressure_mix+deltaP/2;
    *p_mix_left = *pressure_mix-deltaP/2;
    *p_mix = *pressure_mix;
    MixedOut_Root_Function(p_mix_right,val_Averaged_Flux,val_normal,val_right_func,density_mix);
    MixedOut_Root_Function(p_mix_left,val_Averaged_Flux,val_normal,val_left_func,density_mix);
    MixedOut_Root_Function(p_mix,val_Averaged_Flux,val_normal,val_func,density_mix);
    su2double der_func = (*val_right_func-*val_left_func) / deltaP;
    deltaP = -*val_func/der_func;
    resdl = deltaP/val_init_pressure;
    *pressure_mix += relax_factor*(deltaP);
    
    iter += 1;
    if ( abs(resdl) <= toll ) {
      break;
    }
    
  }
  
  MixedOut_Root_Function(pressure_mix,val_Averaged_Flux,val_normal,val_func,density_mix);
  
  /*--- Free locally allocated memory ---*/
  delete val_func;
  delete val_right_func;
  delete val_left_func;
  delete p_mix;
  delete p_mix_right;
  delete p_mix_left;
  
}

void CEulerSolver::MixedOut_Root_Function(su2double *pressure, su2double *val_Averaged_Flux, su2double *val_normal, su2double *valfunc, su2double *density) {
  
  su2double velnormal, velsq;
  
  su2double *vel;
  vel = new su2double[nDim];
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  *valfunc = 0.0;
  *density = 0.0;
  
  velnormal = 0.0;
  velsq = 0.0;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    vel[iDim]  = (val_Averaged_Flux[iDim+1] - (*pressure)*val_normal[iDim]) / val_Averaged_Flux[0];
    velnormal += val_normal[iDim]*vel[iDim];
    velsq += vel[iDim]*vel[iDim];
  }
  *density = val_Averaged_Flux[0] / velnormal;
  if (*density <= 0){
#ifdef HAVE_MPI
		if(size > 1 && rank == MASTER_NODE)cout << " desnity in mixedout routine negative : " << endl;
		else cout << " desnity in mixedout routine negative : " << endl;

#else
  	cout << " desnity in mixedout routine negative : " << endl;
#endif
  }
  FluidModel->SetTDState_Prho(*pressure, *density);
  su2double enthalpy = FluidModel->GetStaticEnergy() + (*pressure)/(*density);
  *valfunc = val_Averaged_Flux[nDim+1]/val_Averaged_Flux[0] - enthalpy - velsq/2;
  if (*valfunc!=*valfunc){
#ifdef HAVE_MPI
		if(size > 1 && rank == MASTER_NODE) cout << " mixedout root func gives nan: " << endl;
		else cout << " mixedout root func gives nan: " << endl;

#else
  cout << " mixedout root func gives nan: " << endl;
#endif
  }
  /*--- Free locally allocated memory ---*/
  delete [] vel;
  
}

void CEulerSolver::Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> >& c4k,signed long& nboundaryvertex) {
  /* Implementation of Fuorier Transformations for non-regfelcting BC will come soon */
}

void CEulerSolver::Boundary_Fourier(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short val_Marker, vector<std::complex<su2double> >& c2k,vector<std::complex<su2double> >& c3k,signed long& nboundaryvertex) {
  /* Implementation of Fuorier Transformations for non-regfelcting BC will come soon */
}

void CEulerSolver::BC_NonReflecting(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, kVar;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double  Area, UnitNormal[3];
  
  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
  su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  su2double Pressure_e;
  su2double ProjVelocity_i;
  su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_b;
  su2double *gridVel;
  su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  su2double *Normal;
  
  Normal = new su2double[nDim];
  
  Velocity_i = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  
  
  Lambda_i = new su2double[nVar];
  
  u_b = new su2double[nVar];
  dw = new su2double[nVar];
  
  S_boundary = new su2double[8];
  
  P_Tensor = new su2double*[nVar];
  invP_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new su2double[nVar];
    invP_Tensor[iVar] = new su2double[nVar];
  }
  
  /*--- new declarations ---*/
  std::vector<std::complex<su2double> > c4k ;//    std::complex<su2double> c3k[nVertex-OddEven]=0;
  std::vector<std::complex<su2double> > c2k ;//    std::complex<su2double> c3k[nVertex-OddEven]=0;
  std::vector<std::complex<su2double> > c3k ;//    std::complex<su2double> c3k[nVertex-OddEven]=0;
  
  su2double  deltaDensity, deltaPressure, AvgMach, deltaTangVelocity, deltaNormalVelocity, cc,rhoc,c1j,c2j,c3j,c4j,
  avg_c1, avg_c2, avg_c3, avg_c4,TangVelocity, NormalVelocity, GilesBeta, c4js, dc4js, *delta_c, **R_Matrix, *deltaprim;
  
  
  delta_c = new su2double[nVar];
  deltaprim = new su2double[nVar];
  R_Matrix= new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    R_Matrix[iVar] = new su2double[nVar];
  }
  
  
//  Mixing_Process(geometry, solver_container,  config, val_marker);
  
  cc = AveragedSoundSpeed[val_marker]*AveragedSoundSpeed[val_marker];
  rhoc = AveragedSoundSpeed[val_marker]*AveragedDensity[val_marker];
  AvgMach = AveragedMach[val_marker];
  
  conv_numerics->GetRMatrix(AveragedSoundSpeed[val_marker], AveragedDensity[val_marker], AveragedNormal[val_marker], R_Matrix);
  
  //  Boundary_Fourier(geometry, solver_container, config, val_marker, c4k, nboundaryvertex);
  //  Boundary_Fourier(geometry, solver_container, config, val_marker, c2k,c3k,nboundaryvertex);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    V_boundary= GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Compute the internal state u_i ---*/
      Velocity2_i = 0;
      for (iDim = 0; iDim < nDim; iDim++)
      {
        Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
        Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }
      
      
      Density_i = node[iPoint]->GetDensity();
      
      Energy_i = node[iPoint]->GetEnergy();
      StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
      
      FluidModel->SetTDState_rhoe(Density_i, StaticEnergy_i);
      
      Pressure_i = FluidModel->GetPressure();
      Enthalpy_i = Energy_i + Pressure_i/Density_i;
      
      SoundSpeed_i = FluidModel->GetSoundSpeed();
      
      Kappa_i = FluidModel->GetdPde_rho() / Density_i;
      Chi_i = FluidModel->GetdPdrho_e() - Kappa_i * StaticEnergy_i;
      
      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      
      
      switch(config->GetKind_Data_NRBC(Marker_Tag))
      {
          
          //TODO(turbo), generilize for 3D case
          //TODO(turbo), generilize for Inlet and Outlet in for backflow treatment
          //TODO(turbo), implement not uniform inlet and radial equilibrium for the outlet
          
        case MIXING_IN:
          
          /* --- Compute jump of primitive variable  --- */
          deltaDensity = ExtAveragedDensity[val_marker] - AveragedDensity[val_marker];
          deltaPressure = ExtAveragedPressure[val_marker] - AveragedPressure[val_marker];
          NormalVelocity= UnitNormal[0]*Velocity_i[0] + UnitNormal[1]*Velocity_i[1];
          deltaTangVelocity= ExtAveragedTangVelocity[val_marker]+AveragedTangVelocity[val_marker];
          deltaNormalVelocity= ExtAveragedNormalVelocity[val_marker]+AveragedNormalVelocity[val_marker];
          
          /* --- Compute characteristic jumps  --- */
          avg_c1= -cc*deltaDensity +deltaPressure;
          avg_c2= (rhoc*deltaTangVelocity);
          avg_c3= (rhoc*deltaNormalVelocity +deltaPressure);
          c4j= -rhoc*(-NormalVelocity +AveragedNormalVelocity[val_marker]) +(Pressure_i - AveragedPressure[val_marker]);
          
          /* --- Impose Inlet BC  --- */
          delta_c[0] = avg_c1;
          delta_c[1] = avg_c2;
          delta_c[2] = avg_c3;
          delta_c[3] = c4j;
          break;
          
        case MIXING_OUT:
          
          /* --- Compute jump of primitive variable  --- */
          deltaDensity = Density_i - AveragedDensity[val_marker];
          deltaPressure = Pressure_i - AveragedPressure[val_marker];
          TangVelocity= UnitNormal[0]*Velocity_i[1] - UnitNormal[1]*Velocity_i[0];
          NormalVelocity= UnitNormal[0]*Velocity_i[0] + UnitNormal[1]*Velocity_i[1];
          deltaTangVelocity= TangVelocity - AveragedTangVelocity[val_marker];
          deltaNormalVelocity= NormalVelocity - AveragedNormalVelocity[val_marker];
          
          /* --- Compute characteristic jumps  --- */
          c1j= -cc*deltaDensity +deltaPressure;
          c2j= rhoc*deltaTangVelocity;
          c3j= rhoc*deltaNormalVelocity + deltaPressure;
          avg_c4 = rhoc*(AveragedNormalVelocity[val_marker]+ExtAveragedNormalVelocity[val_marker]) -(AveragedPressure[val_marker]-ExtAveragedPressure[val_marker]);
          
          /* --- implementation of supersonic NRBC ---*/
          if (AvgMach > 1.001){
            if (AveragedTangVelocity[val_marker] >= 0.0){
              GilesBeta= -sqrt(pow(AvgMach,2)-1.0);
            }else{
              GilesBeta= sqrt(pow(AvgMach,2)-1.0);
            }
            c4js= (2.0 * AveragedNormalMach[val_marker])/(GilesBeta - AveragedTangMach[val_marker])*c2j - (GilesBeta+AveragedTangMach[val_marker])/(GilesBeta-AveragedTangMach[val_marker])*c3j;
            dc4js = c4js;
          }else{
            dc4js = 0.0;
          }
          
          /* --- Impose Outlet BC  --- */
          delta_c[0] = c1j;
          delta_c[1] = c2j;
          delta_c[2] = c3j;
          delta_c[3] = avg_c4 + dc4js;
          break;
          
        case STATIC_PRESSURE:
          
          Pressure_e = config->GetNRBC_Var1(Marker_Tag);
          Pressure_e /= config->GetPressure_Ref();
          
          /* --- Compute jump of primitive variable  --- */
          deltaDensity = Density_i - AveragedDensity[val_marker];
          deltaPressure = Pressure_i - AveragedPressure[val_marker];
          TangVelocity= UnitNormal[0]*Velocity_i[1] - UnitNormal[1]*Velocity_i[0];
          NormalVelocity= UnitNormal[0]*Velocity_i[0] + UnitNormal[1]*Velocity_i[1];
          deltaTangVelocity= TangVelocity - AveragedTangVelocity[val_marker];
          deltaNormalVelocity= NormalVelocity - AveragedNormalVelocity[val_marker];
          
          /* --- Compute characteristic jumps  --- */
          c1j= -cc*deltaDensity +deltaPressure;
          c2j= rhoc*deltaTangVelocity;
          c3j=rhoc*deltaNormalVelocity + deltaPressure;
          c4j=-rhoc*deltaNormalVelocity + deltaPressure;
          avg_c4 = -2.0*(AveragedPressure[val_marker]-Pressure_e);
          
          /* --- implementation of supersonic NRBC ---*/
          if (AvgMach > 1.001){
            if (AveragedTangVelocity[val_marker] >= 0.0){
              GilesBeta= -sqrt(pow(AvgMach,2)-1.0);
            }else{
              GilesBeta= sqrt(pow(AvgMach,2)-1.0);
            }
            c4js= (2.0 * AveragedNormalMach[val_marker])/(GilesBeta - AveragedTangMach[val_marker])*c2j - (GilesBeta+AveragedTangMach[val_marker])/(GilesBeta-AveragedTangMach[val_marker])*c3j;
            dc4js = c4js;
          }else{
            dc4js = 0.0;
          }
          
          /* --- Impose Outlet BC  --- */
          delta_c[0] = c1j;
          delta_c[1] = c2j;
          delta_c[2] = c3j;
          delta_c[3] = avg_c4 + dc4js;
          break;
          
        default:
          cout << "Warning! Invalid NRBC input!" << endl;
          exit(EXIT_FAILURE);
          break;
          
      }
      
      /*--- Compute primitive jump from characteristic variables  ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        deltaprim[iVar]=0;
        for (jVar = 0; jVar < nVar; jVar++)
        {
          deltaprim[iVar] +=  R_Matrix[iVar][jVar]*delta_c[jVar];
        }
      }
      
      /*--- Compute P (matrix of right eigenvectors) ---*/
      conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);
      
      /*--- Compute inverse P (matrix of left eigenvectors)---*/
      conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);
      
      /*--- eigenvalues contribution due to grid motion ---*/
      if (grid_movement){
        gridVel = geometry->node[iPoint]->GetGridVel();
        su2double ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel   += gridVel[iDim]*UnitNormal[iDim];
        ProjVelocity_i -= ProjGridVel;
      }
      
      /*--- Flow eigenvalues ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Lambda_i[iDim] = ProjVelocity_i;
      Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
      Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;
      
      //TODO(turbo), provide the under relaxation factor sigma from cfg file
      su2double sigma;
      sigma = 1.0;
      
      /*--- retrieve boundary variables ---*/
      Density_b = AveragedDensity[val_marker] + sigma*deltaprim[0];
      Pressure_b = AveragedPressure[val_marker] + sigma*deltaprim[3];
      switch(config->GetKind_Data_NRBC(Marker_Tag)){
        case MIXING_IN:
          NormalVelocity = AveragedNormalVelocity[val_marker] - sigma*deltaprim[1];
          TangVelocity = AveragedTangVelocity[val_marker] - sigma*deltaprim[2];
          break;
        case MIXING_OUT: case STATIC_PRESSURE:
          NormalVelocity = AveragedNormalVelocity[val_marker] + sigma*deltaprim[1];
          TangVelocity = AveragedTangVelocity[val_marker] + sigma*deltaprim[2];
          break;
        default:
          cout << "Warning! Invalid NRBC input!" << endl;
          exit(EXIT_FAILURE);
          break;
      }
      
      Velocity_b[0] = NormalVelocity*UnitNormal[0] - TangVelocity*UnitNormal[1];
      Velocity_b[1]	= NormalVelocity*UnitNormal[1] + TangVelocity*UnitNormal[0];
      Velocity2_b = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2_b+= Velocity_b[iDim]*Velocity_b[iDim];
      }
      
      FluidModel->SetTDState_Prho(Pressure_b, Density_b);
      Energy_b = FluidModel->GetStaticEnergy() + 0.5*Velocity2_b;
      StaticEnergy_b = FluidModel->GetStaticEnergy();
      Temperature_b= FluidModel->GetTemperature();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;
      Kappa_b = FluidModel->GetdPde_rho() / Density_b;
      Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
      
      /*--- Compute the thermodynamic state in u_b ---*/
      u_b[0]=Density_b;
      u_b[1]=Density_b*Velocity_b[0];
      u_b[2]=Density_b*Velocity_b[1];
      u_b[3]=Energy_b*Density_b;
      
      /*--- Compute the residuals ---*/
      conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);
      
      /*--- Residual contribution due to grid motion ---*/
      if (grid_movement) {
        gridVel = geometry->node[iPoint]->GetGridVel();
        su2double projVelocity = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          projVelocity +=  gridVel[iDim]*Normal[iDim];
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] -= projVelocity *(u_b[iVar]);
      }
      
      if (implicit) {
        /*--- Residual contribution due to grid motion ---*/
        Jacobian_b = new su2double*[nVar];
        DubDu = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++)
        {
          Jacobian_b[iVar] = new su2double[nVar];
          DubDu[iVar] = new su2double[nVar];
        }
        
        /*--- Initialize DubDu to unit matrix---*/
        for (iVar = 0; iVar < nVar; iVar++)
        {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0;
          
          DubDu[iVar][iVar]= 1;
        }
        
        /*--- Compute DubDu -= RNL---*/
        for (iVar=0; iVar<nVar; iVar++)
        {
          for (jVar=0; jVar<nVar; jVar++)
          {
            for (kVar=0; kVar<nVar; kVar++)
            {
              if (Lambda_i[kVar]<0)
                DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
            }
          }
        }
        
        /*--- Compute flux Jacobian in state b ---*/
        conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);
        
        /*--- Jacobian contribution due to grid motion ---*/
        if (grid_movement)
        {
          su2double projVelocity = 0.0;
          gridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            projVelocity +=  gridVel[iDim]*Normal[iDim];
          for (iVar = 0; iVar < nVar; iVar++){
            Residual[iVar] -= projVelocity *(u_b[iVar]);
            Jacobian_b[iVar][iVar] -= projVelocity;
          }
          
        }
        
        /*--- initiate Jacobian_i to zero matrix ---*/
        for (iVar=0; iVar<nVar; iVar++)
          for (jVar=0; jVar<nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        /*--- Compute numerical flux Jacobian at node i ---*/
        
        for (iVar=0; iVar<nVar; iVar++) {
          for (jVar=0; jVar<nVar; jVar++) {
            for (kVar=0; kVar<nVar; kVar++) {
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
            }
          }
          
        }
        
        for (iVar = 0; iVar < nVar; iVar++) {
          delete [] Jacobian_b[iVar];
          delete [] DubDu[iVar];
        }
        delete [] Jacobian_b;
        delete [] DubDu;
      }
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Primitive variables, using the derived quantities ---*/
        V_boundary[0] = Temperature_b;
        for (iDim = 0; iDim < nDim; iDim++)
          V_boundary[iDim+1] = Velocity_b[iDim];
        V_boundary[nDim+1] = Pressure_b;
        V_boundary[nDim+2] = Density_b;
        V_boundary[nDim+3] = Enthalpy_b;
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_boundary[nDim+5] = FluidModel->GetLaminarViscosity();
        V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();
        V_boundary[nDim+7] = FluidModel->GetThermalConductivity();
        V_boundary[nDim+8] = FluidModel->GetCp();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_boundary);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Secondary variables ---*/
        S_domain = node[iPoint]->GetSecondary();
        
        /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
        
        S_boundary[0]= FluidModel->GetdPdrho_e();
        S_boundary[1]= FluidModel->GetdPde_rho();
        
        S_boundary[2]= FluidModel->GetdTdrho_e();
        S_boundary[3]= FluidModel->GetdTde_rho();
        
        /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
        
        S_boundary[4]= FluidModel->Getdmudrho_T();
        S_boundary[5]= FluidModel->GetdmudT_rho();
        
        S_boundary[6]= FluidModel->Getdktdrho_T();
        S_boundary[7]= FluidModel->GetdktdT_rho();
        
        visc_numerics->SetSecondary(S_domain, S_boundary);
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
  delete [] Velocity_b;
  delete [] Velocity_i;
  
  delete [] S_boundary;
  delete [] Lambda_i;
  delete [] u_b;
  delete [] dw;
  
  
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
  
  delete []	delta_c;
  delete []	deltaprim;
  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] R_Matrix[iVar];
  }
  delete [] R_Matrix;
  
  
}

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3];
  su2double *V_inlet, *V_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  su2double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
  su2double Gas_Constant       = config->GetGas_ConstantND();
  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the inlet ---*/
    
    V_inlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/
      
      if (compressible) {
        
        /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/
        
        switch (Kind_Inlet) {
            
            /*--- Total properties have been specified at the inlet. ---*/
            
          case TOTAL_CONDITIONS:
            
            /*--- Retrieve the specified total conditions for this inlet. ---*/
            
            if (gravity) P_Total = config->GetInlet_Ptotal(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
            else P_Total  = config->GetInlet_Ptotal(Marker_Tag);
            T_Total  = config->GetInlet_Ttotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
            
            /*--- Non-dim. the inputs if necessary. ---*/
            
            P_Total /= config->GetPressure_Ref();
            T_Total /= config->GetTemperature_Ref();
            
            /*--- Store primitives and set some variables for clarity. ---*/
            
            Density = V_domain[nDim+2];
            Velocity2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity[iDim] = V_domain[iDim+1];
              Velocity2 += Velocity[iDim]*Velocity[iDim];
            }
            Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
            Pressure    = V_domain[nDim+1];
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
            
            alpha = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              alpha += UnitNormal[iDim]*Flow_Dir[iDim];
            
            /*--- Coefficients in the quadratic equation for the velocity ---*/
            
            aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
            bb = -1.0*Gamma_Minus_One*alpha*Riemann;
            cc =  0.5*Gamma_Minus_One*Riemann*Riemann
            -2.0*SoundSpeed_Total2/Gamma_Minus_One;
            
            /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/
            
            dd = bb*bb - 4.0*aa*cc;
            dd = sqrt(max(0.0, dd));
            Vel_Mag   = (-bb + dd)/(2.0*aa);
            Vel_Mag   = max(0.0, Vel_Mag);
            Velocity2 = Vel_Mag*Vel_Mag;
            
            /*--- Compute speed of sound from total speed of sound eqn. ---*/
            
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
            
            /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
            
            Mach2 = Velocity2/SoundSpeed2;
            Mach2 = min(1.0, Mach2);
            Velocity2   = Mach2*SoundSpeed2;
            Vel_Mag     = sqrt(Velocity2);
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
            
            /*--- Compute new velocity vector at the inlet ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
            
            /*--- Static temperature from the speed of sound relation ---*/
            
            Temperature = SoundSpeed2/(Gamma*Gas_Constant);
            
            /*--- Static pressure using isentropic relation at a point ---*/
            
            Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);
            
            /*--- Density at the inlet from the gas law ---*/
            
            Density = Pressure/(Gas_Constant*Temperature);
            
            /*--- Using pressure, density, & velocity, compute the energy ---*/
            
            Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
            if (tkeNeeded) Energy += GetTke_Inf();
            
            /*--- Primitive variables, using the derived quantities ---*/
            
            V_inlet[0] = Temperature;
            for (iDim = 0; iDim < nDim; iDim++)
              V_inlet[iDim+1] = Velocity[iDim];
            V_inlet[nDim+1] = Pressure;
            V_inlet[nDim+2] = Density;
            V_inlet[nDim+3] = Energy + Pressure/Density;
            
            break;
            
            /*--- Mass flow has been specified at the inlet. ---*/
            
          case MASS_FLOW:
            
            /*--- Retrieve the specified mass flow for the inlet. ---*/
            
            Density  = config->GetInlet_Ttotal(Marker_Tag);
            Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
            
            /*--- Non-dim. the inputs if necessary. ---*/
            
            Density /= config->GetDensity_Ref();
            Vel_Mag /= config->GetVelocity_Ref();
            
            /*--- Get primitives from current inlet state. ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
            Pressure    = node[iPoint]->GetPressure();
            SoundSpeed2 = Gamma*Pressure/V_domain[nDim+2];
            
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
            
            Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Vel_Mag*Vel_Mag;
            if (tkeNeeded) Energy += GetTke_Inf();
            
            /*--- Primitive variables, using the derived quantities ---*/
            
            V_inlet[0] = Pressure / ( Gas_Constant * Density);
            for (iDim = 0; iDim < nDim; iDim++)
              V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
            V_inlet[nDim+1] = Pressure;
            V_inlet[nDim+2] = Density;
            V_inlet[nDim+3] = Energy + Pressure/Density;
            
            break;
        }
      }
      if (incompressible) {
        
        /*--- Retrieve the specified velocity for the inlet. ---*/
        
        Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag)/config->GetVelocity_Ref();
        Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
        
        /*--- Store the velocity in the primitive variable vector ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
        
        /*--- Neumann condition for pressure ---*/
        
        V_inlet[0] = node[iPoint]->GetPressureInc();
        
        /*--- Constant value of density ---*/
        
        V_inlet[nDim+1] = GetDensity_Inf();
        
        /*--- Beta coefficient from the config file ---*/
        
        V_inlet[nDim+2] = config->GetArtComp_Factor();
        
      }
      if (freesurface) {
        
        /*--- Neumann condition for pressure, density, level set, and distance ---*/
        
        V_inlet[0] = node[iPoint]->GetPressureInc();
        V_inlet[nDim+1] = node[iPoint]->GetDensityInc();
        V_inlet[nDim+5] = node[iPoint]->GetLevelSet();
        V_inlet[nDim+6] = node[iPoint]->GetDistance();
        
        /*--- The velocity is computed from the infinity values ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          V_inlet[iDim+1] = GetVelocity_Inf(iDim);
        }
        
        /*--- The y/z velocity is interpolated due to the
         free surface effect on the pressure ---*/
        
        V_inlet[nDim] = node[iPoint]->GetPrimitive(nDim);
        
        /*--- Neumann condition for artifical compresibility factor ---*/
        
        V_inlet[nDim+2] = config->GetArtComp_Factor();
        
      }
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        if (compressible) {
          V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
          V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        }
        if (incompressible || freesurface) {
          V_inlet[nDim+3] = node[iPoint]->GetLaminarViscosityInc();
          V_inlet[nDim+4] = node[iPoint]->GetEddyViscosityInc();
        }
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double LevelSet, Density_Outlet = 0.0, Pressure, P_Exit, Velocity[3],
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
  Area, UnitNormal[3], Height, yCoordRef, yCoord;
  su2double *V_outlet, *V_domain;
  
  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  su2double Gas_Constant     = config->GetGas_ConstantND();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement      = config->GetGrid_Movement();
  su2double FreeSurface_Zero = config->GetFreeSurface_Zero();
  su2double epsilon          = config->GetFreeSurface_Thickness();
  su2double RatioDensity     = config->GetRatioDensity();
  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  su2double PressFreeSurface = GetPressure_Inf();
  su2double Froude           = config->GetFroude();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/
      if (compressible) {
        
        /*--- Retrieve the specified back pressure for this outlet. ---*/
        if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;
        else P_Exit = config->GetOutlet_Pressure(Marker_Tag);
        
        /*--- Non-dim. the inputs if necessary. ---*/
        P_Exit = P_Exit/config->GetPressure_Ref();
        
        /*--- Check whether the flow is supersonic at the exit. The type
         of boundary update depends on this. ---*/
        Density = V_domain[nDim+2];
        Velocity2 = 0.0; Vn = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = V_domain[iDim+1];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
          Vn += Velocity[iDim]*UnitNormal[iDim];
        }
        Pressure   = V_domain[nDim+1];
        SoundSpeed = sqrt(Gamma*Pressure/Density);
        Mach_Exit  = sqrt(Velocity2)/SoundSpeed;
        
        if (Mach_Exit >= 1.0) {
          
          /*--- Supersonic exit flow: there are no incoming characteristics,
           so no boundary condition is necessary. Set outlet state to current
           state so that upwinding handles the direction of propagation. ---*/
          for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];
          
        } else {
          
          /*--- Subsonic exit flow: there is one incoming characteristic,
           therefore one variable can be specified (back pressure) and is used
           to update the conservative variables. Compute the entropy and the
           acoustic Riemann variable. These invariants, as well as the
           tangential velocity components, are extrapolated. Adapted from an
           original implementation in the Stanford University multi-block
           (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
           der Weide, last modified 09-10-2007. ---*/
          
          Entropy = Pressure*pow(1.0/Density, Gamma);
          Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
          
          /*--- Compute the new fictious state at the outlet ---*/
          Density    = pow(P_Exit/Entropy,1.0/Gamma);
          Pressure   = P_Exit;
          SoundSpeed = sqrt(Gamma*P_Exit/Density);
          Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
          Velocity2  = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Energy = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;
          if (tkeNeeded) Energy += GetTke_Inf();
          
          /*--- Conservative variables, using the derived quantities ---*/
          V_outlet[0] = Pressure / ( Gas_Constant * Density);
          for (iDim = 0; iDim < nDim; iDim++)
            V_outlet[iDim+1] = Velocity[iDim];
          V_outlet[nDim+1] = Pressure;
          V_outlet[nDim+2] = Density;
          V_outlet[nDim+3] = Energy + Pressure/Density;
          
        }
      }
      if (incompressible) {
        
        /*--- The pressure is computed from the infinity values ---*/
        if (gravity) {
          yCoordRef = 0.0;
          yCoord = geometry->node[iPoint]->GetCoord(nDim-1);
          V_outlet[0] = GetPressure_Inf() + GetDensity_Inf()*((yCoordRef-yCoord)/(config->GetFroude()*config->GetFroude()));
        }
        else {
          V_outlet[0] = GetPressure_Inf();
        }
        
        /*--- Neumann condition for the velocity ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          V_outlet[iDim+1] = node[Point_Normal]->GetPrimitive(iDim+1);
        }
        
        /*--- Constant value of density ---*/
        V_outlet[nDim+1] = GetDensity_Inf();
        
        /*--- Beta coefficient from the config file ---*/
        V_outlet[nDim+2] = config->GetArtComp_Factor();
        
      }
      if (freesurface) {
        
        /*--- Imposed pressure, density, level set and distance ---*/
        Height = geometry->node[iPoint]->GetCoord(nDim-1);
        LevelSet = Height - FreeSurface_Zero;
        if (LevelSet < -epsilon) Density_Outlet = config->GetDensity_FreeStreamND();
        if (LevelSet > epsilon) Density_Outlet = RatioDensity*config->GetDensity_FreeStreamND();
        V_outlet[0] = PressFreeSurface + Density_Outlet*((FreeSurface_Zero-Height)/(Froude*Froude));
        V_outlet[nDim+1] = Density_Outlet;
        V_outlet[nDim+5] = LevelSet;
        V_outlet[nDim+6] = LevelSet;
        
        /*--- Neumann condition in the interface for the pressure, density and level set and distance ---*/
        if (fabs(LevelSet) <= epsilon) {
          V_outlet[0] = node[Point_Normal]->GetPressureInc();
          V_outlet[nDim+1] = node[Point_Normal]->GetDensityInc();
          V_outlet[nDim+5] = node[Point_Normal]->GetLevelSet();
          V_outlet[nDim+6] = node[Point_Normal]->GetDistance();
        }
        
        /*--- Neumann condition for the velocity ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          V_outlet[iDim+1] = node[Point_Normal]->GetPrimitive(iDim+1);
        }
        
        /*--- Neumann condition for artifical compresibility factor ---*/
        V_outlet[nDim+2] = config->GetArtComp_Factor();
        
      }
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        if (compressible) {
          V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
          V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        }
        if (incompressible || freesurface) {
          V_outlet[nDim+3] = node[iPoint]->GetLaminarViscosityInc();
          V_outlet[nDim+4] = node[iPoint]->GetEddyViscosityInc();
        }
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                       CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain;
  
  su2double Density, Pressure, Temperature, Energy, *Velocity, Velocity2;
  su2double Gas_Constant = config->GetGas_ConstantND();
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous              = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double *Normal = new su2double[nDim];
  
  /*--- Supersonic inlet flow: there are no outgoing characteristics,
   so all flow variables can be imposed at the inlet.
   First, retrieve the specified values for the primitive variables. ---*/
  
  Temperature = config->GetInlet_Temperature(Marker_Tag);
  Pressure    = config->GetInlet_Pressure(Marker_Tag);
  Velocity    = config->GetInlet_Velocity(Marker_Tag);
  
  /*--- Density at the inlet from the gas law ---*/
  
  Density = Pressure/(Gas_Constant*Temperature);
  
  /*--- Non-dim. the inputs if necessary. ---*/
  
  Temperature = Temperature/config->GetTemperature_Ref();
  Pressure    = Pressure/config->GetPressure_Ref();
  Density     = Density/config->GetDensity_Ref();
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = Velocity[iDim]/config->GetVelocity_Ref();
  
  /*--- Compute the energy from the specified state ---*/
  
  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity2 += Velocity[iDim]*Velocity[iDim];
  Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;
  if (tkeNeeded) Energy += GetTke_Inf();
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    
    V_inlet = GetCharacPrimVar(val_marker, iVertex);
    
    /*--- Primitive variables, using the derived quantities ---*/
    
    V_inlet[0] = Temperature;
    for (iDim = 0; iDim < nDim; iDim++)
      V_inlet[iDim+1] = Velocity[iDim];
    V_inlet[nDim+1] = Pressure;
    V_inlet[nDim+2] = Density;
    V_inlet[nDim+3] = Energy + Pressure/Density;
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                        CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous              = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Supersonic outlet flow: there are no ingoing characteristics,
   so all flow variables can should be interpolated from the domain. ---*/
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Primitive variables, using the derived quantities ---*/
      
      V_outlet[0] = V_domain[0];
      for (iDim = 0; iDim < nDim; iDim++)
        V_outlet[iDim+1] = V_domain[iDim+1];
      V_outlet[nDim+1] = V_domain[nDim+1];
      V_outlet[nDim+2] = V_domain[nDim+2];
      V_outlet[nDim+3] = V_domain[nDim+3];
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Pressure, Inflow_Pressure, Velocity[3], Velocity2, Entropy, Target_Inflow_Mach = 0.0, Density, Energy,
  Riemann, Area, UnitNormal[3], Vn, SoundSpeed, Vn_Exit, Inflow_Pressure_inc, Inflow_Pressure_old, Inflow_Mach_old;
  su2double *V_inflow, *V_domain;
  
  su2double DampingFactor = config->GetDamp_Engine_Inflow();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous              = config->GetViscous();
  su2double Gas_Constant = config->GetGas_ConstantND();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Retrieve the specified target fan face mach in the nacelle. ---*/
  
  Target_Inflow_Mach = config->GetInflow_Mach_Target(Marker_Tag);
  
  /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
  
  Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
  Inflow_Mach_old = config->GetInflow_Mach(Marker_Tag);
  
  /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/
  
  Inflow_Pressure_inc = - (1.0 - (Inflow_Mach_old/Target_Inflow_Mach)) * Baseline_Press;
  
  /*--- Estimate the new fan face pressure ---*/
  
  Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    
    V_inflow = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Subsonic nacelle inflow: there is one incoming characteristic,
       therefore one variable can be specified (back pressure) and is used
       to update the conservative variables.
       
       Compute the entropy and the acoustic variable. These
       riemann invariants, as well as the tangential velocity components,
       are extrapolated. ---*/
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Entropy = Pressure*pow(1.0/Density, Gamma);
      Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
      
      /*--- Compute the new fictious state at the outlet ---*/
      
      Density    = pow(Inflow_Pressure/Entropy,1.0/Gamma);
      Pressure   = Inflow_Pressure;
      SoundSpeed = sqrt(Gamma*Inflow_Pressure/Density);
      Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
      Velocity2  = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      
      Energy = Inflow_Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
      if (tkeNeeded) Energy += GetTke_Inf();
      
      /*--- Conservative variables, using the derived quantities ---*/
      
      V_inflow[0] = Pressure / ( Gas_Constant * Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_inflow[iDim+1] = Velocity[iDim];
      V_inflow[nDim+1] = Pressure;
      V_inflow[nDim+2] = Density;
      V_inflow[nDim+3] = Energy + Pressure/Density;
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inflow);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_inflow[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_inflow[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inflow);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Exhaust_Pressure, Exhaust_Temperature, Velocity[3], Velocity2, H_Exhaust, Temperature, Riemann, Area, UnitNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Exhaust2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
  su2double *V_exhaust, *V_domain;
  su2double Gas_Constant = config->GetGas_ConstantND();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  //  su2double Target_Exhaust_Pressure, Exhaust_Pressure_old, Exhaust_Pressure_inc;
  //  su2double DampingFactor = config->GetDamp_Engine_Exhaust();
  //  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Retrieve the specified exhaust pressure in the engine (non-dimensional). ---*/
  
  //  Target_Exhaust_Pressure = config->GetExhaust_Pressure_Target(Marker_Tag) / config->GetPressure_Ref();
  
  /*--- Retrieve the old exhaust pressure in the engine exhaust (this has been computed in a preprocessing). ---*/
  
  //  Exhaust_Pressure_old = config->GetExhaust_Pressure(Marker_Tag);
  
  /*--- Compute the Pressure increment ---*/
  
  //  Exhaust_Pressure_inc = (1.0 - (Exhaust_Pressure_old/Target_Exhaust_Pressure)) * Baseline_Press;
  
  /*--- Estimate the new exhaust pressure ---*/
  
  //  Exhaust_Pressure = (1.0 - DampingFactor) * Exhaust_Pressure_old + DampingFactor * (Exhaust_Pressure_old + Exhaust_Pressure_inc);
  
  /*--- The temperature is given (no iteration is required) ---*/
  
  Exhaust_Temperature  = config->GetExhaust_Temperature_Target(Marker_Tag);
  Exhaust_Temperature /= config->GetTemperature_Ref();
  
  /*--- The pressure is given (no iteration is required) ---*/
  /*--- CHECK: the above iterative process is overwritten on the next line. ---*/
  
  Exhaust_Pressure  = config->GetExhaust_Pressure_Target(Marker_Tag);
  Exhaust_Pressure /= config->GetPressure_Ref();
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the exhaust ---*/
    
    V_exhaust = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
       therefore we can specify all but one state variable at the inlet.
       The outgoing Riemann invariant provides the final piece of info. ---*/
      
      /*--- Store primitives and set some variables for clarity. ---*/
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure    = V_domain[nDim+1];
      H_Exhaust   = (Gamma*Gas_Constant/Gamma_Minus_One)*Exhaust_Temperature;
      SoundSpeed2 = Gamma*Pressure/Density;
      
      /*--- Compute the acoustic Riemann invariant that is extrapolated
       from the domain interior. ---*/
      
      Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
      for (iDim = 0; iDim < nDim; iDim++)
        Riemann += Velocity[iDim]*UnitNormal[iDim];
      
      /*--- Total speed of sound ---*/
      
      SoundSpeed_Exhaust2 = Gamma_Minus_One*(H_Exhaust - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
      
      /*--- The flow direction is defined by the surface normal ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir[iDim] = -UnitNormal[iDim];
      
      /*--- Dot product of normal and flow direction. This should
       be negative due to outward facing boundary normal convention. ---*/
      
      alpha = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        alpha += UnitNormal[iDim]*Flow_Dir[iDim];
      
      /*--- Coefficients in the quadratic equation for the velocity ---*/
      
      aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
      bb = -1.0*Gamma_Minus_One*alpha*Riemann;
      cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Exhaust2/Gamma_Minus_One;
      
      /*--- Solve quadratic equation for velocity magnitude. Value must
       be positive, so the choice of root is clear. ---*/
      
      dd      = bb*bb - 4.0*aa*cc;
      dd      = sqrt(max(0.0, dd));
      Vel_Mag = (-bb + dd)/(2.0*aa);
      
      if (Vel_Mag >= 0.0) {
        
        Velocity2 = Vel_Mag*Vel_Mag;
        
        /*--- Compute speed of sound from total speed of sound eqn. ---*/
        
        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
        Mach2       = Velocity2/SoundSpeed2;
        Velocity2   = Mach2*SoundSpeed2;
        Vel_Mag     = sqrt(Velocity2);
        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
        
        /*--- Compute new velocity vector at the inlet ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
        
        /*--- Static temperature from the speed of sound relation ---*/
        
        Temperature = SoundSpeed2/(Gamma*Gas_Constant);
        
        /*--- Static pressure using isentropic relation at a point ---*/
        
        Pressure = Exhaust_Pressure*pow((Temperature/Exhaust_Temperature), Gamma/Gamma_Minus_One);
        
        /*--- Density at the exhaust from the gas law ---*/
        
        Density = Pressure/(Gas_Constant*Temperature);
        
        /*--- Using pressure, density, & velocity, compute the energy ---*/
        
        Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Primitive variables, using the derived quantities ---*/
        
        V_exhaust[0] = Temperature;
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = Velocity[iDim];
        V_exhaust[nDim+1] = Pressure;
        V_exhaust[nDim+2] = Density;
        V_exhaust[nDim+3] = Energy + Pressure/Density;
        
      }
      
      /*--- The flow goes in the wrong direction ---*/
      
      else {
        
        V_exhaust[0] = V_domain[0];
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = V_domain[iDim+1];
        V_exhaust[nDim+1] = V_domain[nDim+1];
        V_exhaust[nDim+2] = V_domain[nDim+2];
        V_exhaust[nDim+3] = V_domain[nDim+3];
        
      }
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_exhaust);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_exhaust[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_exhaust[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_exhaust);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Engine_Bleed(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Bleed_Pressure, Target_Bleed_MassFlow, Bleed_Pressure_old, Bleed_MassFlow_old, Bleed_Pressure_inc, Bleed_Temperature, Velocity[3], Velocity2, H_Bleed, Temperature, Riemann, Area, UnitNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Bleed2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
  su2double *V_bleed, *V_domain;
  su2double Gas_Constant = config->GetGas_ConstantND();
  
  su2double DampingFactor = config->GetDamp_Engine_Bleed();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  su2double Baseline_Press = 0.25 * config->GetPressure_FreeStreamND();
  
  su2double *Normal = new su2double[nDim];
  
  /*--- Retrieve the specified target bleed mass flow in the engine (non-dimensional). ---*/
  
  Target_Bleed_MassFlow = config->GetBleed_MassFlow_Target(Marker_Tag) / (config->GetDensity_Ref() * config->GetVelocity_Ref());
  
  /*--- Retrieve the old Bleed pressure and mass flow rate number in the engine bleed (this has been computed in a preprocessing). ---*/
  
  Bleed_Pressure_old = config->GetBleed_Pressure(Marker_Tag);
  Bleed_MassFlow_old = config->GetBleed_MassFlow(Marker_Tag);
  
  /*--- Compute the Pressure increment (note that increasing pressure also increases mass flow rate) ---*/
  
  Bleed_Pressure_inc = (1.0 - (Bleed_MassFlow_old/Target_Bleed_MassFlow)) * Baseline_Press;
  
  /*--- Estimate the new bleed pressure ---*/
  
  Bleed_Pressure = (1.0 - DampingFactor)*Bleed_Pressure_old + DampingFactor * (Bleed_Pressure_old + Bleed_Pressure_inc);
  
  /*--- The temperature is given (no iteration is required) ---*/
  
  Bleed_Temperature  = config->GetBleed_Temperature_Target(Marker_Tag);
  Bleed_Temperature /= config->GetTemperature_Ref();
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the exhaust ---*/
    
    V_bleed = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();
      
      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
       therefore we can specify all but one state variable at the inlet.
       The outgoing Riemann invariant provides the final piece of info. ---*/
      
      /*--- Store primitives and set some variables for clarity. ---*/
      
      Density = V_domain[nDim+2];
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Energy = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure = V_domain[nDim+1];
      H_Bleed     = (Gamma*Gas_Constant/Gamma_Minus_One)*Bleed_Temperature;
      SoundSpeed2 = Gamma*Pressure/Density;
      
      /*--- Compute the acoustic Riemann invariant that is extrapolated
       from the domain interior. ---*/
      
      Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
      for (iDim = 0; iDim < nDim; iDim++)
        Riemann += Velocity[iDim]*UnitNormal[iDim];
      
      /*--- Total speed of sound ---*/
      
      SoundSpeed_Bleed2 = Gamma_Minus_One*(H_Bleed - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
      
      /*--- The flow direction is defined by the surface normal ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir[iDim] = -UnitNormal[iDim];
      
      /*--- Dot product of normal and flow direction. This should
       be negative due to outward facing boundary normal convention. ---*/
      
      alpha = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        alpha += UnitNormal[iDim]*Flow_Dir[iDim];
      
      /*--- Coefficients in the quadratic equation for the velocity ---*/
      
      aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
      bb = -1.0*Gamma_Minus_One*alpha*Riemann;
      cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Bleed2/Gamma_Minus_One;
      
      /*--- Solve quadratic equation for velocity magnitude. Value must
       be positive, so the choice of root is clear. ---*/
      
      dd = bb*bb - 4.0*aa*cc;
      dd = sqrt(max(0.0, dd));
      Vel_Mag   = (-bb + dd)/(2.0*aa);
      
      if (Vel_Mag >= 0.0) {
        
        Velocity2 = Vel_Mag*Vel_Mag;
        
        /*--- Compute speed of sound from total speed of sound eqn. ---*/
        
        SoundSpeed2 = SoundSpeed_Bleed2 - 0.5*Gamma_Minus_One*Velocity2;
        
        /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
        
        Mach2 = Velocity2/SoundSpeed2;
        Mach2 = min(1.0, Mach2);
        Velocity2   = Mach2*SoundSpeed2;
        Vel_Mag     = sqrt(Velocity2);
        SoundSpeed2 = SoundSpeed_Bleed2 - 0.5*Gamma_Minus_One*Velocity2;
        
        /*--- Compute new velocity vector at the inlet ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
        
        /*--- Static temperature from the speed of sound relation ---*/
        
        Temperature = SoundSpeed2/(Gamma*Gas_Constant);
        
        /*--- Static pressure using isentropic relation at a point ---*/
        
        Pressure = Bleed_Pressure*pow((Temperature/Bleed_Temperature), Gamma/Gamma_Minus_One);
        
        /*--- Density at the inlet from the gas law ---*/
        
        Density = Pressure/(Gas_Constant*Temperature);
        
        /*--- Using pressure, density, & velocity, compute the energy ---*/
        
        Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Primitive variables, using the derived quantities ---*/
        
        V_bleed[0] = Temperature;
        for (iDim = 0; iDim < nDim; iDim++)
          V_bleed[iDim+1] = Velocity[iDim];
        V_bleed[nDim+1] = Pressure;
        V_bleed[nDim+2] = Density;
        V_bleed[nDim+3] = Energy + Pressure/Density;
        
      }
      
      /*--- The flow goes in the wrong direction ---*/
      
      else {
        
        V_bleed[0] = V_domain[0];
        for (iDim = 0; iDim < nDim; iDim++)
          V_bleed[iDim+1] = V_domain[iDim+1];
        V_bleed[nDim+1] = V_domain[nDim+1];
        V_bleed[nDim+2] = V_domain[nDim+2];
        V_bleed[nDim+3] = V_domain[nDim+3];
        
      }
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_bleed);
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        
        V_bleed[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_bleed[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_bleed);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  delete [] Normal;
  
}

void CEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler residual ---*/
  
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
  
}

void CEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config) {
  
  unsigned long iVertex, iPoint, jPoint;
  unsigned short iDim, iVar, iMarker;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
#ifndef HAVE_MPI
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
          
          if (iPoint != jPoint) {
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
            }
            
            /*--- Set primitive variables ---*/
            
            numerics->SetPrimitive(PrimVar_i, PrimVar_j);
            
            /*--- Set the normal vector ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, Residual);
            if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
          }
          
        }
      }
    }
  }
  
#else
  
  int rank, jProcessor;
  MPI_Status status;
  //MPI_Status send_stat[1], recv_stat[1];
  //MPI_Request send_req[1], recv_req[1];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  bool compute;
  su2double *Buffer_Send_V = new su2double [nPrimVar];
  su2double *Buffer_Receive_V = new su2double [nPrimVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          /*--- We only send the information that belong to other boundary, -1 processor
           means that the boundary condition is not applied ---*/
          
          if (compute) {
            
            if (jProcessor != rank) {
              
              /*--- Copy the primitive variable ---*/
              
              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Send_V[iVar] = node[iPoint]->GetPrimitive(iVar);
              
              SU2_MPI::Bsend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD);
              
              //          SU2_MPI::Isend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD, &send_req[0]);
              
              //          /*--- Wait for this set of non-blocking comm. to complete ---*/
              //
              //          SU2_MPI::Waitall(1, send_req, send_stat);
              
            }
            
          }
          
        }
      }
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          if (compute) {
            
            /*--- We only receive the information that belong to other boundary ---*/
            
            if (jProcessor != rank) {
              
              SU2_MPI::Recv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &status);
              
              //         SU2_MPI::Irecv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &recv_req[0]);
              
              /*--- Wait for the this set of non-blocking recv's to complete ---*/
              
              //          SU2_MPI::Waitall(1, recv_req, recv_stat);
              
            }
            else {
              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Receive_V[iVar] = node[jPoint]->GetPrimitive(iVar);
            }
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = Buffer_Receive_V[iVar];
            }
            
            /*--- Set Conservative Variables ---*/
            
            numerics->SetPrimitive(PrimVar_i, PrimVar_j);
            
            /*--- Set Normal ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, Residual);
            if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
          }
          
        }
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  delete[] Buffer_Send_V;
  delete[] Buffer_Receive_V;
  
#endif
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config) {
  
  unsigned long iVertex, iPoint, jPoint;
  unsigned short iDim, iVar, iMarker;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
#ifndef HAVE_MPI
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
          
          if (iPoint != jPoint) {
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
            }
            
            /*--- Set primitive variables ---*/
            
            numerics->SetPrimitive(PrimVar_i, PrimVar_j);
            
            /*--- Set the normal vector ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, Residual);
            if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
          }
          
        }
      }
    }
  }
  
#else
  
  int rank, jProcessor;
  MPI_Status status;
  //MPI_Status send_stat[1], recv_stat[1];
  //MPI_Request send_req[1], recv_req[1];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  bool compute;
  su2double *Buffer_Send_V    = new su2double [nPrimVar];
  su2double *Buffer_Receive_V = new su2double [nPrimVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          /*--- We only send the information that belong to other boundary, -1 processor
           means that the boundary condition is not applied ---*/
          
          if (compute) {
            
            if (jProcessor != rank) {
              
              /*--- Copy the primitive variable ---*/
              
              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Send_V[iVar] = node[iPoint]->GetPrimitive(iVar);
              
              SU2_MPI::Bsend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD);
              
              //          SU2_MPI::Isend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD, &send_req[0]);
              
              //          /*--- Wait for this set of non-blocking comm. to complete ---*/
              //
              //          SU2_MPI::Waitall(1, send_req, send_stat);
              
            }
            
          }
          
        }
      }
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          if (compute) {
            
            /*--- We only receive the information that belong to other boundary ---*/
            
            if (jProcessor != rank) {
              
              SU2_MPI::Recv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &status);
              
              //         SU2_MPI::Irecv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &recv_req[0]);
              
              /*--- Wait for the this set of non-blocking recv's to complete ---*/
              
              //          SU2_MPI::Waitall(1, recv_req, recv_stat);
              
            }
            else {
              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Receive_V[iVar] = node[jPoint]->GetPrimitive(iVar);
            }
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
              PrimVar_j[iVar] = Buffer_Receive_V[iVar];
            }
            
            /*--- Set Conservative Variables ---*/
            
            numerics->SetPrimitive(PrimVar_i, PrimVar_j);
            
            /*--- Set Normal ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, Residual);
            if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
          }
          
        }
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  delete[] Buffer_Send_V;
  delete[] Buffer_Receive_V;
  
#endif
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CEulerSolver::BC_ActDisk_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config) {
  
  unsigned long iVertex, iPoint, jPoint, Pin = 0, Pout = 0, jProcessor;
  unsigned short iDim, iVar, iMarker;
  int iProcessor;
  su2double *Coord, radius, DeltaP_avg, DeltaP_tip, DeltaT_avg, DeltaT_tip,
  Radial[3] = {0.0,0.0,0.0}, Tangent[3] = {0.0,0.0,0.0}, Normal[3] = {0.0,0.0,0.0},
  UnitRadial[3] = {0.0,0.0,0.0}, UnitTangent[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};
  su2double H_in_ghost, P_in_ghost, Vel_Swirl_out_ghost, T_in_ghost, Rho_in_ghost, sos_in_ghost, Area, Vel_in_ghost;
  su2double H_out_ghost, P_out_ghost, T_out_ghost, Rho_out_ghost, sos_out_ghost, Vel_Normal_in, Rho_in, Vel_Normal_out_ghost, Vel_out_ghost;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  unsigned short nMarker_ActDisk_Inlet = config->GetnMarker_ActDisk_Inlet();
  //  su2double DampingFactor = 0.75;
  
  if (nMarker_ActDisk_Inlet != 0) {
    
    su2double *PrimVar_out = new su2double[nPrimVar];
    su2double *PrimVar_in = new su2double[nPrimVar];
    su2double *MeanPrimVar = new su2double[nPrimVar];
    su2double *PrimVar_out_ghost = new su2double[nPrimVar];
    su2double *PrimVar_in_ghost = new su2double[nPrimVar];
    su2double *ActDisk_Jump = new su2double[nPrimVar];
    su2double *Buffer_Send_V = new su2double [nPrimVar];
    su2double *Buffer_Receive_V = new su2double [nPrimVar];
    
#ifndef HAVE_MPI
    iProcessor = MASTER_NODE;
#else
    MPI_Status status;
    //MPI_Status send_stat[1], recv_stat[1];
    //MPI_Request send_req[1], recv_req[1];
    MPI_Comm_rank(MPI_COMM_WORLD, &iProcessor);
#endif
    
    /*--- Identify the points that should be sended in a MPI implementation.
     Do the send process, by the moment we are sending each node individually, this must be changed---*/
    
#ifdef HAVE_MPI
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            /*--- Find the associate pair to the original node ---*/
            
            jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
            jProcessor = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            /*--- We only send the information that belong to other boundary, using jPoint as the ID for the message  ---*/
            
            if ((int)jProcessor != iProcessor) {
              
              /*--- Copy the primitive variables ---*/
              
              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Send_V[iVar] = node[iPoint]->GetPrimitive(iVar);
              
              SU2_MPI::Bsend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD);
              
              //              SU2_MPI::Isend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD, &send_req[0]);
              //              SU2_MPI::Waitall(1, send_req, send_stat);
              
            }
            
          }
        }
      }
    }
    
#endif
    
    /*--- Evaluate the fluxes, the donor solution has been sended using MPI ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
        
        unsigned short boundary = config->GetMarker_All_KindBC(iMarker);
        su2double *Origin = config->GetActDisk_Origin(config->GetMarker_All_TagBound(iMarker));
        su2double R_root = config->GetActDisk_RootRadius(config->GetMarker_All_TagBound(iMarker));
        su2double R_tip = config->GetActDisk_TipRadius(config->GetMarker_All_TagBound(iMarker));
        su2double Target_Press_Jump = config->GetActDisk_PressJump(config->GetMarker_All_TagBound(iMarker))/ config->GetPressure_Ref();
        su2double Target_Temp_Jump = config->GetActDisk_TempJump(config->GetMarker_All_TagBound(iMarker))/ config->GetTemperature_Ref();
        su2double Omega = config->GetActDisk_Omega(config->GetMarker_All_TagBound(iMarker))*(PI_NUMBER/30.0);
        unsigned short Distribution = config->GetActDisk_Distribution(config->GetMarker_All_TagBound(iMarker));
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            /*--- Find the associate pair to the original node ---*/
            
            jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
            jProcessor = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            /*--- Receive the information, using jPoint as the ID for the message ---*/
            
            if ((int)jProcessor != iProcessor) {
#ifdef HAVE_MPI
              
              SU2_MPI::Recv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &status);
              
              //              SU2_MPI::Irecv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &recv_req[0]);
              //              SU2_MPI::Waitall(1, send_req, send_stat);
              
#endif
            }
            else {
              
              /*--- The point is in the same processor... no MPI required ---*/
              
              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Receive_V[iVar] = node[jPoint]->GetPrimitive(iVar);
              
            }
            
            /*--- Identify the inner and the outer point (based on the normal direction) ---*/
            
            if (boundary == ACTDISK_INLET) {
              
              Pin = iPoint;
              
              for (iVar = 0; iVar < nPrimVar; iVar++) {
                PrimVar_out[iVar] = Buffer_Receive_V[iVar];
                PrimVar_in[iVar] = node[Pin]->GetPrimitive(iVar);
              }
              
            }
            if (boundary == ACTDISK_OUTLET) {
              
              Pout = iPoint;
              
              for (iVar = 0; iVar < nPrimVar; iVar++) {
                PrimVar_out[iVar] = node[Pout]->GetPrimitive(iVar);
                PrimVar_in[iVar] = Buffer_Receive_V[iVar];
              }
              
            }
            
            /*--- Set the jump in the actuator disk ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              ActDisk_Jump[iVar] = 0.0;
            }
            
            /*--- Compute the distance to the center of the rotor ---*/
            
            Coord = geometry->node[iPoint]->GetCoord();
            radius = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              radius += (Coord[iDim]-Origin[iDim])*(Coord[iDim]-Origin[iDim]);
            radius = sqrt(radius);
            
            /*--- Compute the pressure jump ---*/
            
            //            su2double Press_Jump = fabs(PrimVar_in[nDim+1]-PrimVar_in[nDim+1]);
            //            su2double Press_inc = (1.0 - (Press_Jump/Target_Press_Jump)) * 0.01 * config->GetPressure_FreeStreamND();
            //            DeltaP_avg = (1.0 - DampingFactor)*Press_Jump + DampingFactor * (Press_Jump + Press_inc);
            
            DeltaP_avg = Target_Press_Jump;
            if (Distribution == 1) {
              DeltaP_tip = (3.0/2.0)*DeltaP_avg*R_tip*(R_tip*R_tip-R_root*R_root)/(R_tip*R_tip*R_tip-R_root*R_root*R_root);
              ActDisk_Jump[nDim+1] = DeltaP_tip*radius/R_tip;
            }
            else
              ActDisk_Jump[nDim+1] = DeltaP_avg;
            
            /*--- Compute the temperature jump ---*/
            
            //            su2double Temp_Jump = fabs(PrimVar_in[0]-PrimVar_in[0]);
            //            su2double Temp_inc = (1.0 - (Temp_Jump/Target_Temp_Jump)) * 0.01 * config->GetTemperature_FreeStreamND();
            //            DeltaT_avg = (1.0 - DampingFactor)*Temp_Jump + DampingFactor * (Temp_Jump + Temp_inc);
            
            DeltaT_avg = Target_Temp_Jump;
            if (Distribution == 1) {
              DeltaT_tip = (3.0/2.0)*DeltaT_avg*R_tip*(R_tip*R_tip-R_root*R_root)/(R_tip*R_tip*R_tip-R_root*R_root*R_root);
              ActDisk_Jump[0] = DeltaT_tip*radius/R_tip;
            }
            else
              ActDisk_Jump[0] = DeltaT_avg;
            
            /*--- Inner point ---*/
            
            if (boundary == ACTDISK_INLET) {
              
              for (iVar = 0; iVar < nPrimVar; iVar++)
                PrimVar_in_ghost[iVar] = PrimVar_out[iVar] - ActDisk_Jump[iVar];
              
              /*--- Check that this is a meaningful state ---*/
              
              P_in_ghost = PrimVar_in_ghost[nDim+1];
              T_in_ghost = PrimVar_in_ghost[0];
              
              FluidModel->SetTDState_PT(P_in_ghost, T_in_ghost);
              Rho_in_ghost = FluidModel->GetDensity();
              sos_in_ghost = FluidModel->GetSoundSpeed();
              
              /*--- Find unit normal to the actuator disk ---*/
              
              geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
              
              Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
              for (iDim = 0; iDim < nDim; iDim++) { UnitNormal[iDim] = -Normal[iDim]/Area; }
              
              /*--- Impose only normal velocity on the disk inflow ---*/
              
              Vel_in_ghost = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                Vel_in_ghost += PrimVar_in[iDim+1]*UnitNormal[iDim];
              }
              
              /*--- COmpute the enthalpy ---*/
              
              H_in_ghost = FluidModel->GetStaticEnergy() + 0.5*Vel_in_ghost*Vel_in_ghost + P_in_ghost/Rho_in_ghost;
              
              PrimVar_in_ghost[0] = T_in_ghost;
              PrimVar_in_ghost[1] = Vel_in_ghost*UnitNormal[0];
              PrimVar_in_ghost[2] = Vel_in_ghost*UnitNormal[1];
              PrimVar_in_ghost[3] = Vel_in_ghost*UnitNormal[2];
              PrimVar_in_ghost[nDim+1]= P_in_ghost;
              PrimVar_in_ghost[nDim+2]= Rho_in_ghost;
              PrimVar_in_ghost[nDim+3]= H_in_ghost;
              PrimVar_in_ghost[nDim+4]= sos_in_ghost;
              
              numerics->SetPrimitive(PrimVar_in, PrimVar_in_ghost);
              
            }
            
            /*--- Outer point ---*/
            
            if (boundary == ACTDISK_OUTLET) {
              
              for (iVar = 0; iVar < nPrimVar; iVar++)
                PrimVar_out_ghost[iVar] = PrimVar_in[iVar] + ActDisk_Jump[iVar];
              
              /*--- Check that this is a meaningful state ---*/
              
              P_out_ghost = PrimVar_out_ghost[nDim+1];
              T_out_ghost = PrimVar_out_ghost[0];
              
              FluidModel->SetTDState_PT(P_out_ghost, T_out_ghost);
              Rho_out_ghost = FluidModel->GetDensity();
              sos_out_ghost = FluidModel->GetSoundSpeed();
              
              /*--- Find unit normal to the actuator disk ---*/
              
              geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
              
              Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
              for (iDim = 0; iDim < nDim; iDim++) { UnitNormal[iDim] = Normal[iDim]/Area; }
              
              /*--- Find unit radial to the actuator disk ---*/
              
              for (iDim = 0; iDim < nDim; iDim++)
                Radial[iDim] = Coord[iDim]-Origin[iDim];
              
              Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Radial[iDim]*Radial[iDim]; Area = sqrt(Area);
              for (iDim = 0; iDim < nDim; iDim++) { UnitRadial[iDim] = Radial[iDim]/Area; }
              
              /*--- Find unit tangent to the actuator disk ---*/
              
              Tangent[0] = UnitNormal[1]*UnitRadial[2] - UnitNormal[2]*UnitRadial[1];
              Tangent[1] = UnitNormal[2]*UnitRadial[0] - UnitNormal[0]*UnitRadial[2];
              Tangent[2] = UnitNormal[0]*UnitRadial[1] - UnitNormal[1]*UnitRadial[0];
              
              Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Tangent[iDim]*Tangent[iDim]; Area = sqrt(Area);
              for (iDim = 0; iDim < nDim; iDim++) { UnitTangent[iDim] = Tangent[iDim]/Area; }
              
              /*--- Normal velocity to the disk and density at the inlet---*/
              
              Vel_Normal_in = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) {
                Vel_Normal_in += PrimVar_in[iDim+1]*UnitNormal[iDim];
              }
              
              Rho_in = PrimVar_in[nDim+2];
              
              /*--- Conservation of mass to evaluate velocity at the outlet ---*/
              
              Vel_out_ghost = (Rho_in*Vel_Normal_in)/Rho_out_ghost;
              
              /*--- Compute the swirl velocity (check the formula) ---*/
              
              if (Omega != 0.0)
                Vel_Swirl_out_ghost = Omega*radius*(1.0-(1.0-sqrt(2.0*ActDisk_Jump[nDim+1]/(Rho_in*pow(Omega*radius, 2.0)))));
              else
                Vel_Swirl_out_ghost = 0.0;
              
              /*--- Compute normal component of the velocity ---*/
              
              Vel_Normal_out_ghost = sqrt(Vel_out_ghost*Vel_out_ghost-Vel_Swirl_out_ghost*Vel_Swirl_out_ghost);
              
              /*--- Compute total entalphy ---*/
              
              H_out_ghost = FluidModel->GetStaticEnergy() + 0.5*(Vel_Normal_out_ghost+Vel_Swirl_out_ghost)*(Vel_Normal_out_ghost+Vel_Swirl_out_ghost) + P_out_ghost/Rho_out_ghost;
              
              PrimVar_out_ghost[0] = T_out_ghost;
              PrimVar_out_ghost[1] = Vel_Normal_out_ghost*UnitNormal[0] + Vel_Swirl_out_ghost*UnitTangent[0];
              PrimVar_out_ghost[2] = Vel_Normal_out_ghost*UnitNormal[1] + Vel_Swirl_out_ghost*UnitTangent[1];
              PrimVar_out_ghost[3] = Vel_Normal_out_ghost*UnitNormal[2] + Vel_Swirl_out_ghost*UnitTangent[2];
              PrimVar_out_ghost[nDim+1]= P_out_ghost;
              PrimVar_out_ghost[nDim+2]= Rho_out_ghost;
              PrimVar_out_ghost[nDim+3]= H_out_ghost;
              PrimVar_out_ghost[nDim+4]= sos_out_ghost;
              
              numerics->SetPrimitive(PrimVar_out, PrimVar_out_ghost);
              
            }
            
            /*--- Set the normal vector ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, Residual);
            if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
          }
          
        }
        
      }
      
    }
    
    /*--- We are using non-blocking communications ---*/
    
#ifdef HAVE_MPI
    
    MPI_Barrier(MPI_COMM_WORLD);
    
#endif
    
    /*--- Free locally allocated memory ---*/
    
    delete [] PrimVar_out;
    delete [] PrimVar_in;
    delete [] MeanPrimVar;
    delete [] PrimVar_out_ghost;
    delete [] PrimVar_in_ghost;
    delete [] ActDisk_Jump;
    delete[] Buffer_Send_V;
    delete[] Buffer_Receive_V;
    
  }
  
}

void CEulerSolver::BC_Dirichlet(CGeometry *geometry, CSolver **solver_container,
                                CConfig *config, unsigned short val_marker) { }

void CEulerSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

void CEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool FlowEq         = (RunTime_EqSystem == RUNTIME_FLOW_SYS);
  bool AdjEq          = (RunTime_EqSystem == RUNTIME_ADJFLOW_SYS);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }
      if ((incompressible || freesurface) && (FlowEq || AdjEq)) Residual[0] = 0.0;
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        if ((incompressible || freesurface) && (FlowEq || AdjEq)) Jacobian_i[0][0] = 0.0;
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
    
  }
  
  else {
    
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Get indices for nodes i & j plus the face normal ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      Normal = geometry->edge[iEdge]->GetNormal();
      
      /*--- Grid velocities stored at nodes i & j ---*/
      
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      
      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/
      
      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      
      /*--- Compute the GCL component of the source term for node i ---*/
      
      U_time_n = node[iPoint]->GetSolution_time_n();
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      if ((incompressible || freesurface) && (FlowEq || AdjEq)) Residual[0] = 0.0;
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = node[jPoint]->GetSolution_time_n();
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      if ((incompressible || freesurface) && (FlowEq || AdjEq)) Residual[0] = 0.0;
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---	Loop over the boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
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
        
        U_time_n = node[iPoint]->GetSolution_time_n();
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
        if ((incompressible || freesurface) && (FlowEq || AdjEq)) Residual[0] = 0.0;
        LinSysRes.AddBlock(iPoint, Residual);
        
      }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }
      if ((incompressible || freesurface) && (FlowEq || AdjEq)) Residual[0] = 0.0;
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        if ((incompressible || freesurface) && (FlowEq || AdjEq)) Jacobian_i[0][0] = 0.0;
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CEulerSolver::SetFlow_Displacement(CGeometry **flow_geometry, CVolumetricMovement *flow_grid_movement,
                                        CConfig *flow_config, CConfig *fea_config, CGeometry **fea_geometry, CSolver ***fea_solution) {
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint;
  su2double *Coord, VarCoord[3];
  
  //#ifndef HAVE_MPI
  unsigned long iPoint_Donor;
  su2double *CoordDonor, *DisplacementDonor;
  
  for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {
    
    if (flow_config->GetMarker_All_FSIinterface(iMarker) != 0) {
      
      for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
        
        iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
        
        iPoint_Donor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();
        
        Coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();
        
        CoordDonor = fea_geometry[MESH_0]->node[iPoint_Donor]->GetCoord();
        
        /*--- The displacements come from the predicted solution ---*/
        DisplacementDonor = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution_Pred();
        
        for (iDim = 0; iDim < nDim; iDim++)
          
          VarCoord[iDim] = (CoordDonor[iDim]+DisplacementDonor[iDim])-Coord[iDim];
        
        flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
    }
  }
  flow_grid_movement->SetVolume_Deformation(flow_geometry[MESH_0], flow_config, true);
  
  
  //#else
  //
  //  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
  //  su2double *Buffer_Send_Coord = new su2double [nDim];
  //  su2double *Buffer_Receive_Coord = new su2double [nDim];
  //  unsigned long jPoint;
  //
  //  /*--- Do the send process, by the moment we are sending each
  //   node individually, this must be changed ---*/
  //  for (iMarker = 0; iMarker < fea_config->GetnMarker_All(); iMarker++) {
  //    if (fea_config->GetMarker_All_KindBC(iMarker) == LOAD_BOUNDARY) {
  //      for (iVertex = 0; iVertex < fea_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
  //        iPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
  //
  //        if (fea_geometry[MESH_0]->node[iPoint]->GetDomain()) {
  //
  //          /*--- Find the associate pair to the original node (index and processor) ---*/
  //          jPoint = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
  //          jProcessor = fea_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //          /*--- We only send the pressure that belong to other boundary ---*/
  //          if (jProcessor != rank) {
  //            for (iDim = 0; iDim < nDim; iDim++)
  //              Buffer_Send_Coord[iDim] = fea_geometry[MESH_0]->node[iPoint]->GetCoord(iDim);
  //
  //            MPI::COMM_WORLD.Bsend(Buffer_Send_Coord, nDim, MPI::DOUBLE, jProcessor, iPoint);
  //          }
  //
  //        }
  //      }
  //    }
  //  }
  //
  //  /*--- Now the loop is over the fea points ---*/
  //  for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {
  //    if ((flow_config->GetMarker_All_KindBC(iMarker) == EULER_WALL) &&
  //        (flow_config->GetMarker_All_Moving(iMarker) == YES)) {
  //      for (iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
  //        iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
  //        if (flow_geometry[MESH_0]->node[iPoint]->GetDomain()) {
  //
  //          /*--- Find the associate pair to the original node ---*/
  //          jPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
  //          jProcessor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
  //
  //          /*--- We only receive the information that belong to other boundary ---*/
  //          if (jProcessor != rank)
  //            MPI::COMM_WORLD.Recv(Buffer_Receive_Coord, nDim, MPI::DOUBLE, jProcessor, jPoint);
  //          else {
  //            for (iDim = 0; iDim < nDim; iDim++)
  //              Buffer_Send_Coord[iDim] = fea_geometry[MESH_0]->node[jPoint]->GetCoord(iDim);
  //          }
  //
  //          /*--- Store the solution for both points ---*/
  //          Coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();
  //
  //          for (iDim = 0; iDim < nDim; iDim++)
  //            VarCoord[iDim] = Buffer_Send_Coord[iDim]-Coord[iDim];
  //
  //          flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
  //
  //
  //        }
  //      }
  //    }
  //  }
  //  delete[] Buffer_Send_Coord;
  //  delete[] Buffer_Receive_Coord;
  //
  //#endif
  //
}

void CEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter) {
  
  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine, dull_val;
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  string UnstExt, text_line;
  ifstream restart_file;
  
  string restart_filename = config->GetSolution_FlowFileName();
  
  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Modify file name for an unsteady restart ---*/
  
  if (dual_time)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);
  
  /*--- Open the restart file, and throw an error if this fails. ---*/
  
  restart_file.open(restart_filename.data(), ios::in);
  if (restart_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  
  long *Global2Local = NULL;
  Global2Local = new long[geometry[MESH_0]->GetGlobal_nPointDomain()];
  /*--- First, set all indices to a negative value by default ---*/
  for (iPoint = 0; iPoint < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint++) {
    Global2Local[iPoint] = -1;
  }
  
  /*--- Now fill array with the transform values only for local points ---*/
  
  for (iPoint = 0; iPoint < geometry[MESH_0]->GetnPointDomain(); iPoint++) {
    Global2Local[geometry[MESH_0]->node[iPoint]->GetGlobalIndex()] = iPoint;
  }
  
  /*--- Read all lines in the restart file ---*/
  
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  
  /*--- The first line is the header ---*/
  
  getline (restart_file, text_line);
  
  while (getline (restart_file, text_line)) {
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1, as
     initialized above. Otherwise, the local index for this node on the
     current processor will be returned and used to instantiate the vars. ---*/
    
    iPoint_Local = Global2Local[iPoint_Global];
    if (iPoint_Local >= 0) {
      
      if (compressible) {
        if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
      }
      if (incompressible) {
        if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1] >> Solution[0] >> Solution[1] >> Solution[2];
        if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
      }
      if (freesurface) {
        if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
      }
      
      node[iPoint_Local]->SetSolution(Solution);
      
      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/
      
      if (grid_movement) {
        
        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/
        
        if (turb_model == SA || turb_model == SA_NEG) {
          point_line >> dull_val;
        } else if (turb_model == SST) {
          point_line >> dull_val >> dull_val;
        }
        
        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        
        su2double GridVel[3] = {0.0,0.0,0.0};
        if (nDim == 2) point_line >> GridVel[0] >> GridVel[1];
        else point_line >> GridVel[0] >> GridVel[1] >> GridVel[2];
        
        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[MESH_0]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
        
      }
      
    }
    iPoint_Global++;
  }
  
  /*--- Close the restart file ---*/
  
  restart_file.close();
  
  /*--- Free memory needed for the transformation ---*/
  
  delete [] Global2Local;
  
  /*--- MPI solution ---*/
  
  solver[MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  
  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
  }
  
  /*--- Update the geometry for flows on dynamic meshes ---*/
  
  if (grid_movement) {
    
    /*--- Communicate the new coordinates and grid velocities at the halos ---*/
    
    geometry[MESH_0]->Set_MPI_Coord(config);
    geometry[MESH_0]->Set_MPI_GridVel(config);
    
    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/
    
    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);
    
    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/
    
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine], config);
    }
  }
  
  delete [] Coord;
  
}

void CEulerSolver::SetFreeSurface_Distance(CGeometry *geometry, CConfig *config) {
  su2double *coord = NULL, dist2, *iCoord = NULL, *jCoord = NULL, LevelSet_i, LevelSet_j,
  **Coord_LevelSet = NULL, *xCoord = NULL, *yCoord = NULL, *zCoord = NULL, auxCoordx, auxCoordy,
  auxCoordz, FreeSurface, volume, LevelSetDiff_Squared, LevelSetDiff, dist, Min_dist;
  unsigned short iDim;
  unsigned long iPoint, jPoint, iVertex, jVertex, nVertex_LevelSet, iEdge;
  ifstream index_file;
  ofstream LevelSet_file;
  string text_line;
  int rank = MASTER_NODE;
  char cstr[200], buffer[50];
  
  unsigned short nDim = geometry->GetnDim();
  unsigned long iExtIter = config->GetExtIter();
  
#ifndef HAVE_MPI
  
  /*--- Identification of the 0 level set points and coordinates ---*/
  nVertex_LevelSet = 0;
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0); LevelSet_i = node[iPoint]->GetSolution(nDim+1);
    jPoint = geometry->edge[iEdge]->GetNode(1); LevelSet_j = node[jPoint]->GetSolution(nDim+1);
    if (LevelSet_i*LevelSet_j < 0.0) nVertex_LevelSet ++;
  }
  
  /*--- Allocate vector of boundary coordinates ---*/
  Coord_LevelSet = new su2double* [nVertex_LevelSet];
  for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
    Coord_LevelSet[iVertex] = new su2double [nDim];
  
  /*--- Get coordinates of the points of the surface ---*/
  nVertex_LevelSet = 0;
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0); LevelSet_i = node[iPoint]->GetSolution(nDim+1); iCoord = geometry->node[iPoint]->GetCoord();
    jPoint = geometry->edge[iEdge]->GetNode(1); LevelSet_j = node[jPoint]->GetSolution(nDim+1); jCoord = geometry->node[jPoint]->GetCoord();
    if (LevelSet_i*LevelSet_j < 0.0) {
      for (iDim = 0; iDim < nDim; iDim++)
        Coord_LevelSet[nVertex_LevelSet][iDim] = iCoord[iDim]-LevelSet_i*(jCoord[iDim]-iCoord[iDim])/(LevelSet_j-LevelSet_i);
      nVertex_LevelSet++;
    }
  }
  
#else
  
  int nProcessor, iProcessor;
  unsigned long *Buffer_Send_nVertex = NULL, *Buffer_Receive_nVertex = NULL,
  nLocalVertex_LevelSet = 0, nGlobalVertex_LevelSet = 0, MaxLocalVertex_LevelSet = 0, nBuffer;
  su2double *Buffer_Send_Coord = NULL, *Buffer_Receive_Coord = NULL;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  Buffer_Send_nVertex = new unsigned long [1];
  Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  nLocalVertex_LevelSet = 0;
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0); LevelSet_i = node[iPoint]->GetSolution(nDim+1);
    jPoint = geometry->edge[iEdge]->GetNode(1); LevelSet_j = node[jPoint]->GetSolution(nDim+1);
    if (LevelSet_i*LevelSet_j < 0.0) nLocalVertex_LevelSet ++;
  }
  
  Buffer_Send_nVertex[0] = nLocalVertex_LevelSet;
  
  SU2_MPI::Allreduce(&nLocalVertex_LevelSet, &nGlobalVertex_LevelSet, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalVertex_LevelSet, &MaxLocalVertex_LevelSet, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  nBuffer = MaxLocalVertex_LevelSet*nDim;
  Buffer_Send_Coord = new su2double [nBuffer];
  Buffer_Receive_Coord = new su2double [nProcessor*nBuffer];
  
  for (iVertex = 0; iVertex < MaxLocalVertex_LevelSet; iVertex++)
    for (iDim = 0; iDim < nDim; iDim++)
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
  
  nLocalVertex_LevelSet = 0;
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0); LevelSet_i = node[iPoint]->GetSolution(nDim+1); iCoord = geometry->node[iPoint]->GetCoord();
    jPoint = geometry->edge[iEdge]->GetNode(1); LevelSet_j = node[jPoint]->GetSolution(nDim+1); jCoord = geometry->node[jPoint]->GetCoord();
    
    if (LevelSet_i*LevelSet_j < 0.0) {
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[nLocalVertex_LevelSet*nDim+iDim] = iCoord[iDim]-LevelSet_i*(jCoord[iDim]-iCoord[iDim])/(LevelSet_j-LevelSet_i);
      nLocalVertex_LevelSet++;
    }
  }
  
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
  
  /*--- Identification of the 0 level set points and coordinates ---*/
  nVertex_LevelSet = 0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
    nVertex_LevelSet += Buffer_Receive_nVertex[iProcessor];
  
  /*--- Allocate vector of boundary coordinates ---*/
  Coord_LevelSet = new su2double* [nVertex_LevelSet];
  for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
    Coord_LevelSet[iVertex] = new su2double [nDim];
  
  /*--- Set the value of the coordinates at the level set zero ---*/
  nVertex_LevelSet = 0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
      for (iDim = 0; iDim < nDim; iDim++)
        Coord_LevelSet[nVertex_LevelSet][iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_LevelSet+iVertex)*nDim+iDim];
      nVertex_LevelSet++;
    }
  
  delete [] Buffer_Send_Coord;
  delete [] Buffer_Receive_Coord;
  delete [] Buffer_Send_nVertex;
  delete [] Buffer_Receive_nVertex;
  
#endif
  
  /*--- Get coordinates of the points and compute distances to the surface ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Compute the min distance ---*/
    Min_dist = 1E20;
    for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
      
      dist2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist2 += (coord[iDim] - Coord_LevelSet[iVertex][iDim])*(coord[iDim]-Coord_LevelSet[iVertex][iDim]);
      dist = sqrt(dist2);
      if (dist < Min_dist) { Min_dist = dist; }
      
    }
    
    /*--- Compute the sign using the current solution ---*/
    su2double NumberSign = 1.0;
    if (node[iPoint]->GetSolution(0) != 0.0)
      NumberSign = node[iPoint]->GetSolution(nDim+1)/fabs(node[iPoint]->GetSolution(nDim+1));
    
    /*--- Store the value of the Level Set and the Distance (primitive variables) ---*/
    node[iPoint]->SetPrimitive(nDim+5, node[iPoint]->GetSolution(nDim+1));
    node[iPoint]->SetPrimitive(nDim+6, Min_dist*NumberSign);
    
  }
  
  if (config->GetIntIter() == 0) {
    
    /*--- Order the arrays (x Coordinate, y Coordinate, z Coordiante) ---*/
    for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
      for (jVertex = 0; jVertex < nVertex_LevelSet - 1 - iVertex; jVertex++) {
        if (Coord_LevelSet[jVertex][0] > Coord_LevelSet[jVertex+1][0]) {
          auxCoordx = Coord_LevelSet[jVertex][0]; Coord_LevelSet[jVertex][0] = Coord_LevelSet[jVertex+1][0]; Coord_LevelSet[jVertex+1][0] = auxCoordx;
          auxCoordy = Coord_LevelSet[jVertex][1]; Coord_LevelSet[jVertex][1] = Coord_LevelSet[jVertex+1][1]; Coord_LevelSet[jVertex+1][1] = auxCoordy;
          if (nDim == 3) { auxCoordz = Coord_LevelSet[jVertex][2]; Coord_LevelSet[jVertex][2] = Coord_LevelSet[jVertex+1][2]; Coord_LevelSet[jVertex+1][2] = auxCoordz; }
        }
      }
    }
    
    /*--- Get coordinates of the points and compute distances to the surface ---*/
    FreeSurface = 0.0;
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      coord = geometry->node[iPoint]->GetCoord();
      volume = geometry->node[iPoint]->GetVolume();
      LevelSetDiff = (node[iPoint]->GetSolution(nDim+1) - coord[nDim-1]);
      LevelSetDiff_Squared = LevelSetDiff*LevelSetDiff;
      FreeSurface += 0.5*LevelSetDiff_Squared*volume;
      
      node[iPoint]->SetDiffLevelSet(LevelSetDiff);
      
    }
    
    if ((rank == MASTER_NODE) && (iExtIter % config->GetWrt_Sol_Freq_DualTime() == 0)) {
      
      /*--- Write the Level Set distribution, the target level set---*/
      LevelSet_file.precision(15);
      
      /*--- Write file name with extension ---*/
      strcpy (cstr, "LevelSet");
      if (config->GetUnsteady_Simulation()) {
        if ((SU2_TYPE::Int(iExtIter) >= 0) && (SU2_TYPE::Int(iExtIter) < 10)) SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
        if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d.dat", SU2_TYPE::Int(iExtIter));
        if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d.dat", SU2_TYPE::Int(iExtIter));
        if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat", SU2_TYPE::Int(iExtIter));
        if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
      }
      else {
        SPRINTF (buffer, ".dat");
      }
      
      strcat(cstr, buffer);
      
      LevelSet_file.open(cstr, ios::out);
      LevelSet_file << "TITLE = \"SU2 Free surface simulation\"" << endl;
      if (nDim == 2) LevelSet_file << "VARIABLES = \"x coord\",\"y coord\"" << endl;
      if (nDim == 3) LevelSet_file << "VARIABLES = \"x coord\",\"y coord\",\"z coord\"" << endl;
      LevelSet_file << "ZONE T= \"Free Surface\"" << endl;
      
      for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
        if (nDim == 2) LevelSet_file << scientific << Coord_LevelSet[iVertex][0] << ", " << Coord_LevelSet[iVertex][1] << endl;
        if (nDim == 3) LevelSet_file << scientific << Coord_LevelSet[iVertex][0] << ", " << Coord_LevelSet[iVertex][1] << ", " << Coord_LevelSet[iVertex][2] << endl;
      }
      LevelSet_file.close();
      
    }
    
    /*--- Store the value of the free surface coefficient ---*/
    SetTotal_CFreeSurface(FreeSurface);
    
    delete [] xCoord;
    delete [] yCoord;
    if (nDim == 3) delete [] zCoord;
    
  }
  
  /*--- Deallocate vector of boundary coordinates ---*/
  for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
    delete Coord_LevelSet[iVertex];
  delete [] Coord_LevelSet;
  
}

CNSSolver::CNSSolver(void) : CEulerSolver() {
  
  /*--- Basic array initialization ---*/
  
  CDrag_Visc = NULL; CLift_Visc = NULL; CSideForce_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;
  
  /*--- Surface based array initialization ---*/
  
  Surface_CLift_Visc = NULL; Surface_CDrag_Visc = NULL; Surface_CSideForce_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  
  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Visc = NULL; CT_Visc = NULL; CQ_Visc = NULL;
  
}

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CEulerSolver() {
  
  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double Density, Velocity2, Pressure, Temperature, dull_val, StaticEnergy;
  int Unst_RestartIter;
  ifstream restart_file;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool adjoint = config->GetAdjoint();
  string filename = config->GetSolution_FlowFileName();
  
  unsigned short direct_diff = config->GetDirectDiff();
  unsigned short nMarkerTurboPerf = config->Get_nMarkerTurboPerf();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Array initialization ---*/
  
  CDrag_Visc = NULL; CLift_Visc = NULL; CSideForce_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  
  Surface_CLift_Visc = NULL; Surface_CDrag_Visc = NULL; Surface_CSideForce_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  
  CMerit_Visc = NULL;      CT_Visc = NULL;      CQ_Visc = NULL;
  MaxHeatFlux_Visc = NULL; ForceViscous = NULL; MomentViscous = NULL;
  CSkinFriction = NULL;    Cauchy_Serie = NULL; Heat_Visc = NULL;
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp)
   Incompressible flow, primitive variables (P, vx, vy, vz, rho, beta, lamMu, EddyMu),
   FreeSurface Incompressible flow, primitive variables (P, vx, vy, vz, rho, beta, lamMu, EddyMu, LevelSet, Dist),
   ---*/
  
  nDim = geometry->GetnDim();
  
  if (incompressible) { nVar = nDim+1; nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  if (freesurface)    { nVar = nDim+2; nPrimVar = nDim+7; nPrimVarGrad = nDim+6; }
  if (compressible)   { nVar = nDim+2;
    nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
    nSecondaryVar = 8; nSecondaryVarGrad = 2;
  }
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(geometry, config, iMesh);
  
  /*--- Allocate the node variables ---*/
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  
  Primitive   = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new su2double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the Secondary solution ---*/
  
  if (compressible){
    Secondary   = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
    Secondary_i = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
    Secondary_j = new su2double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;
  }
  
  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  
  if (roe_turkel) {
    LowMach_Precontioner = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new su2double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  }
  
  else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  CharacPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Inviscid force definition and coefficient in all the markers ---*/
  
  CPressure = new su2double* [nMarker];
  CPressureTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];
    CPressureTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Heat flux in all the markers ---*/
  
  HeatFlux = new su2double* [nMarker];
  HeatFluxTarget = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatFlux[iMarker] = new su2double [geometry->nVertex[iMarker]];
    HeatFluxTarget[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      HeatFlux[iMarker][iVertex] = 0.0;
      HeatFluxTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Y plus in all the markers ---*/
  
  YPlus = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    YPlus[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      YPlus[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Skin friction in all the markers ---*/
  
  CSkinFriction = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSkinFriction[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CSkinFriction[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Non dimensional coefficients ---*/
  
  ForceInviscid  = new su2double[3];
  MomentInviscid = new su2double[3];
  CDrag_Inv      = new su2double[nMarker];
  CLift_Inv      = new su2double[nMarker];
  CSideForce_Inv = new su2double[nMarker];
  CMx_Inv        = new su2double[nMarker];
  CMy_Inv        = new su2double[nMarker];
  CMz_Inv        = new su2double[nMarker];
  CEff_Inv       = new su2double[nMarker];
  CFx_Inv        = new su2double[nMarker];
  CFy_Inv        = new su2double[nMarker];
  CFz_Inv        = new su2double[nMarker];
  
  ForceViscous     = new su2double[3];
  MomentViscous    = new su2double[3];
  CDrag_Visc       = new su2double[nMarker];
  CLift_Visc       = new su2double[nMarker];
  CSideForce_Visc  = new su2double[nMarker];
  CMx_Visc         = new su2double[nMarker];
  CMy_Visc         = new su2double[nMarker];
  CMz_Visc         = new su2double[nMarker];
  CEff_Visc        = new su2double[nMarker];
  CFx_Visc         = new su2double[nMarker];
  CFy_Visc         = new su2double[nMarker];
  CFz_Visc         = new su2double[nMarker];
  
  Surface_CLift_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce_Inv = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CLift          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];
  
  Surface_CLift_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSideForce_Visc = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Rotational coefficients ---*/
  
  CMerit_Inv = new su2double[nMarker];
  CT_Inv     = new su2double[nMarker];
  CQ_Inv     = new su2double[nMarker];
  
  CMerit_Visc = new su2double[nMarker];
  CT_Visc     = new su2double[nMarker];
  CQ_Visc     = new su2double[nMarker];
  
  /*--- Heat based coefficients ---*/
  
  Heat_Visc        = new su2double[nMarker];
  MaxHeatFlux_Visc = new su2double[nMarker];
  
  /*--- Supersonic coefficients ---*/
  
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];
  
  /*--- Engine simulation ---*/
  
  Inflow_MassFlow     = new su2double[nMarker];
  Inflow_Pressure     = new su2double[nMarker];
  Inflow_Mach         = new su2double[nMarker];
  Inflow_Area         = new su2double[nMarker];
  
  Exhaust_MassFlow    = new su2double[nMarker];
  Exhaust_Pressure    = new su2double[nMarker];
  Exhaust_Temperature = new su2double[nMarker];
  Exhaust_Area        = new su2double[nMarker];
  
  Bleed_MassFlow      = new su2double[nMarker];
  Bleed_Pressure      = new su2double[nMarker];
  Bleed_Temperature   = new su2double[nMarker];
  Bleed_Area          = new su2double[nMarker];
  
  /*--- Init total coefficients ---*/
  
  Total_CDrag   = 0.0;	Total_CLift        = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy          = 0.0;  Total_CMz          = 0.0;
  Total_CEff    = 0.0;	Total_CEquivArea   = 0.0;  Total_CNearFieldOF = 0.0;
  Total_CFx     = 0.0;	Total_CFy          = 0.0;  Total_CFz          = 0.0;
  Total_CT      = 0.0;	Total_CQ           = 0.0;  Total_CMerit       = 0.0;
  Total_MaxHeat = 0.0;  Total_Heat         = 0.0;
  Total_CpDiff  = 0.0;  Total_HeatFluxDiff = 0.0;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Mach_Inf        = config->GetMach();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  Tke_Inf         = config->GetTke_FreeStreamND();
  
  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  
  switch(direct_diff){
    case NO_DERIVATIVE:
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
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }
  
  /*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;
    
    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;
    
    Bleed_MassFlow[iMarker]      = 0.0;
    Bleed_Temperature[iMarker]   = Temperature_Inf;
    Bleed_Pressure[iMarker]      = Pressure_Inf;
    Bleed_Area[iMarker]          = 0.0;
  }
  
  /*--- Initializate quantities for the mixing process*/
  
  AveragedVelocity = new su2double* [nMarker];
  AveragedNormal = new su2double* [nMarker];
  AveragedGridVel = new su2double* [nMarker];
  AveragedFlux = new su2double* [nMarker];
  TotalFlux = new su2double* [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedVelocity[iMarker] = new su2double [nDim];
    AveragedNormal[iMarker] = new su2double [nDim];
    AveragedGridVel[iMarker] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      AveragedVelocity[iMarker][iDim] = 0.0;
      AveragedNormal[iMarker][iDim] = 0.0;
      AveragedGridVel[iMarker][iDim] = 0.0;
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    AveragedFlux[iMarker] = new su2double [nVar];
    TotalFlux[iMarker] = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      AveragedFlux[iMarker][iVar] = 0.0;
      TotalFlux[iMarker][iVar] = 0.0;
    }
  }
  
  AveragedNormalVelocity = new su2double[nMarker];
  AveragedTangVelocity = new su2double[nMarker];
  ExtAveragedNormalVelocity = new su2double[nMarker];
  ExtAveragedTangVelocity = new su2double[nMarker];
  MassFlow= new su2double[nMarker];
  FlowAngle= new su2double[nMarker];
  AveragedEnthalpy  = new su2double[nMarker];
  AveragedPressure  = new su2double[nMarker];
  AveragedTotPressure  = new su2double[nMarker];
  AveragedTotTemperature  = new su2double[nMarker];
  ExtAveragedTotPressure  = new su2double[nMarker];
  ExtAveragedTotTemperature  = new su2double[nMarker];
  ExtAveragedPressure  = new su2double[nMarker];
  AveragedDensity   = new su2double[nMarker];
  ExtAveragedDensity   = new su2double[nMarker];
  AveragedSoundSpeed= new su2double[nMarker];
  AveragedEntropy   = new su2double[nMarker];
  AveragedTangGridVelocity = new su2double[nMarker];
  AveragedMach = new su2double[nMarker];
  AveragedNormalMach = new su2double[nMarker];
  AveragedTangMach = new su2double[nMarker];
  
  
  /*--- Initializate quantities for turboperformace ---*/
  
  TotalStaticEfficiency = new su2double[nMarkerTurboPerf];
  TotalTotalEfficiency = new su2double[nMarkerTurboPerf];
  KineticEnergyLoss= new su2double[nMarkerTurboPerf];
  TotalPressureLoss= new su2double[nMarkerTurboPerf];
  MassFlowIn= new su2double[nMarkerTurboPerf];
  MassFlowOut= new su2double[nMarkerTurboPerf];
  FlowAngleIn= new su2double[nMarkerTurboPerf];
  FlowAngleOut= new su2double[nMarkerTurboPerf];
  EulerianWork= new su2double[nMarkerTurboPerf];
  TotalEnthalpyIn= new su2double[nMarkerTurboPerf];
  PressureRatio= new su2double[nMarkerTurboPerf];
  PressureOut= new su2double[nMarkerTurboPerf];
  EnthalpyOut= new su2double[nMarkerTurboPerf];
  MachIn= new su2double[nMarkerTurboPerf];
  MachOut= new su2double[nMarkerTurboPerf];
  NormalMachIn= new su2double[nMarkerTurboPerf];
  NormalMachOut= new su2double[nMarkerTurboPerf];
  VelocityOutIs= new su2double[nMarkerTurboPerf];
  
  for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
    TotalStaticEfficiency[iMarker]= 0.0;
    TotalTotalEfficiency[iMarker]= 0.0;
    KineticEnergyLoss[iMarker]= 0.0;
    TotalPressureLoss[iMarker]= 0.0;
    MassFlowIn[iMarker]= 0.0;
    MassFlowOut[iMarker]= 0.0;
    FlowAngleIn[iMarker]= 0.0;
    FlowAngleOut[iMarker]= 0.0;
    EulerianWork[iMarker]= 0.0;
    TotalEnthalpyIn[iMarker]= 0.0;
    PressureRatio[iMarker]= 0.0;
    PressureOut[iMarker]= 0.0;
    EnthalpyOut[iMarker]= 0.0;
    MachIn[iMarker]= 0.0;
    MachOut[iMarker]= 0.0;
    NormalMachIn[iMarker]= 0.0;
    NormalMachOut[iMarker]= 0.0;
    VelocityOutIs[iMarker]= 0.0;
  }
  
  
  /*--- Initialize the cauchy critera array for fixed CL mode ---*/
  
  if (config->GetFixed_CL_Mode())
    
    Cauchy_Serie = new su2double [config->GetCauchy_Elems()+1];
  
  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  
  if (!restart || (iMesh != MESH_0)) {
    
    /*--- Restart the solution from the free-stream state ---*/
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
    
  }
  
  else {
    
    /*--- Modify file name for an unsteady restart ---*/
    
    if (dual_time) {
      
      if (adjoint) {
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter()) - 1;
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
      
    }
    
    if (nZone >1)
      filename= config->GetRestart_FlowFileName(filename, iZone);
    
    /*--- Open the restart file, throw an error if this fails. ---*/
    
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    
    long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    
    /*--- First, set all indices to a negative value by default ---*/
    
    for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
      Global2Local[iPoint] = -1;
    
    /*--- Now fill array with the transform values only for local points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    
    /*--- Read all lines in the restart file ---*/
    
    long iPoint_Local;
    unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
    
    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    while (getline (restart_file, text_line)) {
      istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
      
      if (iPoint_Global >= geometry->GetGlobal_nPointDomain()) { sbuf_NotMatching = 1; break; }
      
      iPoint_Local = Global2Local[iPoint_Global];
      
      /*--- Load the solution for this node. Note that the first entry
       on the restart file line is the global index, followed by the
       node coordinates, and then the conservative variables. ---*/
      
      if (iPoint_Local >= 0) {
        if (compressible) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        }
        if (incompressible) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        }
        if (freesurface) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        }
        node[iPoint_Local] = new CNSVariable(Solution, nDim, nVar, config);
        iPoint_Global_Local++;
      }
      iPoint_Global++;
    }
    
    /*--- Detect a wrong solution file ---*/
    
    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
    
#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
        cout << "It could be empty lines at the end of the file." << endl << endl;
      }
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Solution, nDim, nVar, config);
    
    /*--- Close the restart file ---*/
    
    restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    
    delete [] Global2Local;
    
  }
  
  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  
  if (compressible) {
    
    counter_local = 0;
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      Density = node[iPoint]->GetSolution(0);
      
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);
      
      StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;
      
      FluidModel->SetTDState_rhoe(Density, StaticEnergy);
      Pressure= FluidModel->GetPressure();
      Temperature= FluidModel->GetTemperature();
      
      /*--- Use the values at the infinity ---*/
      
      if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
        Solution[0] = Density_Inf;
        for (iDim = 0; iDim < nDim; iDim++)
          Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
        Solution[nDim+1] = Energy_Inf*Density_Inf;
        node[iPoint]->SetSolution(Solution);
        node[iPoint]->SetSolution_Old(Solution);
        counter_local++;
      }
      
    }
    
    /*--- Warning message about non-physical points ---*/
    
    if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
      SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
      counter_global = counter_local;
#endif
      if ((rank == MASTER_NODE) && (counter_global != 0))
        cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
    }
    
  }
  
  /*--- For incompressible solver set the initial values for the density and viscosity,
   unless a freesurface problem, this must be constant during the computation ---*/
  
  if (incompressible || freesurface) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint]->SetDensityInc(Density_Inf);
      node[iPoint]->SetLaminarViscosityInc(Viscosity_Inf);
    }
  }
  
  /*--- Define solver parameters needed for execution of destructor ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) space_centered = true;
  else space_centered = false;
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;
  
  /*--- Perform the MPI communication of the solution ---*/
  
//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);
  
}

CNSSolver::~CNSSolver(void) {
  unsigned short iMarker;
  
  if (CDrag_Visc != NULL)       delete [] CDrag_Visc;
  if (CLift_Visc != NULL)       delete [] CLift_Visc;
  if (CSideForce_Visc != NULL)  delete [] CSideForce_Visc;
  if (CMx_Visc != NULL)         delete [] CMx_Visc;
  if (CMy_Visc != NULL)         delete [] CMy_Visc;
  if (CMz_Visc != NULL)         delete [] CMz_Visc;
  if (CFx_Visc != NULL)         delete [] CFx_Visc;
  if (CFy_Visc != NULL)         delete [] CFy_Visc;
  if (CFz_Visc != NULL)         delete [] CFz_Visc;
  if (CEff_Visc != NULL)        delete [] CEff_Visc;
  if (CMerit_Visc != NULL)      delete [] CMerit_Visc;
  if (CT_Visc != NULL)          delete [] CT_Visc;
  if (CQ_Visc != NULL)          delete [] CQ_Visc;
  if (Heat_Visc != NULL)        delete [] Heat_Visc;
  if (MaxHeatFlux_Visc != NULL) delete [] MaxHeatFlux_Visc;
  if (ForceViscous != NULL)     delete [] ForceViscous;
  if (MomentViscous != NULL)    delete [] MomentViscous;
  
  
  if (Surface_CLift_Visc != NULL)      delete [] Surface_CLift_Visc;
  if (Surface_CDrag_Visc != NULL)      delete [] Surface_CDrag_Visc;
  if (Surface_CSideForce_Visc != NULL) delete [] Surface_CSideForce_Visc;
  if (Surface_CEff_Visc != NULL)       delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)        delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)        delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)        delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)        delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)        delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)        delete [] Surface_CMz_Visc;
  
  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;
  
  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }
  
}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double StrainMag = 0.0, Omega = 0.0, *Vorticity;
  
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned long ExtIter     = config->GetExtIter();
  bool adjoint              = config->GetAdjoint();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center               = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst           = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow         = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_turb         = ((config->GetSpatialOrder_Turb() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_adjflow      = ((config->GetSpatialOrder_AdjFlow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_visc         = config->GetViscous_Limiter_Flow();
  bool freesurface          = (config->GetKind_Regime() == FREESURFACE);
  bool fixed_cl             = config->GetFixed_CL_Mode();
  bool engine               = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineBleed() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk        = ((config->GetnMarker_ActDisk_Inlet() != 0) || (config->GetnMarker_ActDisk_Outlet() != 0));
  
  /*--- Compute the engine properties ---*/
  
  if (engine) { GetEngine_Properties(geometry, config, iMesh, Output); }
  
  /*--- Compute the actuator disk properties ---*/
  
  if (actuator_disk) { GetActuatorDisk_Properties(geometry, config, iMesh, Output); }
  
  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/
  
  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }
  
  /*--- Compute distance function to zero level set ---*/
  
  if (freesurface) { SetFreeSurface_Distance(geometry, config); }
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Artificial dissipation ---*/
  
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Compute gradient of the primitive variables ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
    //	  if (compressible && !ideal_gas) SetSecondary_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
    //	  if (compressible && !ideal_gas) SetSecondary_Gradient_LS(geometry, config);
  }
  
  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/
  
  if ((iMesh == MESH_0) && (limiter_flow || limiter_turb || limiter_adjflow || limiter_visc) && !Output) { SetPrimitive_Limiter(geometry, config);
    //  if (compressible && !ideal_gas) SetSecondary_Limiter(geometry, config);
  }
  
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  
  StrainMag_Max = 0.0, Omega_Max = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    solver_container[FLOW_SOL]->node[iPoint]->SetVorticity(limiter_visc);
    solver_container[FLOW_SOL]->node[iPoint]->SetStrainMag(limiter_visc);
    
    StrainMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    Vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);
    
    StrainMag_Max = max(StrainMag_Max, StrainMag);
    Omega_Max = max(Omega_Max, Omega);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !config->GetDiscrete_Adjoint()) Jacobian.SetValZero();
  
  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    su2double MyOmega_Max = Omega_Max; Omega_Max = 0.0;
    su2double MyStrainMag_Max = StrainMag_Max; StrainMag_Max = 0.0;
    
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    if (iMesh == MESH_0) {
      config->SetNonphysical_Points(ErrorCounter);
      solver_container[FLOW_SOL]->SetStrainMag_Max(StrainMag_Max);
      solver_container[FLOW_SOL]->SetOmega_Max(Omega_Max);
    }
    
  }
  
}

unsigned long CNSSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0;
  unsigned short turb_model = config->GetKind_Turb_Model();
  bool RightSol = true;
  
  bool tkeNeeded            = (turb_model == SST);
  bool compressible         = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible       = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface          = (config->GetKind_Regime() == FREESURFACE);
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Retrieve the value of the kinetic energy (if need it) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
    }
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Incompressible flow, primitive variables nDim+3, (P, vx, vy, vz, rho, beta),
     FreeSurface Incompressible flow, primitive variables nDim+4, (P, vx, vy, vz, rho, beta, dist),
     Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
    
    if (compressible) {
      RightSol = node[iPoint]->SetPrimVar_Compressible(eddy_visc, turb_ke, FluidModel);
      node[iPoint]->SetSecondaryVar_Compressible(FluidModel);
    }
    
    if (incompressible){
      RightSol = node[iPoint]->SetPrimVar_Incompressible(Density_Inf, Viscosity_Inf, eddy_visc, turb_ke, config);
    }
    
    if (freesurface){
      RightSol = node[iPoint]->SetPrimVar_FreeSurface(eddy_visc, turb_ke, config);
    }
    
    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
}

void CNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  
  su2double Mean_BetaInc2, *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time, Local_Delta_Time_Visc, Mean_DensityInc,
  Global_Delta_Time = 1E6, Mean_LaminarVisc = 0.0, Mean_EddyVisc = 0.0, Mean_Density = 0.0, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  unsigned short iDim, iMarker;
  su2double ProjVel, ProjVel_i, ProjVel_j;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetMax_Lambda_Inv(0.0);
    node[iPoint]->SetMax_Lambda_Visc(0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    if (compressible) {
      Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
      Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    }
    if (incompressible || freesurface) {
      Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
      Mean_BetaInc2 = 0.5 * (node[iPoint]->GetBetaInc2() + node[jPoint]->GetBetaInc2());
      Mean_DensityInc = 0.5 * (node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
      Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
    }
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j) ;
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
    /*--- Viscous contribution ---*/
    
    if (compressible) {
      Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
      Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosity() + node[jPoint]->GetEddyViscosity());
      Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));
    }
    if (incompressible || freesurface) {
      Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosityInc() + node[jPoint]->GetLaminarViscosityInc());
      Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosityInc() + node[jPoint]->GetEddyViscosityInc());
      Mean_Density     = 0.5*(node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
    }
    
    Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
    //TODO (REAL_GAS) removing Gamma it cannot work with FLUIDPROP
    Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
    Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      if (compressible) {
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      }
      if (incompressible || freesurface) {
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_BetaInc2 = node[iPoint]->GetBetaInc2();
        Mean_DensityInc = node[iPoint]->GetDensityInc();
        Mean_SoundSpeed = sqrt(Mean_ProjVel*Mean_ProjVel + (Mean_BetaInc2/Mean_DensityInc)*Area*Area);
      }
      
      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
      /*--- Viscous contribution ---*/
      
      if (compressible) {
        Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
        Mean_EddyVisc    = node[iPoint]->GetEddyViscosity();
        Mean_Density     = node[iPoint]->GetSolution(0);
      }
      if (incompressible || freesurface) {
        Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosityInc() + node[jPoint]->GetLaminarViscosityInc());
        Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosityInc() + node[jPoint]->GetEddyViscosityInc());
        Mean_Density     = 0.5*(node[iPoint]->GetDensityInc() + node[jPoint]->GetDensityInc());
      }
      
      Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
      Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
      Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
      Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;
    
    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive and secondary variables ---*/
    
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    numerics->SetSecondary(node[iPoint]->GetSecondary(), node[jPoint]->GetSecondary());
    
    /*--- Gradient and limiters ---*/
    
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
    numerics->SetPrimVarLimiter(node[iPoint]->GetLimiter_Primitive(), node[jPoint]->GetLimiter_Primitive());
    
    /*--- Turbulent kinetic energy ---*/
    
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                     solver_container[TURB_SOL]->node[jPoint]->GetSolution(0));
    
    /*--- Compute and update residual ---*/
    
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.SubtractBlock(iPoint, Res_Visc);
    LinSysRes.AddBlock(jPoint, Res_Visc);
    
    /*--- Implicit part ---*/
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    }
    
  }
  
}

void CNSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  su2double Viscosity = 0.0, div_vel, *Normal, MomentDist[3] = {0.0, 0.0, 0.0}, WallDist[3] = {0.0, 0.0, 0.0},
  *Coord, *Coord_Normal, Area, WallShearStress, TauNormal, factor, RefTemp, RefVel2,
  RefDensity, GradTemperature, Density = 0.0, WallDistMod, FrictionVel,
  Mach2Vel, Mach_Motion, UnitNormal[3] = {0.0, 0.0, 0.0}, TauElem[3] = {0.0, 0.0, 0.0}, TauTangent[3] = {0.0, 0.0, 0.0},
  Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Force[3] = {0.0, 0.0, 0.0}, Cp, thermal_conductivity, MaxNorm = 8.0,
  Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Grad_Temp[3] = {0.0, 0.0, 0.0},
  delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
#ifdef HAVE_MPI
  su2double MyAllBound_CDrag_Visc, MyAllBound_CLift_Visc, MyAllBound_CSideForce_Visc, MyAllBound_CMx_Visc, MyAllBound_CMy_Visc, MyAllBound_CMz_Visc, MyAllBound_CFx_Visc, MyAllBound_CFy_Visc, MyAllBound_CFz_Visc, MyAllBound_CT_Visc, MyAllBound_CQ_Visc, MyAllBound_HeatFlux_Visc, MyAllBound_MaxHeatFlux_Visc, *MySurface_CLift_Visc = NULL, *MySurface_CDrag_Visc = NULL, *MySurface_CSideForce_Visc = NULL, *MySurface_CEff_Visc = NULL, *MySurface_CFx_Visc = NULL, *MySurface_CFy_Visc = NULL, *MySurface_CFz_Visc = NULL, *MySurface_CMx_Visc = NULL, *MySurface_CMy_Visc = NULL, *MySurface_CMz_Visc = NULL;
#endif
  
  string Marker_Tag, Monitoring_Tag;
  
  su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefAreaCoeff    = config->GetRefAreaCoeff();
  su2double RefLengthMoment = config->GetRefLengthMoment();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  su2double *Origin         = config->GetRefOriginMoment(0);
  bool grid_movement        = config->GetGrid_Movement();
  bool compressible         = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible       = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface          = (config->GetKind_Regime() == FREESURFACE);
  su2double Prandtl_Lam     = config->GetPrandtl_Lam();
  
  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  RefTemp    = Temperature_Inf;
  RefDensity = Density_Inf;
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  } else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  /*--- Variables initialization ---*/
  
  AllBound_CDrag_Visc = 0.0;    AllBound_CLift_Visc = 0.0;       AllBound_CSideForce_Visc = 0.0;
  AllBound_CMx_Visc = 0.0;      AllBound_CMy_Visc = 0.0;         AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc = 0.0;      AllBound_CFy_Visc = 0.0;         AllBound_CFz_Visc = 0.0;
  AllBound_CT_Visc = 0.0;       AllBound_CQ_Visc = 0.0;          AllBound_CMerit_Visc = 0.0;
  AllBound_HeatFlux_Visc = 0.0; AllBound_MaxHeatFlux_Visc = 0.0; AllBound_CEff_Visc = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift_Visc[iMarker_Monitoring]      = 0.0; Surface_CDrag_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CSideForce_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0; Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0; Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0; Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
  }
  
  /*--- Loop over the Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CDrag_Visc[iMarker] = 0.0; CLift_Visc[iMarker] = 0.0;       CSideForce_Visc[iMarker] = 0.0;
      CMx_Visc[iMarker] = 0.0;   CMy_Visc[iMarker] = 0.0;         CMz_Visc[iMarker] = 0.0;
      CFx_Visc[iMarker] = 0.0;   CFy_Visc[iMarker] = 0.0;         CFz_Visc[iMarker] = 0.0;
      CT_Visc[iMarker] = 0.0;    CQ_Visc[iMarker] = 0.0;          CMerit_Visc[iMarker] = 0.0;
      Heat_Visc[iMarker] = 0.0;  MaxHeatFlux_Visc[iMarker] = 0.0; CEff_Visc[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
      MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
        
        Coord = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPointNormal]->GetCoord();
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
          Grad_Temp[iDim] = node[iPoint]->GetGradient_Primitive(0, iDim);
        }
        
        if (compressible) {
          Viscosity = node[iPoint]->GetLaminarViscosity();
          Density = node[iPoint]->GetDensity();
        }
        if (incompressible || freesurface) {
          Viscosity = node[iPoint]->GetLaminarViscosityInc();
          Density = node[iPoint]->GetDensityInc();
        }
        
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
          MomentDist[iDim] = Coord[iDim] - Origin[iDim];
        }
        
        /*--- Evaluate Tau ---*/
        
        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*delta[iDim][jDim];
          }
        }
        
        /*--- Project Tau in each surface element ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++) {
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
          }
        }
        
        /*--- Compute wall shear stress (using the stress tensor) ---*/
        
        TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++) TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
        WallShearStress = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallShearStress += TauTangent[iDim]*TauTangent[iDim];
        WallShearStress = sqrt(WallShearStress);
        
        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);
        
        /*--- Compute wall skin friction coefficient, and heat flux on the wall ---*/
        
        CSkinFriction[iMarker][iVertex] = WallShearStress / (0.5*RefDensity*RefVel2);
        
        /*--- Compute y+ and non-dimensional velocity ---*/
        
        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);
        
        /*--- Compute total and maximum heat flux on the wall ---*/
        
        if (compressible) {
          
          GradTemperature = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GradTemperature -= Grad_Temp[iDim]*UnitNormal[iDim];
          
          Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
          thermal_conductivity = Cp * Viscosity/Prandtl_Lam;
          HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature;
          Heat_Visc[iMarker] += HeatFlux[iMarker][iVertex]*Area;
          MaxHeatFlux_Visc[iMarker] += pow(HeatFlux[iMarker][iVertex], MaxNorm);
          
        }
        
        /*--- Note that y+, and heat are computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {
          
          /*--- Force computation ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim]*Area*factor;
            ForceViscous[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (iDim == 3) {
            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLengthMoment;
            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLengthMoment;
          
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if (Monitoring == YES) {
        if (nDim == 2) {
          CDrag_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CLift_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker]        = CLift_Visc[iMarker] / (CDrag_Visc[iMarker]+EPS);
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CT_Visc[iMarker]          = -CFx_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker]+EPS);
          MaxHeatFlux_Visc[iMarker] = pow(MaxHeatFlux_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          CDrag_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          CLift_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSideForce_Visc[iMarker]  = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker]        = CLift_Visc[iMarker]/(CDrag_Visc[iMarker] + EPS);
          CMx_Visc[iMarker]         = MomentViscous[0];
          CMy_Visc[iMarker]         = MomentViscous[1];
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CFz_Visc[iMarker]         = ForceViscous[2];
          CT_Visc[iMarker]          = -CFz_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker] + EPS);
          MaxHeatFlux_Visc[iMarker] = pow(MaxHeatFlux_Visc[iMarker], 1.0/MaxNorm);
        }
        
        AllBound_CDrag_Visc       += CDrag_Visc[iMarker];
        AllBound_CLift_Visc       += CLift_Visc[iMarker];
        AllBound_CSideForce_Visc  += CSideForce_Visc[iMarker];
        AllBound_CMx_Visc         += CMx_Visc[iMarker];
        AllBound_CMy_Visc         += CMy_Visc[iMarker];
        AllBound_CMz_Visc         += CMz_Visc[iMarker];
        AllBound_CFx_Visc         += CFx_Visc[iMarker];
        AllBound_CFy_Visc         += CFy_Visc[iMarker];
        AllBound_CFz_Visc         += CFz_Visc[iMarker];
        AllBound_CT_Visc          += CT_Visc[iMarker];
        AllBound_CQ_Visc          += CQ_Visc[iMarker];
        AllBound_HeatFlux_Visc    += Heat_Visc[iMarker];
        AllBound_MaxHeatFlux_Visc += pow(MaxHeatFlux_Visc[iMarker], MaxNorm);
        
        /*--- Compute the coefficients per surface ---*/
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CLift_Visc[iMarker_Monitoring]      += CLift_Visc[iMarker];
            Surface_CDrag_Visc[iMarker_Monitoring]      += CDrag_Visc[iMarker];
            Surface_CSideForce_Visc[iMarker_Monitoring] += CSideForce_Visc[iMarker];
            Surface_CEff_Visc[iMarker_Monitoring]       += CEff_Visc[iMarker];
            Surface_CFx_Visc[iMarker_Monitoring]        += CFx_Visc[iMarker];
            Surface_CFy_Visc[iMarker_Monitoring]        += CFy_Visc[iMarker];
            Surface_CFz_Visc[iMarker_Monitoring]        += CFz_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]        += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]        += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]        += CMz_Visc[iMarker];
          }
        }
        
      }
      
    }
  }
  
  /*--- Update some global coeffients ---*/
  
  AllBound_CEff_Visc = AllBound_CLift_Visc / (AllBound_CDrag_Visc + EPS);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  AllBound_MaxHeatFlux_Visc = pow(AllBound_MaxHeatFlux_Visc, 1.0/MaxNorm);
  
  
#ifdef HAVE_MPI
  
  /*--- Add AllBound information using all the nodes ---*/
  
  MyAllBound_CDrag_Visc        = AllBound_CDrag_Visc;                      AllBound_CDrag_Visc = 0.0;
  MyAllBound_CLift_Visc        = AllBound_CLift_Visc;                      AllBound_CLift_Visc = 0.0;
  MyAllBound_CSideForce_Visc   = AllBound_CSideForce_Visc;                 AllBound_CSideForce_Visc = 0.0;
  AllBound_CEff_Visc = 0.0;
  MyAllBound_CMx_Visc          = AllBound_CMx_Visc;                        AllBound_CMx_Visc = 0.0;
  MyAllBound_CMy_Visc          = AllBound_CMy_Visc;                        AllBound_CMy_Visc = 0.0;
  MyAllBound_CMz_Visc          = AllBound_CMz_Visc;                        AllBound_CMz_Visc = 0.0;
  MyAllBound_CFx_Visc          = AllBound_CFx_Visc;                        AllBound_CFx_Visc = 0.0;
  MyAllBound_CFy_Visc          = AllBound_CFy_Visc;                        AllBound_CFy_Visc = 0.0;
  MyAllBound_CFz_Visc          = AllBound_CFz_Visc;                        AllBound_CFz_Visc = 0.0;
  MyAllBound_CT_Visc           = AllBound_CT_Visc;                         AllBound_CT_Visc = 0.0;
  MyAllBound_CQ_Visc           = AllBound_CQ_Visc;                         AllBound_CQ_Visc = 0.0;
  AllBound_CMerit_Visc = 0.0;
  MyAllBound_HeatFlux_Visc     = AllBound_HeatFlux_Visc;                   AllBound_HeatFlux_Visc = 0.0;
  MyAllBound_MaxHeatFlux_Visc  = pow(AllBound_MaxHeatFlux_Visc, MaxNorm);  AllBound_MaxHeatFlux_Visc = 0.0;
  
  SU2_MPI::Allreduce(&MyAllBound_CDrag_Visc, &AllBound_CDrag_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CLift_Visc, &AllBound_CLift_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CSideForce_Visc, &AllBound_CSideForce_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Visc = AllBound_CLift_Visc / (AllBound_CDrag_Visc + EPS);
  SU2_MPI::Allreduce(&MyAllBound_CMx_Visc, &AllBound_CMx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMy_Visc, &AllBound_CMy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CMz_Visc, &AllBound_CMz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFx_Visc, &AllBound_CFx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFy_Visc, &AllBound_CFy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CFz_Visc, &AllBound_CFz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CT_Visc, &AllBound_CT_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_CQ_Visc, &AllBound_CQ_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  SU2_MPI::Allreduce(&MyAllBound_HeatFlux_Visc, &AllBound_HeatFlux_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyAllBound_MaxHeatFlux_Visc, &AllBound_MaxHeatFlux_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_MaxHeatFlux_Visc = pow(AllBound_MaxHeatFlux_Visc, 1.0/MaxNorm);
  
  /*--- Add the forces on the surfaces using all the nodes ---*/
  
  MySurface_CLift_Visc      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CDrag_Visc      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSideForce_Visc = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    
    MySurface_CLift_Visc[iMarker_Monitoring]      = Surface_CLift_Visc[iMarker_Monitoring];
    MySurface_CDrag_Visc[iMarker_Monitoring]      = Surface_CDrag_Visc[iMarker_Monitoring];
    MySurface_CSideForce_Visc[iMarker_Monitoring] = Surface_CSideForce_Visc[iMarker_Monitoring];
    MySurface_CEff_Visc[iMarker_Monitoring]       = Surface_CEff_Visc[iMarker_Monitoring];
    MySurface_CFx_Visc[iMarker_Monitoring]        = Surface_CFx_Visc[iMarker_Monitoring];
    MySurface_CFy_Visc[iMarker_Monitoring]        = Surface_CFy_Visc[iMarker_Monitoring];
    MySurface_CFz_Visc[iMarker_Monitoring]        = Surface_CFz_Visc[iMarker_Monitoring];
    MySurface_CMx_Visc[iMarker_Monitoring]        = Surface_CMx_Visc[iMarker_Monitoring];
    MySurface_CMy_Visc[iMarker_Monitoring]        = Surface_CMy_Visc[iMarker_Monitoring];
    MySurface_CMz_Visc[iMarker_Monitoring]        = Surface_CMz_Visc[iMarker_Monitoring];
    
    Surface_CLift_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CDrag_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CSideForce_Visc[iMarker_Monitoring] = 0.0;
    Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
  }
  
  SU2_MPI::Allreduce(MySurface_CLift_Visc, Surface_CLift_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CDrag_Visc, Surface_CDrag_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CSideForce_Visc, Surface_CSideForce_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
    Surface_CEff_Visc[iMarker_Monitoring] = Surface_CLift_Visc[iMarker_Monitoring] / (Surface_CDrag_Visc[iMarker_Monitoring] + EPS);
  SU2_MPI::Allreduce(MySurface_CFx_Visc, Surface_CFx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFy_Visc, Surface_CFy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CFz_Visc, Surface_CFz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMx_Visc, Surface_CMx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMy_Visc, Surface_CMy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(MySurface_CMz_Visc, Surface_CMz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  delete [] MySurface_CLift_Visc; delete [] MySurface_CDrag_Visc; delete [] MySurface_CSideForce_Visc;
  delete [] MySurface_CEff_Visc;  delete [] MySurface_CFx_Visc;   delete [] MySurface_CFy_Visc;
  delete [] MySurface_CFz_Visc;   delete [] MySurface_CMx_Visc;   delete [] MySurface_CMy_Visc;
  delete [] MySurface_CMz_Visc;
  
#endif
  
  /*--- Update the total coefficients (note that all the nodes have the same value)---*/
  
  Total_CDrag       += AllBound_CDrag_Visc;
  Total_CLift       += AllBound_CLift_Visc;
  Total_CSideForce  += AllBound_CSideForce_Visc;
  Total_CEff        = Total_CLift / (Total_CDrag + EPS);
  Total_CMx         += AllBound_CMx_Visc;
  Total_CMy         += AllBound_CMy_Visc;
  Total_CMz         += AllBound_CMz_Visc;
  Total_CFx         += AllBound_CFx_Visc;
  Total_CFy         += AllBound_CFy_Visc;
  Total_CFz         += AllBound_CFz_Visc;
  Total_CT          += AllBound_CT_Visc;
  Total_CQ          += AllBound_CQ_Visc;
  Total_CMerit       = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  Total_Heat        = AllBound_HeatFlux_Visc;
  Total_MaxHeat     = AllBound_MaxHeatFlux_Visc;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift[iMarker_Monitoring]      += Surface_CLift_Visc[iMarker_Monitoring];
    Surface_CDrag[iMarker_Monitoring]      += Surface_CDrag_Visc[iMarker_Monitoring];
    Surface_CSideForce[iMarker_Monitoring] += Surface_CSideForce_Visc[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CLift[iMarker_Monitoring] / (Surface_CDrag[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        += Surface_CFx_Visc[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        += Surface_CFy_Visc[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        += Surface_CFz_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        += Surface_CMz_Visc[iMarker_Monitoring];
  }
  
}

void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, jDim, iVar, jVar;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  su2double Wall_HeatFlux, dist_ij, *Coord_i, *Coord_j, theta2;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, *Normal, Area, Pressure = 0.0;
  su2double total_viscosity, div_vel, Density, tau_vel[3] = {0.0, 0.0, 0.0}, UnitNormal[3] = {0.0, 0.0, 0.0};
  su2double laminar_viscosity = 0.0, eddy_viscosity = 0.0, Grad_Vel[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}},
  tau[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Identify the boundary by string name ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Get the specified wall heat flux from config ---*/
  
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);
  
  /*--- Loop over all of the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Compute dual-grid area and boundary normal ---*/
      
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }
      
      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      
      if (compressible) node[iPoint]->SetVelocity_Old(Vector);
      if (incompressible || freesurface) node[iPoint]->SetVelocityInc_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      
      if (compressible) Res_Visc[nDim+1] = Wall_HeatFlux * Area;
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        
        Density  = node[iPoint]->GetSolution(0);
        if (compressible) {
          Pressure = node[iPoint]->GetPressure();
          laminar_viscosity = node[iPoint]->GetLaminarViscosity();
          eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        }
        if (incompressible || freesurface) {
          Pressure = node[iPoint]->GetPressureInc();
          laminar_viscosity = node[iPoint]->GetLaminarViscosityInc();
          eddy_viscosity    = node[iPoint]->GetEddyViscosityInc();
        }
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
        }
        
        /*--- Divergence of the velocity ---*/
        
        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim]+Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
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
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
          /*--- Now the Jacobian contribution related to the shear stress ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          /*--- Compute closest normal neighbor ---*/
          
          Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
          
          /*--- Get coordinates of i & nearest normal and compute distance ---*/
          
          Coord_i = geometry->node[iPoint]->GetCoord();
          Coord_j = geometry->node[Point_Normal]->GetCoord();
          
          dist_ij = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
          dist_ij = sqrt(dist_ij);
          
          theta2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            theta2 += UnitNormal[iDim]*UnitNormal[iDim];
          
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
          } else {
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
          
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
          
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

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij, theta2;
  su2double Twall, dTdn, dTdrho, thermal_conductivity;
  su2double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  su2double ProjGridVel, *GridVel, GridVel2, Pressure = 0.0, Density, Vel2;
  su2double total_viscosity, div_vel, tau_vel[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};
  su2double laminar_viscosity, eddy_viscosity, Grad_Vel[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
  tau[3][3] = {{0.0, 0.0, 0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}}, delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Identify the boundary ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Retrieve the specified wall temperature ---*/
  
  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();
  
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
      
      if (grid_movement) {
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
      
      if (compressible) node[iPoint]->SetVelocity_Old(Vector);
      if (incompressible || freesurface) node[iPoint]->SetVelocityInc_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Compute the normal gradient in temperature using Twall ---*/
      
      dTdn = -(node[Point_Normal]->GetPrimitive(0) - Twall)/dist_ij;
      
      /*--- Get transport coefficients ---*/
      
      laminar_viscosity    = node[iPoint]->GetLaminarViscosity();
      eddy_viscosity       = node[iPoint]->GetEddyViscosity();
      thermal_conductivity = Cp * ( laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);
      
      // work in progress on real-gases...
      //thermal_conductivity = node[iPoint]->GetThermalConductivity();
      //Cp = node[iPoint]->GetSpecificHeatCp();
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
        
        Density = node[iPoint]->GetPrimitive(nDim+2);
        Vel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Vel2 += node[iPoint]->GetPrimitive(iDim+1) * node[iPoint]->GetPrimitive(iDim+1);
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
        
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        
        Density  = node[iPoint]->GetSolution(0);
        if (compressible) {
          Pressure = node[iPoint]->GetPressure();
          laminar_viscosity = node[iPoint]->GetLaminarViscosity();
          eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        }
        if (incompressible || freesurface) {
          Pressure = node[iPoint]->GetPressureInc();
          laminar_viscosity = node[iPoint]->GetLaminarViscosityInc();
          eddy_viscosity    = node[iPoint]->GetEddyViscosityInc();
        }
        
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
        }
        
        /*--- Divergence of the velocity ---*/
        
        div_vel = 0.0; for (iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
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
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
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
          
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
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
