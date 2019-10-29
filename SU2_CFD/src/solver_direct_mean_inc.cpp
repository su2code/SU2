/*!
 * \file solution_direct_mean_inc.cpp
 * \brief Main subroutines for solving incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include "../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../include/variables/CIncEulerVariable.hpp"
#include "../include/variables/CIncNSVariable.hpp"

CIncEulerSolver::CIncEulerSolver(void) : CSolver() {
  /*--- Basic array initialization ---*/

  CD_Inv  = NULL; CL_Inv  = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt  = NULL; CL_Mnt  = NULL; CSF_Mnt = NULL;  CEff_Mnt = NULL;
  CMx_Mnt = NULL; CMy_Mnt = NULL; CMz_Mnt = NULL;
  CFx_Mnt = NULL; CFy_Mnt = NULL; CFz_Mnt = NULL;
  CoPx_Mnt = NULL; CoPy_Mnt = NULL; CoPz_Mnt = NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL; MomentMomentum = NULL;

  /*--- Surface based array initialization ---*/

  Surface_CL_Inv  = NULL; Surface_CD_Inv  = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt  = NULL; Surface_CD_Mnt  = NULL; Surface_CSF_Mnt = NULL; Surface_CEff_Mnt = NULL;
  Surface_CFx_Mnt = NULL; Surface_CFy_Mnt = NULL; Surface_CFz_Mnt = NULL;
  Surface_CMx_Mnt = NULL; Surface_CMy_Mnt = NULL; Surface_CMz_Mnt = NULL;

  Surface_CL  = NULL; Surface_CD  = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;
  
  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;
  
  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  Smatrix = NULL; Cvector = NULL;
  Preconditioner = NULL;

  FluidModel = NULL;

  SlidingState     = NULL;
  SlidingStateNodes = NULL;
  
  nodes = nullptr;
}

CIncEulerSolver::CIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  string filename = config->GetSolution_FileName();
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  bool fsi     = config->GetFSI_Simulation();
  bool multizone = config->GetMultizone_Problem();
  string filename_ = config->GetSolution_FileName();

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  unsigned short direct_diff = config->GetDirectDiff();

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone, ".dat");

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }

    /*--- Read and store the restart metadata. ---*/

//    Read_SU2_Restart_Metadata(geometry, config, false, filename_);
    
  }

  /*--- Basic array initialization ---*/

  CD_Inv  = NULL; CL_Inv  = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt  = NULL; CL_Mnt  = NULL; CSF_Mnt = NULL; CEff_Mnt = NULL;
  CMx_Mnt = NULL; CMy_Mnt = NULL; CMz_Mnt = NULL;
  CFx_Mnt = NULL; CFy_Mnt = NULL; CFz_Mnt = NULL;
  CoPx_Mnt= NULL;   CoPy_Mnt= NULL;   CoPz_Mnt= NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL;  MomentMomentum = NULL;

  /*--- Surface based array initialization ---*/

  Surface_CL_Inv  = NULL; Surface_CD_Inv  = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt  = NULL; Surface_CD_Mnt  = NULL; Surface_CSF_Mnt = NULL; Surface_CEff_Mnt= NULL;
  Surface_CFx_Mnt = NULL; Surface_CFy_Mnt = NULL; Surface_CFz_Mnt = NULL;
  Surface_CMx_Mnt = NULL; Surface_CMy_Mnt = NULL; Surface_CMz_Mnt = NULL;

  Surface_CL  = NULL; Surface_CD  = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Rotorcraft simulation array initialization ---*/

  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;

  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  Smatrix = NULL; Cvector = NULL;
  Preconditioner = NULL;

  /*--- Fluid model pointer initialization ---*/

  FluidModel = NULL;

  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure.
   * Incompressible flow, primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/
  
  nDim = geometry->GetnDim();
  
  nVar = nDim+2; nPrimVar = nDim+9; nPrimVarGrad = nDim+4;

  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nPrimVarGrad;
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
 
  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) 
    nVertex[iMarker] = geometry->nVertex[iMarker];
  
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(config, iMesh);
  
  /*--- Check if we are executing a verification case. If so, the
   VerificationSolution object will be instantiated for a particular
   option from the available library of verification solutions. Note
   that this is done after SetNondim(), as problem-specific initial
   parameters are needed by the solution constructors. ---*/
  
  SetVerificationSolution(nDim, nVar, config);
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]     = 0.0;
  Residual_RMS  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_Max  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  Res_Conv      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]     = 0.0;
  Res_Visc      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]     = 0.0;
  Res_Sour      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]     = 0.0;
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  
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
  
  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }

  Preconditioner = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar ++)
    Preconditioner[iVar] = new su2double[nVar];

  /*--- Initialize the solution and right-hand side vectors for storing
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
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
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
    
    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];
    
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
  
  /*--- Store the value of the Total Pressure at the inlet BC ---*/
  
  Inlet_Ttotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_Ttotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_Ttotal[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Total Temperature at the inlet BC ---*/
  
  Inlet_Ptotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_Ptotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_Ptotal[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Flow direction at the inlet BC ---*/
  
  Inlet_FlowDir = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_FlowDir[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_FlowDir[iMarker][iVertex] = new su2double [nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        Inlet_FlowDir[iMarker][iVertex][iDim] = 0;
      }
    }
  }
  
  /*--- Non-dimensional coefficients ---*/

  ForceInviscid  = new su2double[nDim];
  MomentInviscid = new su2double[3];
  CD_Inv         = new su2double[nMarker];
  CL_Inv         = new su2double[nMarker];
  CSF_Inv        = new su2double[nMarker];
  CMx_Inv        = new su2double[nMarker];
  CMy_Inv        = new su2double[nMarker];
  CMz_Inv        = new su2double[nMarker];
  CEff_Inv       = new su2double[nMarker];
  CFx_Inv        = new su2double[nMarker];
  CFy_Inv        = new su2double[nMarker];
  CFz_Inv        = new su2double[nMarker];
  CoPx_Inv       = new su2double[nMarker];
  CoPy_Inv       = new su2double[nMarker];
  CoPz_Inv       = new su2double[nMarker];

  ForceMomentum  = new su2double[nDim];
  MomentMomentum = new su2double[3];
  CD_Mnt         = new su2double[nMarker];
  CL_Mnt         = new su2double[nMarker];
  CSF_Mnt        = new su2double[nMarker];
  CMx_Mnt        = new su2double[nMarker];
  CMy_Mnt        = new su2double[nMarker];
  CMz_Mnt        = new su2double[nMarker];
  CEff_Mnt       = new su2double[nMarker];
  CFx_Mnt        = new su2double[nMarker];
  CFy_Mnt        = new su2double[nMarker];
  CFz_Mnt        = new su2double[nMarker];
  CoPx_Mnt       = new su2double[nMarker];
  CoPy_Mnt       = new su2double[nMarker];
  CoPz_Mnt       = new su2double[nMarker];

  Surface_CL_Inv   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv  = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL_Mnt   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Mnt   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Mnt  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Mnt = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Mnt  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Mnt  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Mnt  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Mnt  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Mnt  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Mnt  = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz  = new su2double[config->GetnMarker_Monitoring()];

  /*--- Rotorcraft coefficients ---*/

  CT_Inv           = new su2double[nMarker];
  CQ_Inv           = new su2double[nMarker];
  CMerit_Inv       = new su2double[nMarker];

  CT_Mnt           = new su2double[nMarker];
  CQ_Mnt           = new su2double[nMarker];
  CMerit_Mnt       = new su2double[nMarker];

  /*--- Init total coefficients ---*/

  Total_CD       = 0.0;    Total_CL           = 0.0;    Total_CSF            = 0.0;
  Total_CMx      = 0.0;    Total_CMy          = 0.0;    Total_CMz            = 0.0;
  Total_CoPx     = 0.0;    Total_CoPy         = 0.0;    Total_CoPz           = 0.0;
  Total_CEff     = 0.0;
  Total_CFx      = 0.0;    Total_CFy          = 0.0;    Total_CFz            = 0.0;
  Total_CT       = 0.0;    Total_CQ           = 0.0;    Total_CMerit         = 0.0;
  Total_MaxHeat  = 0.0;    Total_Heat         = 0.0;    Total_ComboObj       = 0.0;
  Total_CpDiff   = 0.0;    Total_HeatFluxDiff = 0.0;    Total_Custom_ObjFunc = 0.0;
  
  /*--- Read farfield conditions ---*/

  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();

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
  
  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++){
    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }

    }
  }

  /*--- Only initialize when there is a Marker_Fluid_Load defined
   *--- (this avoids overhead in all other cases while a more permanent structure is being developed) ---*/
  if((config->GetnMarker_Fluid_Load() > 0) && (MGLevel == MESH_0)){

    InitVertexTractionContainer(geometry, config);

    if (config->GetDiscrete_Adjoint())
      InitVertexTractionAdjointContainer(geometry, config);

  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CIncEulerVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (fsi || multizone){
    Residual_BGS      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
    Residual_Max_BGS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }
  }

  /*--- Define solver parameters needed for execution of destructor ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Communicate and store volume and the number of neighbors for
   any dual CVs that lie on on periodic markers. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
  }
  SetImplicitPeriodic(euler_implicit);
  if (iMesh == MESH_0) SetRotatePeriodic(true);
  
  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "INC.FLOW";

}

CIncEulerSolver::~CIncEulerSolver(void) {

  unsigned short iMarker, iVar;
  unsigned long iVertex;

  /*--- Array deallocation ---*/

  if (CD_Inv  != NULL)  delete [] CD_Inv;
  if (CL_Inv  != NULL)  delete [] CL_Inv;
  if (CSF_Inv != NULL)  delete [] CSF_Inv;
  if (CMx_Inv != NULL)  delete [] CMx_Inv;
  if (CMy_Inv != NULL)  delete [] CMy_Inv;
  if (CMz_Inv != NULL)  delete [] CMz_Inv;
  if (CFx_Inv != NULL)  delete [] CFx_Inv;
  if (CFy_Inv != NULL)  delete [] CFy_Inv;
  if (CFz_Inv != NULL)  delete [] CFz_Inv;
  if (CoPx_Inv != NULL) delete [] CoPx_Inv;
  if (CoPy_Inv != NULL) delete [] CoPy_Inv;
  if (CoPz_Inv != NULL) delete [] CoPz_Inv;

  if (Surface_CL_Inv   != NULL) delete [] Surface_CL_Inv;
  if (Surface_CD_Inv   != NULL) delete [] Surface_CD_Inv;
  if (Surface_CSF_Inv  != NULL) delete [] Surface_CSF_Inv;
  if (Surface_CEff_Inv != NULL) delete [] Surface_CEff_Inv;
  if (Surface_CFx_Inv  != NULL) delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv  != NULL) delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv  != NULL) delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv  != NULL) delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv  != NULL) delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv  != NULL) delete [] Surface_CMz_Inv;

  if (CD_Mnt  != NULL)  delete [] CD_Mnt;
  if (CL_Mnt  != NULL)  delete [] CL_Mnt;
  if (CSF_Mnt != NULL)  delete [] CSF_Mnt;
  if (CMx_Mnt != NULL)  delete [] CMx_Mnt;
  if (CMy_Mnt != NULL)  delete [] CMy_Mnt;
  if (CMz_Mnt != NULL)  delete [] CMz_Mnt;
  if (CFx_Mnt != NULL)  delete [] CFx_Mnt;
  if (CFy_Mnt != NULL)  delete [] CFy_Mnt;
  if (CFz_Mnt != NULL)  delete [] CFz_Mnt;
  if (CoPx_Mnt != NULL) delete [] CoPx_Mnt;
  if (CoPy_Mnt != NULL) delete [] CoPy_Mnt;
  if (CoPz_Mnt != NULL) delete [] CoPz_Mnt;

  if (Surface_CL_Mnt   != NULL) delete [] Surface_CL_Mnt;
  if (Surface_CD_Mnt   != NULL) delete [] Surface_CD_Mnt;
  if (Surface_CSF_Mnt  != NULL) delete [] Surface_CSF_Mnt;
  if (Surface_CEff_Mnt != NULL) delete [] Surface_CEff_Mnt;
  if (Surface_CFx_Mnt  != NULL) delete [] Surface_CFx_Mnt;
  if (Surface_CFy_Mnt  != NULL) delete [] Surface_CFy_Mnt;
  if (Surface_CFz_Mnt  != NULL) delete [] Surface_CFz_Mnt;
  if (Surface_CMx_Mnt  != NULL) delete [] Surface_CMx_Mnt;
  if (Surface_CMy_Mnt  != NULL) delete [] Surface_CMy_Mnt;
  if (Surface_CMz_Mnt  != NULL) delete [] Surface_CMz_Mnt;

  if (Surface_CL   != NULL) delete [] Surface_CL;
  if (Surface_CD   != NULL) delete [] Surface_CD;
  if (Surface_CSF  != NULL) delete [] Surface_CSF;
  if (Surface_CEff != NULL) delete [] Surface_CEff;
  if (Surface_CFx  != NULL) delete [] Surface_CFx;
  if (Surface_CFy  != NULL) delete [] Surface_CFy;
  if (Surface_CFz  != NULL) delete [] Surface_CFz;
  if (Surface_CMx  != NULL) delete [] Surface_CMx;
  if (Surface_CMy  != NULL) delete [] Surface_CMy;
  if (Surface_CMz  != NULL) delete [] Surface_CMz;
  
  if (CEff_Inv   != NULL) delete [] CEff_Inv;
  if (CMerit_Inv != NULL) delete [] CMerit_Inv;
  if (CT_Inv     != NULL) delete [] CT_Inv;
  if (CQ_Inv     != NULL) delete [] CQ_Inv;

  if (CEff_Mnt   != NULL) delete [] CEff_Mnt;
  if (CMerit_Mnt != NULL) delete [] CMerit_Mnt;
  if (CT_Mnt     != NULL) delete [] CT_Mnt;
  if (CQ_Mnt     != NULL) delete [] CQ_Mnt;

  if (ForceInviscid  != NULL) delete [] ForceInviscid;
  if (MomentInviscid != NULL) delete [] MomentInviscid;
  if (ForceMomentum  != NULL) delete [] ForceMomentum;
  if (MomentMomentum != NULL) delete [] MomentMomentum;

  if (Primitive   != NULL) delete [] Primitive;
  if (Primitive_i != NULL) delete [] Primitive_i;
  if (Primitive_j != NULL) delete [] Primitive_j;

  if (Preconditioner != NULL) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete [] Preconditioner[iVar];
    delete [] Preconditioner;
  }

  if (CPressure != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CPressure[iMarker];
    delete [] CPressure;
  }
  
  if (CPressureTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CPressureTarget[iMarker];
    delete [] CPressureTarget;
  }

  if (CharacPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex<nVertex[iMarker]; iVertex++)
        delete [] CharacPrimVar[iMarker][iVertex];
      delete [] CharacPrimVar[iMarker];
    }
    delete [] CharacPrimVar;
  }

  if (SlidingState != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }
  
  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
        if (SlidingStateNodes[iMarker] != NULL)
            delete [] SlidingStateNodes[iMarker];  
    }
    delete [] SlidingStateNodes;
  }

  if (Inlet_Ttotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Inlet_Ttotal[iMarker] != NULL)
        delete [] Inlet_Ttotal[iMarker];
    delete [] Inlet_Ttotal;
  }
  
  if (Inlet_Ptotal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Inlet_Ptotal[iMarker] != NULL)
        delete [] Inlet_Ptotal[iMarker];
    delete [] Inlet_Ptotal;
  }
  
  if (Inlet_FlowDir != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Inlet_FlowDir[iMarker] != NULL) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          delete [] Inlet_FlowDir[iMarker][iVertex];
        delete [] Inlet_FlowDir[iMarker];
      }
    }
    delete [] Inlet_FlowDir;
  }

  if (HeatFlux != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] HeatFlux[iMarker];
    }
    delete [] HeatFlux;
  }
  
  if (HeatFluxTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] HeatFluxTarget[iMarker];
    }
    delete [] HeatFluxTarget;
  }
  
  if (YPlus != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] YPlus[iMarker];
    }
    delete [] YPlus;
  }
  
  if (FluidModel != NULL) delete FluidModel;
  
  if (nodes != nullptr) delete nodes;
}

void CIncEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {
  
  su2double Temperature_FreeStream = 0.0,  ModVel_FreeStream = 0.0,Energy_FreeStream = 0.0,
  ModVel_FreeStreamND = 0.0, Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Pressure_Thermodynamic = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Temperature_Ref = 0.0, Velocity_Ref = 0.0, Time_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Heat_Flux_Ref = 0.0, Energy_Ref= 0.0, Pressure_FreeStreamND = 0.0, Pressure_ThermodynamicND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0, Specific_Heat_CpND = 0.0, Specific_Heat_CvND = 0.0, Thermal_Expansion_CoeffND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;
  
  unsigned short iDim, iVar;
  
  /*--- Local variables ---*/
  
  su2double Mach     = config->GetMach();
  su2double Reynolds = config->GetReynolds();
  
  bool unsteady      = (config->GetTime_Marching() != NO);
  bool viscous       = config->GetViscous();
  bool turbulent     = ((config->GetKind_Solver() == INC_RANS) ||
                        (config->GetKind_Solver() == DISC_ADJ_INC_RANS));
  bool tkeNeeded     = ((turbulent) && ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)));
  bool energy        = config->GetEnergy_Equation();
  bool boussinesq    = (config->GetKind_DensityModel() == BOUSSINESQ);

  /*--- Compute dimensional free-stream values. ---*/

  Density_FreeStream     = config->GetInc_Density_Init();     config->SetDensity_FreeStream(Density_FreeStream);
  Temperature_FreeStream = config->GetInc_Temperature_Init(); config->SetTemperature_FreeStream(Temperature_FreeStream);
  Pressure_FreeStream    = 0.0; config->SetPressure_FreeStream(Pressure_FreeStream);

  ModVel_FreeStream   = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ModVel_FreeStream += config->GetInc_Velocity_Init()[iDim]*config->GetInc_Velocity_Init()[iDim];
    config->SetVelocity_FreeStream(config->GetInc_Velocity_Init()[iDim],iDim);
  }
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Depending on the density model chosen, select a fluid model. ---*/

  switch (config->GetKind_FluidModel()) {

    case CONSTANT_DENSITY:

      FluidModel = new CConstantDensity(Density_FreeStream, config->GetSpecific_Heat_Cp());
      FluidModel->SetTDState_T(Temperature_FreeStream);
      break;

    case INC_IDEAL_GAS:

      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT/(config->GetMolecular_Weight()/1000.0));
      Pressure_Thermodynamic = Density_FreeStream*Temperature_FreeStream*config->GetGas_Constant();
      FluidModel = new CIncIdealGas(config->GetSpecific_Heat_Cp(), config->GetGas_Constant(), Pressure_Thermodynamic);
      FluidModel->SetTDState_T(Temperature_FreeStream);
      Pressure_Thermodynamic = FluidModel->GetPressure();
      config->SetPressure_Thermodynamic(Pressure_Thermodynamic);
      break;
      
    case INC_IDEAL_GAS_POLY:
      
      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT/(config->GetMolecular_Weight()/1000.0));
      Pressure_Thermodynamic = Density_FreeStream*Temperature_FreeStream*config->GetGas_Constant();
      FluidModel = new CIncIdealGasPolynomial(config->GetGas_Constant(), Pressure_Thermodynamic);
      if (viscous) {
        /*--- Variable Cp model via polynomial. ---*/
        for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++)
          config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(iVar), iVar);
        FluidModel->SetCpModel(config);
      }
      FluidModel->SetTDState_T(Temperature_FreeStream);
      Pressure_Thermodynamic = FluidModel->GetPressure();
      config->SetPressure_Thermodynamic(Pressure_Thermodynamic);
      break;

    default:

      SU2_MPI::Error("Fluid model not implemented for incompressible solver.", CURRENT_FUNCTION);
      break;
  }

  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
      To accomplish this, simply set the non-dimensional coefficients to the
      dimensional ones. This will be overruled later.---*/

    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());
    config->SetMu_ConstantND(config->GetMu_Constant());
    
    for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(iVar), iVar);

    /*--- Use the fluid model to compute the dimensional viscosity/conductivity. ---*/

    FluidModel->SetLaminarViscosityModel(config);
    Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
    config->SetViscosity_FreeStream(Viscosity_FreeStream);

    Reynolds = Density_FreeStream*ModVel_FreeStream/Viscosity_FreeStream; config->SetReynolds(Reynolds);

    /*--- Turbulence kinetic energy ---*/

    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }

  /*--- The non-dim. scheme for incompressible flows uses the following ref. values:
     Reference length      = 1 m (fixed by default, grid in meters)
     Reference density     = liquid density or freestream (input)
     Reference velocity    = liquid velocity or freestream (input)
     Reference temperature = liquid temperature or freestream (input)
     Reference pressure    = Reference density * Reference velocity * Reference velocity
     Reference viscosity   = Reference Density * Reference velocity * Reference length
     This is the same non-dim. scheme as in the compressible solver.
     Note that the Re and Re Length are not used as part of initialization. ---*/

  if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
    Density_Ref     = 1.0;
    Velocity_Ref    = 1.0;
    Temperature_Ref = 1.0;
    Pressure_Ref    = 1.0;
  }
  else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
    Density_Ref     = Density_FreeStream;
    Velocity_Ref    = ModVel_FreeStream;
    Temperature_Ref = Temperature_FreeStream;
    Pressure_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref;
  } 
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    Density_Ref     = config->GetInc_Density_Ref();
    Velocity_Ref    = config->GetInc_Velocity_Ref();
    Temperature_Ref = config->GetInc_Temperature_Ref();
    Pressure_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref;
  }
  config->SetDensity_Ref(Density_Ref);
  config->SetVelocity_Ref(Velocity_Ref);
  config->SetTemperature_Ref(Temperature_Ref);
  config->SetPressure_Ref(Pressure_Ref);

  /*--- More derived reference values ---*/
  
  Length_Ref       = 1.0;                                                config->SetLength_Ref(Length_Ref);
  Time_Ref         = Length_Ref/Velocity_Ref;                            config->SetTime_Ref(Time_Ref);
  Omega_Ref        = Velocity_Ref/Length_Ref;                            config->SetOmega_Ref(Omega_Ref);
  Force_Ref        = Velocity_Ref*Velocity_Ref/Length_Ref;               config->SetForce_Ref(Force_Ref);
  Heat_Flux_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref*Velocity_Ref; config->SetHeat_Flux_Ref(Heat_Flux_Ref);
  Gas_Constant_Ref = Velocity_Ref*Velocity_Ref/Temperature_Ref;          config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref    = Density_Ref*Velocity_Ref*Length_Ref;                config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref = Viscosity_Ref*Gas_Constant_Ref;                     config->SetConductivity_Ref(Conductivity_Ref);

  /*--- Get the freestream energy. Only useful if energy equation is active. ---*/

  Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
  config->SetEnergy_FreeStream(Energy_FreeStream);
  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute Mach number ---*/

  if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
    Mach = ModVel_FreeStream / sqrt(config->GetBulk_Modulus()/Density_FreeStream);
  } else {
    Mach = 0.0;
  }
  config->SetMach(Mach);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  
  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref(); config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Pressure_ThermodynamicND = Pressure_Thermodynamic/config->GetPressure_Ref(); config->SetPressure_ThermodynamicND(Pressure_ThermodynamicND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();   config->SetDensity_FreeStreamND(Density_FreeStreamND);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  Gas_ConstantND      = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);
  Specific_Heat_CpND  = config->GetSpecific_Heat_Cp()/Gas_Constant_Ref; config->SetSpecific_Heat_CpND(Specific_Heat_CpND);
  
  /*--- We assume that Cp = Cv for our incompressible fluids. ---*/
  Specific_Heat_CvND  = config->GetSpecific_Heat_Cp()/Gas_Constant_Ref; config->SetSpecific_Heat_CvND(Specific_Heat_CvND);
  
  Thermal_Expansion_CoeffND = config->GetThermal_Expansion_Coeff()*config->GetTemperature_Ref(); config->SetThermal_Expansion_CoeffND(Thermal_Expansion_CoeffND);

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
 
  /*--- Delete the original (dimensional) FluidModel object. No fluid is used for inscompressible cases. ---*/
  
  delete FluidModel;

  switch (config->GetKind_FluidModel()) {
      
    case CONSTANT_DENSITY:
      FluidModel = new CConstantDensity(Density_FreeStreamND, Specific_Heat_CpND);
      FluidModel->SetTDState_T(Temperature_FreeStreamND);
      break;

    case INC_IDEAL_GAS:
      FluidModel = new CIncIdealGas(Specific_Heat_CpND, Gas_ConstantND, Pressure_ThermodynamicND);
      FluidModel->SetTDState_T(Temperature_FreeStreamND);
      break;
      
    case INC_IDEAL_GAS_POLY:
      FluidModel = new CIncIdealGasPolynomial(Gas_ConstantND, Pressure_ThermodynamicND);
      if (viscous) {
        /*--- Variable Cp model via polynomial. ---*/
        config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(0)/Gas_Constant_Ref, 0);
        for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
          config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Gas_Constant_Ref, iVar);
        FluidModel->SetCpModel(config);
      }
      FluidModel->SetTDState_T(Temperature_FreeStreamND);
      break;
      
  }
  
  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
  
  if (viscous) {
    
    /*--- Constant viscosity model ---*/

    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);
    
    /*--- Sutherland's model ---*/
    
    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());
    
    /*--- Viscosity model via polynomial. ---*/

    config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(0)/Viscosity_Ref, 0);
    for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Viscosity_Ref, iVar);
    
    /*--- Constant thermal conductivity model ---*/

    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);
    
    /*--- Conductivity model via polynomial. ---*/

    config->SetKt_PolyCoeffND(config->GetKt_PolyCoeff(0)/Conductivity_Ref, 0);
    for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetKt_PolyCoeffND(config->GetKt_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Conductivity_Ref, iVar);
    
    /*--- Set up the transport property models. ---*/

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

    if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
      cout << "Incompressible flow: rho_ref, vel_ref, temp_ref, p_ref" << endl;
      cout << "are set to 1.0 in order to perform a dimensional calculation." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    }
    else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
      cout << "Incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
      cout << "are based on the initial values, p_ref = rho_ref*vel_ref^2." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    } 
    else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      cout << "Incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
      cout << "are user-provided reference values, p_ref = rho_ref*vel_ref^2." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using reference values." << endl;
    }
    cout << "The reference area for force coeffs. is " << config->GetRefArea() << " m^2." << endl;
    cout << "The reference length for force coeffs. is " << config->GetRefLength() << " m." << endl;

    cout << "The pressure is decomposed into thermodynamic and dynamic components." << endl;
    cout << "The initial value of the dynamic pressure is 0." << endl;

    cout << "Mach number: "<< config->GetMach();
    if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
      cout << ", computed using the Bulk modulus." << endl;
    } else {
      cout << ", computed using fluid speed of sound." << endl;
    }

    cout << "For external flows, the initial state is imposed at the far-field." << endl;
    cout << "Angle of attack (deg): "<< config->GetAoA() << ", computed using the initial velocity." << endl;
    cout << "Side slip angle (deg): "<< config->GetAoS() << ", computed using the initial velocity." << endl;

    if (viscous) { 
      cout << "Reynolds number per meter: " << config->GetReynolds() << ", computed using initial values."<< endl;
      cout << "Reynolds number is a byproduct of inputs only (not used internally)." << endl;
    }
    cout << "SI units only. The grid should be dimensional (meters)." << endl;
    
    switch (config->GetKind_DensityModel()) {
      
      case CONSTANT:
        if (energy) cout << "Energy equation is active and decoupled." << endl;
        else cout << "No energy equation." << endl;
        break;

      case BOUSSINESQ:
        if (energy) cout << "Energy equation is active and coupled through Boussinesq approx." << endl;
        break;

      case VARIABLE:
        if (energy) cout << "Energy equation is active and coupled for variable density." << endl;
        break;

    }
    
    stringstream NonDimTableOut, ModelTableOut;
    stringstream Unit;  
    
    cout << endl;
    PrintingToolbox::CTablePrinter ModelTable(&ModelTableOut);
    ModelTableOut <<"-- Models:"<< endl;

    ModelTable.AddColumn("Viscosity Model", 25);
    ModelTable.AddColumn("Conductivity Model", 26);
    ModelTable.AddColumn("Fluid Model", 25);
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
    
    if (viscous){
      
      switch(config->GetKind_ViscosityModel()){
      case CONSTANT_VISCOSITY:
        ModelTable << "CONSTANT_VISCOSITY";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << config->GetMu_Constant() << config->GetMu_Constant()/config->GetMu_ConstantND() << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case SUTHERLAND:
        ModelTable << "SUTHERLAND";        
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Ref. Viscosity" <<  config->GetMu_Ref() <<  config->GetViscosity_Ref() << Unit.str() << config->GetMu_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Temp." << config->GetMu_Temperature_Ref() <<  config->GetTemperature_Ref() << Unit.str() << config->GetMu_Temperature_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Const." << config->GetMu_S() << config->GetTemperature_Ref() << Unit.str() << config->GetMu_SND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;
        
      case POLYNOMIAL_VISCOSITY:
        ModelTable << "POLYNOMIAL_VISCOSITY";
        for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
          stringstream ss;
          ss << iVar;
          if (config->GetMu_PolyCoeff(iVar) != 0.0)
            NonDimTable << "Mu(T) Poly. Coeff. " + ss.str()  << config->GetMu_PolyCoeff(iVar) << config->GetMu_PolyCoeff(iVar)/config->GetMu_PolyCoeffND(iVar) << "-" << config->GetMu_PolyCoeffND(iVar);
        }
        Unit.str("");
        NonDimTable.PrintFooter();
        break;
      }

      switch(config->GetKind_ConductivityModel()){
      case CONSTANT_PRANDTL:
        ModelTable << "CONSTANT_PRANDTL";
        NonDimTable << "Prandtl (Lam.)"  << "-" << "-" << "-" << config->GetPrandtl_Lam();         
        Unit.str("");
        NonDimTable << "Prandtl (Turb.)" << "-" << "-" << "-" << config->GetPrandtl_Turb();         
        Unit.str("");
        NonDimTable.PrintFooter();
        break;
        
      case CONSTANT_CONDUCTIVITY:
        ModelTable << "CONSTANT_CONDUCTIVITY";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << config->GetKt_Constant() << config->GetKt_Constant()/config->GetKt_ConstantND() << Unit.str() << config->GetKt_ConstantND();         
        Unit.str("");
        NonDimTable.PrintFooter();
        break;
        
      case POLYNOMIAL_CONDUCTIVITY:
        ModelTable << "POLYNOMIAL_CONDUCTIVITY";
        for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
          stringstream ss;
          ss << iVar;
          if (config->GetKt_PolyCoeff(iVar) != 0.0)
            NonDimTable << "Kt(T) Poly. Coeff. " + ss.str()  << config->GetKt_PolyCoeff(iVar) << config->GetKt_PolyCoeff(iVar)/config->GetKt_PolyCoeffND(iVar) << "-" << config->GetKt_PolyCoeffND(iVar);
        }
        Unit.str("");
        NonDimTable.PrintFooter();  
        break;
      }
    } else {
      ModelTable << "-" << "-";
    }
    
    switch (config->GetKind_FluidModel()){
    case CONSTANT_DENSITY:
      ModelTable << "CONSTANT_DENSITY";
      if (energy){
        Unit << "N.m/kg.K";
        NonDimTable << "Spec. Heat (Cp)" << config->GetSpecific_Heat_Cp() << config->GetSpecific_Heat_Cp()/config->GetSpecific_Heat_CpND() << Unit.str() << config->GetSpecific_Heat_CpND();         
        Unit.str("");
      }
      if (boussinesq){
        Unit << "K^-1";
        NonDimTable << "Thermal Exp." << config->GetThermal_Expansion_Coeff() << config->GetThermal_Expansion_Coeff()/config->GetThermal_Expansion_CoeffND() << Unit.str() <<  config->GetThermal_Expansion_CoeffND();         
        Unit.str("");
      }
      Unit << "Pa";
      NonDimTable << "Bulk Modulus" << config->GetBulk_Modulus() << 1.0 << Unit.str() <<  config->GetBulk_Modulus();         
      Unit.str("");
      NonDimTable.PrintFooter();
      break;
      
    case INC_IDEAL_GAS:
      ModelTable << "INC_IDEAL_GAS";      
      Unit << "N.m/kg.K";
      NonDimTable << "Spec. Heat (Cp)" << config->GetSpecific_Heat_Cp() << config->GetSpecific_Heat_Cp()/config->GetSpecific_Heat_CpND() << Unit.str() << config->GetSpecific_Heat_CpND();         
      Unit.str("");
      Unit << "g/mol";
      NonDimTable << "Molecular weight" << config->GetMolecular_Weight()<< 1.0 << Unit.str() << config->GetMolecular_Weight();         
      Unit.str("");
      Unit << "N.m/kg.K";
      NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();         
      Unit.str("");
      Unit << "Pa";
      NonDimTable << "Therm. Pressure" << config->GetPressure_Thermodynamic() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_ThermodynamicND();         
      Unit.str("");
      NonDimTable.PrintFooter();
      break;
      
    case INC_IDEAL_GAS_POLY:
      ModelTable << "INC_IDEAL_GAS_POLY";             
      Unit.str("");
      Unit << "g/mol";
      NonDimTable << "Molecular weight" << config->GetMolecular_Weight()<< 1.0 << Unit.str() << config->GetMolecular_Weight();         
      Unit.str("");
      Unit << "N.m/kg.K";
      NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();         
      Unit.str("");
      Unit << "Pa";
      NonDimTable << "Therm. Pressure" << config->GetPressure_Thermodynamic() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_ThermodynamicND();         
      Unit.str("");
      for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
        stringstream ss;
        ss << iVar;
        if (config->GetCp_PolyCoeff(iVar) != 0.0)
          NonDimTable << "Cp(T) Poly. Coeff. " + ss.str()  << config->GetCp_PolyCoeff(iVar) << config->GetCp_PolyCoeff(iVar)/config->GetCp_PolyCoeffND(iVar) << "-" << config->GetCp_PolyCoeffND(iVar);
      }
      Unit.str("");
      NonDimTable.PrintFooter();
      break;
      
    }

    
    NonDimTableOut <<"-- Initial and free-stream conditions:"<< endl;
    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Dynamic Pressure" << config->GetPressure_FreeStream() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Total Pressure" << config->GetPressure_FreeStream() + 0.5*config->GetDensity_FreeStream()*config->GetModVel_FreeStream()*config->GetModVel_FreeStream() 
                << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND() + 0.5*config->GetDensity_FreeStreamND()*config->GetModVel_FreeStreamND()*config->GetModVel_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "kg/m^3";
    else if (config->GetSystemMeasurements() == US) Unit << "slug/ft^3";
    NonDimTable << "Density" << config->GetDensity_FreeStream() << config->GetDensity_Ref() << Unit.str() << config->GetDensity_FreeStreamND();
    Unit.str("");
    if (energy){
      if      (config->GetSystemMeasurements() == SI) Unit << "K";
      else if (config->GetSystemMeasurements() == US) Unit << "R";
      NonDimTable << "Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
      Unit.str("");
    }
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
      if      (config->GetSystemMeasurements() == SI) Unit << "W/m^2.K";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf/ft.s.R";
      NonDimTable << "Conductivity" << "-" << config->GetConductivity_Ref() << Unit.str() << "-";
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

void CIncEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {
  
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, *Solution_Fine, *Solution;
    
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool rans      = ((config->GetKind_Solver() == INC_RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_INC_RANS));
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  
  /*--- Check if a verification solution is to be computed. ---*/
  if ((VerificationSolution) && (TimeIter == 0) && !restart) {
    
    /*--- Loop over the multigrid levels. ---*/
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      
      /*--- Loop over all grid points. ---*/
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        
        /* Set the pointers to the coordinates and solution of this DOF. */
        const su2double *coor = geometry[iMesh]->node[iPoint]->GetCoord();
        su2double *solDOF     = solver_container[iMesh][FLOW_SOL]->GetNodes()->GetSolution(iPoint);
        
        /* Set the solution in this DOF to the initial condition provided by
           the verification solution class. This can be the exact solution,
           but this is not necessary. */
        VerificationSolution->GetInitialCondition(coor, solDOF);
      }
    }
  }
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  
  if (restart && (TimeIter == 0)) {
    
    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][FLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint,Solution);
      }
      solver_container[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver_container[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
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
            Solution_Fine = solver_container[iMesh-1][TURB_SOL]->GetNodes()->GetSolution(Point_Fine);
            for (iVar = 0; iVar < nVar_Turb; iVar++) {
              Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
            }
          }
          solver_container[iMesh][TURB_SOL]->GetNodes()->SetSolution(iPoint,Solution);
        }
        solver_container[iMesh][TURB_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION_EDDY);
        solver_container[iMesh][TURB_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION_EDDY);
        solver_container[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver_container[iMesh], config, iMesh);
      }
      delete [] Solution;
    }
    
  }
  
  /*--- The value of the solution for the first iteration of the dual time ---*/
  
  if (dual_time && (TimeIter == 0 || (restart && (long)TimeIter == (long)config->GetRestart_Iter()))) {
    
    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/
    
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
      solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n1();
      if (rans) {
        solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
        solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1();
      }
    }
    
    if ((restart && (long)TimeIter == (long)config->GetRestart_Iter()) &&
        (config->GetTime_Marching() == DT_STEPPING_2ND)) {
      
      /*--- Load an additional restart file for a 2nd-order restart ---*/
      
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetRestart_Iter()-1), true);
      
      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetRestart_Iter()-1), false);
      
      /*--- Push back this new solution to time level N. ---*/
      
      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
        if (rans) {
          solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
        }
      }
    }
  }
}

void CIncEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long ErrorCounter = 0;

  unsigned long InnerIter = config->GetInnerIter();
  bool cont_adjoint     = config->GetContinuous_Adjoint();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
  bool limiter          = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool center           = ((config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED));
  bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  bool outlet           = ((config->GetnMarker_Outlet() != 0));

  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Upwind second order reconstruction ---*/
  
  if ((muscl && !center) && (iMesh == MESH_0) && !Output) {
    
    /*--- Gradient computation ---*/
    
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
    }
    
    /*--- Limiter computation ---*/
    
    if ((limiter) && (iMesh == MESH_0) && !Output && !van_albada) {
      SetPrimitive_Limiter(geometry, config);
    }
    
  }
  
  /*--- Artificial dissipation ---*/
  
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Update the beta value based on the maximum velocity. ---*/

  SetBeta_Parameter(geometry, solver_container, config, iMesh);
  
  /*--- Compute properties needed for mass flow BCs. ---*/
  
  if (outlet) GetOutlet_Properties(geometry, config, iMesh, Output);

  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !disc_adjoint) Jacobian.SetValZero();

  /*--- Error message ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CIncEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) { }

unsigned long CIncEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  bool physical = true;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the non-physical points vector ---*/
    
    nodes->SetNon_Physical(iPoint,false);
    
    /*--- Incompressible flow, primitive variables ---*/
    
    physical = nodes->SetPrimVar(iPoint,FluidModel);

    /*--- Record any non-physical points. ---*/

    if (!physical) { nodes->SetNon_Physical(iPoint,true); ErrorCounter++; }
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }

  return ErrorCounter;
}

void CIncEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {
  
  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0,
  Mean_BetaInc2, Lambda, Local_Delta_Time,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
  
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool time_steping  = config->GetTime_Marching() == TIME_STEPPING;
  bool dual_time     = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetMax_Lambda_Inv(iPoint,0.0);
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    
    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);
    
    /*--- Mean Values ---*/

    Mean_ProjVel    = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
    Mean_BetaInc2   = 0.5 * (nodes->GetBetaInc2(iPoint)      + nodes->GetBetaInc2(jPoint));
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

    /*--- Adjustment for grid movement ---*/
    
    if (dynamic_grid) {
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
    if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(iPoint,Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(jPoint,Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      
      /*--- Mean Values ---*/

      Mean_ProjVel    = nodes->GetProjVel(iPoint,Normal);
      Mean_BetaInc2   = nodes->GetBetaInc2(iPoint);
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

      /*--- Adjustment for grid movement ---*/
      
      if (dynamic_grid) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        nodes->AddMax_Lambda_Inv(iPoint,Lambda);
      }
      
    }
    }
  }
  
  /*--- Local time-stepping: each element uses their own speed for steady state 
   simulations or for pseudo time steps in a dual time simulation. ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time  = config->GetCFL(iMesh)*Vol / nodes->GetMax_Lambda_Inv(iPoint);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time    = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time    = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      nodes->SetDelta_Time(iPoint,Local_Delta_Time);
    }
    else {
      nodes->SetDelta_Time(iPoint,0.0);
    }
    
  }
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
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
  
  /*--- For time-accurate simulations use the minimum delta time of the whole mesh (global) ---*/
  
  if (time_steping) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      
      /*--- Sets the regular CFL equal to the unsteady CFL ---*/
      
      config->SetCFL(iMesh,config->GetUnst_CFL());
      
      /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
       it computes the time step based on the unsteady CFL ---*/
      
      if (config->GetCFL(iMesh) == 0.0){
        nodes->SetDelta_Time(iPoint,config->GetDelta_UnstTime());
      } else {
        nodes->SetDelta_Time(iPoint,Global_Delta_Time);
      }
    }
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
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), nodes->GetDelta_Time(iPoint));
        nodes->SetDelta_Time(iPoint,Local_Delta_Time);
      }
    }
  
}

void CIncEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool jst_scheme  = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, set normal vectors, and number of neighbors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
    /*--- Set primitive variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(jPoint));
    
    /*--- Set the largest convective eigenvalue ---*/
    
    numerics->SetLambda(nodes->GetLambda(iPoint), nodes->GetLambda(jPoint));
    
    /*--- Set undivided laplacian and pressure-based sensor ---*/
    
    if (jst_scheme) {
      numerics->SetUndivided_Laplacian(nodes->GetUndivided_Laplacian(iPoint), nodes->GetUndivided_Laplacian(jPoint));
      numerics->SetSensor(nodes->GetSensor(iPoint), nodes->GetSensor(jPoint));
    }
    
    /*--- Grid movement ---*/
    
    if (dynamic_grid) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- Compute residuals, and Jacobians ---*/

    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    
    /*--- Update convective and artificial dissipation residuals ---*/

    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Store implicit contributions from the residual calculation. ---*/
    
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
  }
  
}

void CIncEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j,
  *V_i, *V_j, *S_i, *S_j, *Limiter_i = NULL, *Limiter_j = NULL, Non_Physical = 1.0;
  
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;
  
  unsigned long InnerIter = config->GetInnerIter();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  bool limiter          = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;

  /*--- Loop over all the edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Grid movement ---*/
    
    if (dynamic_grid)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    /*--- Get primitive variables ---*/
    
    V_i = nodes->GetPrimitive(iPoint); V_j = nodes->GetPrimitive(jPoint);
    S_i = nodes->GetSecondary(iPoint); S_j = nodes->GetSecondary(jPoint);

    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (muscl) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      Gradient_i = nodes->GetGradient_Primitive(iPoint);
      Gradient_j = nodes->GetGradient_Primitive(jPoint);
      if (limiter) {
        Limiter_i = nodes->GetLimiter_Primitive(iPoint);
        Limiter_j = nodes->GetLimiter_Primitive(jPoint);
      }
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        Non_Physical = nodes->GetNon_Physical(iPoint)*nodes->GetNon_Physical(jPoint);
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim]*Non_Physical;
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim]*Non_Physical;
        }
        if (limiter) {
          if (van_albada){
            Limiter_i[iVar] = (V_j[iVar]-V_i[iVar])*(2.0*Project_Grad_i + V_j[iVar]-V_i[iVar])/(4*Project_Grad_i*Project_Grad_i+(V_j[iVar]-V_i[iVar])*(V_j[iVar]-V_i[iVar])+EPS);
            Limiter_j[iVar] = (V_j[iVar]-V_i[iVar])*(-2.0*Project_Grad_j + V_j[iVar]-V_i[iVar])/(4*Project_Grad_j*Project_Grad_j+(V_j[iVar]-V_i[iVar])*(V_j[iVar]-V_i[iVar])+EPS);
          }
          Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }

      for (iVar = nPrimVarGrad; iVar < nPrimVar; iVar++) {
        Primitive_i[iVar] = V_i[iVar];
        Primitive_j[iVar] = V_j[iVar];
      }

      numerics->SetPrimitive(Primitive_i, Primitive_j);
      
    } else {
      
      /*--- Set conservative variables without reconstruction ---*/
      
      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);
      
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
  }
  
  /*--- Warning message about non-physical reconstructions. ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);
  }
  
}

void CIncEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  bool axisymmetric   = config->GetAxisymmetric();
  bool body_force     = config->GetBody_Force();
  bool boussinesq     = (config->GetKind_DensityModel() == BOUSSINESQ);
  bool viscous        = config->GetViscous();


  /*--- Initialize the source residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

  if (body_force) {

    /*--- Loop over all points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/

      numerics->SetConservative(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Set incompressible density  ---*/
      
      numerics->SetDensity(nodes->GetDensity(iPoint),
                           nodes->GetDensity(iPoint));

      /*--- Load the volume of the dual mesh cell ---*/

      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute the rotating frame source residual ---*/

      numerics->ComputeResidual(Residual, config);

      /*--- Add the source residual to the total ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }

  if (boussinesq) {

    /*--- Loop over all points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/

      numerics->SetConservative(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Set incompressible density  ---*/
      
      numerics->SetDensity(nodes->GetDensity(iPoint),
                           nodes->GetDensity(iPoint));

      /*--- Load the volume of the dual mesh cell ---*/

      numerics->SetVolume(geometry->node[iPoint]->GetVolume());

      /*--- Compute the rotating frame source residual ---*/

      numerics->ComputeResidual(Residual, config);

      /*--- Add the source residual to the total ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }

  if (rotating_frame) {

    /*--- Loop over all points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the conservative variables ---*/
      
      numerics->SetConservative(nodes->GetSolution(iPoint), NULL);
      
      /*--- Set incompressible density ---*/

      numerics->SetDensity(nodes->GetDensity(iPoint), 0.0);

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

    /*--- For viscous problems, we need an additional gradient. ---*/

    if (viscous) {

      su2double AuxVar, Total_Viscosity, yCoord, yVelocity;

      for (iPoint = 0; iPoint < nPoint; iPoint++) {

        yCoord          = geometry->node[iPoint]->GetCoord(1);
        yVelocity       = nodes->GetVelocity(iPoint,1);
        Total_Viscosity = (nodes->GetLaminarViscosity(iPoint) +
                           nodes->GetEddyViscosity(iPoint));

        if (yCoord > EPS) {
          AuxVar = Total_Viscosity*yVelocity/yCoord;
        } else {
          AuxVar = 0.0;
        }

        /*--- Set the auxilairy variable for this node. ---*/

        nodes->SetAuxVar(iPoint,AuxVar);

      }

      /*--- Compute the auxiliary variable gradient with GG or WLS. ---*/

      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
        SetAuxVar_Gradient_GG(geometry, config);
      }
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
        SetAuxVar_Gradient_LS(geometry, config);
      }
      
    }
    
    /*--- loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Conservative variables w/o reconstruction ---*/

      numerics->SetPrimitive(nodes->GetPrimitive(iPoint), NULL);

      /*--- Set incompressible density  ---*/
      
      numerics->SetDensity(nodes->GetDensity(iPoint),
                           nodes->GetDensity(iPoint));

      /*--- Set control volume ---*/
      
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                         geometry->node[iPoint]->GetCoord());

      /*--- If viscous, we need gradients for extra terms. ---*/

      if (viscous) {

        /*--- Gradient of the primitive variables ---*/

        numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), NULL);

        /*--- Load the aux variable gradient that we already computed. ---*/

        numerics->SetAuxVarGrad(nodes->GetAuxVarGradient(iPoint), NULL);
        
      }

      /*--- Compute Source term Residual ---*/
      
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Check if a verification solution is to be computed. ---*/
  
  if (VerificationSolution) {
    if ( VerificationSolution->IsManufacturedSolution() ) {
      
      /*--- Get the physical time. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching()) time = config->GetPhysicalTime();
      
      /*--- Loop over points ---*/
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
        /*--- Get control volume size. ---*/
        su2double Volume = geometry->node[iPoint]->GetVolume();
        
        /*--- Get the current point coordinates. ---*/
        const su2double *coor = geometry->node[iPoint]->GetCoord();
        
        /*--- Get the MMS source term. ---*/
        vector<su2double> sourceMan(nVar,0.0);
        VerificationSolution->GetMMSSourceTerm(coor, time, sourceMan.data());
        
        /*--- Compute the residual for this control volume. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] = sourceMan[iVar]*Volume;
        }
        
        /*--- Subtract Residual ---*/
        LinSysRes.SubtractBlock(iPoint, Residual);
        
      }
    }
  }
  
}

void CIncEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */
  
  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = nodes->GetSolution(iPoint) */
  /* Get the volume (or area in 2D) associated with iPoint = nodes->GetVolume(iPoint) */
  /* Get the vector of geometric coordinates of point iPoint = nodes->GetCoord(iPoint) */
  
}

void CIncEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
  
  su2double *Normal, Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0,
  Mean_BetaInc2, Lambda, ProjVel, ProjVel_i, ProjVel_j, *GridVel, *GridVel_i, *GridVel_j;
  
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    nodes->SetLambda(iPoint,0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);
    
    /*--- Mean Values ---*/

    Mean_ProjVel    = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
    Mean_BetaInc2   = 0.5 * (nodes->GetBetaInc2(iPoint)      + nodes->GetBetaInc2(jPoint));
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

    /*--- Adjustment for grid movement ---*/
    
    if (dynamic_grid) {
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
    if (geometry->node[iPoint]->GetDomain()) nodes->AddLambda(iPoint,Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddLambda(jPoint,Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel    = nodes->GetProjVel(iPoint,Normal);
      Mean_BetaInc2   = nodes->GetBetaInc2(iPoint);
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);
      
      /*--- Adjustment for grid movement ---*/
      
      if (dynamic_grid) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        nodes->AddLambda(iPoint,Lambda);
      }
      
    }
    }
  }
  
  /*--- Correct the eigenvalue values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
  }
  
  /*--- MPI parallelization ---*/
  
  InitiateComms(geometry, config, MAX_EIGENVALUE);
  CompleteComms(geometry, config, MAX_EIGENVALUE);
  
}

void CIncEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, jPoint, iEdge;
  su2double *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;
  
  Diff = new su2double[nVar];
  
  nodes->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Solution differences ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = nodes->GetSolution(iPoint,iVar) - nodes->GetSolution(jPoint,iVar);
    
    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) nodes->SubtractUnd_Lapl(iPoint,Diff);
      if (geometry->node[jPoint]->GetDomain()) nodes->AddUnd_Lapl(jPoint,Diff);
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) nodes->SubtractUnd_Lapl(iPoint,Diff);
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) nodes->AddUnd_Lapl(jPoint,Diff);
    
  }
  
  /*--- Correct the Laplacian values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
  }
  
  /*--- MPI parallelization ---*/
  
  InitiateComms(geometry, config, UNDIVIDED_LAPLACIAN);
  CompleteComms(geometry, config, UNDIVIDED_LAPLACIAN);
  
  delete [] Diff;
  
}

void CIncEulerSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  su2double Pressure_i = 0.0, Pressure_j = 0.0;
  bool boundary_i, boundary_j;
  
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

    Pressure_i = nodes->GetDensity(iPoint);
    Pressure_j = nodes->GetDensity(jPoint);

    boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both on the boundary ---*/
    
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      
      if (geometry->node[iPoint]->GetDomain()) {
        iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i);
        jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j);
      }
      
      if (geometry->node[jPoint]->GetDomain()) {
        iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j);
        jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j);
      }
      
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) {
        iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i);
        jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j);
      }
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) {
        iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j);
        jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j);
      }
    
  }
  
  /*--- Correct the sensor values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
  }
  
  /*--- Set pressure switch for each point ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetSensor(iPoint,fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
  
  /*--- MPI parallelization ---*/
  
  InitiateComms(geometry, config, SENSOR);
  CompleteComms(geometry, config, SENSOR);
  
}

void CIncEulerSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double Pressure = 0.0, *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord,
  factor, RefVel2 = 0.0, RefDensity = 0.0, RefPressure,
  Force[3] = {0.0,0.0,0.0};
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0}, MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;

  bool axisymmetric = config->GetAxisymmetric();

  string Marker_Tag, Monitoring_Tag;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Inv, MyAllBound_CL_Inv, MyAllBound_CSF_Inv, MyAllBound_CMx_Inv, MyAllBound_CMy_Inv, MyAllBound_CMz_Inv, MyAllBound_CoPx_Inv, MyAllBound_CoPy_Inv, MyAllBound_CoPz_Inv, MyAllBound_CFx_Inv, MyAllBound_CFy_Inv, MyAllBound_CFz_Inv, MyAllBound_CT_Inv, MyAllBound_CQ_Inv, *MySurface_CL_Inv = NULL, *MySurface_CD_Inv = NULL, *MySurface_CSF_Inv = NULL, *MySurface_CEff_Inv = NULL, *MySurface_CFx_Inv = NULL, *MySurface_CFy_Inv = NULL, *MySurface_CFz_Inv = NULL, *MySurface_CMx_Inv = NULL, *MySurface_CMy_Inv = NULL, *MySurface_CMz_Inv = NULL;
#endif

  su2double Alpha     = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta      = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea   = config->GetRefArea();
  su2double RefLength = config->GetRefLength();

  su2double *Origin = NULL;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }

  /*--- Evaluate reference values for non-dimensionalization.
   For dimensional or non-dim based on initial values, use
   the far-field state (inf). For a custom non-dim based
   on user-provided reference values, use the ref values
   to compute the forces. ---*/

  if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) || 
      (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
    RefDensity  = Density_Inf;
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    RefDensity = config->GetInc_Density_Ref();
    RefVel2    = config->GetInc_Velocity_Ref()*config->GetInc_Velocity_Ref();
  }

  /*--- Reference pressure is always the far-field value. ---*/

  RefPressure = Pressure_Inf;

  /*--- Compute factor for force coefficients. ---*/

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*-- Variables initialization ---*/

  Total_CD   = 0.0; Total_CL  = 0.0; Total_CSF = 0.0; Total_CEff = 0.0;
  Total_CMx  = 0.0; Total_CMy = 0.0; Total_CMz = 0.0;
  Total_CoPx = 0.0; Total_CoPy = 0.0;  Total_CoPz = 0.0;
  Total_CFx  = 0.0; Total_CFy = 0.0; Total_CFz = 0.0;
  Total_CT   = 0.0; Total_CQ  = 0.0; Total_CMerit = 0.0;
  Total_Heat = 0.0; Total_MaxHeat = 0.0;

  AllBound_CD_Inv   = 0.0; AllBound_CL_Inv  = 0.0;  AllBound_CSF_Inv    = 0.0;
  AllBound_CMx_Inv  = 0.0; AllBound_CMy_Inv = 0.0;  AllBound_CMz_Inv    = 0.0;
  AllBound_CoPx_Inv = 0.0; AllBound_CoPy_Inv = 0.0; AllBound_CoPz_Inv = 0.0;
  AllBound_CFx_Inv  = 0.0; AllBound_CFy_Inv = 0.0;  AllBound_CFz_Inv    = 0.0;
  AllBound_CT_Inv   = 0.0; AllBound_CQ_Inv  = 0.0;  AllBound_CMerit_Inv = 0.0;
  AllBound_CEff_Inv = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Inv[iMarker_Monitoring]  = 0.0; Surface_CD_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring] = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring] = 0.0; Surface_CFy_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring] = 0.0; Surface_CMx_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring] = 0.0; Surface_CMz_Inv[iMarker_Monitoring]  = 0.0;
    
    Surface_CL[iMarker_Monitoring]  = 0.0; Surface_CD[iMarker_Monitoring]   = 0.0;
    Surface_CSF[iMarker_Monitoring] = 0.0; Surface_CEff[iMarker_Monitoring] = 0.0;
    Surface_CFx[iMarker_Monitoring] = 0.0; Surface_CFy[iMarker_Monitoring]  = 0.0;
    Surface_CFz[iMarker_Monitoring] = 0.0; Surface_CMx[iMarker_Monitoring]  = 0.0;
    Surface_CMy[iMarker_Monitoring] = 0.0; Surface_CMz[iMarker_Monitoring]  = 0.0;
  }

  /*--- Loop over the Euler and Navier-Stokes markers ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Boundary   = config->GetMarker_All_KindBC(iMarker);
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

    if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
        (Boundary == ISOTHERMAL) || (Boundary == NEARFIELD_BOUNDARY) ||
        (Boundary == CHT_WALL_INTERFACE) ||
        (Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {

      /*--- Forces initialization at each Marker ---*/

      CD_Inv[iMarker]   = 0.0; CL_Inv[iMarker]  = 0.0;  CSF_Inv[iMarker]    = 0.0;
      CMx_Inv[iMarker]  = 0.0; CMy_Inv[iMarker] = 0.0;  CMz_Inv[iMarker]    = 0.0;
      CoPx_Inv[iMarker] = 0.0; CoPy_Inv[iMarker] = 0.0; CoPz_Inv[iMarker] = 0.0;
      CFx_Inv[iMarker]  = 0.0; CFy_Inv[iMarker] = 0.0;  CFz_Inv[iMarker]    = 0.0;
      CT_Inv[iMarker]   = 0.0; CQ_Inv[iMarker]  = 0.0;  CMerit_Inv[iMarker] = 0.0;
      CEff_Inv[iMarker] = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
      MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
      MomentX_Force[0] = 0.0; MomentX_Force[1] = 0.0; MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0; MomentY_Force[1] = 0.0; MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0; MomentZ_Force[1] = 0.0; MomentZ_Force[2] = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        Pressure = nodes->GetPressure(iPoint);

        CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefArea;

        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();

          for (iDim = 0; iDim < nDim; iDim++) {
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/

          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf) * Normal[iDim] * factor * AxiFactor;
            ForceInviscid[iDim] += Force[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (nDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1]  += (-Force[1]*Coord[2]);
            MomentX_Force[2]  += (Force[2]*Coord[1]);

            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2]  += (-Force[2]*Coord[0]);
            MomentY_Force[0]  += (Force[0]*Coord[2]);
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0]  += (-Force[0]*Coord[1]);
          MomentZ_Force[1]  += (Force[1]*Coord[0]);
        }

      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if (Monitoring == YES) {

        if (Boundary != NEARFIELD_BOUNDARY) {
          if (nDim == 2) {
            CD_Inv[iMarker]  =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            CL_Inv[iMarker]  = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
            CEff_Inv[iMarker]   = CL_Inv[iMarker] / (CD_Inv[iMarker]+EPS);
            CMz_Inv[iMarker]    = MomentInviscid[2];
            CoPx_Inv[iMarker]   = MomentZ_Force[1];
            CoPy_Inv[iMarker]   = -MomentZ_Force[0];
            CFx_Inv[iMarker]    = ForceInviscid[0];
            CFy_Inv[iMarker]    = ForceInviscid[1];
            CT_Inv[iMarker]     = -CFx_Inv[iMarker];
            CQ_Inv[iMarker]     = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker] = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }
          if (nDim == 3) {
            CD_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
            CL_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
            CSF_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            CEff_Inv[iMarker]       = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
            CMx_Inv[iMarker]        = MomentInviscid[0];
            CMy_Inv[iMarker]        = MomentInviscid[1];
            CMz_Inv[iMarker]        = MomentInviscid[2];
            CoPx_Inv[iMarker]    = -MomentY_Force[0];
            CoPz_Inv[iMarker]    = MomentY_Force[2];
            CFx_Inv[iMarker]        = ForceInviscid[0];
            CFy_Inv[iMarker]        = ForceInviscid[1];
            CFz_Inv[iMarker]        = ForceInviscid[2];
            CT_Inv[iMarker]         = -CFz_Inv[iMarker];
            CQ_Inv[iMarker]         = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker]     = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
          }

          AllBound_CD_Inv     += CD_Inv[iMarker];
          AllBound_CL_Inv     += CL_Inv[iMarker];
          AllBound_CSF_Inv    += CSF_Inv[iMarker];
          AllBound_CEff_Inv    = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
          AllBound_CMx_Inv    += CMx_Inv[iMarker];
          AllBound_CMy_Inv    += CMy_Inv[iMarker];
          AllBound_CMz_Inv    += CMz_Inv[iMarker];
          AllBound_CoPx_Inv   += CoPx_Inv[iMarker];
          AllBound_CoPy_Inv   += CoPy_Inv[iMarker];
          AllBound_CoPz_Inv   += CoPz_Inv[iMarker];
          AllBound_CFx_Inv    += CFx_Inv[iMarker];
          AllBound_CFy_Inv    += CFy_Inv[iMarker];
          AllBound_CFz_Inv    += CFz_Inv[iMarker];
          AllBound_CT_Inv     += CT_Inv[iMarker];
          AllBound_CQ_Inv     += CQ_Inv[iMarker];
          AllBound_CMerit_Inv  = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);

          /*--- Compute the coefficients per surface ---*/

          for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
            Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            if (Marker_Tag == Monitoring_Tag) {
              Surface_CL_Inv[iMarker_Monitoring]   += CL_Inv[iMarker];
              Surface_CD_Inv[iMarker_Monitoring]   += CD_Inv[iMarker];
              Surface_CSF_Inv[iMarker_Monitoring]  += CSF_Inv[iMarker];
              Surface_CEff_Inv[iMarker_Monitoring]  = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
              Surface_CFx_Inv[iMarker_Monitoring]  += CFx_Inv[iMarker];
              Surface_CFy_Inv[iMarker_Monitoring]  += CFy_Inv[iMarker];
              Surface_CFz_Inv[iMarker_Monitoring]  += CFz_Inv[iMarker];
              Surface_CMx_Inv[iMarker_Monitoring]  += CMx_Inv[iMarker];
              Surface_CMy_Inv[iMarker_Monitoring]  += CMy_Inv[iMarker];
              Surface_CMz_Inv[iMarker_Monitoring]  += CMz_Inv[iMarker];
            }
          }

        }

      }

    }
  }

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  MyAllBound_CD_Inv        = AllBound_CD_Inv;        AllBound_CD_Inv = 0.0;
  MyAllBound_CL_Inv        = AllBound_CL_Inv;        AllBound_CL_Inv = 0.0;
  MyAllBound_CSF_Inv   = AllBound_CSF_Inv;   AllBound_CSF_Inv = 0.0;
  AllBound_CEff_Inv = 0.0;
  MyAllBound_CMx_Inv          = AllBound_CMx_Inv;          AllBound_CMx_Inv = 0.0;
  MyAllBound_CMy_Inv          = AllBound_CMy_Inv;          AllBound_CMy_Inv = 0.0;
  MyAllBound_CMz_Inv          = AllBound_CMz_Inv;          AllBound_CMz_Inv = 0.0;
  MyAllBound_CoPx_Inv          = AllBound_CoPx_Inv;          AllBound_CoPx_Inv = 0.0;
  MyAllBound_CoPy_Inv          = AllBound_CoPy_Inv;          AllBound_CoPy_Inv = 0.0;
  MyAllBound_CoPz_Inv          = AllBound_CoPz_Inv;          AllBound_CoPz_Inv = 0.0;
  MyAllBound_CFx_Inv          = AllBound_CFx_Inv;          AllBound_CFx_Inv = 0.0;
  MyAllBound_CFy_Inv          = AllBound_CFy_Inv;          AllBound_CFy_Inv = 0.0;
  MyAllBound_CFz_Inv          = AllBound_CFz_Inv;          AllBound_CFz_Inv = 0.0;
  MyAllBound_CT_Inv           = AllBound_CT_Inv;           AllBound_CT_Inv = 0.0;
  MyAllBound_CQ_Inv           = AllBound_CQ_Inv;           AllBound_CQ_Inv = 0.0;
  AllBound_CMerit_Inv = 0.0;

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(&MyAllBound_CD_Inv, &AllBound_CD_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CL_Inv, &AllBound_CL_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CSF_Inv, &AllBound_CSF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_CEff_Inv = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
    SU2_MPI::Allreduce(&MyAllBound_CMx_Inv, &AllBound_CMx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMy_Inv, &AllBound_CMy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMz_Inv, &AllBound_CMz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPx_Inv, &AllBound_CoPx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPy_Inv, &AllBound_CoPy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPz_Inv, &AllBound_CoPz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFx_Inv, &AllBound_CFx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFy_Inv, &AllBound_CFy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFz_Inv, &AllBound_CFz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CT_Inv, &AllBound_CT_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CQ_Inv, &AllBound_CQ_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_CMerit_Inv = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
  }
  
  /*--- Add the forces on the surfaces using all the nodes ---*/

  MySurface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Inv = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CL_Inv[iMarker_Monitoring]      = Surface_CL_Inv[iMarker_Monitoring];
    MySurface_CD_Inv[iMarker_Monitoring]      = Surface_CD_Inv[iMarker_Monitoring];
    MySurface_CSF_Inv[iMarker_Monitoring] = Surface_CSF_Inv[iMarker_Monitoring];
    MySurface_CEff_Inv[iMarker_Monitoring]       = Surface_CEff_Inv[iMarker_Monitoring];
    MySurface_CFx_Inv[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    MySurface_CFy_Inv[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    MySurface_CFz_Inv[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    MySurface_CMx_Inv[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    MySurface_CMy_Inv[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    MySurface_CMz_Inv[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];

    Surface_CL_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CD_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0;
    Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
  }

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(MySurface_CL_Inv, Surface_CL_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CD_Inv, Surface_CD_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CSF_Inv, Surface_CSF_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
      Surface_CEff_Inv[iMarker_Monitoring] = Surface_CL_Inv[iMarker_Monitoring] / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
    SU2_MPI::Allreduce(MySurface_CFx_Inv, Surface_CFx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CFy_Inv, Surface_CFy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CFz_Inv, Surface_CFz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMx_Inv, Surface_CMx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMy_Inv, Surface_CMy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMz_Inv, Surface_CMz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  
  delete [] MySurface_CL_Inv; delete [] MySurface_CD_Inv; delete [] MySurface_CSF_Inv;
  delete [] MySurface_CEff_Inv;  delete [] MySurface_CFx_Inv;   delete [] MySurface_CFy_Inv;
  delete [] MySurface_CFz_Inv;   delete [] MySurface_CMx_Inv;   delete [] MySurface_CMy_Inv;
  delete [] MySurface_CMz_Inv;

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/

  Total_CD            = AllBound_CD_Inv;
  Total_CL            = AllBound_CL_Inv;
  Total_CSF           = AllBound_CSF_Inv;
  Total_CEff          = Total_CL / (Total_CD + EPS);
  Total_CMx           = AllBound_CMx_Inv;
  Total_CMy           = AllBound_CMy_Inv;
  Total_CMz           = AllBound_CMz_Inv;
  Total_CoPx          = AllBound_CoPx_Inv;
  Total_CoPy          = AllBound_CoPy_Inv;
  Total_CoPz          = AllBound_CoPz_Inv;
  Total_CFx           = AllBound_CFx_Inv;
  Total_CFy           = AllBound_CFy_Inv;
  Total_CFz           = AllBound_CFz_Inv;
  Total_CT            = AllBound_CT_Inv;
  Total_CQ            = AllBound_CQ_Inv;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]      = Surface_CL_Inv[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      = Surface_CD_Inv[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] = Surface_CSF_Inv[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CL_Inv[iMarker_Monitoring] / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        = Surface_CFx_Inv[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        = Surface_CFy_Inv[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        = Surface_CFz_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        = Surface_CMz_Inv[iMarker_Monitoring];
  }

}

void CIncEulerSolver::Momentum_Forces(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord, Area,
  factor, RefVel2 = 0.0, RefDensity = 0.0,
  Force[3] = {0.0,0.0,0.0}, Velocity[3], MassFlow, Density;
  string Marker_Tag, Monitoring_Tag;
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0}, MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Mnt, MyAllBound_CL_Mnt, MyAllBound_CSF_Mnt,
  MyAllBound_CMx_Mnt, MyAllBound_CMy_Mnt, MyAllBound_CMz_Mnt,
  MyAllBound_CoPx_Mnt, MyAllBound_CoPy_Mnt, MyAllBound_CoPz_Mnt,
  MyAllBound_CFx_Mnt, MyAllBound_CFy_Mnt, MyAllBound_CFz_Mnt, MyAllBound_CT_Mnt,
  MyAllBound_CQ_Mnt,
  *MySurface_CL_Mnt = NULL, *MySurface_CD_Mnt = NULL, *MySurface_CSF_Mnt = NULL,
  *MySurface_CEff_Mnt = NULL, *MySurface_CFx_Mnt = NULL, *MySurface_CFy_Mnt = NULL,
  *MySurface_CFz_Mnt = NULL,
  *MySurface_CMx_Mnt = NULL, *MySurface_CMy_Mnt = NULL,  *MySurface_CMz_Mnt = NULL;
#endif

  su2double Alpha     = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta      = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea   = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double *Origin = NULL;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }
  bool axisymmetric          = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dimensional or non-dim based on initial values, use
   the far-field state (inf). For a custom non-dim based
   on user-provided reference values, use the ref values
   to compute the forces. ---*/

  if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) || 
      (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
    RefDensity  = Density_Inf;
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    RefDensity = config->GetInc_Density_Ref();
    RefVel2    = config->GetInc_Velocity_Ref()*config->GetInc_Velocity_Ref();
  }

  /*--- Compute factor for force coefficients. ---*/

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*-- Variables initialization ---*/

  AllBound_CD_Mnt = 0.0;        AllBound_CL_Mnt = 0.0; AllBound_CSF_Mnt = 0.0;
  AllBound_CMx_Mnt = 0.0;          AllBound_CMy_Mnt = 0.0;   AllBound_CMz_Mnt = 0.0;
  AllBound_CoPx_Mnt = 0.0;          AllBound_CoPy_Mnt = 0.0;   AllBound_CoPz_Mnt = 0.0;
  AllBound_CFx_Mnt = 0.0;          AllBound_CFy_Mnt = 0.0;   AllBound_CFz_Mnt = 0.0;
  AllBound_CT_Mnt = 0.0;           AllBound_CQ_Mnt = 0.0;    AllBound_CMerit_Mnt = 0.0;
  AllBound_CEff_Mnt = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Mnt[iMarker_Monitoring]      = 0.0; Surface_CD_Mnt[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring] = 0.0; Surface_CEff_Mnt[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Mnt[iMarker_Monitoring]        = 0.0; Surface_CFy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Mnt[iMarker_Monitoring]        = 0.0; Surface_CMy_Mnt[iMarker_Monitoring]        = 0.0; Surface_CMz_Mnt[iMarker_Monitoring]        = 0.0;
  }

  /*--- Loop over the Inlet / Outlet Markers  ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Boundary   = config->GetMarker_All_KindBC(iMarker);
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

    if ((Boundary == INLET_FLOW) || (Boundary == OUTLET_FLOW) ||
        (Boundary == ACTDISK_INLET) || (Boundary == ACTDISK_OUTLET)||
        (Boundary == ENGINE_INFLOW) || (Boundary == ENGINE_EXHAUST)) {

      /*--- Forces initialization at each Marker ---*/

      CD_Mnt[iMarker] = 0.0;        CL_Mnt[iMarker] = 0.0; CSF_Mnt[iMarker] = 0.0;
      CMx_Mnt[iMarker] = 0.0;          CMy_Mnt[iMarker] = 0.0;   CMz_Mnt[iMarker] = 0.0;
      CFx_Mnt[iMarker] = 0.0;          CFy_Mnt[iMarker] = 0.0;   CFz_Mnt[iMarker] = 0.0;
      CoPx_Mnt[iMarker] = 0.0;         CoPy_Mnt[iMarker] = 0.0;  CoPz_Mnt[iMarker] = 0.0;
      CT_Mnt[iMarker] = 0.0;           CQ_Mnt[iMarker] = 0.0;    CMerit_Mnt[iMarker] = 0.0;
      CEff_Mnt[iMarker] = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) ForceMomentum[iDim] = 0.0;
      MomentMomentum[0] = 0.0; MomentMomentum[1] = 0.0; MomentMomentum[2] = 0.0;
      MomentX_Force[0] = 0.0; MomentX_Force[1] = 0.0; MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0; MomentY_Force[1] = 0.0; MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0; MomentZ_Force[1] = 0.0; MomentZ_Force[2] = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          Density   = nodes->GetDensity(iPoint);

          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);

          MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim]   = nodes->GetVelocity(iPoint,iDim);
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
            MassFlow -= Normal[iDim]*Velocity[iDim]*Density;
          }

          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/

          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = MassFlow * Velocity[iDim] * factor * AxiFactor;
            ForceMomentum[iDim] += Force[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (iDim == 3) {
            MomentMomentum[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLength;
            MomentX_Force[1]  += (-Force[1]*Coord[2]);
            MomentX_Force[2]  += (Force[2]*Coord[1]);

            MomentMomentum[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLength;
            MomentY_Force[2]  += (-Force[2]*Coord[0]);
            MomentY_Force[0]  += (Force[0]*Coord[2]);
          }
          MomentMomentum[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLength;
          MomentZ_Force[0]  += (-Force[0]*Coord[1]);
          MomentZ_Force[1]  += (Force[1]*Coord[0]);

        }

      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if (Monitoring == YES) {

        if (nDim == 2) {
          CD_Mnt[iMarker]  =  ForceMomentum[0]*cos(Alpha) + ForceMomentum[1]*sin(Alpha);
          CL_Mnt[iMarker]  = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[1]*cos(Alpha);
          CEff_Mnt[iMarker]   = CL_Mnt[iMarker] / (CD_Mnt[iMarker]+EPS);
          CMz_Mnt[iMarker]    = MomentInviscid[2];
          CFx_Mnt[iMarker]    = ForceMomentum[0];
          CFy_Mnt[iMarker]    = ForceMomentum[1];
          CoPx_Mnt[iMarker]   = MomentZ_Force[1];
          CoPy_Mnt[iMarker]   = -MomentZ_Force[0];
          CT_Mnt[iMarker]     = -CFx_Mnt[iMarker];
          CQ_Mnt[iMarker]     = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker] = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }
        if (nDim == 3) {
          CD_Mnt[iMarker]      =  ForceMomentum[0]*cos(Alpha)*cos(Beta) + ForceMomentum[1]*sin(Beta) + ForceMomentum[2]*sin(Alpha)*cos(Beta);
          CL_Mnt[iMarker]      = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[2]*cos(Alpha);
          CSF_Mnt[iMarker] = -ForceMomentum[0]*sin(Beta)*cos(Alpha) + ForceMomentum[1]*cos(Beta) - ForceMomentum[2]*sin(Beta)*sin(Alpha);
          CEff_Mnt[iMarker]       = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
          CMx_Mnt[iMarker]        = MomentInviscid[0];
          CMy_Mnt[iMarker]        = MomentInviscid[1];
          CMz_Mnt[iMarker]        = MomentInviscid[2];
          CFx_Mnt[iMarker]        = ForceMomentum[0];
          CFy_Mnt[iMarker]        = ForceMomentum[1];
          CFz_Mnt[iMarker]        = ForceMomentum[2];
          CoPx_Mnt[iMarker]       = -MomentY_Force[0];
          CoPz_Mnt[iMarker]       =  MomentY_Force[2];
          CT_Mnt[iMarker]         = -CFz_Mnt[iMarker];
          CQ_Mnt[iMarker]         = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker]     = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }

        AllBound_CD_Mnt        += CD_Mnt[iMarker];
        AllBound_CL_Mnt        += CL_Mnt[iMarker];
        AllBound_CSF_Mnt   += CSF_Mnt[iMarker];
        AllBound_CEff_Mnt          = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
        AllBound_CMx_Mnt          += CMx_Mnt[iMarker];
        AllBound_CMy_Mnt          += CMy_Mnt[iMarker];
        AllBound_CMz_Mnt          += CMz_Mnt[iMarker];
        AllBound_CFx_Mnt          += CFx_Mnt[iMarker];
        AllBound_CFy_Mnt          += CFy_Mnt[iMarker];
        AllBound_CFz_Mnt          += CFz_Mnt[iMarker];
        AllBound_CoPx_Mnt         += CoPx_Mnt[iMarker];
        AllBound_CoPy_Mnt         += CoPy_Mnt[iMarker];
        AllBound_CoPz_Mnt         += CoPz_Mnt[iMarker];
        AllBound_CT_Mnt           += CT_Mnt[iMarker];
        AllBound_CQ_Mnt           += CQ_Mnt[iMarker];
        AllBound_CMerit_Mnt        += AllBound_CT_Mnt / (AllBound_CQ_Mnt + EPS);

        /*--- Compute the coefficients per surface ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Mnt[iMarker_Monitoring]      += CL_Mnt[iMarker];
            Surface_CD_Mnt[iMarker_Monitoring]      += CD_Mnt[iMarker];
            Surface_CSF_Mnt[iMarker_Monitoring] += CSF_Mnt[iMarker];
            Surface_CEff_Mnt[iMarker_Monitoring]        = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
            Surface_CFx_Mnt[iMarker_Monitoring]        += CFx_Mnt[iMarker];
            Surface_CFy_Mnt[iMarker_Monitoring]        += CFy_Mnt[iMarker];
            Surface_CFz_Mnt[iMarker_Monitoring]        += CFz_Mnt[iMarker];
            Surface_CMx_Mnt[iMarker_Monitoring]        += CMx_Mnt[iMarker];
            Surface_CMy_Mnt[iMarker_Monitoring]        += CMy_Mnt[iMarker];
            Surface_CMz_Mnt[iMarker_Monitoring]        += CMz_Mnt[iMarker];
          }
        }

      }


    }
  }

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  MyAllBound_CD_Mnt        = AllBound_CD_Mnt;        AllBound_CD_Mnt = 0.0;
  MyAllBound_CL_Mnt        = AllBound_CL_Mnt;        AllBound_CL_Mnt = 0.0;
  MyAllBound_CSF_Mnt   = AllBound_CSF_Mnt;   AllBound_CSF_Mnt = 0.0;
  AllBound_CEff_Mnt = 0.0;
  MyAllBound_CMx_Mnt          = AllBound_CMx_Mnt;          AllBound_CMx_Mnt = 0.0;
  MyAllBound_CMy_Mnt          = AllBound_CMy_Mnt;          AllBound_CMy_Mnt = 0.0;
  MyAllBound_CMz_Mnt          = AllBound_CMz_Mnt;          AllBound_CMz_Mnt = 0.0;
  MyAllBound_CFx_Mnt          = AllBound_CFx_Mnt;          AllBound_CFx_Mnt = 0.0;
  MyAllBound_CFy_Mnt          = AllBound_CFy_Mnt;          AllBound_CFy_Mnt = 0.0;
  MyAllBound_CFz_Mnt          = AllBound_CFz_Mnt;          AllBound_CFz_Mnt = 0.0;
  MyAllBound_CoPx_Mnt         = AllBound_CoPx_Mnt;         AllBound_CoPx_Mnt = 0.0;
  MyAllBound_CoPy_Mnt         = AllBound_CoPy_Mnt;         AllBound_CoPy_Mnt = 0.0;
  MyAllBound_CoPz_Mnt         = AllBound_CoPz_Mnt;         AllBound_CoPz_Mnt = 0.0;
  MyAllBound_CT_Mnt           = AllBound_CT_Mnt;           AllBound_CT_Mnt = 0.0;
  MyAllBound_CQ_Mnt           = AllBound_CQ_Mnt;           AllBound_CQ_Mnt = 0.0;
  AllBound_CMerit_Mnt = 0.0;

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(&MyAllBound_CD_Mnt, &AllBound_CD_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CL_Mnt, &AllBound_CL_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CSF_Mnt, &AllBound_CSF_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_CEff_Mnt = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
    SU2_MPI::Allreduce(&MyAllBound_CMx_Mnt, &AllBound_CMx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMy_Mnt, &AllBound_CMy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMz_Mnt, &AllBound_CMz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFx_Mnt, &AllBound_CFx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFy_Mnt, &AllBound_CFy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFz_Mnt, &AllBound_CFz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPx_Mnt, &AllBound_CoPx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPy_Mnt, &AllBound_CoPy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPz_Mnt, &AllBound_CoPz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CT_Mnt, &AllBound_CT_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CQ_Mnt, &AllBound_CQ_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_CMerit_Mnt = AllBound_CT_Mnt / (AllBound_CQ_Mnt + EPS);
  }
  
  /*--- Add the forces on the surfaces using all the nodes ---*/

  MySurface_CL_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Mnt = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CL_Mnt[iMarker_Monitoring]      = Surface_CL_Mnt[iMarker_Monitoring];
    MySurface_CD_Mnt[iMarker_Monitoring]      = Surface_CD_Mnt[iMarker_Monitoring];
    MySurface_CSF_Mnt[iMarker_Monitoring] = Surface_CSF_Mnt[iMarker_Monitoring];
    MySurface_CEff_Mnt[iMarker_Monitoring]       = Surface_CEff_Mnt[iMarker_Monitoring];
    MySurface_CFx_Mnt[iMarker_Monitoring]        = Surface_CFx_Mnt[iMarker_Monitoring];
    MySurface_CFy_Mnt[iMarker_Monitoring]        = Surface_CFy_Mnt[iMarker_Monitoring];
    MySurface_CFz_Mnt[iMarker_Monitoring]        = Surface_CFz_Mnt[iMarker_Monitoring];
    MySurface_CMx_Mnt[iMarker_Monitoring]        = Surface_CMx_Mnt[iMarker_Monitoring];
    MySurface_CMy_Mnt[iMarker_Monitoring]        = Surface_CMy_Mnt[iMarker_Monitoring];
    MySurface_CMz_Mnt[iMarker_Monitoring]        = Surface_CMz_Mnt[iMarker_Monitoring];

    Surface_CL_Mnt[iMarker_Monitoring]      = 0.0;
    Surface_CD_Mnt[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring] = 0.0;
    Surface_CEff_Mnt[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Mnt[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Mnt[iMarker_Monitoring]        = 0.0;
  }

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(MySurface_CL_Mnt, Surface_CL_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CD_Mnt, Surface_CD_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CSF_Mnt, Surface_CSF_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
      Surface_CEff_Mnt[iMarker_Monitoring] = Surface_CL_Mnt[iMarker_Monitoring] / (Surface_CD_Mnt[iMarker_Monitoring] + EPS);
    SU2_MPI::Allreduce(MySurface_CFx_Mnt, Surface_CFx_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CFy_Mnt, Surface_CFy_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CFz_Mnt, Surface_CFz_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMx_Mnt, Surface_CMx_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMy_Mnt, Surface_CMy_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMz_Mnt, Surface_CMz_Mnt, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  
  delete [] MySurface_CL_Mnt; delete [] MySurface_CD_Mnt; delete [] MySurface_CSF_Mnt;
  delete [] MySurface_CEff_Mnt;  delete [] MySurface_CFx_Mnt;   delete [] MySurface_CFy_Mnt;
  delete [] MySurface_CFz_Mnt;
  delete [] MySurface_CMx_Mnt;   delete [] MySurface_CMy_Mnt;  delete [] MySurface_CMz_Mnt;

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/

  Total_CD            += AllBound_CD_Mnt;
  Total_CL            += AllBound_CL_Mnt;
  Total_CSF           += AllBound_CSF_Mnt;
  Total_CEff          = Total_CL / (Total_CD + EPS);
  Total_CMx           += AllBound_CMx_Mnt;
  Total_CMy           += AllBound_CMy_Mnt;
  Total_CMz           += AllBound_CMz_Mnt;
  Total_CFx           += AllBound_CFx_Mnt;
  Total_CFy           += AllBound_CFy_Mnt;
  Total_CFz           += AllBound_CFz_Mnt;
  Total_CoPx          += AllBound_CoPx_Mnt;
  Total_CoPy          += AllBound_CoPy_Mnt;
  Total_CoPz          += AllBound_CoPz_Mnt;
  Total_CT            += AllBound_CT_Mnt;
  Total_CQ            += AllBound_CQ_Mnt;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]   += Surface_CL_Mnt[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]   += Surface_CD_Mnt[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring]  += Surface_CSF_Mnt[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring] += Surface_CL_Mnt[iMarker_Monitoring] / (Surface_CD_Mnt[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]  += Surface_CFx_Mnt[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]  += Surface_CFy_Mnt[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]  += Surface_CFz_Mnt[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]  += Surface_CMx_Mnt[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]  += Surface_CMy_Mnt[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]  += Surface_CMz_Mnt[iMarker_Monitoring];
  }

}

void CIncEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar, jVar;
  unsigned long iPoint;
  
  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetContinuous_Adjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = (geometry->node[iPoint]->GetVolume() +
           geometry->node[iPoint]->GetPeriodicVolume());
    Delta = nodes->GetDelta_Time(iPoint) / Vol;

    Res_TruncError = nodes->GetResTruncError(iPoint);
    Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        Res = 0.0;
        for (jVar = 0; jVar < nVar; jVar ++ )
          Res += Preconditioner[iVar][jVar]*(Residual[jVar] + Res_TruncError[jVar]);
        nodes->AddSolution(iPoint,iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
  }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
  /*--- For verification cases, compute the global error metrics. ---*/
  
  ComputeVerificationError(geometry, config);
  
}

void CIncEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar, jVar;
  unsigned long iPoint;
  
  bool adjoint = config->GetContinuous_Adjoint();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = (geometry->node[iPoint]->GetVolume() +
           geometry->node[iPoint]->GetPeriodicVolume());
    Delta = nodes->GetDelta_Time(iPoint) / Vol;
    
    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    local_Residual = LinSysRes.GetBlock(iPoint);


    if (!adjoint) {
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        Res = 0.0;
        for (jVar = 0; jVar < nVar; jVar ++ )
          Res += Preconditioner[iVar][jVar]*(local_Residual[jVar] + local_Res_TruncError[jVar]);
        nodes->AddSolution(iPoint,iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
      }
    }
  }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
  /*--- For verification cases, compute the global error metrics. ---*/
  
  ComputeVerificationError(geometry, config);
  
}

void CIncEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar, jVar;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Delta, *local_Res_TruncError, Vol;
  
  bool adjoint = config->GetContinuous_Adjoint();
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/
    
    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    
    /*--- Read the volume ---*/
    
    Vol = (geometry->node[iPoint]->GetVolume() +
           geometry->node[iPoint]->GetPeriodicVolume());
    
    /*--- Apply the preconditioner and add to the diagonal. ---*/
    
    if (nodes->GetDelta_Time(iPoint) != 0.0) {
      Delta = Vol / nodes->GetDelta_Time(iPoint);
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        for (jVar = 0; jVar < nVar; jVar ++ ) {
          Preconditioner[iVar][jVar] = Delta*Preconditioner[iVar][jVar];
        }
      }
      Jacobian.AddBlock(iPoint, iPoint, Preconditioner);
    } else {
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
  
  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
  /*--- Set the residual --- */
  
  valResidual = System.GetResidual();
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        nodes->AddSolution(iPoint,iVar, config->GetRelaxation_Factor_Flow()*LinSysSol[iPoint*nVar+iVar]);
      }
    }
  }
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
  /*--- For verification cases, compute the global error metrics. ---*/

  ComputeVerificationError(geometry, config);
  
}

void CIncEulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal, Vol;
  
  /*--- Incompressible flow, primitive variables nDim+4, (P, vx, vy, vz, T, rho, beta) ---*/
  
  PrimVar_Vertex = new su2double [nPrimVarGrad];
  PrimVar_i = new su2double [nPrimVarGrad];
  PrimVar_j = new su2double [nPrimVarGrad];
  
  /*--- Set Gradient_Primitive to zero ---*/
  nodes->SetGradient_PrimitiveZero();
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
      PrimVar_j[iVar] = nodes->GetPrimitive(jPoint,iVar);
    }
    
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          nodes->AddGradient_Primitive(iPoint,iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          nodes->SubtractGradient_Primitive(jPoint,iVar, iDim, Partial_Res);
      }
    }
  }
  
  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
       (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          PrimVar_Vertex[iVar] = nodes->GetPrimitive(iPoint,iVar);
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
            nodes->SubtractGradient_Primitive(iPoint,iVar, iDim, Partial_Res);
          }
      }
    }
    }
  }
  
  /*--- Correct the gradient values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRIM_GG);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRIM_GG);
  }
  
  /*--- Update gradient value ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Get the volume, which may include periodic components. ---*/
    
    Vol = (geometry->node[iPoint]->GetVolume() +
           geometry->node[iPoint]->GetPeriodicVolume());
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = nodes->GetGradient_Primitive(iPoint,iVar, iDim)/Vol;
        nodes->SetGradient_Primitive(iPoint,iVar, iDim, Partial_Gradient);
      }
    }
  }
  
  /*--- Communicate the gradient values via MPI. ---*/
  
  InitiateComms(geometry, config, PRIMITIVE_GRADIENT);
  CompleteComms(geometry, config, PRIMITIVE_GRADIENT);
  
  delete [] PrimVar_Vertex;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j;
  su2double r11, r12, r13, r22, r23, r23_a, r23_b, r33, weight;
  su2double z11, z12, z13, z22, z23, z33, detR2;
  bool singular;
  
  /*--- Clear Rmatrix, which could eventually be computed once
   and stored for static meshes, as well as the prim gradient. ---*/

  nodes->SetRmatrixZero();
  nodes->SetGradient_PrimitiveZero();

  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get primitives from CVariable ---*/
    
    PrimVar_i = nodes->GetPrimitive(iPoint);
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;
    
    /*--- Clear Rmatrix, which could eventually be computed once
     and stored for static meshes, as well as the prim gradient. ---*/
    
    AD::StartPreacc();
    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    AD::SetPreaccIn(Coord_i, nDim);
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      PrimVar_j = nodes->GetPrimitive(jPoint);
      
      AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);
      
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (weight != 0.0) {
        
        nodes->AddRmatrix(iPoint,0, 0, (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight);
        nodes->AddRmatrix(iPoint,0, 1, (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight);
        nodes->AddRmatrix(iPoint,1, 1, (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight);
        
        if (nDim == 3) {
          nodes->AddRmatrix(iPoint,0, 2, (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight);
          nodes->AddRmatrix(iPoint,1, 2, (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight);
          nodes->AddRmatrix(iPoint,2, 1, (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight);
          nodes->AddRmatrix(iPoint,2, 2, (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight);
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          for (iDim = 0; iDim < nDim; iDim++) {
            nodes->AddGradient_Primitive(iPoint,iVar,iDim, (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight);
          }
        }
        
      }
    }
    AD::SetPreaccOut(nodes->GetRmatrix(iPoint), nDim, nDim);
    AD::SetPreaccOut(nodes->GetGradient_Primitive(iPoint), nPrimVarGrad, nDim);
    AD::EndPreacc();
  }
  
  /*--- Correct the gradient values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRIM_LS);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRIM_LS);
  }
  
  /*--- Second loop over points of the grid to compute final gradient ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set the value of the singular ---*/
    
    singular = false;
    
    /*--- Entries of upper triangular matrix R ---*/
    
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
    
    r11 = nodes->GetRmatrix(iPoint,0,0);
    r12 = nodes->GetRmatrix(iPoint,0,1);
    r22 = nodes->GetRmatrix(iPoint,1,1);
    
    AD::StartPreacc();
    AD::SetPreaccIn(r11);
    AD::SetPreaccIn(r12);
    AD::SetPreaccIn(r22);
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    
    if (nDim == 3) {
      r13   = nodes->GetRmatrix(iPoint,0,2);
      r23_a = nodes->GetRmatrix(iPoint,1,2);
      r23_b = nodes->GetRmatrix(iPoint,2,1);
      r33   = nodes->GetRmatrix(iPoint,2,2);
      
      AD::SetPreaccIn(r13);
      AD::SetPreaccIn(r23_a);
      AD::SetPreaccIn(r23_b);
      AD::SetPreaccIn(r33);
      
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
    
    AD::SetPreaccOut(Smatrix, nDim, nDim);
    AD::EndPreacc();
    
    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Cvector[iVar][iDim] = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          Cvector[iVar][iDim] += Smatrix[iDim][jDim]*nodes->GetGradient_Primitive(iPoint,iVar, jDim);
        }
      }
    }
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetGradient_Primitive(iPoint,iVar, iDim, Cvector[iVar][iDim]);
      }
    }
    
  }
  
  /*--- Communicate the gradient values via MPI. ---*/
  
  InitiateComms(geometry, config, PRIMITIVE_GRADIENT);
  CompleteComms(geometry, config, PRIMITIVE_GRADIENT);
  
}

void CIncEulerSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j,
  *Primitive, *Primitive_i, *Primitive_j,
  *LocalMinPrimitive = NULL, *LocalMaxPrimitive = NULL,
  *GlobalMinPrimitive = NULL, *GlobalMaxPrimitive = NULL,
  dave, LimK, eps2, eps1, dm, dp, du, y, limiter;
  
#ifdef CODI_REVERSE_TYPE
  bool TapeActive = false;

  if (config->GetDiscrete_Adjoint() && config->GetFrozen_Limiter_Disc()) {
    /*--- If limiters are frozen do not record the computation ---*/
    TapeActive = AD::globalTape.isActive();
    AD::StopRecording();
  }
#endif
  
  dave = config->GetRefElemLength();
  LimK = config->GetVenkat_LimiterCoeff();

  if (config->GetKind_SlopeLimit_Flow() == NO_LIMITER) {
   
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        nodes->SetLimiter_Primitive(iPoint,iVar, 1.0);
      }
    }
    
  }
  
  else {
    
    /*--- Initialize solution max and solution min and the limiter in the entire domain --*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        nodes->SetSolution_Max(iPoint,iVar, -EPS);
        nodes->SetSolution_Min(iPoint,iVar, EPS);
        nodes->SetLimiter_Primitive(iPoint,iVar, 2.0);
      }
    }
    
    /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      
      /*--- Get the primitive variables ---*/
      
      Primitive_i = nodes->GetPrimitive(iPoint);
      Primitive_j = nodes->GetPrimitive(jPoint);
      
      /*--- Compute the maximum, and minimum values for nodes i & j ---*/
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        du = (Primitive_j[iVar] - Primitive_i[iVar]);
        nodes->SetSolution_Min(iPoint, iVar, min(nodes->GetSolution_Min(iPoint, iVar), du));
        nodes->SetSolution_Max(iPoint, iVar, max(nodes->GetSolution_Max(iPoint, iVar), du));
        nodes->SetSolution_Min(jPoint, iVar, min(nodes->GetSolution_Min(jPoint, iVar), -du));
        nodes->SetSolution_Max(jPoint, iVar, max(nodes->GetSolution_Max(jPoint, iVar), -du));
      }
      
    }
    
    /*--- Correct the limiter values across any periodic boundaries. ---*/

    for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
      InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_PRIM_1);
      CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_PRIM_1);
    }
    
  }
  
  
  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = nodes->GetGradient_Primitive(iPoint);
      Gradient_j = nodes->GetGradient_Primitive(jPoint);
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        AD::SetPreaccIn(nodes->GetSolution_Max(iPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(iPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Max(jPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(jPoint,iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = nodes->GetSolution_Max(iPoint,iVar);
          else dp = nodes->GetSolution_Min(iPoint,iVar);
          limiter = dp/dm;
        }
        
        if (limiter < nodes->GetLimiter_Primitive(iPoint,iVar)) {
          nodes->SetLimiter_Primitive(iPoint,iVar, limiter);
          AD::SetPreaccOut(nodes->GetLimiter_Primitive(iPoint)[iVar]);
        }
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = nodes->GetSolution_Max(jPoint,iVar);
          else dp = nodes->GetSolution_Min(jPoint,iVar);
          limiter = dp/dm;
        }
        
        if (limiter < nodes->GetLimiter_Primitive(jPoint,iVar)) {
          nodes->SetLimiter_Primitive(jPoint,iVar, limiter);
          AD::SetPreaccOut(nodes->GetLimiter_Primitive(jPoint)[iVar]);
        }
        
      }
      
      AD::EndPreacc();
      
    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        y =  nodes->GetLimiter_Primitive(iPoint,iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        nodes->SetLimiter_Primitive(iPoint,iVar, limiter);
      }
    }
    
  }
  
  /*--- Venkatakrishnan limiter ---*/
  
  if ((config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) ||
      (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG)) {
    
    if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
      
      /*--- Allocate memory for the max and min primitive value --*/
      
      LocalMinPrimitive = new su2double [nPrimVarGrad]; GlobalMinPrimitive = new su2double [nPrimVarGrad];
      LocalMaxPrimitive = new su2double [nPrimVarGrad]; GlobalMaxPrimitive = new su2double [nPrimVarGrad];
      
      /*--- Compute the max value and min value of the solution ---*/
      
      Primitive = nodes->GetPrimitive(0);
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        LocalMinPrimitive[iVar] = Primitive[iVar];
        LocalMaxPrimitive[iVar] = Primitive[iVar];
      }
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Get the primitive variables ---*/
        
        Primitive = nodes->GetPrimitive(iPoint);

        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          LocalMinPrimitive[iVar] = min (LocalMinPrimitive[iVar], Primitive[iVar]);
          LocalMaxPrimitive[iVar] = max (LocalMaxPrimitive[iVar], Primitive[iVar]);
        }
        
      }

#ifdef HAVE_MPI
      SU2_MPI::Allreduce(LocalMinPrimitive, GlobalMinPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(LocalMaxPrimitive, GlobalMaxPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        GlobalMinPrimitive[iVar] = LocalMinPrimitive[iVar];
        GlobalMaxPrimitive[iVar] = LocalMaxPrimitive[iVar];
      }
#endif
    }
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = nodes->GetGradient_Primitive(iPoint);
      Gradient_j = nodes->GetGradient_Primitive(jPoint);
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          
        AD::StartPreacc();
        AD::SetPreaccIn(Gradient_i[iVar], nDim);
        AD::SetPreaccIn(Gradient_j[iVar], nDim);
        AD::SetPreaccIn(Coord_i, nDim);
        AD::SetPreaccIn(Coord_j, nDim);
        AD::SetPreaccIn(nodes->GetSolution_Max(iPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(iPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Max(jPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(jPoint,iVar));
        
        if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
          AD::SetPreaccIn(GlobalMaxPrimitive[iVar]);
          AD::SetPreaccIn(GlobalMinPrimitive[iVar]);
          eps1 = LimK * (GlobalMaxPrimitive[iVar] - GlobalMinPrimitive[iVar]);
          eps2 = eps1*eps1;
        }
        else {
          eps1 = LimK*dave;
          eps2 = eps1*eps1*eps1;
        }
        
        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = nodes->GetSolution_Max(iPoint,iVar);
        else dp = nodes->GetSolution_Min(iPoint,iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < nodes->GetLimiter_Primitive(iPoint,iVar)){
          nodes->SetLimiter_Primitive(iPoint,iVar, limiter);
          AD::SetPreaccOut(nodes->GetLimiter_Primitive(iPoint)[iVar]);
        }
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = nodes->GetSolution_Max(jPoint,iVar);
        else dp = nodes->GetSolution_Min(jPoint,iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < nodes->GetLimiter_Primitive(jPoint,iVar)){
          nodes->SetLimiter_Primitive(jPoint,iVar, limiter);
          AD::SetPreaccOut(nodes->GetLimiter_Primitive(jPoint)[iVar]);
        }
        
        AD::EndPreacc();
      }
    }
    
    if (LocalMinPrimitive  != NULL) delete [] LocalMinPrimitive;
    if (LocalMaxPrimitive  != NULL) delete [] LocalMaxPrimitive;
    if (GlobalMinPrimitive != NULL) delete [] GlobalMinPrimitive;
    if (GlobalMaxPrimitive != NULL) delete [] GlobalMaxPrimitive;

  }
  
  /*--- Correct the limiter values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_PRIM_2);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_PRIM_2);
  }
  
  /*--- Limiter MPI ---*/
  
  InitiateComms(geometry, config, PRIMITIVE_LIMITER);
  CompleteComms(geometry, config, PRIMITIVE_LIMITER);

#ifdef CODI_REVERSE_TYPE
  if (TapeActive) AD::StartRecording();
#endif
}

void CIncEulerSolver::SetInletAtVertex(su2double *val_inlet,
                                       unsigned short iMarker,
                                       unsigned long iVertex) {
  
  /*--- Alias positions within inlet file for readability ---*/
  
  unsigned short T_position       = nDim;
  unsigned short P_position       = nDim+1;
  unsigned short FlowDir_position = nDim+2;
  
  /*--- Check that the norm of the flow unit vector is actually 1 ---*/
  
  su2double norm = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    norm += pow(val_inlet[FlowDir_position + iDim], 2);
  }
  norm = sqrt(norm);
  
  /*--- The tolerance here needs to be loose.  When adding a very
   * small number (1e-10 or smaller) to a number close to 1.0, floating
   * point roundoff errors can occur. ---*/
  
  if (abs(norm - 1.0) > 1e-6) {
    ostringstream error_msg;
    error_msg << "ERROR: Found these values in columns ";
    error_msg << FlowDir_position << " - ";
    error_msg << FlowDir_position + nDim - 1 << endl;
    error_msg << std::scientific;
    error_msg << "  [" << val_inlet[FlowDir_position];
    error_msg << ", " << val_inlet[FlowDir_position + 1];
    if (nDim == 3) error_msg << ", " << val_inlet[FlowDir_position + 2];
    error_msg << "]" << endl;
    error_msg << "  These values should be components of a unit vector for direction," << endl;
    error_msg << "  but their magnitude is: " << norm << endl;
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  }
  
  /*--- Store the values in our inlet data structures. ---*/
  
  Inlet_Ttotal[iMarker][iVertex] = val_inlet[T_position];
  Inlet_Ptotal[iMarker][iVertex] = val_inlet[P_position];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Inlet_FlowDir[iMarker][iVertex][iDim] =  val_inlet[FlowDir_position + iDim];
  }
  
}

su2double CIncEulerSolver::GetInletAtVertex(su2double *val_inlet,
                                            unsigned long val_inlet_point,
                                            unsigned short val_kind_marker,
                                            string val_marker,
                                            CGeometry *geometry,
                                            CConfig *config) {
  
  /*--- Local variables ---*/
  
  unsigned short iMarker, iDim;
  unsigned long iPoint, iVertex;
  su2double Area = 0.0;
  su2double Normal[3] = {0.0,0.0,0.0};
  
  /*--- Alias positions within inlet file for readability ---*/
  
    unsigned short T_position       = nDim;
    unsigned short P_position       = nDim+1;
    unsigned short FlowDir_position = nDim+2;
  
  if (val_kind_marker == INLET_FLOW) {
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++){
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (iPoint == val_inlet_point) {
            
            /*-- Compute boundary face area for this vertex. ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);
            
            /*--- Access and store the inlet variables for this vertex. ---*/
            
            val_inlet[T_position] = Inlet_Ttotal[iMarker][iVertex];
            val_inlet[P_position] = Inlet_Ptotal[iMarker][iVertex];
            for (iDim = 0; iDim < nDim; iDim++) {
              val_inlet[FlowDir_position + iDim] = Inlet_FlowDir[iMarker][iVertex][iDim];
            }
            
            /*--- Exit once we find the point. ---*/
            
            return Area;
            
          }
        }
      }
    }
  }
  
  /*--- If we don't find a match, then the child point is not on the
   current inlet boundary marker. Return zero area so this point does
   not contribute to the restriction operator and continue. ---*/
  
  return Area;
  
}

void CIncEulerSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {
  
  if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
    
    string Marker_Tag   = config->GetMarker_All_TagBound(iMarker);
    su2double p_total   = config->GetInlet_Ptotal(Marker_Tag);
    su2double t_total   = config->GetInlet_Ttotal(Marker_Tag);
    su2double* flow_dir = config->GetInlet_FlowDir(Marker_Tag);
    
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_Ttotal[iMarker][iVertex] = t_total;
      Inlet_Ptotal[iMarker][iVertex] = p_total;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Inlet_FlowDir[iMarker][iVertex][iDim] = flow_dir[iDim];
    }
    
  } else {
    
    /*--- For now, non-inlets just get set to zero. In the future, we
     can do more customization for other boundary types here. ---*/
    
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_Ttotal[iMarker][iVertex] = 0.0;
      Inlet_Ptotal[iMarker][iVertex] = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Inlet_FlowDir[iMarker][iVertex][iDim] = 0.0;
    }
  }
  
}

void CIncEulerSolver::Evaluate_ObjFunc(CConfig *config) {

  unsigned short iMarker_Monitoring, Kind_ObjFunc;
  su2double Weight_ObjFunc;

  Total_ComboObj = 0.0;
  
  /*--- Loop over all monitored markers, add to the 'combo' objective ---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {

    Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
    Kind_ObjFunc = config->GetKind_ObjFunc(iMarker_Monitoring);
   
    switch(Kind_ObjFunc) {
      case DRAG_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CD[iMarker_Monitoring]);
        break;
      case LIFT_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CL[iMarker_Monitoring]);
        break;
      case SIDEFORCE_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CSF[iMarker_Monitoring]);
        break;
      case EFFICIENCY:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CEff[iMarker_Monitoring]);
        break;
      case MOMENT_X_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMx[iMarker_Monitoring]);
        break;
      case MOMENT_Y_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMy[iMarker_Monitoring]);
        break;
      case MOMENT_Z_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMz[iMarker_Monitoring]);
        break;
      case FORCE_X_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Surface_CFx[iMarker_Monitoring];
        break;
      case FORCE_Y_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Surface_CFy[iMarker_Monitoring];
        break;
      case FORCE_Z_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*Surface_CFz[iMarker_Monitoring];
        break;
      case TOTAL_HEATFLUX:
        Total_ComboObj+=Weight_ObjFunc*Surface_HF_Visc[iMarker_Monitoring];
        break;
      case MAXIMUM_HEATFLUX:
        Total_ComboObj+=Weight_ObjFunc*Surface_MaxHF_Visc[iMarker_Monitoring];
        break;
      default:
        break;

    }
  }
  
  /*--- The following are not per-surface, and so to avoid that they are
   double-counted when multiple surfaces are specified, they have been
   placed outside of the loop above. In addition, multi-objective mode is
   also disabled for these objective functions (error thrown at start). ---*/
  
  Weight_ObjFunc = config->GetWeight_ObjFunc(0);
  Kind_ObjFunc   = config->GetKind_ObjFunc(0);
  
  switch(Kind_ObjFunc) {
    case INVERSE_DESIGN_PRESSURE:
      Total_ComboObj+=Weight_ObjFunc*Total_CpDiff;
      break;
    case INVERSE_DESIGN_HEATFLUX:
      Total_ComboObj+=Weight_ObjFunc*Total_HeatFluxDiff;
      break;
    case THRUST_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Total_CT;
      break;
    case TORQUE_COEFFICIENT:
      Total_ComboObj+=Weight_ObjFunc*Total_CQ;
      break;
    case FIGURE_OF_MERIT:
      Total_ComboObj+=Weight_ObjFunc*Total_CMerit;
      break;
    case SURFACE_TOTAL_PRESSURE:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_TotalPressure(0);
      break;
    case SURFACE_STATIC_PRESSURE:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_Pressure(0);
      break;
    case SURFACE_MASSFLOW:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_MassFlow(0);
      break;
    case SURFACE_UNIFORMITY:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_Uniformity(0);
      break;
    case SURFACE_SECONDARY:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_SecondaryStrength(0);
      break;
    case SURFACE_MOM_DISTORTION:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_MomentumDistortion(0);
      break;
    case SURFACE_SECOND_OVER_UNIFORM:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_SecondOverUniform(0);
      break;
    case SURFACE_PRESSURE_DROP:
      Total_ComboObj+=Weight_ObjFunc*config->GetSurface_PressureDrop(0);
      break;
    case CUSTOM_OBJFUNC:
      Total_ComboObj+=Weight_ObjFunc*Total_Custom_ObjFunc;
      break;
    default:
      break;
  }
  
}

void CIncEulerSolver::SetBeta_Parameter(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh) {
  
  su2double epsilon2  = config->GetBeta_Factor();
  su2double epsilon2_default = 4.1;
  su2double maxVel2 = 0.0;
  su2double Beta = 1.0;

  unsigned long iPoint;

  /*--- For now, only the finest mesh level stores the Beta for all levels. ---*/
  
  if (iMesh == MESH_0) {
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      
      /*--- Store the local maximum of the squared velocity in the field. ---*/
      
      if (nodes->GetVelocity2(iPoint) > maxVel2)
        maxVel2 = nodes->GetVelocity2(iPoint);
      
    }
    
    /*--- Communicate the max globally to give a conservative estimate. ---*/
    
#ifdef HAVE_MPI
    su2double myMaxVel2 = maxVel2; maxVel2 = 0.0;
    SU2_MPI::Allreduce(&myMaxVel2, &maxVel2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    Beta = max(1e-10,maxVel2);
    config->SetMax_Vel2(Beta);
    
  }

  /*--- Allow an override if user supplies a large epsilon^2. ---*/

  epsilon2 = max(epsilon2_default,epsilon2);

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    nodes->SetBetaInc2(iPoint,epsilon2*config->GetMax_Vel2());

}

void CIncEulerSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) {

  unsigned short iDim, jDim;

  su2double  BetaInc2, Density, dRhodT, Temperature, oneOverCp, Cp;
  su2double  Velocity[3] = {0.0,0.0,0.0};

  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool energy           = config->GetEnergy_Equation();

  /*--- Access the primitive variables at this node. ---*/

  Density     = nodes->GetDensity(iPoint);
  BetaInc2    = nodes->GetBetaInc2(iPoint);
  Cp          = nodes->GetSpecificHeatCp(iPoint);
  oneOverCp   = 1.0/Cp;
  Temperature = nodes->GetTemperature(iPoint);

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);

  /*--- We need the derivative of the equation of state to build the
   preconditioning matrix. For now, the only option is the ideal gas
   law, but in the future, dRhodT should be in the fluid model. ---*/

  if (variable_density) {
    dRhodT = -Density/Temperature;
  } else {
    dRhodT = 0.0;
  }

  /*--- Calculating the inverse of the preconditioning matrix
   that multiplies the time derivative during time integration. ---*/

  if (implicit) {

    /*--- For implicit calculations, we multiply the preconditioner
     by the cell volume over the time step and add to the Jac diagonal. ---*/

    Preconditioner[0][0] = 1.0/BetaInc2;
    for (iDim = 0; iDim < nDim; iDim++)
      Preconditioner[iDim+1][0] = Velocity[iDim]/BetaInc2;

    if (energy) Preconditioner[nDim+1][0] = Cp*Temperature/BetaInc2;
    else        Preconditioner[nDim+1][0] = 0.0;

    for (jDim = 0; jDim < nDim; jDim++) {
      Preconditioner[0][jDim+1] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        if (iDim == jDim) Preconditioner[iDim+1][jDim+1] = Density;
        else Preconditioner[iDim+1][jDim+1] = 0.0;
      }
      Preconditioner[nDim+1][jDim+1] = 0.0;
    }

    Preconditioner[0][nDim+1] = dRhodT;
    for (iDim = 0; iDim < nDim; iDim++)
      Preconditioner[iDim+1][nDim+1] = Velocity[iDim]*dRhodT;

    if (energy) Preconditioner[nDim+1][nDim+1] = Cp*(dRhodT*Temperature + Density);
    else        Preconditioner[nDim+1][nDim+1] = 1.0;

  } else {

    /*--- For explicit calculations, we move the residual to the
     right-hand side and pre-multiply by the preconditioner inverse.
     Therefore, we build inv(Precon) here and multiply by the residual
     later in the R-K and Euler Explicit time integration schemes. ---*/

    Preconditioner[0][0] = Temperature*BetaInc2*dRhodT/Density + BetaInc2;
    for (iDim = 0; iDim < nDim; iDim ++)
      Preconditioner[iDim+1][0] = -1.0*Velocity[iDim]/Density;

    if (energy) Preconditioner[nDim+1][0] = -1.0*Temperature/Density;
    else        Preconditioner[nDim+1][0] = 0.0;


    for (jDim = 0; jDim < nDim; jDim++) {
      Preconditioner[0][jDim+1] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        if (iDim == jDim) Preconditioner[iDim+1][jDim+1] = 1.0/Density;
        else Preconditioner[iDim+1][jDim+1] = 0.0;
      }
      Preconditioner[nDim+1][jDim+1] = 0.0;
    }

    Preconditioner[0][nDim+1] = -1.0*BetaInc2*dRhodT*oneOverCp/Density;
    for (iDim = 0; iDim < nDim; iDim ++)
      Preconditioner[iDim+1][nDim+1] = 0.0;

    if (energy) Preconditioner[nDim+1][nDim+1] = oneOverCp/Density;
    else        Preconditioner[nDim+1][nDim+1] = 0.0;
    
  }
  
}

void CIncEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  
  su2double *V_infty, *V_domain;
  
  bool implicit      = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool viscous       = config->GetViscous();
  
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
      
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Recompute and store the velocity in the primitive variable vector. ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = GetVelocity_Inf(iDim);

      /*--- Far-field pressure set to static pressure (0.0). ---*/

      V_infty[0] = GetPressure_Inf();

      /*--- Dirichlet condition for temperature at far-field (if energy is active). ---*/

      V_infty[nDim+1] = GetTemperature_Inf();

      /*--- Store the density.  ---*/

      V_infty[nDim+2] = GetDensity_Inf();

      /*--- Beta coefficient stored at the node ---*/

      V_infty[nDim+3] = nodes->GetBetaInc2(iPoint);

      /*--- Cp is needed for Temperature equation. ---*/

      V_infty[nDim+7] = nodes->GetSpecificHeatCp(iPoint);

      /*--- Set various quantities in the numerics class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Convective Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous residual contribution ---*/
      
      if (viscous) {
        
        /*--- Set transport properties at infinity. ---*/

        V_infty[nDim+4] = nodes->GetLaminarViscosity(iPoint);
        V_infty[nDim+5] = nodes->GetEddyViscosity(iPoint);
        V_infty[nDim+6] = nodes->GetThermalConductivity(iPoint);

        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));
        
        /*--- Turbulent kinetic energy ---*/
        
        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
        
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

void CIncEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  unsigned long Point_Normal;
  su2double *Flow_Dir, Flow_Dir_Mag, Vel_Mag, Area, P_total, P_domain, Vn;
  su2double *V_inlet, *V_domain;
  su2double UnitFlowDir[3] = {0.0,0.0,0.0};
  su2double dV[3] = {0.0,0.0,0.0};
  su2double Damping = config->GetInc_Inlet_Damping();

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous       = config->GetViscous();

  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);

  unsigned short Kind_Inlet = config->GetKind_Inc_Inlet(Marker_Tag);

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
    
      /*--- Both types of inlets may use the prescribed flow direction.
       Ensure that the flow direction is a unit vector. ---*/
      
      Flow_Dir = Inlet_FlowDir[val_marker][iVertex];
      Flow_Dir_Mag = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir_Mag += Flow_Dir[iDim]*Flow_Dir[iDim];
      Flow_Dir_Mag = sqrt(Flow_Dir_Mag);
      
      /*--- Store the unit flow direction vector. ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitFlowDir[iDim] = Flow_Dir[iDim]/Flow_Dir_Mag;

      /*--- Retrieve solution at this boundary node. ---*/
      
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Neumann condition for dynamic pressure ---*/
      
      V_inlet[0] = nodes->GetPressure(iPoint);
      
      /*--- The velocity is either prescribed or computed from total pressure. ---*/

      switch (Kind_Inlet) {
          
          /*--- Velocity and temperature (if required) been specified at the inlet. ---*/
          
        case VELOCITY_INLET:
          
          /*--- Retrieve the specified velocity and temperature for the inlet. ---*/
          
          Vel_Mag  = Inlet_Ptotal[val_marker][iVertex]/config->GetVelocity_Ref();
          
          /*--- Store the velocity in the primitive variable vector. ---*/
          
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*UnitFlowDir[iDim];
          
          /*--- Dirichlet condition for temperature (if energy is active) ---*/
          
          V_inlet[nDim+1] = Inlet_Ttotal[val_marker][iVertex]/config->GetTemperature_Ref();
          
          break;
          
          /*--- Stagnation pressure has been specified at the inlet. ---*/
          
        case PRESSURE_INLET:
          
          /*--- Retrieve the specified total pressure for the inlet. ---*/
          
          P_total = Inlet_Ptotal[val_marker][iVertex]/config->GetPressure_Ref();
          
          /*--- Store the current static pressure for clarity. ---*/
          
          P_domain = nodes->GetPressure(iPoint);
          
          /*--- Check for back flow through the inlet. ---*/
          
          Vn = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Vn += V_domain[iDim+1]*(-1.0*Normal[iDim]/Area);
          }
          
          /*--- If the local static pressure is larger than the specified
           total pressure or the velocity is directed upstream, we have a
           back flow situation. The specified total pressure should be used
           as a static pressure condition and the velocity from the domain
           is used for the BC. ---*/
          
          if ((P_domain > P_total) || (Vn < 0.0)) {
            
            /*--- Back flow: use the prescribed P_total as static pressure. ---*/
            
            V_inlet[0] = Inlet_Ptotal[val_marker][iVertex]/config->GetPressure_Ref();
            
            /*--- Neumann condition for velocity. ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              V_inlet[iDim+1] = V_domain[iDim+1];
            
            /*--- Neumann condition for the temperature. ---*/
            
            V_inlet[nDim+1] = nodes->GetTemperature(iPoint);
            
          } else {
            
            /*--- Update the velocity magnitude using the total pressure. ---*/
            
            Vel_Mag = sqrt((P_total - P_domain)/(0.5*nodes->GetDensity(iPoint)));
            
            /*--- If requested, use the local boundary normal (negative),
             instead of the prescribed flow direction in the config. ---*/
            
            if (config->GetInc_Inlet_UseNormal()) {
              for (iDim = 0; iDim < nDim; iDim++)
                UnitFlowDir[iDim] = -Normal[iDim]/Area;
            }
            
            /*--- Compute the delta change in velocity in each direction. ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              dV[iDim] = Vel_Mag*UnitFlowDir[iDim] - V_domain[iDim+1];
            
            /*--- Update the velocity in the primitive variable vector.
             Note we use damping here to improve stability/convergence. ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              V_inlet[iDim+1] = V_domain[iDim+1] + Damping*dV[iDim];
            
            /*--- Dirichlet condition for temperature (if energy is active) ---*/
            
            V_inlet[nDim+1] = Inlet_Ttotal[val_marker][iVertex]/config->GetTemperature_Ref();
            
          }
          
          break;
          
      }

      /*--- Access density at the node. This is either constant by
        construction, or will be set fixed implicitly by the temperature
        and equation of state. ---*/

      V_inlet[nDim+2] = nodes->GetDensity(iPoint);

      /*--- Beta coefficient from the config file ---*/

      V_inlet[nDim+3] = nodes->GetBetaInc2(iPoint);

      /*--- Cp is needed for Temperature equation. ---*/

      V_inlet[nDim+7] = nodes->GetSpecificHeatCp(iPoint);

      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution, commented out because serious convergence problems ---*/

      if (viscous) {
        
        /*--- Set transport properties at the inlet ---*/
        
        V_inlet[nDim+4] = nodes->GetLaminarViscosity(iPoint);
        V_inlet[nDim+5] = nodes->GetEddyViscosity(iPoint);
        V_inlet[nDim+6] = nodes->GetThermalConductivity(iPoint);

        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));
        
        /*--- Turbulent kinetic energy ---*/
        
        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
        
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

void CIncEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Area;
  su2double *V_outlet, *V_domain, P_Outlet = 0.0, P_domain;
  su2double mDot_Target, mDot_Old, dP, Density_Avg, Area_Outlet;
  su2double Damping = config->GetInc_Outlet_Damping();

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous       = config->GetViscous();
  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);

  su2double *Normal = new su2double[nDim];
  
  unsigned short Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);
  
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
       
      /*--- Current solution at this boundary node ---*/
      
      V_domain = nodes->GetPrimitive(iPoint);
      
      /*--- Store the current static pressure for clarity. ---*/
      
      P_domain = nodes->GetPressure(iPoint);
      
      /*--- Compute a boundary value for the pressure depending on whether
       we are prescribing a back pressure or a mass flow target. ---*/
      
      switch (Kind_Outlet) {
          
          /*--- Velocity and temperature (if required) been specified at the inlet. ---*/
          
        case PRESSURE_OUTLET:
          
          /*--- Retrieve the specified back pressure for this outlet. ---*/
          
          P_Outlet = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();
          
          /*--- The pressure is prescribed at the outlet. ---*/
          
          V_outlet[0] = P_Outlet;
          
          /*--- Neumann condition for the velocity. ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            V_outlet[iDim+1] = nodes->GetVelocity(iPoint,iDim);
          }
          
          break;
          
          /*--- A mass flow target has been specified for the outlet. ---*/
          
        case MASS_FLOW_OUTLET:
          
          /*--- Retrieve the specified target mass flow at the outlet. ---*/
          
          mDot_Target = config->GetOutlet_Pressure(Marker_Tag)/(config->GetDensity_Ref() * config->GetVelocity_Ref());

          /*--- Retrieve the old mass flow, density, and area of the outlet,
           which has been computed in a preprocessing step. These values
           were stored in non-dim. form in the config container. ---*/
          
          mDot_Old    = config->GetOutlet_MassFlow(Marker_Tag);
          Density_Avg = config->GetOutlet_Density(Marker_Tag);
          Area_Outlet = config->GetOutlet_Area(Marker_Tag);

          /*--- Compute the pressure increment based on the difference
           between the current and target mass flow. Note that increasing
           pressure decreases flow speed. ---*/
          
          dP = 0.5*Density_Avg*(mDot_Old*mDot_Old - mDot_Target*mDot_Target)/((Density_Avg*Area_Outlet)*(Density_Avg*Area_Outlet));
          
          /*--- Update the new outlet pressure. Note that we use damping
           here to improve stability/convergence. ---*/
          
          P_Outlet = P_domain + Damping*dP;

          /*--- The pressure is prescribed at the outlet. ---*/
          
          V_outlet[0] = P_Outlet;
          
          /*--- Neumann condition for the velocity ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            V_outlet[iDim+1] = nodes->GetVelocity(iPoint,iDim);
          }
          
          break;
          
      }
      
      /*--- Neumann condition for the temperature. ---*/

      V_outlet[nDim+1] = nodes->GetTemperature(iPoint);

      /*--- Access density at the interior node. This is either constant by
        construction, or will be set fixed implicitly by the temperature
        and equation of state. ---*/
      
      V_outlet[nDim+2] = nodes->GetDensity(iPoint);

      /*--- Beta coefficient from the config file ---*/
      
      V_outlet[nDim+3] = nodes->GetBetaInc2(iPoint);

      /*--- Cp is needed for Temperature equation. ---*/

      V_outlet[nDim+7] = nodes->GetSpecificHeatCp(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      
      /*--- Viscous contribution, commented out because serious convergence problems ---*/

      if (viscous) {

        /*--- Set transport properties at the outlet. ---*/

        V_outlet[nDim+4] = nodes->GetLaminarViscosity(iPoint);
        V_outlet[nDim+5] = nodes->GetEddyViscosity(iPoint);
        V_outlet[nDim+6] = nodes->GetThermalConductivity(iPoint);

        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));
        
        /*--- Turbulent kinetic energy ---*/
        
        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
        
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


void CIncEulerSolver::BC_Euler_Wall(CGeometry      *geometry,
                                    CSolver        **solver_container,
                                    CNumerics      *conv_numerics,
                                    CNumerics      *visc_numerics,
                                    CConfig        *config,
                                    unsigned short val_marker) {

  /*--- Call the equivalent symmetry plane boundary condition. ---*/
  BC_Sym_Plane(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}


void CIncEulerSolver::BC_Sym_Plane(CGeometry      *geometry,
                                   CSolver        **solver_container,
                                   CNumerics      *conv_numerics,
                                   CNumerics      *visc_numerics,
                                   CConfig        *config,
                                   unsigned short val_marker) {
  
  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT),
       viscous  = config->GetViscous();
  
  /*--- Allocation of variables necessary for convective fluxes. ---*/
  su2double Area, ProjVelocity_i,
            *V_reflected, 
            *V_domain,
            *Normal     = new su2double[nDim],
            *UnitNormal = new su2double[nDim];

  /*--- Allocation of variables necessary for viscous fluxes. ---*/
  su2double ProjGradient, ProjNormVelGrad, ProjTangVelGrad, TangentialNorm,
            *Tangential  = new su2double[nDim],
            *GradNormVel = new su2double[nDim],
            *GradTangVel = new su2double[nDim];

  /*--- Allocation of primitive gradient arrays for viscous fluxes. ---*/
  su2double **Grad_Reflected = new su2double*[nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Grad_Reflected[iVar] = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker. ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    if (iVertex == 0 || 
        geometry->bound_is_straight[val_marker] != true) {

      /*----------------------------------------------------------------------------------------------*/
      /*--- Preprocessing:                                                                         ---*/
      /*--- Compute the unit normal and (in case of viscous flow) a corresponding unit tangential  ---*/
      /*--- to that normal. On a straight(2D)/plane(3D) boundary these two vectors are constant.   ---*/
      /*--- This circumstance is checked in gemoetry->ComputeSurf_Straightness(...) and stored     ---*/
      /*--- such that the recomputation does not occur for each node. On true symmetry planes, the ---*/
      /*--- normal is constant but this routines is used for Symmetry, Euler-Wall in inviscid flow ---*/
      /*--- and Euler Wall in viscous flow as well. In the latter curvy boundaries are likely to   ---*/
      /*--- happen. In doubt, the conditional above which checks straightness can be thrown out    ---*/
      /*--- such that the recomputation is done for each node (which comes with a tiny performance ---*/
      /*--- penalty).                                                                              ---*/
      /*----------------------------------------------------------------------------------------------*/

      /*--- Normal vector for a random vertex (zero) on this marker (negate for outward convention). ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal); 
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];

      /*--- Compute unit normal, to be used for unit tangential, projected velocity and velocity 
            component gradients. ---*/
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Preprocessing: Compute unit tangential, the direction is arbitrary as long as 
            t*n=0 && |t|_2 = 1 ---*/
      if (viscous) {
        switch( nDim ) {
          case 2: {
            Tangential[0] = -UnitNormal[1];
            Tangential[1] =  UnitNormal[0];
            break;
          }
          case 3: {
            /*--- n = ai + bj + ck, if |b| > |c| ---*/
            if( abs(UnitNormal[1]) > abs(UnitNormal[2])) {
              /*--- t = bi + (c-a)j - bk  ---*/
              Tangential[0] = UnitNormal[1];
              Tangential[1] = UnitNormal[2] - UnitNormal[0];
              Tangential[2] = -UnitNormal[1];
            } else {
              /*--- t = ci - cj + (b-a)k  ---*/
              Tangential[0] = UnitNormal[2];
              Tangential[1] = -UnitNormal[2];
              Tangential[2] = UnitNormal[1] - UnitNormal[0];
            }
            /*--- Make it a unit vector. ---*/
            TangentialNorm = sqrt(pow(Tangential[0],2) + pow(Tangential[1],2) + pow(Tangential[2],2));
            Tangential[0] = Tangential[0] / TangentialNorm;
            Tangential[1] = Tangential[1] / TangentialNorm;
            Tangential[2] = Tangential[2] / TangentialNorm;
            break;
          }
        }// switch
      }//if viscous
    }//if bound_is_straight
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*-------------------------------------------------------------------------------*/
      /*--- Step 1: For the convective fluxes, create a reflected state of the      ---*/
      /*---         Primitive variables by copying all interior values to the       ---*/
      /*---         reflected. Only the velocity is mirrored along the symmetry     ---*/
      /*---         axis. Based on the Upwind_Residual routine.                     ---*/
      /*-------------------------------------------------------------------------------*/

      /*--- Allocate the reflected state at the symmetry boundary. ---*/
      V_reflected = GetCharacPrimVar(val_marker, iVertex);

      /*--- Grid movement ---*/
      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Normal vector for this vertex (negate for outward convention). ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Get current solution at this boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);
      
      /*--- Set the reflected state based on the boundary node. Scalars are copied and 
            the velocity is mirrored along the symmetry boundary, i.e. the velocity in 
            normal direction is substracted twice. ---*/
      for(iVar = 0; iVar < nPrimVar; iVar++)
        V_reflected[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Compute velocity in normal direction (ProjVelcity_i=(v*n)) und substract twice from
            velocity in normal direction: v_r = v - 2 (v*n)n ---*/
      ProjVelocity_i = nodes->GetProjVel(iPoint,UnitNormal);
      
      for (iDim = 0; iDim < nDim; iDim++)
        V_reflected[iDim+1] = nodes->GetVelocity(iPoint,iDim) - 2.0 * ProjVelocity_i*UnitNormal[iDim];
      
      /*--- Set Primitive and Secondary for numerics class. ---*/
      conv_numerics->SetPrimitive(V_domain, V_reflected);
      conv_numerics->SetSecondary(nodes->GetSecondary(iPoint), nodes->GetSecondary(iPoint));

      /*--- Compute the residual using an upwind scheme. ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration. ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      
      if (viscous) {

        /*-------------------------------------------------------------------------------*/
        /*--- Step 2: The viscous fluxes of the Navier-Stokes equations depend on the ---*/
        /*---         Primitive variables and their gradients. The viscous numerics   ---*/
        /*---         container is filled just as the convective numerics container,  ---*/
        /*---         but the primitive gradients of the reflected state have to be   ---*/
        /*---         determined additionally such that symmetry at the boundary is   ---*/
        /*---         enforced. Based on the Viscous_Residual routine.                ---*/
        /*-------------------------------------------------------------------------------*/

        /*--- Set the normal vector and the coordinates. ---*/
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), 
                                geometry->node[iPoint]->GetCoord());
        visc_numerics->SetNormal(Normal);
        
        /*--- Set the primitive and Secondary variables. ---*/
        visc_numerics->SetPrimitive(V_domain, V_reflected);
        visc_numerics->SetSecondary(nodes->GetSecondary(iPoint), nodes->GetSecondary(iPoint));
        
        /*--- For viscous Fluxes also the gradients of the primitives need to be determined.
              1. The gradients of scalars are mirrored along the sym plane just as velocity for the primitives
              2. The gradients of the velocity components need more attention, i.e. the gradient of the
                 normal velocity in tangential direction is mirrored and the gradient of the tangential velocity in 
                 normal direction is mirrored. ---*/

        /*--- Get gradients of primitives of boundary cell ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Grad_Reflected[iVar][iDim] = nodes->GetGradient_Primitive(iPoint,iVar, iDim);
        
        /*--- Reflect the gradients for all scalars including the velocity components.
              The gradients of the velocity components are set later with the 
              correct values: grad(V)_r = grad(V) - 2 [grad(V)*n]n, V beeing any primitive ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          if(iVar == 0 || iVar > nDim) { // Exclude velocity component gradients

            /*--- Compute projected part of the gradient in a dot product ---*/
            ProjGradient = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              ProjGradient += Grad_Reflected[iVar][iDim]*UnitNormal[iDim];
              
            for (iDim = 0; iDim < nDim; iDim++)           
              Grad_Reflected[iVar][iDim] = Grad_Reflected[iVar][iDim] - 2.0 * ProjGradient*UnitNormal[iDim];
          }
        }
        
        /*--- Compute gradients of normal and tangential velocity:
              grad(v*n) = grad(v_x) n_x + grad(v_y) n_y (+ grad(v_z) n_z)
              grad(v*t) = grad(v_x) t_x + grad(v_y) t_y (+ grad(v_z) t_z) ---*/
        for (iVar = 0; iVar < nDim; iVar++) { // counts gradient components
          GradNormVel[iVar] = 0.0;
          GradTangVel[iVar] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) { // counts sum with unit normal/tangential
            GradNormVel[iVar] += Grad_Reflected[iDim+1][iVar] * UnitNormal[iDim];
            GradTangVel[iVar] += Grad_Reflected[iDim+1][iVar] * Tangential[iDim];
          }
        }

        /*--- Refelect gradients in tangential and normal direction by substracting the normal/tangential
              component twice, just as done with velocity above.
              grad(v*n)_r = grad(v*n) - 2 {grad([v*n])*t}t
              grad(v*t)_r = grad(v*t) - 2 {grad([v*t])*n}n ---*/
        ProjNormVelGrad = 0.0;
        ProjTangVelGrad = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjNormVelGrad += GradNormVel[iDim]*Tangential[iDim]; //grad([v*n])*t
          ProjTangVelGrad += GradTangVel[iDim]*UnitNormal[iDim]; //grad([v*t])*n
        }
        
        for (iDim = 0; iDim < nDim; iDim++) {
          GradNormVel[iDim] = GradNormVel[iDim] - 2.0 * ProjNormVelGrad * Tangential[iDim];
          GradTangVel[iDim] = GradTangVel[iDim] - 2.0 * ProjTangVelGrad * UnitNormal[iDim];
        }
        
        /*--- Transfer reflected gradients back into the Cartesian Coordinate system:
              grad(v_x)_r = grad(v*n)_r n_x + grad(v*t)_r t_x
              grad(v_y)_r = grad(v*n)_r n_y + grad(v*t)_r t_y
              ( grad(v_z)_r = grad(v*n)_r n_z + grad(v*t)_r t_z ) ---*/
        for (iVar = 0; iVar < nDim; iVar++) // loops over the velocity component gradients
          for (iDim = 0; iDim < nDim; iDim++) // loops over the entries of the above
            Grad_Reflected[iVar+1][iDim] = GradNormVel[iDim]*UnitNormal[iVar] + GradTangVel[iDim]*Tangential[iVar];
        
        /*--- Set the primitive gradients of the boundary and reflected state. ---*/
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), Grad_Reflected);
        
        /*--- Turbulent kinetic energy. ---*/
        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
        
        /*--- Compute and update residual. Note that the viscous shear stress tensor is computed in the 
              following routine based upon the velocity-component gradients. ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration. ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }//if viscous
    }//if GetDomain
  }//for iVertex
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] UnitNormal;
  delete [] Tangential;
  delete [] GradNormVel;
  delete [] GradTangVel;

  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Grad_Reflected[iVar];
  delete [] Grad_Reflected;
}


void CIncEulerSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config) {
  
  unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
  unsigned short iDim, iVar, iMarker, nDonorVertex;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous       = config->GetViscous();
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  su2double *tmp_residual = new su2double[nVar];
  
  su2double weight;
   
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          nDonorVertex = GetnSlidingStates(iMarker, iVertex);
          
          /*--- Initialize Residual, this will serve to accumulate the average ---*/

          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = 0.0;

          /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/

          for (jVertex = 0; jVertex < nDonorVertex; jVertex++){

            Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
              PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);
            }
            
            /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/

            weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

            /*--- Set primitive variables ---*/

            conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );
          
            /*--- Set the normal vector ---*/
 
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) 
              Normal[iDim] = -Normal[iDim];

            conv_numerics->SetNormal(Normal);

            if (dynamic_grid)
              conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the convective residual using an upwind scheme ---*/

            conv_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);

            /*--- Accumulate the residuals to compute the average ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] += weight*tmp_residual[iVar];

          }

          /*--- Add Residuals and Jacobians ---*/
  
          LinSysRes.AddBlock(iPoint, Residual);
          if (implicit) 
            Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          if (viscous) {
            
            /*--- Initialize Residual, this will serve to accumulate the average ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              Residual[iVar] = 0.0;
              
            /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/
            
            for (jVertex = 0; jVertex < nDonorVertex; jVertex++){
              PrimVar_j[nDim+5] = GetSlidingState(iMarker, iVertex, nDim+5, jVertex); 
              PrimVar_j[nDim+6] = GetSlidingState(iMarker, iVertex, nDim+6, jVertex); 

              /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/
              
              weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);
              
              /*--- Set the normal vector and the coordinates ---*/

              visc_numerics->SetNormal(Normal);
              visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

              /*--- Primitive variables, and gradient ---*/

              visc_numerics->SetPrimitive(PrimVar_i, PrimVar_j);
              visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));

              /*--- Turbulent kinetic energy ---*/

              if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
                visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0), solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

              /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
              
              visc_numerics->SetTauWall(-1.0, -1.0);

              /*--- Compute and update residual ---*/

              visc_numerics->ComputeResidual(tmp_residual, Jacobian_i, Jacobian_j, config);
              
              /*--- Accumulate the residuals to compute the average ---*/
              
              for (iVar = 0; iVar < nVar; iVar++)
                Residual[iVar] += weight*tmp_residual[iVar];
            }
          
            LinSysRes.SubtractBlock(iPoint, Residual);
 
            /*--- Jacobian contribution for implicit integration ---*/

            if (implicit)
              Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
            
          }
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] tmp_residual;
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *numerics, CConfig *config) {
  
  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated correctly during the communications. For implicit calculations,
   the Jacobians and linear system are also correctly adjusted here. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }
  
}

void CIncEulerSolver::BC_Custom(CGeometry      *geometry,
                                CSolver        **solver_container,
                                CNumerics      *conv_numerics,
                                CNumerics      *visc_numerics,
                                CConfig        *config,
                                unsigned short val_marker) {
  
  /* Check for a verification solution. */
  
  if (VerificationSolution) {
    
    unsigned short iVar;
    unsigned long iVertex, iPoint, total_index;
    
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    
    /*--- Get the physical time. ---*/
    
    su2double time = 0.0;
    if (config->GetTime_Marching()) time = config->GetPhysicalTime();
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
      
      /*--- Get the point index for the current node. ---*/
      
      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
      
      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
      
      if (geometry->node[iPoint]->GetDomain()) {
        
        /*--- Get the coordinates for the current node. ---*/
        
        const su2double *coor = geometry->node[iPoint]->GetCoord();
        
        /*--- Get the conservative state from the verification solution. ---*/
        
        VerificationSolution->GetBCState(coor, time, Solution);
        
        /*--- For verification cases, we will apply a strong Dirichlet
         condition by setting the solution values at the boundary nodes
         directly and setting the residual to zero at those nodes. ---*/
        
        nodes->SetSolution_Old(iPoint,Solution);
        nodes->SetSolution(iPoint,Solution);
        nodes->SetRes_TruncErrorZero(iPoint);
        LinSysRes.SetBlock_Zero(iPoint);
        
        /*--- Adjust rows of the Jacobian (includes 1 in the diagonal) ---*/
        
        if (implicit){
          for (iVar = 0; iVar < nVar; iVar++) {
            total_index = iPoint*nVar+iVar;
            Jacobian.DeleteValsRowi(total_index);
          }
        }
        
      }
    }
    
  } else {
    
    /* The user must specify the custom BC's here. */
    SU2_MPI::Error("Implement customized boundary conditions here.", CURRENT_FUNCTION);
    
  }
  
}

void CIncEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double Density, Cp;
  su2double *V_time_nM1, *V_time_n, *V_time_nP1;
  su2double U_time_nM1[5], U_time_n[5], U_time_nP1[5];
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool energy   = config->GetEnergy_Equation();
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!dynamic_grid) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Initialize the Residual / Jacobian container to zero. ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
        if (implicit) {
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
        }
      }
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. These are actually
       the primitive values, but we will convert to conservatives. ---*/
      
      V_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      V_time_n   = nodes->GetSolution_time_n(iPoint);
      V_time_nP1 = nodes->GetSolution(iPoint);
      
      /*--- Access the density and Cp at this node (constant for now). ---*/
      
      Density     = nodes->GetDensity(iPoint);
      Cp          = nodes->GetSpecificHeatCp(iPoint);
      
      /*--- Compute the conservative variable vector for all time levels. ---*/
      
      U_time_nM1[0] = Density;
      U_time_n[0]   = Density;
      U_time_nP1[0] = Density;
      
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_nM1[iDim+1] = Density*V_time_nM1[iDim+1];
        U_time_n[iDim+1]   = Density*V_time_n[iDim+1];
        U_time_nP1[iDim+1] = Density*V_time_nP1[iDim+1];
      }
      
      U_time_nM1[nDim+1] = Density*Cp*V_time_nM1[nDim+1];
      U_time_n[nDim+1]   = Density*Cp*V_time_n[nDim+1];
      U_time_nP1[nDim+1] = Density*Cp*V_time_nP1[nDim+1];
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order). Note that for an
       incompressible problem, the pressure equation does not have a
       contribution, as the time derivative should always be zero. ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }
      
      if (!energy) Residual[nDim+1] = 0.0;
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      if (implicit) {

        SetPreconditioner(config, iPoint);
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_i[iVar][jVar] = Preconditioner[iVar][jVar];
          }
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            if (config->GetTime_Marching() == DT_STEPPING_1ST)
              Jacobian_i[iVar][jVar] *= Volume_nP1 / TimeStep;
            if (config->GetTime_Marching() == DT_STEPPING_2ND)
              Jacobian_i[iVar][jVar] *= (Volume_nP1*3.0)/(2.0*TimeStep);
          }
        }

        if (!energy) {
            for (iVar = 0; iVar < nVar; iVar++) {
              Jacobian_i[iVar][nDim+1] = 0.0;
              Jacobian_i[nDim+1][iVar] = 0.0;
            }
        }
        
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
      
      /*--- Initialize the Residual / Jacobian container to zero. ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
      
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
      
      V_time_n = nodes->GetSolution_time_n(iPoint);
      
      /*--- Access the density and Cp at this node (constant for now). ---*/
      
      Density     = nodes->GetDensity(iPoint);
      Cp          = nodes->GetSpecificHeatCp(iPoint);
      
      /*--- Compute the conservative variable vector for all time levels. ---*/
      
      U_time_n[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_n[iDim+1] = Density*V_time_n[iDim+1];
      }
      U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];
      
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      
      if (!energy) Residual[nDim+1] = 0.0;
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      V_time_n = nodes->GetSolution_time_n(jPoint);
      
      U_time_n[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_n[iDim+1] = Density*V_time_n[iDim+1];
      }
      U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];
      
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      
      if (!energy) Residual[nDim+1] = 0.0;
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---  Loop over the boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- Initialize the Residual / Jacobian container to zero. ---*/
        
        for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
        
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
        
        V_time_n = nodes->GetSolution_time_n(iPoint);
        
        /*--- Access the density and Cp at this node (constant for now). ---*/
        
        Density     = nodes->GetDensity(iPoint);
        Cp          = nodes->GetSpecificHeatCp(iPoint);
        
        U_time_n[0] = Density;
        for (iDim = 0; iDim < nDim; iDim++) {
          U_time_n[iDim+1] = Density*V_time_n[iDim+1];
        }
        U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];
        
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
        
        if (!energy) Residual[nDim+1] = 0.0;
        LinSysRes.AddBlock(iPoint, Residual);
        
      }
      }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Initialize the Residual / Jacobian container to zero. ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      V_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      V_time_n   = nodes->GetSolution_time_n(iPoint);
      V_time_nP1 = nodes->GetSolution(iPoint);
      
      /*--- Access the density and Cp at this node (constant for now). ---*/
      
      Density     = nodes->GetDensity(iPoint);
      Cp          = nodes->GetSpecificHeatCp(iPoint);
      
      /*--- Compute the conservative variable vector for all time levels. ---*/
      
      U_time_nM1[0] = Density;
      U_time_n[0]   = Density;
      U_time_nP1[0] = Density;
      
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_nM1[iDim+1] = Density*V_time_nM1[iDim+1];
        U_time_n[iDim+1]   = Density*V_time_n[iDim+1];
        U_time_nP1[iDim+1] = Density*V_time_nP1[iDim+1];
      }
      
      U_time_nM1[nDim+1] = Density*Cp*V_time_nM1[nDim+1];
      U_time_n[nDim+1]   = Density*Cp*V_time_n[nDim+1];
      U_time_nP1[nDim+1] = Density*Cp*V_time_nP1[nDim+1];
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      if (!energy) Residual[nDim+1] = 0.0;
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        SetPreconditioner(config, iPoint);
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_i[iVar][jVar] = Preconditioner[iVar][jVar];
          }
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            if (config->GetTime_Marching() == DT_STEPPING_1ST)
              Jacobian_i[iVar][jVar] *= Volume_nP1 / TimeStep;
            if (config->GetTime_Marching() == DT_STEPPING_2ND)
              Jacobian_i[iVar][jVar] *= (Volume_nP1*3.0)/(2.0*TimeStep);
          }
        }

        if (!energy) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Jacobian_i[iVar][nDim+1] = 0.0;
            Jacobian_i[nDim+1][iVar] = 0.0;
          }
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CIncEulerSolver::GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
  
  unsigned short iDim, iMarker;
  unsigned long iVertex, iPoint;
  su2double *V_outlet = NULL, Velocity[3], MassFlow,
  Velocity2, Density, Area, AxiFactor;
  unsigned short iMarker_Outlet, nMarker_Outlet;
  string Inlet_TagBound, Outlet_TagBound;
  
  bool axisymmetric = config->GetAxisymmetric();

  bool write_heads = ((((config->GetInnerIter() % (config->GetWrt_Con_Freq()*40)) == 0)
                       && (config->GetInnerIter()!= 0))
                      || (config->GetInnerIter() == 1));
  
  /*--- Get the number of outlet markers and check for any mass flow BCs. ---*/
  
  nMarker_Outlet = config->GetnMarker_Outlet();
  bool Evaluate_BC = false;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
    Outlet_TagBound = config->GetMarker_Outlet_TagBound(iMarker_Outlet);
    if (config->GetKind_Inc_Outlet(Outlet_TagBound) == MASS_FLOW_OUTLET)
      Evaluate_BC = true;
  }
  
  /*--- If we have a massflow outlet BC, then we need to compute and
   communicate the total massflow, density, and area through each outlet
   boundary, so that it can be used in the iterative procedure to update
   the back pressure until we converge to the desired mass flow. This
   routine is called only once per iteration as a preprocessing and the
   values for all outlets are stored and retrieved later in the BC_Outlet
   routines. ---*/
  
  if (Evaluate_BC) {
    
    su2double *Outlet_MassFlow = new su2double[config->GetnMarker_All()];
    su2double *Outlet_Density  = new su2double[config->GetnMarker_All()];
    su2double *Outlet_Area     = new su2double[config->GetnMarker_All()];
    
    /*--- Comute MassFlow, average temp, press, etc. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      Outlet_MassFlow[iMarker] = 0.0;
      Outlet_Density[iMarker]  = 0.0;
      Outlet_Area[iMarker]     = 0.0;
      
      if ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) ) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (geometry->node[iPoint]->GetDomain()) {
            
            V_outlet = nodes->GetPrimitive(iPoint);
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
            
            if (axisymmetric) {
              if (geometry->node[iPoint]->GetCoord(1) != 0.0)
                AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
              else
                AxiFactor = 1.0;
            } else {
              AxiFactor = 1.0;
            }
            
            Density      = V_outlet[nDim+2];
            
            Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;
            
            for (iDim = 0; iDim < nDim; iDim++) {
              Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
              Velocity[iDim] = V_outlet[iDim+1];
              Velocity2 += Velocity[iDim] * Velocity[iDim];
              MassFlow += Vector[iDim] * AxiFactor * Density * Velocity[iDim];
            }
            Area = sqrt (Area);
            
            Outlet_MassFlow[iMarker] += MassFlow;
            Outlet_Density[iMarker]  += Density*Area;
            Outlet_Area[iMarker]     += Area;
            
          }
        }
      }
    }
    
    /*--- Copy to the appropriate structure ---*/
    
    su2double *Outlet_MassFlow_Local = new su2double[nMarker_Outlet];
    su2double *Outlet_Density_Local  = new su2double[nMarker_Outlet];
    su2double *Outlet_Area_Local     = new su2double[nMarker_Outlet];
    
    su2double *Outlet_MassFlow_Total = new su2double[nMarker_Outlet];
    su2double *Outlet_Density_Total  = new su2double[nMarker_Outlet];
    su2double *Outlet_Area_Total     = new su2double[nMarker_Outlet];
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Local[iMarker_Outlet] = 0.0;
      Outlet_Density_Local[iMarker_Outlet]  = 0.0;
      Outlet_Area_Local[iMarker_Outlet]     = 0.0;
      
      Outlet_MassFlow_Total[iMarker_Outlet] = 0.0;
      Outlet_Density_Total[iMarker_Outlet]  = 0.0;
      Outlet_Area_Total[iMarker_Outlet]     = 0.0;
    }
    
    /*--- Copy the values to the local array for MPI ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW)) {
        for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
          Outlet_TagBound = config->GetMarker_Outlet_TagBound(iMarker_Outlet);
          if (config->GetMarker_All_TagBound(iMarker) == Outlet_TagBound) {
            Outlet_MassFlow_Local[iMarker_Outlet] += Outlet_MassFlow[iMarker];
            Outlet_Density_Local[iMarker_Outlet]  += Outlet_Density[iMarker];
            Outlet_Area_Local[iMarker_Outlet]     += Outlet_Area[iMarker];
          }
        }
      }
    }
    
    /*--- All the ranks to compute the total value ---*/
    
#ifdef HAVE_MPI
    
    SU2_MPI::Allreduce(Outlet_MassFlow_Local, Outlet_MassFlow_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Density_Local, Outlet_Density_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Area_Local, Outlet_Area_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
#else
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Total[iMarker_Outlet] = Outlet_MassFlow_Local[iMarker_Outlet];
      Outlet_Density_Total[iMarker_Outlet]  = Outlet_Density_Local[iMarker_Outlet];
      Outlet_Area_Total[iMarker_Outlet]     = Outlet_Area_Local[iMarker_Outlet];
    }
    
#endif
    
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      if (Outlet_Area_Total[iMarker_Outlet] != 0.0) {
        Outlet_Density_Total[iMarker_Outlet] /= Outlet_Area_Total[iMarker_Outlet];
      }
      else {
        Outlet_Density_Total[iMarker_Outlet] = 0.0;
      }
      
      if (iMesh == MESH_0) {
        config->SetOutlet_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
        config->SetOutlet_Density(iMarker_Outlet, Outlet_Density_Total[iMarker_Outlet]);
        config->SetOutlet_Area(iMarker_Outlet, Outlet_Area_Total[iMarker_Outlet]);
      }
    }
    
    /*--- Screen output using the values already stored in the config container ---*/
    
    if ((rank == MASTER_NODE) && (iMesh == MESH_0) ) {
      
      cout.precision(5);
      cout.setf(ios::fixed, ios::floatfield);
      
      if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
        cout << endl   << "---------------------------- Outlet properties --------------------------" << endl;
      }
      
      for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
        Outlet_TagBound = config->GetMarker_Outlet_TagBound(iMarker_Outlet);
        if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
          
          /*--- Geometry defintion ---*/
          
          cout <<"Outlet surface: " << Outlet_TagBound << "." << endl;
          
          if ((nDim ==3) || axisymmetric) {
            cout <<"Area (m^2): " << config->GetOutlet_Area(Outlet_TagBound) << endl;
          }
          if (nDim == 2) {
            cout <<"Length (m): " << config->GetOutlet_Area(Outlet_TagBound) << "." << endl;
          }
          
          cout << setprecision(5) << "Outlet Avg. Density (kg/m^3): " <<  config->GetOutlet_Density(Outlet_TagBound) * config->GetDensity_Ref() << endl;
          su2double Outlet_mDot = fabs(config->GetOutlet_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
          cout << "Outlet mass flow (kg/s): "; cout << setprecision(5) << Outlet_mDot;
          
        }
      }
      
      if (write_heads && Output && !config->GetDiscrete_Adjoint()) {cout << endl;
        cout << "-------------------------------------------------------------------------" << endl << endl;
      }
      
      cout.unsetf(ios_base::floatfield);
      
    }
    
    delete [] Outlet_MassFlow_Local;
    delete [] Outlet_Density_Local;
    delete [] Outlet_Area_Local;
    
    delete [] Outlet_MassFlow_Total;
    delete [] Outlet_Density_Total;
    delete [] Outlet_Area_Total;
    
    delete [] Outlet_MassFlow;
    delete [] Outlet_Density;
    delete [] Outlet_Area;
    
  }
  
}

void CIncEulerSolver::ComputeVerificationError(CGeometry *geometry,
                                               CConfig   *config) {

  /*--- The errors only need to be computed on the finest grid. ---*/
  if(MGLevel != MESH_0) return;
  
  /*--- If this is a verification case, we can compute the global
   error metrics by using the difference between the local error
   and the known solution at each DOF. This is then collected into
   RMS (L2) and maximum (Linf) global error norms. From these
   global measures, one can compute the order of accuracy. ---*/
  
  bool write_heads = ((((config->GetInnerIter() % (config->GetWrt_Con_Freq()*40)) == 0)
                       && (config->GetInnerIter()!= 0))
                      || (config->GetInnerIter() == 1));
  if( !write_heads ) return;
  
  /*--- Check if there actually is an exact solution for this
        verification case, if computed at all. ---*/
  if (VerificationSolution) {
    if (VerificationSolution->ExactSolutionKnown()) {
    
      /*--- Get the physical time if necessary. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching()) time = config->GetPhysicalTime();
    
      /*--- Reset the global error measures to zero. ---*/
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        VerificationSolution->SetError_RMS(iVar, 0.0);
        VerificationSolution->SetError_Max(iVar, 0.0, 0);
      }
    
      /*--- Loop over all owned points. ---*/
      for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
        /* Set the pointers to the coordinates and solution of this DOF. */
        const su2double *coor = geometry->node[iPoint]->GetCoord();
        su2double *solDOF     = nodes->GetSolution(iPoint);
      
        /* Get local error from the verification solution class. */
        vector<su2double> error(nVar,0.0);
        VerificationSolution->GetLocalError(coor, time, solDOF, error.data());
      
        /* Increment the global error measures */
        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
          VerificationSolution->AddError_RMS(iVar, error[iVar]*error[iVar]);
          VerificationSolution->AddError_Max(iVar, fabs(error[iVar]),
                                             geometry->node[iPoint]->GetGlobalIndex(),
                                             geometry->node[iPoint]->GetCoord());
        }
      }
    
      /* Finalize the calculation of the global error measures. */
      VerificationSolution->SetVerificationError(geometry->GetGlobal_nPointDomain(), config);
    
      /*--- Screen output of the error metrics. This can be improved
       once the new output classes are in place. ---*/
    
      if ((rank == MASTER_NODE) && (geometry->GetMGLevel() == MESH_0)) {
      
        cout.precision(6);
        cout.setf(ios::scientific, ios::floatfield);
      
        if (!config->GetDiscrete_Adjoint()) {
        
          cout << endl   << "------------------------ Global Error Analysis --------------------------" << endl;
        
          cout << setw(20) << "RMS Error [P]: " << setw(12) << VerificationSolution->GetError_RMS(0) << "     | ";
          cout << setw(20) << "Max Error [P]: " << setw(12) << VerificationSolution->GetError_Max(0);
          cout << endl;
        
          cout << setw(20) << "RMS Error [U]: " << setw(12) << VerificationSolution->GetError_RMS(1) << "     | ";
          cout << setw(20) << "Max Error [U]: " << setw(12) << VerificationSolution->GetError_Max(1);
          cout << endl;
        
          cout << setw(20) << "RMS Error [V]: " << setw(12) << VerificationSolution->GetError_RMS(2) << "     | ";
          cout << setw(20) << "Max Error [V]: " << setw(12) << VerificationSolution->GetError_Max(2);
          cout << endl;
        
          if (nDim == 3) {
            cout << setw(20) << "RMS Error [W]: " << setw(12) << VerificationSolution->GetError_RMS(3) << "     | ";
            cout << setw(20) << "Max Error [W]: " << setw(12) << VerificationSolution->GetError_Max(3);
            cout << endl;
          }
        
          if (config->GetEnergy_Equation()) {
            cout << setw(20) << "RMS Error [T]: " << setw(12) << VerificationSolution->GetError_RMS(nDim+1) << "     | ";
            cout << setw(20) << "Max Error [T]: " << setw(12) << VerificationSolution->GetError_Max(nDim+1);
            cout << endl;
          }
        
          cout << "-------------------------------------------------------------------------" << endl << endl;
          cout.unsetf(ios_base::floatfield);
        }
      }
    }
    
  }
  
}

void CIncEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {
  
  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool static_fsi = ((config->GetTime_Marching() == STEADY) &&
                     (config->GetFSI_Simulation()));
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool steady_restart = config->GetSteadyRestart();
  bool turbulent     = (config->GetKind_Solver() == INC_RANS) || (config->GetKind_Solver() == DISC_ADJ_INC_RANS);
  
  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;
  
  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Store the number of variables for the turbulence model
   (that could appear in the restart file before the grid velocities). ---*/
  unsigned short turbVars = 0;
  if (turbulent){
    if ((turb_model == SST) || (turb_model == SST_SUST)) turbVars = 2;
    else turbVars = 1;
  }
  
  /*--- Adjust the number of solution variables in the restart. We always
   carry a space in nVar for the energy equation in the solver, but we only
   write it to the restart if it is active. Therefore, we must reduce nVar
   here if energy is inactive so that the restart is read correctly. ---*/
  
  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();
  
  unsigned short nVar_Restart = nVar;
  if ((!energy) && (!weakly_coupled_heat)) nVar_Restart--;
  Solution[nVar-1] = GetTemperature_Inf();
  
  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar_Restart; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);
      iPoint_Global_Local++;

      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/

      if (dynamic_grid && val_update_geo) {

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/

        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar_Restart + turbVars;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[MESH_0]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
      }
      

      /*--- For static FSI problems, grid_movement is 0 but we need to read in and store the
       grid coordinates for each node (but not the grid velocities, as there are none). ---*/

      if (static_fsi && val_update_geo) {
       /*--- Rewind the index to retrieve the Coords. ---*/
        index = counter*Restart_Vars[1];
        for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim];}

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
        }
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
      
    }
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }
  
  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/
  
  solver[MESH_0][FLOW_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][FLOW_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);
  
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][FLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint,Solution);
    }
    solver[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }
  
  /*--- Update the geometry for flows on dynamic meshes ---*/
  
  if (dynamic_grid && val_update_geo) {
    
    /*--- Communicate the new coordinates and grid velocities at the halos ---*/
    
    geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
    geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);
    
    geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, GRID_VELOCITY);
    geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, GRID_VELOCITY);
    
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
  
  /*--- Update the geometry for flows on static FSI problems with moving meshes ---*/
  
  if (static_fsi && val_update_geo) {
    
    /*--- Communicate the new coordinates and grid velocities at the halos ---*/
    
    geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
    geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);
    
    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/
    
    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);
    geometry[MESH_0]->SetMaxLength(config);
    
    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/
    
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      geometry[iMesh]->SetMaxLength(config);
    }
  }
  
  /*--- Update the old geometry (coordinates n and n-1) in dual time-stepping strategy ---*/
  if (dual_time && config->GetGrid_Movement() && (config->GetKind_GridMovement() != RIGID_MOTION))
    Restart_OldGeometry(geometry[MESH_0], config);

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;
  
}

void CIncEulerSolver::SetFreeStream_Solution(CConfig *config){

  unsigned long iPoint;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    nodes->SetSolution(iPoint,0, Pressure_Inf);
    for (iDim = 0; iDim < nDim; iDim++){
      nodes->SetSolution(iPoint,iDim+1, Velocity_Inf[iDim]);
    }
    nodes->SetSolution(iPoint,nDim+1, Temperature_Inf);
  }
}

CIncNSSolver::CIncNSSolver(void) : CIncEulerSolver() {
  
  /*--- Basic array initialization ---*/
  
  CD_Visc = NULL; CL_Visc = NULL; CSF_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  CoPx_Visc = NULL;   CoPy_Visc = NULL;   CoPz_Visc = NULL;

  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;
  
  /*--- Surface based array initialization ---*/
  
  Surface_CL_Visc = NULL; Surface_CD_Visc = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  Surface_HF_Visc = NULL; Surface_MaxHF_Visc = NULL;

  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Visc = NULL; CT_Visc = NULL; CQ_Visc = NULL;
  
  SlidingState      = NULL;
  SlidingStateNodes = NULL;
  
}

CIncNSSolver::CIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CIncEulerSolver() {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  string filename_ = config->GetSolution_FileName();

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  unsigned short direct_diff = config->GetDirectDiff();

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone, ".dat");

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }

    /*--- Read and store the restart metadata. ---*/

//    Read_SU2_Restart_Metadata(geometry, config, false, filename_);
    
  }

  /*--- Array initialization ---*/
  
  CD_Visc = NULL; CL_Visc = NULL; CSF_Visc = NULL; CEff_Visc = NULL;
  CMx_Visc = NULL;   CMy_Visc = NULL;   CMz_Visc = NULL;
  CFx_Visc = NULL;   CFy_Visc = NULL;   CFz_Visc = NULL;
  CoPx_Visc = NULL;   CoPy_Visc = NULL;   CoPz_Visc = NULL;

  Surface_CL_Visc = NULL; Surface_CD_Visc = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL;   Surface_CFy_Visc = NULL;   Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL;   Surface_CMy_Visc = NULL;   Surface_CMz_Visc = NULL;
  Surface_HF_Visc = NULL; Surface_MaxHF_Visc = NULL;

  CMerit_Visc = NULL;      CT_Visc = NULL;      CQ_Visc = NULL;
  MaxHF_Visc = NULL; ForceViscous = NULL; MomentViscous = NULL;
  CSkinFriction = NULL; HF_Visc = NULL;
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   * Incompressible flow, primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) --- */

  nDim = geometry->GetnDim();
  
  nVar = nDim+2; nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nPrimVarGrad;
  
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
 
  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];
 
  /*--- Fluid model intialization. ---*/

  FluidModel = NULL;
  
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  SetNondimensionalization(config, iMesh);
  
  /*--- Check if we are executing a verification case. If so, the
   VerificationSolution object will be instantiated for a particular
   option from the available library of verification solutions. Note
   that this is done after SetNondim(), as problem-specific initial
   parameters are needed by the solution constructors. ---*/
  
  SetVerificationSolution(nDim, nVar, config);
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
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
  
  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }

  Preconditioner = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar ++)
    Preconditioner[iVar] = new su2double[nVar];

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
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
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
    
    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];
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

  /*--- Store the values of the temperature and the heat flux density at the boundaries,
   used for coupling with a solid donor cell ---*/
  unsigned short nHeatConjugateVar = 4;

  HeatConjugateVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatConjugateVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      HeatConjugateVar[iMarker][iVertex] = new su2double [nHeatConjugateVar];
      for (iVar = 1; iVar < nHeatConjugateVar ; iVar++) {
        HeatConjugateVar[iMarker][iVertex][iVar] = 0.0;
      }
      HeatConjugateVar[iMarker][iVertex][0] = config->GetTemperature_FreeStreamND();
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
  
  CSkinFriction = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSkinFriction[iMarker] = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      CSkinFriction[iMarker][iDim] = new su2double[geometry->nVertex[iMarker]];
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        CSkinFriction[iMarker][iDim][iVertex] = 0.0;
      }
    }
  }
  
  /*--- Store the value of the Total Pressure at the inlet BC ---*/
  
  Inlet_Ttotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_Ttotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_Ttotal[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Total Temperature at the inlet BC ---*/
  
  Inlet_Ptotal = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_Ptotal[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_Ptotal[iMarker][iVertex] = 0;
    }
  }
  
  /*--- Store the value of the Flow direction at the inlet BC ---*/
  
  Inlet_FlowDir = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_FlowDir[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      Inlet_FlowDir[iMarker][iVertex] = new su2double [nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        Inlet_FlowDir[iMarker][iVertex][iDim] = 0;
      }
    }
  }
  
  /*--- Non dimensional coefficients ---*/
  
  ForceInviscid  = new su2double[3];
  MomentInviscid = new su2double[3];
  CD_Inv      = new su2double[nMarker];
  CL_Inv      = new su2double[nMarker];
  CSF_Inv = new su2double[nMarker];
  CMx_Inv        = new su2double[nMarker];
  CMy_Inv        = new su2double[nMarker];
  CMz_Inv        = new su2double[nMarker];
  CEff_Inv       = new su2double[nMarker];
  CFx_Inv        = new su2double[nMarker];
  CFy_Inv        = new su2double[nMarker];
  CFz_Inv        = new su2double[nMarker];
  CoPx_Inv        = new su2double[nMarker];
  CoPy_Inv        = new su2double[nMarker];
  CoPz_Inv        = new su2double[nMarker];

  ForceMomentum  = new su2double[3];
  MomentMomentum = new su2double[3];
  CD_Mnt      = new su2double[nMarker];
  CL_Mnt      = new su2double[nMarker];
  CSF_Mnt = new su2double[nMarker];
  CMx_Mnt        = new su2double[nMarker];
  CMy_Mnt        = new su2double[nMarker];
  CMz_Mnt        = new su2double[nMarker];
  CEff_Mnt       = new su2double[nMarker];
  CFx_Mnt        = new su2double[nMarker];
  CFy_Mnt        = new su2double[nMarker];
  CFz_Mnt        = new su2double[nMarker];
  CoPx_Mnt        = new su2double[nMarker];
  CoPy_Mnt        = new su2double[nMarker];
  CoPz_Mnt        = new su2double[nMarker];

  ForceViscous     = new su2double[3];
  MomentViscous    = new su2double[3];
  CD_Visc       = new su2double[nMarker];
  CL_Visc       = new su2double[nMarker];
  CSF_Visc  = new su2double[nMarker];
  CMx_Visc         = new su2double[nMarker];
  CMy_Visc         = new su2double[nMarker];
  CMz_Visc         = new su2double[nMarker];
  CEff_Visc        = new su2double[nMarker];
  CFx_Visc         = new su2double[nMarker];
  CFy_Visc         = new su2double[nMarker];
  CFz_Visc         = new su2double[nMarker];
  CoPx_Visc         = new su2double[nMarker];
  CoPy_Visc         = new su2double[nMarker];
  CoPz_Visc         = new su2double[nMarker];

  Surface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Inv = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Mnt      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Mnt = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Mnt       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Mnt        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Mnt        = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD          = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff           = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy            = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz            = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Visc      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Visc = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  Surface_HF_Visc         = new su2double[config->GetnMarker_Monitoring()];
  Surface_MaxHF_Visc      = new su2double[config->GetnMarker_Monitoring()];
  
  /*--- Rotorcraft coefficients ---*/

  CT_Inv           = new su2double[nMarker];
  CQ_Inv           = new su2double[nMarker];
  CMerit_Inv       = new su2double[nMarker];

  CT_Mnt           = new su2double[nMarker];
  CQ_Mnt           = new su2double[nMarker];
  CMerit_Mnt       = new su2double[nMarker];

  CMerit_Visc      = new su2double[nMarker];
  CT_Visc          = new su2double[nMarker];
  CQ_Visc          = new su2double[nMarker];
  
  /*--- Heat based coefficients ---*/

  HF_Visc    = new su2double[nMarker];
  MaxHF_Visc = new su2double[nMarker];

  /*--- Init total coefficients ---*/

  Total_CD       = 0.0;  Total_CL           = 0.0;  Total_CSF            = 0.0;
  Total_CMx      = 0.0;  Total_CMy          = 0.0;  Total_CMz            = 0.0;
  Total_CoPx     = 0.0;  Total_CoPy         = 0.0;  Total_CoPz           = 0.0;
  Total_CEff     = 0.0;
  Total_CFx      = 0.0;  Total_CFy          = 0.0;  Total_CFz            = 0.0;
  Total_CT       = 0.0;  Total_CQ           = 0.0;  Total_CMerit         = 0.0;
  Total_MaxHeat  = 0.0;  Total_Heat         = 0.0;  Total_ComboObj       = 0.0;
  Total_CpDiff   = 0.0;  Total_HeatFluxDiff = 0.0;  Total_Custom_ObjFunc = 0.0;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
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

  /*--- Initializate quantities for SlidingMesh Interface ---*/
  
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++){

    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;
    
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }

    }
  }

  /*--- Only initialize when there is a Marker_Fluid_Load defined
   *--- (this avoids overhead in all other cases while a more permanent structure is being developed) ---*/
  if((config->GetnMarker_Fluid_Load() > 0) && (MGLevel == MESH_0)){

    InitVertexTractionContainer(geometry, config);

    if (config->GetDiscrete_Adjoint())
      InitVertexTractionAdjointContainer(geometry, config);

  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CIncNSVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (config->GetMultizone_Residual()){
    Residual_BGS      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_Max_BGS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 0.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
    }
  }

  /*--- Define solver parameters needed for execution of destructor ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Communicate and store volume and the number of neighbors for
   any dual CVs that lie on on periodic markers. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
  }
  SetImplicitPeriodic(euler_implicit);
  if (iMesh == MESH_0) SetRotatePeriodic(true);
  
  /*--- Perform the MPI communication of the solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "INC.FLOW";

}

CIncNSSolver::~CIncNSSolver(void) {

  unsigned short iMarker, iDim;

  unsigned long iVertex;

  if (CD_Visc != NULL)       delete [] CD_Visc;
  if (CL_Visc != NULL)       delete [] CL_Visc;
  if (CSF_Visc != NULL)  delete [] CSF_Visc;
  if (CMx_Visc != NULL)         delete [] CMx_Visc;
  if (CMy_Visc != NULL)         delete [] CMy_Visc;
  if (CMz_Visc != NULL)         delete [] CMz_Visc;
  if (CoPx_Visc != NULL)        delete [] CoPx_Visc;
  if (CoPy_Visc != NULL)        delete [] CoPy_Visc;
  if (CoPz_Visc != NULL)        delete [] CoPz_Visc;
  if (CFx_Visc != NULL)         delete [] CFx_Visc;
  if (CFy_Visc != NULL)         delete [] CFy_Visc;
  if (CFz_Visc != NULL)         delete [] CFz_Visc;
  if (CEff_Visc != NULL)        delete [] CEff_Visc;
  if (CMerit_Visc != NULL)      delete [] CMerit_Visc;
  if (CT_Visc != NULL)          delete [] CT_Visc;
  if (CQ_Visc != NULL)          delete [] CQ_Visc;
  if (HF_Visc != NULL)        delete [] HF_Visc;
  if (MaxHF_Visc != NULL) delete [] MaxHF_Visc;
  if (ForceViscous != NULL)     delete [] ForceViscous;
  if (MomentViscous != NULL)    delete [] MomentViscous;

  if (Surface_CL_Visc != NULL)      delete [] Surface_CL_Visc;
  if (Surface_CD_Visc != NULL)      delete [] Surface_CD_Visc;
  if (Surface_CSF_Visc != NULL) delete [] Surface_CSF_Visc;
  if (Surface_CEff_Visc != NULL)       delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)        delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)        delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)        delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)        delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)        delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)        delete [] Surface_CMz_Visc;
  if (Surface_HF_Visc != NULL)      delete [] Surface_HF_Visc;
  if (Surface_MaxHF_Visc != NULL)   delete [] Surface_MaxHF_Visc;

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
  
}

void CIncNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double StrainMag = 0.0, Omega = 0.0, *Vorticity;
  
  unsigned long InnerIter     = config->GetInnerIter();
  bool cont_adjoint         = config->GetContinuous_Adjoint();
  bool disc_adjoint         = config->GetDiscrete_Adjoint();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center               = ((config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED));
  bool center_jst           = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow         = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool limiter_turb         = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool limiter_adjflow      = (cont_adjoint && (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()));
  bool van_albada           = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  bool outlet               = ((config->GetnMarker_Outlet() != 0));

  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Compute gradient of the primitive variables ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/

  if ((iMesh == MESH_0) && (limiter_flow || limiter_turb || limiter_adjflow)
      && !Output && !van_albada) { SetPrimitive_Limiter(geometry, config); }
  
  /*--- Artificial dissipation for centered schemes. ---*/
  
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Update the beta value based on the maximum velocity / viscosity. ---*/

  SetBeta_Parameter(geometry, solver_container, config, iMesh);

  /*--- Compute properties needed for mass flow BCs. ---*/
  
  if (outlet) GetOutlet_Properties(geometry, config, iMesh, Output);
  
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  
  solver_container[FLOW_SOL]->GetNodes()->SetVorticity_StrainMag();

  StrainMag_Max = 0.0; Omega_Max = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    StrainMag = solver_container[FLOW_SOL]->GetNodes()->GetStrainMag(iPoint);
    Vorticity = solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint);
    Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);
    
    StrainMag_Max = max(StrainMag_Max, StrainMag);
    Omega_Max = max(Omega_Max, Omega);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !disc_adjoint) Jacobian.SetValZero();

  /*--- Error message ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
    
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

unsigned long CIncNSSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0, DES_LengthScale = 0.0;
  unsigned short turb_model = config->GetKind_Turb_Model();
  bool physical = true;
  
  bool tkeNeeded = ((turb_model == SST) || (turb_model == SST_SUST));
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Retrieve the value of the kinetic energy (if needed) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
      
      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
        DES_LengthScale = solver_container[TURB_SOL]->GetNodes()->GetDES_LengthScale(iPoint);
      }
    }
    
    /*--- Initialize the non-physical points vector ---*/
    
    nodes->SetNon_Physical(iPoint,false);
    
    /*--- Incompressible flow, primitive variables --- */

    physical = static_cast<CIncNSVariable*>(nodes)->SetPrimVar(iPoint,eddy_visc, turb_ke, FluidModel);
    
    /*--- Record any non-physical points. ---*/

    if (!physical) { nodes->SetNon_Physical(iPoint,true); ErrorCounter++; }

    /*--- Set the DES length scale ---*/
    
    nodes->SetDES_LengthScale(iPoint,DES_LengthScale);    
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }

  return ErrorCounter;

}

void CIncNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  
  su2double Mean_BetaInc2, *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time, Local_Delta_Time_Visc,
  Global_Delta_Time = 1E6, Mean_LaminarVisc = 0.0, Mean_EddyVisc = 0.0, Mean_Density = 0.0, Mean_Thermal_Conductivity = 0.0, Mean_Cv = 0.0, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  unsigned short iDim, iMarker;
  su2double ProjVel, ProjVel_i, ProjVel_j;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool energy = config->GetEnergy_Equation();

  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    nodes->SetMax_Lambda_Inv(iPoint,0.0);
    nodes->SetMax_Lambda_Visc(iPoint,0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel    = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
    Mean_BetaInc2   = 0.5 * (nodes->GetBetaInc2(iPoint)       + nodes->GetBetaInc2(jPoint));
    Mean_Density    = 0.5 * (nodes->GetDensity(iPoint)        + nodes->GetDensity(jPoint));
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);
    
    /*--- Adjustment for grid movement ---*/
    
    if (dynamic_grid) {
      su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(iPoint,Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(jPoint,Lambda);
    
    /*--- Viscous contribution ---*/
    
    Mean_LaminarVisc          = 0.5*(nodes->GetLaminarViscosity(iPoint)    + nodes->GetLaminarViscosity(jPoint));
    Mean_EddyVisc             = 0.5*(nodes->GetEddyViscosity(iPoint)       + nodes->GetEddyViscosity(jPoint));
    Mean_Density              = 0.5*(nodes->GetDensity(iPoint)             + nodes->GetDensity(jPoint));
    Mean_Thermal_Conductivity = 0.5*(nodes->GetThermalConductivity(iPoint) + nodes->GetThermalConductivity(jPoint));
    Mean_Cv                   = 0.5*(nodes->GetSpecificHeatCv(iPoint)      + nodes->GetSpecificHeatCv(jPoint));

    Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
    Lambda_2 = 0.0;
    if (energy) Lambda_2 = (1.0/Mean_Cv)*Mean_Thermal_Conductivity;
    Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
    
    if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Visc(iPoint,Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddMax_Lambda_Visc(jPoint,Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel    = nodes->GetProjVel(iPoint,Normal);
      Mean_BetaInc2   = nodes->GetBetaInc2(iPoint);
      Mean_Density    = nodes->GetDensity(iPoint);
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

      /*--- Adjustment for grid movement ---*/
      
      if (dynamic_grid) {
        su2double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        nodes->AddMax_Lambda_Inv(iPoint,Lambda);
      }
      
      /*--- Viscous contribution ---*/

      Mean_LaminarVisc          = nodes->GetLaminarViscosity(iPoint);
      Mean_EddyVisc             = nodes->GetEddyViscosity(iPoint);
      Mean_Density              = nodes->GetDensity(iPoint);
      Mean_Thermal_Conductivity = nodes->GetThermalConductivity(iPoint);
      Mean_Cv                   = nodes->GetSpecificHeatCv(iPoint);

      Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
      Lambda_2 = 0.0;
      if (energy) Lambda_2 = (1.0/Mean_Cv)*Mean_Thermal_Conductivity;
      Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      
      if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Visc(iPoint,Lambda);
      
    }
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / nodes->GetMax_Lambda_Inv(iPoint);
      Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ nodes->GetMax_Lambda_Visc(iPoint);
      Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      nodes->SetDelta_Time(iPoint,Local_Delta_Time);
    }
    else {
      nodes->SetDelta_Time(iPoint,0.0);
    }
    
  }
  
  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetComm_Level() == COMM_FULL) {
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
  if (config->GetTime_Marching() == TIME_STEPPING) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      nodes->SetDelta_Time(iPoint,Global_Delta_Time);
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
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), nodes->GetDelta_Time(iPoint));
        nodes->SetDelta_Time(iPoint,Local_Delta_Time);
      }
    }
  
}

void CIncNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive and secondary variables ---*/
    
    numerics->SetPrimitive(nodes->GetPrimitive(iPoint),
                           nodes->GetPrimitive(jPoint));
    
    /*--- Gradient and limiters ---*/
    
    numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                 nodes->GetGradient_Primitive(jPoint));
    
    /*--- Turbulent kinetic energy ---*/
    
    if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
      numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                     solver_container[TURB_SOL]->GetNodes()->GetSolution(jPoint,0));
    
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

void CIncNSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  su2double Viscosity = 0.0, div_vel, *Normal, MomentDist[3] = {0.0, 0.0, 0.0}, WallDist[3] = {0.0, 0.0, 0.0},
  *Coord, *Coord_Normal, Area, WallShearStress, TauNormal, factor, RefVel2 = 0.0,
  RefDensity = 0.0, Density = 0.0, WallDistMod, FrictionVel, UnitNormal[3] = {0.0, 0.0, 0.0}, TauElem[3] = {0.0, 0.0, 0.0}, TauTangent[3] = {0.0, 0.0, 0.0},
  Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}}, Force[3] = {0.0, 0.0, 0.0},
  Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}},
  delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}},
  Grad_Temp[3] = {0.0, 0.0, 0.0}, GradTemperature, thermal_conductivity, MaxNorm = 8.0;
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0}, MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Visc, MyAllBound_CL_Visc, MyAllBound_CSF_Visc, MyAllBound_CMx_Visc, MyAllBound_CMy_Visc, MyAllBound_CMz_Visc, MyAllBound_CoPx_Visc, MyAllBound_CoPy_Visc, MyAllBound_CoPz_Visc, MyAllBound_CFx_Visc, MyAllBound_CFy_Visc, MyAllBound_CFz_Visc, MyAllBound_CT_Visc, MyAllBound_CQ_Visc, MyAllBound_HF_Visc, MyAllBound_MaxHF_Visc, *MySurface_CL_Visc = NULL, *MySurface_CD_Visc = NULL, *MySurface_CSF_Visc = NULL, *MySurface_CEff_Visc = NULL, *MySurface_CFx_Visc = NULL, *MySurface_CFy_Visc = NULL, *MySurface_CFz_Visc = NULL, *MySurface_CMx_Visc = NULL, *MySurface_CMy_Visc = NULL, *MySurface_CMz_Visc = NULL, *MySurface_HF_Visc = NULL, *MySurface_MaxHF_Visc = NULL;
#endif

  string Marker_Tag, Monitoring_Tag;

  su2double Alpha       = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta        = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea     = config->GetRefArea();
  su2double RefLength   = config->GetRefLength();
  su2double RefHeatFlux = config->GetHeat_Flux_Ref();
  su2double *Origin = NULL;

  if (config->GetnMarker_Monitoring() != 0) { Origin = config->GetRefOriginMoment(0); }

  bool axisymmetric = config->GetAxisymmetric();
  bool energy       = config->GetEnergy_Equation();

  /*--- Evaluate reference values for non-dimensionalization.
   For dimensional or non-dim based on initial values, use
   the far-field state (inf). For a custom non-dim based
   on user-provided reference values, use the ref values
   to compute the forces. ---*/

  if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) || 
      (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
    RefDensity  = Density_Inf;
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    RefDensity = config->GetInc_Density_Ref();
    RefVel2    = config->GetInc_Velocity_Ref()*config->GetInc_Velocity_Ref();
  }

  /*--- Compute factor for force coefficients. ---*/

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*--- Variables initialization ---*/

  AllBound_CD_Visc = 0.0;    AllBound_CL_Visc = 0.0;       AllBound_CSF_Visc = 0.0;
  AllBound_CMx_Visc = 0.0;      AllBound_CMy_Visc = 0.0;         AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc = 0.0;      AllBound_CFy_Visc = 0.0;         AllBound_CFz_Visc = 0.0;
  AllBound_CoPx_Visc = 0.0;      AllBound_CoPy_Visc = 0.0;         AllBound_CoPz_Visc = 0.0;
  AllBound_CT_Visc = 0.0;       AllBound_CQ_Visc = 0.0;          AllBound_CMerit_Visc = 0.0;
  AllBound_HF_Visc = 0.0; AllBound_MaxHF_Visc = 0.0; AllBound_CEff_Visc = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Visc[iMarker_Monitoring]      = 0.0; Surface_CD_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0; Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0; Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0; Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]              = 0.0; Surface_MaxHF_Visc[iMarker_Monitoring]           = 0.0;
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

    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL) || (Boundary == CHT_WALL_INTERFACE)) {

      /*--- Forces initialization at each Marker ---*/

      CD_Visc[iMarker] = 0.0; CL_Visc[iMarker] = 0.0;       CSF_Visc[iMarker] = 0.0;
      CMx_Visc[iMarker] = 0.0;   CMy_Visc[iMarker] = 0.0;         CMz_Visc[iMarker] = 0.0;
      CFx_Visc[iMarker] = 0.0;   CFy_Visc[iMarker] = 0.0;         CFz_Visc[iMarker] = 0.0;
      CoPx_Visc[iMarker] = 0.0;  CoPy_Visc[iMarker] = 0.0;       CoPz_Visc[iMarker] = 0.0;
      CT_Visc[iMarker] = 0.0;    CQ_Visc[iMarker] = 0.0;          CMerit_Visc[iMarker] = 0.0;
      HF_Visc[iMarker] = 0.0;  MaxHF_Visc[iMarker] = 0.0; CEff_Visc[iMarker] = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
      MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;
      MomentX_Force[0] = 0.0; MomentX_Force[1] = 0.0; MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0; MomentY_Force[1] = 0.0; MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0; MomentZ_Force[1] = 0.0; MomentZ_Force[2] = 0.0;

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
          Grad_Temp[iDim] = nodes->GetGradient_Primitive(iPoint,nDim+1, iDim);
        }

        Viscosity = nodes->GetLaminarViscosity(iPoint);
        Density = nodes->GetDensity(iPoint);

        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
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

        /*--- Compute wall shear stress (using the stress tensor). Compute wall skin friction coefficient, and heat flux on the wall ---*/

        TauNormal = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          TauNormal += TauElem[iDim] * UnitNormal[iDim];

        WallShearStress = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
          CSkinFriction[iMarker][iDim][iVertex] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);
          WallShearStress += TauTangent[iDim] * TauTangent[iDim];
        }
        WallShearStress = sqrt(WallShearStress);

        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);

        /*--- Compute y+ and non-dimensional velocity ---*/

        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);

        /*--- Compute total and maximum heat flux on the wall ---*/

        GradTemperature = 0.0;
        if (energy) {
        for (iDim = 0; iDim < nDim; iDim++)
          GradTemperature -= Grad_Temp[iDim]*UnitNormal[iDim];
        }
        
        thermal_conductivity       = nodes->GetThermalConductivity(iPoint);
        HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature*RefHeatFlux;

        /*--- Note that y+, and heat are computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {

          /*--- Axisymmetric simulations ---*/

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Force computation ---*/

          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim] * Area * factor * AxiFactor;
            ForceViscous[iDim] += Force[iDim];
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (iDim == 3) {
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
          
          HF_Visc[iMarker]          += HeatFlux[iMarker][iVertex]*Area;
          MaxHF_Visc[iMarker]       += pow(HeatFlux[iMarker][iVertex], MaxNorm);
          
        }

      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if (Monitoring == YES) {
        if (nDim == 2) {
          CD_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CL_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker]        = CL_Visc[iMarker] / (CD_Visc[iMarker]+EPS);
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CoPx_Visc[iMarker]        = MomentZ_Force[1];
          CoPy_Visc[iMarker]        = -MomentZ_Force[0];
          CT_Visc[iMarker]          = -CFx_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker]+EPS);
          MaxHF_Visc[iMarker] = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          CD_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          CL_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSF_Visc[iMarker]  = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker]        = CL_Visc[iMarker]/(CD_Visc[iMarker] + EPS);
          CMx_Visc[iMarker]         = MomentViscous[0];
          CMy_Visc[iMarker]         = MomentViscous[1];
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CFz_Visc[iMarker]         = ForceViscous[2];
          CoPx_Visc[iMarker]        =  -MomentY_Force[0];
          CoPz_Visc[iMarker]        = MomentY_Force[2];
          CT_Visc[iMarker]          = -CFz_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker] / (CQ_Visc[iMarker] + EPS);
          MaxHF_Visc[iMarker] = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }

        AllBound_CD_Visc       += CD_Visc[iMarker];
        AllBound_CL_Visc       += CL_Visc[iMarker];
        AllBound_CSF_Visc  += CSF_Visc[iMarker];
        AllBound_CMx_Visc         += CMx_Visc[iMarker];
        AllBound_CMy_Visc         += CMy_Visc[iMarker];
        AllBound_CMz_Visc         += CMz_Visc[iMarker];
        AllBound_CFx_Visc         += CFx_Visc[iMarker];
        AllBound_CFy_Visc         += CFy_Visc[iMarker];
        AllBound_CFz_Visc         += CFz_Visc[iMarker];
        AllBound_CoPx_Visc        += CoPx_Visc[iMarker];
        AllBound_CoPy_Visc        += CoPy_Visc[iMarker];
        AllBound_CoPz_Visc        += CoPz_Visc[iMarker];
        AllBound_CT_Visc          += CT_Visc[iMarker];
        AllBound_CQ_Visc          += CQ_Visc[iMarker];
        AllBound_HF_Visc    += HF_Visc[iMarker];
        AllBound_MaxHF_Visc += pow(MaxHF_Visc[iMarker], MaxNorm);

        /*--- Compute the coefficients per surface ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Visc[iMarker_Monitoring]      += CL_Visc[iMarker];
            Surface_CD_Visc[iMarker_Monitoring]      += CD_Visc[iMarker];
            Surface_CSF_Visc[iMarker_Monitoring] += CSF_Visc[iMarker];
            Surface_CEff_Visc[iMarker_Monitoring]       += CEff_Visc[iMarker];
            Surface_CFx_Visc[iMarker_Monitoring]        += CFx_Visc[iMarker];
            Surface_CFy_Visc[iMarker_Monitoring]        += CFy_Visc[iMarker];
            Surface_CFz_Visc[iMarker_Monitoring]        += CFz_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]        += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]        += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]        += CMz_Visc[iMarker];
            Surface_HF_Visc[iMarker_Monitoring]         += HF_Visc[iMarker];
            Surface_MaxHF_Visc[iMarker_Monitoring]      += pow(MaxHF_Visc[iMarker],MaxNorm);
          }
        }

      }

    }
  }

  /*--- Update some global coeffients ---*/

  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);


#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  MyAllBound_CD_Visc        = AllBound_CD_Visc;                      AllBound_CD_Visc = 0.0;
  MyAllBound_CL_Visc        = AllBound_CL_Visc;                      AllBound_CL_Visc = 0.0;
  MyAllBound_CSF_Visc   = AllBound_CSF_Visc;                 AllBound_CSF_Visc = 0.0;
  AllBound_CEff_Visc = 0.0;
  MyAllBound_CMx_Visc          = AllBound_CMx_Visc;                        AllBound_CMx_Visc = 0.0;
  MyAllBound_CMy_Visc          = AllBound_CMy_Visc;                        AllBound_CMy_Visc = 0.0;
  MyAllBound_CMz_Visc          = AllBound_CMz_Visc;                        AllBound_CMz_Visc = 0.0;
  MyAllBound_CFx_Visc          = AllBound_CFx_Visc;                        AllBound_CFx_Visc = 0.0;
  MyAllBound_CFy_Visc          = AllBound_CFy_Visc;                        AllBound_CFy_Visc = 0.0;
  MyAllBound_CFz_Visc          = AllBound_CFz_Visc;                        AllBound_CFz_Visc = 0.0;
  MyAllBound_CoPx_Visc          = AllBound_CoPx_Visc;                        AllBound_CoPx_Visc = 0.0;
  MyAllBound_CoPy_Visc          = AllBound_CoPy_Visc;                        AllBound_CoPy_Visc = 0.0;
  MyAllBound_CoPz_Visc          = AllBound_CoPz_Visc;                        AllBound_CoPz_Visc = 0.0;
  MyAllBound_CT_Visc           = AllBound_CT_Visc;                         AllBound_CT_Visc = 0.0;
  MyAllBound_CQ_Visc           = AllBound_CQ_Visc;                         AllBound_CQ_Visc = 0.0;
  AllBound_CMerit_Visc = 0.0;
  MyAllBound_HF_Visc     = AllBound_HF_Visc;                   AllBound_HF_Visc = 0.0;
  MyAllBound_MaxHF_Visc  = pow(AllBound_MaxHF_Visc, MaxNorm);  AllBound_MaxHF_Visc = 0.0;

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(&MyAllBound_CD_Visc, &AllBound_CD_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CL_Visc, &AllBound_CL_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CSF_Visc, &AllBound_CSF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
    SU2_MPI::Allreduce(&MyAllBound_CMx_Visc, &AllBound_CMx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMy_Visc, &AllBound_CMy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMz_Visc, &AllBound_CMz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFx_Visc, &AllBound_CFx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFy_Visc, &AllBound_CFy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFz_Visc, &AllBound_CFz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPx_Visc, &AllBound_CoPx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPy_Visc, &AllBound_CoPy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CoPz_Visc, &AllBound_CoPz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CT_Visc, &AllBound_CT_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CQ_Visc, &AllBound_CQ_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
    SU2_MPI::Allreduce(&MyAllBound_HF_Visc, &AllBound_HF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_MaxHF_Visc, &AllBound_MaxHF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);
  }
  
  /*--- Add the forces on the surfaces using all the nodes ---*/

  MySurface_CL_Visc      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Visc      = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Visc = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Visc       = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Visc        = new su2double[config->GetnMarker_Monitoring()];
  MySurface_HF_Visc         = new su2double[config->GetnMarker_Monitoring()];
  MySurface_MaxHF_Visc      = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {

    MySurface_CL_Visc[iMarker_Monitoring]      = Surface_CL_Visc[iMarker_Monitoring];
    MySurface_CD_Visc[iMarker_Monitoring]      = Surface_CD_Visc[iMarker_Monitoring];
    MySurface_CSF_Visc[iMarker_Monitoring] = Surface_CSF_Visc[iMarker_Monitoring];
    MySurface_CEff_Visc[iMarker_Monitoring]       = Surface_CEff_Visc[iMarker_Monitoring];
    MySurface_CFx_Visc[iMarker_Monitoring]        = Surface_CFx_Visc[iMarker_Monitoring];
    MySurface_CFy_Visc[iMarker_Monitoring]        = Surface_CFy_Visc[iMarker_Monitoring];
    MySurface_CFz_Visc[iMarker_Monitoring]        = Surface_CFz_Visc[iMarker_Monitoring];
    MySurface_CMx_Visc[iMarker_Monitoring]        = Surface_CMx_Visc[iMarker_Monitoring];
    MySurface_CMy_Visc[iMarker_Monitoring]        = Surface_CMy_Visc[iMarker_Monitoring];
    MySurface_CMz_Visc[iMarker_Monitoring]        = Surface_CMz_Visc[iMarker_Monitoring];
    MySurface_HF_Visc[iMarker_Monitoring]         = Surface_HF_Visc[iMarker_Monitoring];
    MySurface_MaxHF_Visc[iMarker_Monitoring]      = Surface_MaxHF_Visc[iMarker_Monitoring];

    Surface_CL_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CD_Visc[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring] = 0.0;
    Surface_CEff_Visc[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]        = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]        = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]         = 0.0;
    Surface_MaxHF_Visc[iMarker_Monitoring]      = 0.0;
  }

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(MySurface_CL_Visc, Surface_CL_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CD_Visc, Surface_CD_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CSF_Visc, Surface_CSF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
      Surface_CEff_Visc[iMarker_Monitoring] = Surface_CL_Visc[iMarker_Monitoring] / (Surface_CD_Visc[iMarker_Monitoring] + EPS);
    SU2_MPI::Allreduce(MySurface_CFx_Visc, Surface_CFx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CFy_Visc, Surface_CFy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CFz_Visc, Surface_CFz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMx_Visc, Surface_CMx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMy_Visc, Surface_CMy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_CMz_Visc, Surface_CMz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_HF_Visc, Surface_HF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MySurface_MaxHF_Visc, Surface_MaxHF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  
  delete [] MySurface_CL_Visc; delete [] MySurface_CD_Visc; delete [] MySurface_CSF_Visc;
  delete [] MySurface_CEff_Visc;  delete [] MySurface_CFx_Visc;   delete [] MySurface_CFy_Visc;
  delete [] MySurface_CFz_Visc;   delete [] MySurface_CMx_Visc;   delete [] MySurface_CMy_Visc;
  delete [] MySurface_CMz_Visc; delete [] MySurface_HF_Visc; delete [] MySurface_MaxHF_Visc;

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value)---*/

  Total_CD          += AllBound_CD_Visc;
  Total_CL          += AllBound_CL_Visc;
  Total_CSF         += AllBound_CSF_Visc;
  Total_CEff        = Total_CL / (Total_CD + EPS);
  Total_CMx         += AllBound_CMx_Visc;
  Total_CMy         += AllBound_CMy_Visc;
  Total_CMz         += AllBound_CMz_Visc;
  Total_CFx         += AllBound_CFx_Visc;
  Total_CFy         += AllBound_CFy_Visc;
  Total_CFz         += AllBound_CFz_Visc;
  Total_CoPx        += AllBound_CoPx_Visc;
  Total_CoPy        += AllBound_CoPy_Visc;
  Total_CoPz        += AllBound_CoPz_Visc;
  Total_CT          += AllBound_CT_Visc;
  Total_CQ          += AllBound_CQ_Visc;
  Total_CMerit      = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  Total_Heat        = AllBound_HF_Visc;
  Total_MaxHeat     = AllBound_MaxHF_Visc;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]      += Surface_CL_Visc[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]      += Surface_CD_Visc[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring] += Surface_CSF_Visc[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       = Surface_CL[iMarker_Monitoring] / (Surface_CD[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        += Surface_CFx_Visc[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        += Surface_CFy_Visc[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        += Surface_CFz_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        += Surface_CMz_Visc[iMarker_Monitoring];
  }

}

void CIncNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar;// Wall_Function;
  unsigned long iVertex, iPoint, total_index;
  
  su2double *GridVel, *Normal, Area, Wall_HeatFlux;

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool energy        = config->GetEnergy_Equation();

  /*--- Identify the boundary by string name ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Get the specified wall heat flux from config ---*/
  
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag)/config->GetHeat_Flux_Ref();

//  /*--- Get wall function treatment from config. ---*/
//
//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != NO_WALL_FUNCTION) {
//    SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
//  }

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

      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/
      
      if (dynamic_grid) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      
      nodes->SetVelocity_Old(iPoint,Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);
      
      if (energy) {

        /*--- Apply a weak boundary condition for the energy equation.
        Compute the residual due to the prescribed heat flux. ---*/

        Res_Visc[nDim+1] = Wall_HeatFlux*Area;

        /*--- Viscous contribution to the residual at the wall ---*/

        LinSysRes.SubtractBlock(iPoint, Res_Visc);

      }

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

void CIncNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar, Wall_Function;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  su2double *GridVel;
  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij;
  su2double Twall, dTdn;
  su2double thermal_conductivity;

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool energy        = config->GetEnergy_Equation();

  /*--- Identify the boundary by string name ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Retrieve the specified wall temperature ---*/
  
  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  /*--- Get wall function treatment from config. ---*/

  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  if (Wall_Function != NO_WALL_FUNCTION) {
    SU2_MPI::Error("Wall function treatment not implemented yet.", CURRENT_FUNCTION);
  }

  /*--- Loop over all of the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/
      
      if (dynamic_grid) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      
      nodes->SetVelocity_Old(iPoint,Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);
      
      if (energy) {

        /*--- Compute dual grid area and boundary normal ---*/
        
        Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
        
        Area = 0.0; 
        for (iDim = 0; iDim < nDim; iDim++) 
          Area += Normal[iDim]*Normal[iDim]; 
        Area = sqrt (Area);
        
        /*--- Compute closest normal neighbor ---*/
        
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
        
        /*--- Get coordinates of i & nearest normal and compute distance ---*/
        
        Coord_i = geometry->node[iPoint]->GetCoord();
        Coord_j = geometry->node[Point_Normal]->GetCoord();
        dist_ij = 0;
        for (iDim = 0; iDim < nDim; iDim++)
          dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
        dist_ij = sqrt(dist_ij);

        /*--- Compute the normal gradient in temperature using Twall ---*/
        
        dTdn = -(nodes->GetTemperature(Point_Normal) - Twall)/dist_ij;
        
        /*--- Get thermal conductivity ---*/

        thermal_conductivity = nodes->GetThermalConductivity(iPoint);

        /*--- Apply a weak boundary condition for the energy equation.
        Compute the residual due to the prescribed heat flux. ---*/

        Res_Visc[nDim+1] = thermal_conductivity*dTdn*Area;

        /*--- Jacobian contribution for temperature equation. ---*/
    
        if (implicit) {
          su2double Edge_Vector[3];
          su2double dist_ij_2 = 0, proj_vector_ij = 0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
            dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
            proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
          }
          if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
          else proj_vector_ij = proj_vector_ij/dist_ij_2;

          Jacobian_i[nDim+1][nDim+1] = -thermal_conductivity*proj_vector_ij;

          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }

        /*--- Viscous contribution to the residual at the wall ---*/
        
        LinSysRes.SubtractBlock(iPoint, Res_Visc);

      }

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


void CIncNSSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iVar, jVar, iDim, Wall_Function;
  unsigned long iVertex, iPoint, total_index;

  su2double *GridVel, Tconjugate;
  su2double Temperature_Ref = config->GetTemperature_Ref();

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool energy        = config->GetEnergy_Equation();

  /*--- Identify the boundary ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall function treatment.---*/

  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  if(Wall_Function != NO_WALL_FUNCTION) {
      SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
  }

  /*--- Loop over boundary points ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Initialize the convective & viscous residuals to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/

      if (dynamic_grid) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }

      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/

      nodes->SetVelocity_Old(iPoint,Vector);

      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      nodes->SetVel_ResTruncError_Zero(iPoint);

      if (energy) {

        Tconjugate = GetConjugateHeatVariable(val_marker, iVertex, 0)/Temperature_Ref;

        /*--- Strong imposition of the temperature on the fluid zone. ---*/

        LinSysRes.SetBlock_Zero(iPoint, nDim+1);
        nodes->SetSolution_Old(iPoint, nDim+1, Tconjugate);
        nodes->SetEnergy_ResTruncError_Zero(iPoint);

      }

      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
        if(energy) {
          total_index = iPoint*nVar+nDim+1;
          Jacobian.DeleteValsRowi(total_index);
        }
      }

    }
  }
}
