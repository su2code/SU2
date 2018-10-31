/*!
 * \file solution_direct_mean_inc.cpp
 * \brief Main subroutines for solving incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/
  
  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  
  FluidModel = NULL;
  
  DonorPrimVar = NULL; DonorGlobalIndex = NULL; DonorJacobian = NULL;
 
}

CIncEulerSolver::CIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
  string filename = config->GetSolution_FlowFileName();
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  bool fsi     = config->GetFSI_Simulation();
  string filename_ = config->GetSolution_FlowFileName();

  unsigned short direct_diff = config->GetDirectDiff();

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/

    Read_SU2_Restart_Metadata(geometry, config, false, filename_);
    
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

  /*--- Fixed CL mode initialization (cauchy criteria) ---*/

  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  
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
  
  /*--- Allocate the node variables ---*/
  
  node = new CVariable*[nPoint];
  
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
  
  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/
  
  DonorPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      if (rans) {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar*2];
        for (iVar = 0; iVar < nPrimVar*2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar*2];
        for (iVar = 0; iVar < nPrimVar*2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the Jacobian from the donor for the periodic
   boundary condition (point implicit). ---*/
  
  DonorJacobian = new su2double*** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorJacobian[iMarker] = new su2double** [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorJacobian[iMarker][iVertex] = new su2double* [nPrimVarGrad*nDim];
      for (iVar = 0; iVar < nPrimVarGrad*nDim; iVar++) {
        DonorJacobian[iMarker][iVertex][iVar] = new su2double [nPrimVarGrad*nDim];
        for (unsigned short jVar = 0; jVar < nPrimVarGrad*nDim; jVar++) {
          DonorJacobian[iMarker][iVertex][iVar][jVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
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
  AoA_Prev       = 0.0;
  Total_CL_Prev  = 0.0;    Total_CD_Prev      = 0.0;
  Total_CMx_Prev = 0.0;    Total_CMy_Prev     = 0.0;     Total_CMz_Prev      = 0.0;

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
  
  /*--- Initialize the cauchy critera array for fixed CL mode ---*/

  if (config->GetFixed_CL_Mode())
    Cauchy_Serie = new su2double [config->GetCauchy_Elems()+1];

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CIncEulerVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, config);

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (fsi){
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

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Store volume for periodic if necessary (MG too). ---*/
  
  if (config->GetnMarker_Periodic() != 0) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetAuxVar(geometry->node[iPoint]->GetVolume());
  }
  
  /*--- Perform the MPI communication of the solution ---*/

  Set_MPI_Solution(geometry, config);

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

  if (iPoint_UndLapl != NULL) delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL) delete [] jPoint_UndLapl;

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

  if (DonorPrimVar != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] DonorPrimVar[iMarker][iVertex];
      delete [] DonorPrimVar[iMarker];
    }
    delete [] DonorPrimVar;
  }
  
  if (DonorJacobian != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        for (iVar = 0; iVar < nPrimVarGrad*nDim; iVar++)
          delete [] DonorJacobian[iMarker][iVertex][iVar];
        delete [] DonorJacobian[iMarker][iVertex];
      }
      delete [] DonorJacobian[iMarker];
    }
    delete [] DonorJacobian;
  }
  
  if (DonorGlobalIndex != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] DonorGlobalIndex[iMarker];
    delete [] DonorGlobalIndex;
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
  
  if (nVertex!=NULL) delete [] nVertex;

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
  
  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;
  
  if (FluidModel != NULL) delete FluidModel;
}

void CIncEulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta,
  phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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
      Buffer_Send_U    = new su2double[nBufferS_Vector];
      
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

void CIncEulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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

void CIncEulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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

void CIncEulerSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR, *Buffer_Receive_Neighbor = NULL, *Buffer_Send_Neighbor = NULL;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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

void CIncEulerSolver::Set_MPI_Sensor(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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

void CIncEulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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

void CIncEulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nVar];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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

void CIncEulerSolver::Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
  su2double **Gradient = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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

void CIncEulerSolver::Set_MPI_Primitive_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
  su2double *Limiter = new su2double [nPrimVarGrad];
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
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
  
  unsigned short iDim;
  
  /*--- Local variables ---*/
  
  su2double Mach     = config->GetMach();
  su2double Reynolds = config->GetReynolds();
  
  bool unsteady      = (config->GetUnsteady_Simulation() != NO);
  bool viscous       = config->GetViscous();
  bool grid_movement = config->GetGrid_Movement();
  bool turbulent     = ((config->GetKind_Solver() == RANS) ||
                        (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool tkeNeeded     = ((turbulent) && (config->GetKind_Turb_Model() == SST));
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
      
  }
  
  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
  
  if (viscous) {
    
    /*--- Constant viscosity model ---*/

    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);
    
    /*--- Sutherland's model ---*/    
    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());
    
    /*--- Constant thermal conductivity model ---*/

    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);
    
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
      if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    }
    else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
      cout << "Incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
      cout << "are based on the initial values, p_ref = rho_ref*vel_ref^2." << endl;
      if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    } 
    else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      cout << "Incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
      cout << "are user-provided reference values, p_ref = rho_ref*vel_ref^2." << endl;
      if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
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

    cout <<"-- Input conditions:"<< endl;
    
    switch (config->GetKind_FluidModel()) {
      
      case CONSTANT_DENSITY:
        cout << "Fluid Model: CONSTANT_DENSITY "<< endl;
        if (energy) {
          cout << "Specific heat at constant pressure (Cp): " << config->GetSpecific_Heat_Cp() << " N.m/kg.K." << endl;
        }
        if (boussinesq) cout << "Thermal expansion coefficient: " << config->GetThermal_Expansion_Coeff() << " K^-1." << endl;
        cout << "Thermodynamic pressure not required." << endl;
        break;
        
      case INC_IDEAL_GAS:
        cout << "Fluid Model: INC_IDEAL_GAS "<< endl;
        cout << "Variable density incompressible flow using ideal gas law." << endl;
        cout << "Density is a function of temperature (constant thermodynamic pressure)." << endl;
        cout << "Specific heat at constant pressure (Cp): " << config->GetSpecific_Heat_Cp() << " N.m/kg.K." << endl;
        cout << "Molecular weight : "<< config->GetMolecular_Weight() << " g/mol" << endl;
        cout << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
        cout << "Thermodynamic pressure: " << config->GetPressure_Thermodynamic();
        if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
        break;
      
    }
    if (viscous) {
      switch (config->GetKind_ViscosityModel()) {

        case CONSTANT_VISCOSITY:
          cout << "Viscosity Model: CONSTANT_VISCOSITY  "<< endl;
          cout << "Constant Laminar Viscosity: " << config->GetMu_Constant();
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          break;

        case SUTHERLAND:
          cout << "Viscosity Model: SUTHERLAND "<< endl;
          cout << "Ref. Laminar Viscosity: " << config->GetMu_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
          cout << "Ref. Temperature: " << config->GetMu_Temperature_Ref();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Sutherland Constant: "<< config->GetMu_S();
          if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
          else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
          cout << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< endl;
          cout << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< endl;
          cout << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< endl;
          break;

      }

      if (energy) {
        switch (config->GetKind_ConductivityModel()) {

          case CONSTANT_PRANDTL:
            cout << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
            cout << "Prandtl (Laminar): " << config->GetPrandtl_Lam()<< endl;
            cout << "Prandtl (Turbulent): " << config->GetPrandtl_Turb()<< endl;
            break;

          case CONSTANT_CONDUCTIVITY:
            cout << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
            cout << "Molecular Conductivity: " << config->GetKt_Constant()<< " W/m^2.K." << endl;
            cout << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< endl;
            break;

        }
      }

    }
    
    if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
      cout << "Bulk modulus: " << config->GetBulk_Modulus();
      if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    }

    cout << "Initial dynamic pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    cout << "Initial total pressure: " << config->GetPressure_FreeStream() + 0.5*config->GetDensity_FreeStream()*config->GetModVel_FreeStream()*config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    if (energy) {
      cout << "Initial temperature: " << config->GetTemperature_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " R." << endl;
    }

    cout << "Initial density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;
    
    if (nDim == 2) {
      cout << "Initial velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      cout << "Initial velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      cout << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) cout << " m/s. ";
    else if (config->GetSystemMeasurements() == US) cout << " ft/s. ";
    
    cout << "Magnitude: "  << config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;

    if (viscous) {
      cout << "Initial laminar viscosity: " << config->GetViscosity_FreeStream();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      if (turbulent) {
        cout << "Initial turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;
        cout << "Initial specific dissipation: " << config->GetOmega_FreeStream();
        if (config->GetSystemMeasurements() == SI) cout << " 1/s." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " 1/s." << endl;
      }
    }
    
    if (unsteady) { cout << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << endl; }
    
    /*--- Print out reference values. ---*/
    
    cout <<"-- Reference values:"<< endl;

    if (config->GetKind_FluidModel() != CONSTANT_DENSITY) {
      cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
    } else {
      if (energy) {
        cout << "Reference specific heat: " << config->GetGas_Constant_Ref();
        if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
        else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
      }
    }

    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    
    if (energy) {
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

    cout << "Reference length: " << config->GetLength_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " in." << endl;

    if (viscous) {
      cout << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
    }
    
    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;
    
    /*--- Print out resulting non-dim values here. ---*/
    
    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;
    if (viscous) {
      cout << "Reynolds number (per m): " << config->GetReynolds() << endl;
    }
    
    if (config->GetKind_FluidModel() != CONSTANT_DENSITY) {
      cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
      cout << "Initial thermodynamic pressure (non-dim): " << config->GetPressure_ThermodynamicND() << endl;
    } else {
      if (energy) {
        cout << "Specific heat at constant pressure (non-dim): " << config->GetSpecific_Heat_CpND() << endl;
        if (boussinesq) cout << "Thermal expansion coefficient (non-dim.): " << config->GetThermal_Expansion_CoeffND() << " K^-1." << endl;
      }
    }
    
    if (energy) cout << "Initial temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;
    cout << "Initial pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;
    cout << "Initial density (non-dim): " << config->GetDensity_FreeStreamND() << endl;
    
    if (nDim == 2) {
      cout << "Initial velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Initial velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "   << config->GetModVel_FreeStreamND() << endl;
    
    if (viscous) {
      cout << "Initial viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << endl;
      if (turbulent) {
        cout << "Initial turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << endl;
        cout << "Initial specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << endl;
      }
    }
    
    if (unsteady) {
      cout << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << endl;
      cout << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << endl;
    }
    
    cout << endl;
    
  }
  
}

void CIncEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, *Solution_Fine, *Solution;
    
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool rans      = ((config->GetKind_Solver() == RANS) ||
                    (config->GetKind_Solver() == ADJ_RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  if ((ExtIter == 0) && config->GetTaylorGreen()) {

    /* Write a message that the solution is initialized for the Taylor-Green vortex
     test case. */
    if(rank == MASTER_NODE) {
      cout << endl;
      cout << "Warning: Solution is initialized for the Taylor-Green vortex test case!!!" << endl;
      cout << endl << flush;
    }

    /* The initial conditions are set for the Taylor-Green vortex case, which
     is a DNS case that features vortex breakdown into turbulence. These
     particular settings are for the typical Re = 1600 case (M = 0.08) with
     an initial temperature of 300 K. Note that this condition works in both
     2D and 3D. */

    const su2double tgvLength    = 1.0;     // Taylor-Green length scale.
    const su2double tgvVelocity  = 1.0;     // Taylor-Green velocity.
    const su2double tgvDensity   = 1.0;     // Taylor-Green density.
    const su2double tgvViscosity = 0.1;   // Taylor-Green pressure.

    /* Useful coefficient in which Gamma is present. */
    const su2double ovGm1    = 1.0/Gamma_Minus_One;

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {

        /* Set the pointers to the coordinates and solution of this DOF. */
        const su2double *coor = geometry[iMesh]->node[iPoint]->GetCoord();
        su2double *solDOF     = solver_container[iMesh][FLOW_SOL]->node[iPoint]->GetSolution();

        su2double coorZ = 0.0;
        if (nDim == 3) coorZ = coor[2];

        /* Compute the primitive variables. */
        su2double rho = tgvDensity;
        su2double u   =  tgvVelocity * (sin(coor[0]/tgvLength)*
                                        cos(coor[1]/tgvLength)*
                                        cos(coorZ  /tgvLength));
        su2double v   = -tgvVelocity * (cos(coor[0]/tgvLength)*
                                        sin(coor[1]/tgvLength)*
                                        cos(coorZ  /tgvLength));
        su2double factorA = 1.0;
        su2double factorB = cos(2.0*coor[0]/tgvLength) + cos(2.0*coor[1]/tgvLength);
        su2double p   = (tgvDensity/4.0)*factorA*factorB;

        /* Compute the conservative variables. Note that both 2D and 3D
         cases are treated correctly. */
        solDOF[0]      = p;
        solDOF[1]      = u;
        solDOF[2]      = v;
        solDOF[3]      = 0.0;
        solDOF[nVar-1] = 0.0;
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
      
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1), true);
      
      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter()-1), false);
      
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

void CIncEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long ErrorCounter = 0;

  unsigned long ExtIter = config->GetExtIter();
  bool cont_adjoint     = config->GetContinuous_Adjoint();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
  bool limiter          = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool center           = ((config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED));
  bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
  bool fixed_cl         = config->GetFixed_CL_Mode();
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;

  /*--- Store the original volume for periodic cells on the boundaries,
   since this will be increased as we complete the CVs during our
   updates. ---*/
  
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Reset to the original volume. ---*/
        geometry->node[iPoint]->SetVolume(node[iPoint]->GetAuxVar());
        
      }
    }
  }
  
  /*--- Update the angle of attack at the far-field for fixed CL calculations (only direct problem). ---*/
  
  if ((fixed_cl) && (!disc_adjoint) && (!cont_adjoint)) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }

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

  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !disc_adjoint) Jacobian.SetValZero();

  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
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
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Incompressible flow, primitive variables ---*/
    
    physical = node[iPoint]->SetPrimVar(FluidModel);

    /*--- Record any non-physical points. ---*/

    if (!physical) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
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
  bool grid_movement = config->GetGrid_Movement();
  bool time_steping  = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool dual_time     = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
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
    
    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);
    
    /*--- Mean Values ---*/

    Mean_ProjVel    = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_BetaInc2   = 0.5 * (node[iPoint]->GetBetaInc2()      + node[jPoint]->GetBetaInc2());
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

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

      Mean_ProjVel    = node[iPoint]->GetProjVel(Normal);
      Mean_BetaInc2   = node[iPoint]->GetBetaInc2();
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

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
  }
  
  /*--- Local time-stepping: each element uses their own speed for steady state 
   simulations or for pseudo time steps in a dual time simulation. ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time  = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time    = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time    = max(Max_Delta_Time, Local_Delta_Time);
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
        node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
      } else {
        node[iPoint]->SetDelta_Time(Global_Delta_Time);
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
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CIncEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool jst_scheme  = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
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
    
    /*--- Set undivided laplacian and pressure-based sensor ---*/
    
    if (jst_scheme) {
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
  
  unsigned long ExtIter = config->GetExtIter();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter          = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool grid_movement    = config->GetGrid_Movement();
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;

  /*--- Loop over all the edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Grid movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    /*--- Get primitive variables ---*/
    
    V_i = node[iPoint]->GetPrimitive(); V_j = node[jPoint]->GetPrimitive();
    S_i = node[iPoint]->GetSecondary(); S_j = node[jPoint]->GetSecondary();

    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (muscl) {
      
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
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
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

      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());

      /*--- Set incompressible density  ---*/
      
      numerics->SetDensity(node[iPoint]->GetDensity(),
                           node[iPoint]->GetDensity());

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

      numerics->SetConservative(node[iPoint]->GetSolution(),
                                node[iPoint]->GetSolution());

      /*--- Set incompressible density  ---*/
      
      numerics->SetDensity(node[iPoint]->GetDensity(),
                           node[iPoint]->GetDensity());

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

    /*--- For viscous problems, we need an additional gradient. ---*/

    if (viscous) {

      su2double AuxVar, Total_Viscosity, yCoord, yVelocity;

      for (iPoint = 0; iPoint < nPoint; iPoint++) {

        yCoord          = geometry->node[iPoint]->GetCoord(1);
        yVelocity       = node[iPoint]->GetVelocity(1);
        Total_Viscosity = (node[iPoint]->GetLaminarViscosity() +
                           node[iPoint]->GetEddyViscosity());

        if (yCoord > EPS) {
          AuxVar = Total_Viscosity*yVelocity/yCoord;
        } else {
          AuxVar = 0.0;
        }

        /*--- Set the auxilairy variable for this node. ---*/

        node[iPoint]->SetAuxVar(AuxVar);

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

      numerics->SetPrimitive(node[iPoint]->GetPrimitive(), NULL);

      /*--- Set incompressible density  ---*/
      
      numerics->SetDensity(node[iPoint]->GetDensity(),
                           node[iPoint]->GetDensity());

      /*--- Set control volume ---*/
      
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                         geometry->node[iPoint]->GetCoord());

      /*--- If viscous, we need gradients for extra terms. ---*/

      if (viscous) {

        /*--- Gradient of the primitive variables ---*/

        numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), NULL);

        /*--- Load the aux variable gradient that we already computed. ---*/

        numerics->SetAuxVarGrad(node[iPoint]->GetAuxVarGradient(), NULL);
        
      }

      /*--- Compute Source term Residual ---*/
      
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}

void CIncEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */
  
  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = node[iPoint]->GetSolution() */
  /* Get the volume (or area in 2D) associated with iPoint = node[iPoint]->GetVolume() */
  /* Get the vector of geometric coordinates of point iPoint = node[iPoint]->GetCoord() */
  
}

void CIncEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
  
  su2double *Normal, Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0,
  Mean_BetaInc2, Lambda, ProjVel, ProjVel_i, ProjVel_j, *GridVel, *GridVel_i, *GridVel_j;
  
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

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
    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);
    
    /*--- Mean Values ---*/

    Mean_ProjVel    = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_BetaInc2   = 0.5 * (node[iPoint]->GetBetaInc2()      + node[jPoint]->GetBetaInc2());
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

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
      
      Mean_ProjVel    = node[iPoint]->GetProjVel(Normal);
      Mean_BetaInc2   = node[iPoint]->GetBetaInc2();
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);
      
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
  }
  
  /*--- Correct the gradient values for any periodic boundaries. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    BC_Periodic_Eigenvalue(geometry, config, iPeriodic);
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_MaxEigenvalue(geometry, config);
  
}

void CIncEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint, jPoint, iEdge;
  su2double *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;
  
  Diff = new su2double[nVar];
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Solution differences ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
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
  
  /*--- Correct the undivied Laplacian values for any periodic boundaries. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    BC_Periodic_Laplacian(geometry, config, iPeriodic);
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Undivided_Laplacian(geometry, config);
  
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

    Pressure_i = node[iPoint]->GetDensity();
    Pressure_j = node[jPoint]->GetDensity();

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
  
  /*--- Set pressure switch for each point ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetSensor(fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
  
  /*--- Correct the sensor values for any periodic boundaries. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    BC_Periodic_Sensor(geometry, config, iPeriodic);
  }
  
  /*--- MPI parallelization ---*/
  
  Set_MPI_Sensor(geometry, config);
  
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

        Pressure = node[iPoint]->GetPressure();

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
          Density   = node[iPoint]->GetDensity();

          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);

          MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim]   = node[iPoint]->GetVelocity(iDim);
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
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        Res = 0.0;
        for (jVar = 0; jVar < nVar; jVar ++ )
          Res += Preconditioner[iVar][jVar]*(Residual[jVar] + Res_TruncError[jVar]);
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
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);


    if (!adjoint) {
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        Res = 0.0;
        for (jVar = 0; jVar < nVar; jVar ++ )
          Res += Preconditioner[iVar][jVar]*(local_Residual[jVar] + local_Res_TruncError[jVar]);
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
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Apply the preconditioner and add to the diagonal. ---*/
    
    if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
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

void CIncEulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal;
  
  /*--- Incompressible flow, primitive variables nDim+4, (P, vx, vy, vz, T, rho, beta) ---*/
  
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
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
       (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
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
  
  /*--- Correct the gradient values for any periodic boundaries. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    BC_Periodic_GG(geometry, config, iPeriodic);
  }
  
  Set_MPI_Primitive_Gradient(geometry, config);
  
}

void CIncEulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;
  
  /*--- Incompressible flow, primitive variables nDim+4, (P, vx, vy, vz, T, rho, beta) ---*/
  
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
        Cvector[iVar][iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
    
    AD::StartPreacc();
    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    AD::SetPreaccIn(Coord_i, nDim);
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      PrimVar_j = node[jPoint]->GetPrimitive();
      
      AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

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
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
        
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
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        }
        
        node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
      }
    }
    
    AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    AD::EndPreacc();
  }
  
  /*--- Correct the gradient values for any periodic boundaries. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    BC_Periodic_LS(geometry, config, iPeriodic);
  }
  
  Set_MPI_Primitive_Gradient(geometry, config);
  
}

void CIncEulerSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j,
  *Primitive, *Primitive_i, *Primitive_j, *LocalMinPrimitive, *LocalMaxPrimitive,
  *GlobalMinPrimitive, *GlobalMaxPrimitive,
  dave, LimK, eps2, eps1, dm, dp, du, y, limiter;
  
  dave = config->GetRefElemLength();
  LimK = config->GetVenkat_LimiterCoeff();

  if (config->GetKind_SlopeLimit_Flow() == NO_LIMITER) {
   
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        node[iPoint]->SetLimiter_Primitive(iVar, 1.0);
      }
    }
    
  }
  
  else {
    
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
    
    /*--- Correct min/max values for periodic boundaries. ---*/
    
    for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
      BC_Periodic_Limiter1(geometry, config, iPeriodic);
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
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

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
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }
        
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
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }
        
      }
      
      AD::EndPreacc();
      
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
  
  if ((config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) ||
      (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG)) {
    
    /*--- Allocate memory for the max and min primitive value --*/
    
    LocalMinPrimitive = new su2double [nPrimVarGrad]; GlobalMinPrimitive = new su2double [nPrimVarGrad];
    LocalMaxPrimitive = new su2double [nPrimVarGrad]; GlobalMaxPrimitive = new su2double [nPrimVarGrad];
    
    /*--- Compute the max value and min value of the solution ---*/
    
    Primitive = node[0]->GetPrimitive();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      LocalMinPrimitive[iVar] = Primitive[iVar];
      LocalMaxPrimitive[iVar] = Primitive[iVar];
    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Get the primitive variables ---*/
      
      Primitive = node[iPoint]->GetPrimitive();

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
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Primitive();
      Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();

      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        
        if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
          eps1 = LimK * (GlobalMaxPrimitive[iVar] - GlobalMinPrimitive[iVar]);
          eps2 = eps1*eps1;
        }
        else {
          eps1 = LimK*dave;
          eps2 = eps1*eps1*eps1;
        }
        
        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)){
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)){
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }
        
      }

      AD::EndPreacc();
      
    }
    
    delete [] LocalMinPrimitive; delete [] GlobalMinPrimitive;
    delete [] LocalMaxPrimitive; delete [] GlobalMaxPrimitive;

  }
  
  /*--- Correct the limiter for periodic boundaries to make sure the minimum is chosen. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    BC_Periodic_Limiter2(geometry, config, iPeriodic);
  }
  
  /*--- Limiter MPI ---*/
  
  Set_MPI_Primitive_Limiter(geometry, config);
  
}

void CIncEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output) {
  
  su2double Target_CL = 0.0, AoA = 0.0, Vel_Infty[3] = {0.0,0.0,0.0},
  AoA_inc = 0.0, Vel_Infty_Mag, Old_AoA;
  
  unsigned short iDim;
  
  unsigned long Iter_Fixed_CL = config->GetIter_Fixed_CL();
  unsigned long Update_Alpha = config->GetUpdate_Alpha();
  
  unsigned long ExtIter       = config->GetExtIter();
  bool write_heads = ((ExtIter % Iter_Fixed_CL == 0) && (ExtIter != 0));
  su2double Beta                 = config->GetAoS()*PI_NUMBER/180.0;
  su2double dCL_dAlpha           = config->GetdCL_dAlpha()*180.0/PI_NUMBER;
  bool Update_AoA             = false;
  
  if (ExtIter == 0) AoA_Counter = 0;
  
  /*--- Only the fine mesh level should check the convergence criteria ---*/
  
  if ((iMesh == MESH_0) && Output) {
    
    /*--- Initialize the update flag to false ---*/
    
    Update_AoA = false;
    
    /*--- Reevaluate Angle of Attack at a fixed number of iterations ---*/
    
    if ((ExtIter % Iter_Fixed_CL == 0) && (ExtIter != 0)) {
      AoA_Counter++;
      if ((AoA_Counter <= Update_Alpha)) Update_AoA = true;
      Update_AoA = true;
    }
    
    /*--- Store the update boolean for use on other mesh levels in the MG ---*/
    
    config->SetUpdate_AoA(Update_AoA);
    
  } else {
    Update_AoA = config->GetUpdate_AoA();
  }
  
  if (Update_AoA && Output) {
    
    /*--- Retrieve the specified target CL value. ---*/
    
    Target_CL = config->GetTarget_CL();
    
    /*--- Retrieve the old AoA (radians) ---*/
    
    AoA_old = config->GetAoA()*PI_NUMBER/180.0;
    
    /*--- Estimate the increment in AoA based on dCL_dAlpha (radians) ---*/
    
    AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - Total_CL);
    
    /*--- Compute a new value for AoA on the fine mesh only (radians)---*/
    
    if (iMesh == MESH_0) AoA = AoA_old + AoA_inc;
    else { AoA = config->GetAoA()*PI_NUMBER/180.0; }
    
    /*--- Only the fine mesh stores the updated values for AoA in config ---*/
    
    if (iMesh == MESH_0) {
      config->SetAoA(AoA*180.0/PI_NUMBER);
    }
    
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
    
    /*--- Only the fine mesh stores the updated values for velocity in config ---*/
    
    if (iMesh == MESH_0) {
      for (iDim = 0; iDim < nDim; iDim++)
        config->SetVelocity_FreeStreamND(Vel_Infty[iDim], iDim);
    }
    
    /*--- Output some information to the console with the headers ---*/
    
    if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && !config->GetDiscrete_Adjoint()) {
      Old_AoA = config->GetAoA() - AoA_inc*(180.0/PI_NUMBER);
      
      cout.precision(7);
      cout.setf(ios::fixed, ios::floatfield);
      cout << endl << "----------------------------- Fixed CL Mode -----------------------------" << endl;
      cout << "CL: " << Total_CL;
      cout << " (target: " << config->GetTarget_CL() <<")." << endl;
      cout.precision(4);
      cout << "Previous AoA: " << Old_AoA << " deg";
      cout << ", new AoA: " << config->GetAoA() << " deg." << endl;

      cout << "-------------------------------------------------------------------------" << endl << endl;
    }
    
  }
  
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
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCD_dCL()*(Surface_CL[iMarker_Monitoring]);
        if (config->GetFixed_CM_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCD_dCMy()*(Surface_CMy[iMarker_Monitoring]);
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
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMx_dCL()*(Surface_CL[iMarker_Monitoring]);
        break;
      case MOMENT_Y_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMy[iMarker_Monitoring]);
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMy_dCL()*(Surface_CL[iMarker_Monitoring]);
        break;
      case MOMENT_Z_COEFFICIENT:
        Total_ComboObj+=Weight_ObjFunc*(Surface_CMz[iMarker_Monitoring]);
        if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMz_dCL()*(Surface_CL[iMarker_Monitoring]);
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
      
      if (node[iPoint]->GetVelocity2() > maxVel2)
        maxVel2 = node[iPoint]->GetVelocity2();
      
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
    node[iPoint]->SetBetaInc2(epsilon2*config->GetMax_Vel2());

}

void CIncEulerSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) {

  unsigned short iDim, jDim;

  su2double  BetaInc2, Density, dRhodT, Temperature, oneOverCp, Cp;
  su2double  Velocity[3] = {0.0,0.0,0.0};

  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool energy           = config->GetEnergy_Equation();

  /*--- Access the primitive variables at this node. ---*/

  Density     = node[iPoint]->GetDensity();
  BetaInc2    = node[iPoint]->GetBetaInc2();
  Cp          = node[iPoint]->GetSpecificHeatCp();
  oneOverCp   = 1.0/Cp;
  Temperature = node[iPoint]->GetTemperature();

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = node[iPoint]->GetVelocity(iDim);

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

void CIncEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar;
  unsigned long iPoint, iVertex;

  su2double Density = 0.0, Pressure = 0.0, *Normal = NULL, Area, *NormalArea, turb_ke;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool tkeNeeded = (((config->GetKind_Solver() == RANS ) ||
                     (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  
  Normal     = new su2double[nDim];
  NormalArea = new su2double[nDim];

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
      }

      /*--- Compute the residual ---*/

      Pressure = node[iPoint]->GetPressure();
      Density  = node[iPoint]->GetDensity();

      Residual[0] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual[iDim+1] = Pressure*NormalArea[iDim];
      Residual[nDim+1] = 0.0;

      /*--- Add the Reynolds stress tensor contribution ---*/

      if (tkeNeeded) {
        turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] += (2.0/3.0)*Density*turb_ke*NormalArea[iDim];
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
        
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_i[iDim+1][0] = -Normal[iDim];
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      }
    }
  }
  
  delete [] Normal;
  delete [] NormalArea;
  
}

void CIncEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  
  su2double *V_infty, *V_domain;
  
  bool implicit      = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool grid_movement = config->GetGrid_Movement();
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
      
      V_domain = node[iPoint]->GetPrimitive();

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

      V_infty[nDim+3] = node[iPoint]->GetBetaInc2();

      /*--- Cp is needed for Temperature equation. ---*/

      V_infty[nDim+7] = node[iPoint]->GetSpecificHeatCp();

      /*--- Set various quantities in the numerics class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      if (grid_movement)
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

        V_infty[nDim+4] = node[iPoint]->GetLaminarViscosity();
        V_infty[nDim+5] = node[iPoint]->GetEddyViscosity();
        V_infty[nDim+6] = node[iPoint]->GetThermalConductivity();

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

void CIncEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  unsigned long Point_Normal;
  su2double *Flow_Dir, Flow_Dir_Mag, Vel_Mag, Area;
  su2double *V_inlet, *V_domain;
  su2double UnitFlowDir[3] = {0.0,0.0,0.0};
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
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
      
      V_domain = node[iPoint]->GetPrimitive();

      /*--- The velocity is either prescribed or computed from total pressure. ---*/

      switch (Kind_Inlet) {

        /*--- Velocity and temperature (if required) been specified at the inlet. ---*/

        case VELOCITY_INLET:

        /*--- Retrieve the specified velocity and temperature for the inlet. ---*/
          
        Vel_Mag  = Inlet_Ptotal[val_marker][iVertex]/config->GetVelocity_Ref();

        /*--- Store the velocity in the primitive variable vector. ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim+1] = Vel_Mag*UnitFlowDir[iDim];

        break;

      }
      
      /*--- Neumann condition for dynamic pressure ---*/

      V_inlet[0] = node[iPoint]->GetPressure();

      /*--- Dirichlet condition for temperature (if energy is active) ---*/
      
      V_inlet[nDim+1] = Inlet_Ttotal[val_marker][iVertex]/config->GetTemperature_Ref();

      /*--- Access density at the node. This is either constant by
        construction, or will be set fixed implicitly by the temperature
        and equation of state. ---*/

      V_inlet[nDim+2] = node[iPoint]->GetDensity();

      /*--- Beta coefficient from the config file ---*/

      V_inlet[nDim+3] = node[iPoint]->GetBetaInc2();

      /*--- Cp is needed for Temperature equation. ---*/

      V_inlet[nDim+7] = node[iPoint]->GetSpecificHeatCp();

      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
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
        
        V_inlet[nDim+4] = node[iPoint]->GetLaminarViscosity();
        V_inlet[nDim+5] = node[iPoint]->GetEddyViscosity();
        V_inlet[nDim+6] = node[iPoint]->GetThermalConductivity();

        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                          node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                              solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
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
  su2double *V_outlet, *V_domain, P_Outlet;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool viscous       = config->GetViscous();
  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);

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
       
      /*--- Current solution at this boundary node ---*/
      
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Retrieve the specified back pressure for this outlet. ---*/

      P_Outlet = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();

      /*--- The pressure is prescribed at the outlet. ---*/
      
      V_outlet[0] = P_Outlet;

      /*--- Neumann condition for the velocity ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {
        V_outlet[iDim+1] = node[iPoint]->GetPrimitive(iDim+1);
      }

      /*--- Neumann condition for the temperature. ---*/

      V_outlet[nDim+1] = node[iPoint]->GetTemperature();

      /*--- Access density at the interior node. This is either constant by
        construction, or will be set fixed implicitly by the temperature
        and equation of state. ---*/
      
      V_outlet[nDim+2] = node[iPoint]->GetDensity();

      /*--- Beta coefficient from the config file ---*/
      
      V_outlet[nDim+3] = node[iPoint]->GetBetaInc2();

      /*--- Cp is needed for Temperature equation. ---*/

      V_outlet[nDim+7] = node[iPoint]->GetSpecificHeatCp();

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
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

        V_outlet[nDim+4] = node[iPoint]->GetLaminarViscosity();
        V_outlet[nDim+5] = node[iPoint]->GetEddyViscosity();
        V_outlet[nDim+6] = node[iPoint]->GetThermalConductivity();

        /*--- Set the normal vector and the coordinates ---*/
        
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                          node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                              solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
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

void CIncEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler wall residual method. ---*/
  
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
  
}

void CIncEulerSolver::BC_Periodic_GG(CGeometry *geometry, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, jVar, iMarker, jMarker, iPeriodic, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
  jacMatrix[15][5] ,
  rotJacob[15][5];
  
  unsigned long buff_size = nPrimVarGrad*nDim + 5;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [buff_size];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                for (jVar = 0; jVar < nDim; jVar++) {
                  jacMatrix[iVar][jVar]  = node[iPoint]->GetAuxVar()* node[iPoint]->GetGradient_Primitive(iVar, jVar);
                  rotJacob[iVar][jVar]  = node[iPoint]->GetAuxVar()* node[iPoint]->GetGradient_Primitive(iVar, jVar);
                  
                }
              }
              
              /*--- Rotate the Jacobian momentum components. ---*/
              
//              for (unsigned short iVar = 0; iVar < nPrimVarGrad; iVar++) {
//
//                if (nDim == 2) {
//                  rotJacob[iVar][0] = (rotMatrix[0][0]*jacMatrix[iVar][0] +
//                                       rotMatrix[0][1]*jacMatrix[iVar][1]);
//                  rotJacob[iVar][1] = (rotMatrix[1][0]*jacMatrix[iVar][0] +
//                                       rotMatrix[1][1]*jacMatrix[iVar][1]);
//                } else {
//
//                  rotJacob[iVar][0] = (rotMatrix[0][0]*jacMatrix[iVar][0] +
//                                       rotMatrix[0][1]*jacMatrix[iVar][1] +
//                                       rotMatrix[0][2]*jacMatrix[iVar][2]);
//                  rotJacob[iVar][1] = (rotMatrix[1][0]*jacMatrix[iVar][0] +
//                                       rotMatrix[1][1]*jacMatrix[iVar][1] +
//                                       rotMatrix[1][2]*jacMatrix[iVar][2]);
//                  rotJacob[iVar][2] = (rotMatrix[2][0]*jacMatrix[iVar][0] +
//                                       rotMatrix[2][1]*jacMatrix[iVar][1] +
//                                       rotMatrix[2][2]*jacMatrix[iVar][2]);
//                }
//              }
              
              int ii = 0;
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                for (jVar = 0; jVar < nDim; jVar++) {
                  Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+ii] = rotJacob[iVar][jVar];
                  ii++;
                }
              }
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+3)]  = node[iPoint]->GetAuxVar();
              
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+2)]);
        
        iPrimVar[0] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nPrimVarGrad*nDim+3];
        SetDonorPrimVar(iMarker, iVertex, 0, iPrimVar[0]);
        
        int ii = 0;
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          for (jVar = 0; jVar < nDim; jVar++) {
            SetDonorJacobian(iMarker, iVertex, iVar, jVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+ii]);
            ii++;
          }
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+2)]);
        
        iPrimVar[0] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nPrimVarGrad*nDim+3];
        SetDonorPrimVar(iMarker, iVertex, 0, iPrimVar[0]);
        
        int ii = 0;
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          for (jVar = 0; jVar < nDim; jVar++) {
            SetDonorJacobian(iMarker, iVertex, iVar, jVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+ii]);
            ii++;
          }
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            /*--- Update the volume and time step using the donor info. ---*/
            
            
            su2double Vol = geometry->node[iPoint]->GetVolume() + GetDonorPrimVar(iMarker, iVertex, 0);
            
            
            /*--- Access the Jacobian from the donor. ---*/
            
            su2double Grad[20][10];
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (jVar = 0; jVar < nDim; jVar++) {
                Grad[iVar][jVar] = (GetDonorJacobian(iMarker, iVertex, iVar, jVar) +
                                    node[iPoint]->GetAuxVar()* node[iPoint]->GetGradient_Primitive(iVar, jVar))/Vol;
              }
            }
            
            // Store grad full
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (jVar = 0; jVar < nDim; jVar++) {
                node[iPoint]->SetGradient_Primitive(iVar,jVar,Grad[iVar][jVar]);
              }
            }
            
            
          }
        }
      }
    }
  }
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Periodic_LS(CGeometry *geometry, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, jPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iDim, jDim, iNeigh, iVar, jVar, iMarker, jMarker, iPeriodic, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
  jacMatrix[15][5] ,
  rotJacob[15][5];
  
  unsigned long buff_size = nPrimVarGrad*nDim + 13;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [buff_size];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12,
              r13, r22, r23, r23_a, r23_b, r33, weight;
              
              /*--- Get coordinates ---*/
              
              Coord_i = geometry->node[iPoint]->GetCoord();
              
              /*--- Get primitives from CVariable ---*/
              
              PrimVar_i = node[iPoint]->GetPrimitive();
              
              /*--- Inizialization of variables ---*/
              
              for (iVar = 0; iVar < nPrimVarGrad; iVar++)
                for (iDim = 0; iDim < nDim; iDim++)
                  Cvector[iVar][iDim] = 0.0;
              
              r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
              r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
              
              AD::StartPreacc();
              AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
              AD::SetPreaccIn(Coord_i, nDim);
              
              for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
                jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
                
                /*--- avoid halos and boundary points so that we don't duplicate edges ---*/
                
                if (geometry->node[jPoint]->GetDomain() && (!geometry->node[jPoint]->GetBoundary())) {
                  
                  Coord_j = geometry->node[jPoint]->GetCoord();
                  
                  PrimVar_j = node[jPoint]->GetPrimitive();
                  
                  AD::SetPreaccIn(Coord_j, nDim);
                  AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);
                  
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
                        Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
                    
                  }
                }
              }
              
              /*--- Store the Cvector for this point. ---*/
              
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                for (jVar = 0; jVar < nDim; jVar++) {
                  jacMatrix[iVar][jVar]  = Cvector[iVar][jVar];
                  rotJacob[iVar][jVar]   = Cvector[iVar][jVar];
                  
                }
              }
              
              /*--- Rotate the Jacobian momentum components. ---*/
              
              for (unsigned short iVar = 0; iVar < nPrimVarGrad; iVar++) {
                
                if (nDim == 2) {
                  rotJacob[iVar][0] = (rotMatrix[0][0]*jacMatrix[iVar][0] +
                                       rotMatrix[0][1]*jacMatrix[iVar][1]);
                  rotJacob[iVar][1] = (rotMatrix[1][0]*jacMatrix[iVar][0] +
                                       rotMatrix[1][1]*jacMatrix[iVar][1]);
                } else {
                  
                  rotJacob[iVar][0] = (rotMatrix[0][0]*jacMatrix[iVar][0] +
                                       rotMatrix[0][1]*jacMatrix[iVar][1] +
                                       rotMatrix[0][2]*jacMatrix[iVar][2]);
                  rotJacob[iVar][1] = (rotMatrix[1][0]*jacMatrix[iVar][0] +
                                       rotMatrix[1][1]*jacMatrix[iVar][1] +
                                       rotMatrix[1][2]*jacMatrix[iVar][2]);
                  rotJacob[iVar][2] = (rotMatrix[2][0]*jacMatrix[iVar][0] +
                                       rotMatrix[2][1]*jacMatrix[iVar][1] +
                                       rotMatrix[2][2]*jacMatrix[iVar][2]);
                }
              }
              
              int ii = 0;
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                for (jVar = 0; jVar < nDim; jVar++) {
                  Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+ii] = rotJacob[iVar][jVar];
                  ii++;
                }
              }
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+3)]  = node[iPoint]->GetAuxVar();
              
              // send the other matrix terms
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+4)]  = r11;
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+5)]  = r12;
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+6)]  = r22;
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+7)]  = r13;
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+8)]  = r23_a;
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+9)]  = r23_b;
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVarGrad*nDim+10)]  = r33;
              
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+2)]);
        
        for (iVar = 0; iVar < 7; iVar++) {
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nPrimVarGrad*nDim+4+iVar];
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        }
        
        int ii = 0;
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          for (jVar = 0; jVar < nDim; jVar++) {
            SetDonorJacobian(iMarker, iVertex, iVar, jVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+ii]);
            ii++;
          }
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVarGrad*nDim+2)]);
        
        for (iVar = 0; iVar < 7; iVar++) {
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nPrimVarGrad*nDim+4+iVar];
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        }
        
        int ii = 0;
        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          for (jVar = 0; jVar < nDim; jVar++) {
            SetDonorJacobian(iMarker, iVertex, iVar, jVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+ii]);
            ii++;
          }
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
            r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
            bool singular;
            
            /*--- Set the value of the singular ---*/
            singular = false;
            
            /*--- Get coordinates ---*/
            
            Coord_i = geometry->node[iPoint]->GetCoord();
            
            /*--- Get primitives from CVariable ---*/
            
            PrimVar_i = node[iPoint]->GetPrimitive();
            
            /*--- Inizialization of variables ---*/
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++)
              for (iDim = 0; iDim < nDim; iDim++)
                Cvector[iVar][iDim] = 0.0;
            
            r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
            r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
            
            AD::StartPreacc();
            AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
            AD::SetPreaccIn(Coord_i, nDim);
            
            for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
              
              /*--- avoid halos and boundary points so that we don't duplicate edges ---*/
              
              Coord_j = geometry->node[jPoint]->GetCoord();
              
              PrimVar_j = node[jPoint]->GetPrimitive();
              
              AD::SetPreaccIn(Coord_j, nDim);
              AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);
              
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
                    Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
                
              }
            }
            
            // add the Cvector for this point from donor before computing final grad
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (jVar = 0; jVar < nDim; jVar++) {
                Cvector[iVar][jVar] += GetDonorJacobian(iMarker, iVertex, iVar, jVar);
              }
            }
            
            // get the other matrix entries
            
            r11 += GetDonorPrimVar(iMarker, iVertex, 0);
            r12 += GetDonorPrimVar(iMarker, iVertex, 1);
            r22 += GetDonorPrimVar(iMarker, iVertex, 2);
            
            if (nDim == 3) {
              r13 += GetDonorPrimVar(iMarker, iVertex, 3);
              r23_a += GetDonorPrimVar(iMarker, iVertex, 4);
              r23_b += GetDonorPrimVar(iMarker, iVertex, 5);
              r33 += GetDonorPrimVar(iMarker, iVertex, 6);
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
                  product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
                }
                
                node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
              }
            }
            
          }
        }
      }
    }
  }
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Periodic_Limiter1(CGeometry *geometry, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker, iPeriodic, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  
  unsigned long buff_size = nPrimVar*2 + 5;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [buff_size];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              for (iVar = 0; iVar < nPrimVar; iVar++) {
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetSolution_Min(iVar);
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+nPrimVar+iVar] = node[iPoint]->GetSolution_Max(iVar);
              }
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+3)]  = node[iPoint]->GetAuxVar();
              
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        for (iVar = 0; iVar < nPrimVar; iVar++) {
          SetDonorPrimVar(iMarker, iVertex, iVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar]);
          SetDonorPrimVar(iMarker, iVertex, iVar+nPrimVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+nPrimVar+iVar]);
          
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        for (iVar = 0; iVar < nPrimVar; iVar++) {
          SetDonorPrimVar(iMarker, iVertex, iVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar]);
          SetDonorPrimVar(iMarker, iVertex, iVar+nPrimVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+nPrimVar+iVar]);
          
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            // set the proper min and max for this node
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), GetDonorPrimVar(iMarker, iVertex, iVar)));
              node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), GetDonorPrimVar(iMarker, iVertex, iVar+nPrimVar)));
            }
            
            
          }
        }
      }
    }
  }
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Periodic_Limiter2(CGeometry *geometry, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker, iPeriodic, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  
  unsigned long buff_size = nPrimVar*2 + 5;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [buff_size];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              for (iVar = 0; iVar < nPrimVar; iVar++) {
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+iVar] = node[iPoint]->GetLimiter_Primitive(iVar);
              }
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+3)]  = node[iPoint]->GetAuxVar();
              
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        for (iVar = 0; iVar < nPrimVar; iVar++) {
          SetDonorPrimVar(iMarker, iVertex, iVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar]);
          
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        for (iVar = 0; iVar < nPrimVar; iVar++) {
          SetDonorPrimVar(iMarker, iVertex, iVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar]);
          
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            // set the proper min for the limiter
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              node[iPoint]->SetLimiter_Primitive(iVar, min(node[iPoint]->GetLimiter_Primitive(iVar), GetDonorPrimVar(iMarker, iVertex, iVar)));
            }
            
            
          }
        }
      }
    }
  }
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Periodic_Eigenvalue(CGeometry *geometry, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal, iNeighbor, Neighbor_Point, nNeighbor;
  unsigned short iMarker, jMarker, iPeriodic, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  
  unsigned long buff_size = nPrimVar*2 + 5;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [buff_size];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+0] = node[iPoint]->GetLambda();
              
              /*--- Communicate number of neighbors that do not lie on the BC. ---*/
              
              nNeighbor = 0;
              for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
                Neighbor_Point = geometry->node[iPoint]->GetPoint(iNeighbor);
                
                /*--- Check if this neighbor lies on the surface. If not,
                 increment the  ---*/
                if (geometry->node[iNeighbor]->GetBoundary()) nNeighbor++;
                
              }
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+1] = su2double(nNeighbor);
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+3)]  = node[iPoint]->GetAuxVar();
              
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        SetDonorPrimVar(iMarker, iVertex, 0, Buffer_Receive_PrimVar[iPoint*(buff_size)+0]);
        SetDonorPrimVar(iMarker, iVertex, 1, Buffer_Receive_PrimVar[iPoint*(buff_size)+1]);
        
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        SetDonorPrimVar(iMarker, iVertex, 0, Buffer_Receive_PrimVar[iPoint*(buff_size)+0]);
        SetDonorPrimVar(iMarker, iVertex, 1, Buffer_Receive_PrimVar[iPoint*(buff_size)+1]);
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            // Add Lambda
            
            node[iPoint]->AddLambda(GetDonorPrimVar(iMarker, iVertex, 0));
            
            // Add to neighbor count
            
            //geometry->node[iPoint]->SetnNeighbor(geometry->node[iPoint]->GetnNeighbor() + SU2_TYPE::Int(GetDonorPrimVar(iMarker, iVertex, 1)));
            
          }
        }
      }
    }
  }
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Periodic_Laplacian(CGeometry *geometry, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, jPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, iMarker, jMarker, iPeriodic, iNeigh, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  
  unsigned long buff_size = nPrimVar*2 + 5;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double *Diff = new su2double[nVar];
  su2double *Und_Lapl = new su2double[nVar];
  
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [buff_size];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            for (iVar = 0; iVar< nVar; iVar++)
              Und_Lapl[iVar] = 0.0;
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              bool boundary_i, boundary_j;
              su2double Pressure_i, Pressure_j;
              
              for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
                jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
                
                /*--- avoid halos and boundary points so that we don't duplicate edges ---*/
                
                if (geometry->node[jPoint]->GetDomain() && (!geometry->node[jPoint]->GetBoundary())) {
                  
                  /*--- Solution differences ---*/
                  
                  for (iVar = 0; iVar < nVar; iVar++)
                    Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);

                  boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
                  boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
                  
                  /*--- Both points inside the domain, or both in the boundary ---*/
                  
                  if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
                    if (geometry->node[iPoint]->GetDomain()) {
                      for (iVar = 0; iVar< nVar; iVar++)
                        Und_Lapl[iVar] -= Diff[iVar];
                    }
                  }
                  
                  /*--- iPoint inside the domain, jPoint on the boundary ---*/
                  
                  if (!boundary_i && boundary_j)
                    if (geometry->node[iPoint]->GetDomain()){
                      for (iVar = 0; iVar< nVar; iVar++)
                        Und_Lapl[iVar] -= Diff[iVar];
                    }
                  
                }
              }
              
              for (iVar = 0; iVar < nVar; iVar++)
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+iVar] = Und_Lapl[iVar];
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+3)]  = node[iPoint]->GetAuxVar();
              
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        for (iVar = 0; iVar < nVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar]);
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        for (iVar = 0; iVar < nVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar]);
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            // adjust undivided Laplacian
            
            for (iVar = 0; iVar < nVar; iVar++)
              Diff[iVar] = GetDonorPrimVar(iMarker, iVertex, iVar);
            
            // they were subtracted during accumulation, so now just add
            
            node[iPoint]->AddUnd_Lapl(Diff);
            
          }
        }
      }
    }
  }
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
  delete [] Diff;
  delete [] Und_Lapl;
  
}

void CIncEulerSolver::BC_Periodic_Sensor(CGeometry *geometry, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, jPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iMarker, jMarker, iPeriodic, iNeigh, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  
  unsigned long buff_size = nPrimVar*2 + 5;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [buff_size];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              bool boundary_i, boundary_j;
              su2double Sensor_i = 0.0, Sensor_j = 0.0, Pressure_i, Pressure_j;
              
              for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
                jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
                
                /*--- avoid halos and boundary points so that we don't duplicate edges ---*/
                
                if (geometry->node[jPoint]->GetDomain() && (!geometry->node[jPoint]->GetBoundary())) {
                  
                  
                  Pressure_i = node[iPoint]->GetDensity();
                  Pressure_j = node[jPoint]->GetDensity();
                  
                  boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
                  boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
                  
                  /*--- Both points inside the domain, or both on the boundary ---*/
                  
                  if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
                    if (geometry->node[iPoint]->GetDomain()) { Sensor_i += (Pressure_j - Pressure_i); Sensor_j += (Pressure_i + Pressure_j); }
                    
                  }
                  
                  /*--- iPoint inside the domain, jPoint on the boundary ---*/
                  
                  if (!boundary_i && boundary_j)
                    if (geometry->node[iPoint]->GetDomain()) { Sensor_i += (Pressure_j - Pressure_i); Sensor_j += (Pressure_i + Pressure_j); }
                  
                }
                
              }
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+0] = Sensor_i;
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+1] = Sensor_j;
              
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nPrimVar*2+3)]  = node[iPoint]->GetAuxVar();
              
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        SetDonorPrimVar(iMarker, iVertex, 0, Buffer_Receive_PrimVar[iPoint*(buff_size)+0]);
        SetDonorPrimVar(iMarker, iVertex, 1, Buffer_Receive_PrimVar[iPoint*(buff_size)+1]);
        
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nPrimVar*2+2)]);
        
        SetDonorPrimVar(iMarker, iVertex, 0, Buffer_Receive_PrimVar[iPoint*(buff_size)+0]);
        SetDonorPrimVar(iMarker, iVertex, 1, Buffer_Receive_PrimVar[iPoint*(buff_size)+1]);
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            bool boundary_i, boundary_j;
            su2double Sensor_i = 0.0, Sensor_j = 0.0, Pressure_i, Pressure_j;
            
            for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
              
              /*--- avoid halos and boundary points so that we don't duplicate edges ---*/
              
              if (geometry->node[jPoint]->GetDomain()) {
                
                Pressure_i = node[iPoint]->GetDensity();
                Pressure_j = node[jPoint]->GetDensity();
                
                boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
                boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
                
                /*--- Both points inside the domain, or both on the boundary ---*/
                
                if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
                  if (geometry->node[iPoint]->GetDomain()) { Sensor_i += (Pressure_j - Pressure_i); Sensor_j += (Pressure_i + Pressure_j); }
                  
                }
                
                /*--- iPoint inside the domain, jPoint on the boundary ---*/
                
                if (!boundary_i && boundary_j)
                  if (geometry->node[iPoint]->GetDomain()) { Sensor_i += (Pressure_j - Pressure_i); Sensor_j += (Pressure_i + Pressure_j); }
                
              }
              
            }
            
            /*--- Add in the missing values from the periodic neighbor. ---*/
            
            Sensor_i += GetDonorPrimVar(iMarker, iVertex, 0);
            Sensor_j += GetDonorPrimVar(iMarker, iVertex, 1);
            
            /*--- Set the final sensor. ---*/
            
            node[iPoint]->SetSensor(fabs(Sensor_i) / Sensor_j);
            
          }
        }
      }
    }
  }
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *numerics, CConfig *config, unsigned short val_periodic) {
  
  unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
  Buffer_Send_nPointTotal = 0, iGlobalIndex, iGlobal;
  unsigned short iVar, jVar, iMarker, jMarker, iPeriodic, nPeriodic = 0;
  long nDomain = 0, iDomain, jDomain;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
  jacMatrix[5][5] = {{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0}},
  rotJacob[5][5] = {{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0}};
  
  unsigned long buff_size = nVar + nVar*nVar + 5;
  
  /*--- Evaluate the number of periodic boundary conditions ---*/
  
  nPeriodic = config->GetnMarker_Periodic();
  
#ifdef HAVE_MPI
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  SU2_MPI::Status status, status_;
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  su2double        *Buffer_Send_PrimVar          = NULL;
  su2double        *iPrimVar          = new su2double [nPrimVar];
  
  unsigned long *nPointTotal_s = new unsigned long[size];
  unsigned long *nPointTotal_r = new unsigned long[size];
  
  unsigned long Buffer_Size_PrimVar          = 0;
  unsigned long PointTotal_Counter = 0;
  
  /*--- Allocate the memory that we only need if we have MPI support ---*/
  
  su2double        *Buffer_Receive_PrimVar          = NULL;
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  
  /*--- This loop gets the array sizes of points for each
   rank to send to each other rank. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Loop over the markers to perform the dimensionalizaton
     of the domain variables ---*/
    
    Buffer_Send_nPointTotal = 0;
    
    /*--- Loop over all of the markers and count the number of each
     type of point and element that needs to be sent. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              Buffer_Send_nPointTotal++;
            }
          }
        }
      }
    }
    
    /*--- Store the counts on a partition by partition basis. ---*/
    
    nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
    
    /*--- Total counts for allocating send buffers below ---*/
    
    Buffer_Size_PrimVar          += nPointTotal_s[iDomain]*(buff_size);
    
  }
  
  /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
  
  Buffer_Send_PrimVar          = new su2double[Buffer_Size_PrimVar];
  
  /*--- Now that we know the sizes of the point, we can
   allocate and send the information in large chunks to all processors. ---*/
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- A rank does not communicate with itself through MPI ---*/
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the counts to iDomain with non-blocking sends ---*/
      
      SU2_MPI::Bsend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD);
      
#endif
      
    } else {
      
      /*--- If iDomain = rank, we simply copy values into place in memory ---*/
      
      nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
      
    }
    
    /*--- Receive the counts. All processors are sending their counters to
     iDomain up above, so only iDomain needs to perform the recv here from
     all other ranks. ---*/
    
    if (rank == iDomain) {
      
      for (jDomain = 0; jDomain < size; jDomain++) {
        
        /*--- A rank does not communicate with itself through MPI ---*/
        
        if (rank != jDomain) {
          
#ifdef HAVE_MPI
          
          /*--- Recv the data by probing for the current sender, jDomain,
           first and then receiving the values from it. ---*/
          
          SU2_MPI::Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain, rank, MPI_COMM_WORLD, &status);
          
#endif
          
        }
      }
      
    }
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Initialize the counters for the larger send buffers (by domain) ---*/
  
  PointTotal_Counter  = 0;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Set the value of the interior geometry. Initialize counters. ---*/
    
    iPointTotal = 0;
    
    /*--- Load up the actual values into the buffers for sending. ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        
        iPeriodic = config->GetMarker_All_PerBound(iMarker);
        if ((iPeriodic == val_periodic) ||
            (iPeriodic == val_periodic + nPeriodic/2)) {
          
          /*--- Retrieve the supplied periodic information. ---*/
          
          center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
          angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
          trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
          
          /*--- Store (center+trans) as it is constant and will be added on. ---*/
          
          translation[0] = center[0] + trans[0];
          translation[1] = center[1] + trans[1];
          translation[2] = center[2] + trans[2];
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          theta = angles[0];
          phi   = angles[1];
          psi   = angles[2];
          
          cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
          sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
            
            if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
              
              iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
              jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
              jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
              
              for (iVar = 0; iVar < nVar; iVar++) {
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+iVar] = LinSysRes.GetBlock(iPoint, iVar);
              }
              
              /*--- Rotate the momentum components. ---*/
              if (nDim == 2) {
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+1] = (rotMatrix[0][0]*LinSysRes.GetBlock(iPoint, 1) +
                                                                                       rotMatrix[0][1]*LinSysRes.GetBlock(iPoint, 2));
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+2] = (rotMatrix[1][0]*LinSysRes.GetBlock(iPoint, 1) +
                                                                                       rotMatrix[1][1]*LinSysRes.GetBlock(iPoint, 2));
              }
              else {
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+1] = (rotMatrix[0][0]*LinSysRes.GetBlock(iPoint, 1) +
                                                                                       rotMatrix[0][1]*LinSysRes.GetBlock(iPoint, 2) +
                                                                                       rotMatrix[0][2]*LinSysRes.GetBlock(iPoint, 3));
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+2] = (rotMatrix[1][0]*LinSysRes.GetBlock(iPoint, 1) +
                                                                                       rotMatrix[1][1]*LinSysRes.GetBlock(iPoint, 2) +
                                                                                       rotMatrix[1][2]*LinSysRes.GetBlock(iPoint, 3));
                Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+3] = (rotMatrix[2][0]*LinSysRes.GetBlock(iPoint, 1) +
                                                                                       rotMatrix[2][1]*LinSysRes.GetBlock(iPoint, 2) +
                                                                                       rotMatrix[2][2]*LinSysRes.GetBlock(iPoint, 3));
              }
              
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nVar+0)]  = su2double(iGlobalIndex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nVar+1)]  = su2double(jVertex);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nVar+2)]  = su2double(jMarker);
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nVar+3)]  = geometry->node[iPoint]->GetVolume();
              Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nVar+4)]  = node[iPoint]->GetDelta_Time();
              
              if (implicit) {
                
                for (iVar = 0; iVar < nVar; iVar++) {
                  for (jVar = 0; jVar < nVar; jVar++) {
                    jacMatrix[iVar][jVar] = Jacobian.GetBlock(iPoint, iPoint, iVar, jVar);
                    rotJacob[iVar][jVar]  = Jacobian.GetBlock(iPoint, iPoint, iVar, jVar);
                  }
                }
                
                //Rotate the Jacobian momentum components
                //              for (unsigned short jVar = 1; jVar <= nDim; jVar++) {
                //              //  for (unsigned short jVar = 0; jVar < nVar; jVar++) {
                //
                //                if (nDim == 2) {
                //                  rotJacob[1][jVar] = (rotMatrix[0][0]*jacMatrix[1][jVar] +
                //                                       rotMatrix[0][1]*jacMatrix[2][jVar]);
                //                  rotJacob[2][jVar] = (rotMatrix[1][0]*jacMatrix[1][jVar] +
                //                                       rotMatrix[1][1]*jacMatrix[2][jVar]);
                //                } else {
                //
                //                rotJacob[1][jVar] = (rotMatrix[0][0]*jacMatrix[1][jVar] +
                //                                     rotMatrix[0][1]*jacMatrix[2][jVar] +
                //                                     rotMatrix[0][2]*jacMatrix[3][jVar]);
                //                rotJacob[2][jVar] = (rotMatrix[1][0]*jacMatrix[1][jVar] +
                //                                     rotMatrix[1][1]*jacMatrix[2][jVar] +
                //                                     rotMatrix[1][2]*jacMatrix[3][jVar]);
                //                rotJacob[3][jVar] = (rotMatrix[2][0]*jacMatrix[1][jVar] +
                //                                     rotMatrix[2][1]*jacMatrix[2][jVar] +
                //                                     rotMatrix[2][2]*jacMatrix[3][jVar]);
                //              }
                //              }
                
                int ii = 0;
                for (iVar = 0; iVar < nVar; iVar++) {
                  for (jVar = 0; jVar < nVar; jVar++) {
                    Buffer_Send_PrimVar[(buff_size)*(PointTotal_Counter+iPointTotal)+(nVar+5)+ii] = rotJacob[iVar][jVar];
                    ii++;
                  }
                }
                
              }
              
              iPointTotal++;
              
            }
            
          }
          
        }
        
      }
      
    }
    
    /*--- Send the buffers with the geometrical information ---*/
    
    if (iDomain != rank) {
      
#ifdef HAVE_MPI
      
      /*--- Communicate the coordinates, global index, colors, and element
       date to iDomain with non-blocking sends. ---*/
      
      SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)],
                     nPointTotal_s[iDomain]*(buff_size), MPI_DOUBLE, iDomain,
                     iDomain,  MPI_COMM_WORLD);
      
#endif
      
    }
    
    else {
      
      /*--- Allocate local memory for the local recv of the elements ---*/
      
      Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(buff_size)];
      
      for (iter = 0; iter < nPointTotal_s[iDomain]*(buff_size); iter++)
        Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(buff_size)+iter];
      
      /*--- Recv the point data from ourselves (same procedure as above) ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+2)]);
        for (iVar = 0; iVar < nVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar];
        
        iPrimVar[nVar+3] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nVar+3];
        iPrimVar[nVar+4] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nVar+4];
        
        if (implicit){
          int ii = 0;
          for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_i[iVar][jVar] = Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+5)+ii];
              SetDonorJacobian(iMarker, iVertex, iVar, jVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+5)+ii]);
              ii++;
            }
          }
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
    }
    
    /*--- Increment the counters for the send buffers (iDomain loop) ---*/
    
    PointTotal_Counter += iPointTotal;
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- The next section begins the recv of all data for the interior
   points/elements in the mesh. First, create the domain structures for
   the points on this rank. First, we recv all of the point data ---*/
  
  for (iDomain = 0; iDomain < size; iDomain++) {
    
    if (rank != iDomain) {
      
#ifdef HAVE_MPI
      
      /*--- Allocate the receive buffer vector. Send the colors so that we
       know whether what we recv is an owned or halo node. ---*/
      
      Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(buff_size)];
      
      /*--- Receive the buffers with the coords, global index, and colors ---*/
      
      SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(buff_size) , MPI_DOUBLE,
                    iDomain, rank, MPI_COMM_WORLD, &status_);
      
      /*--- Loop over all of the points that we have recv'd and store the
       coords, global index vertex and markers ---*/
      
      for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
        
        iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+0)]);
        iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+1)]);
        iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+2)]);
        for (iVar = 0; iVar < nPrimVar; iVar++)
          iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(buff_size)+iVar];
        
        iPrimVar[nVar+3] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nVar+3];
        iPrimVar[nVar+4] = Buffer_Receive_PrimVar[iPoint*(buff_size)+nVar+4];
        
        if (implicit){
          
          int ii = 0;
          for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
              Jacobian_i[iVar][jVar] = Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+5)+ii];
              SetDonorJacobian(iMarker, iVertex, iVar, jVar, Buffer_Receive_PrimVar[iPoint*(buff_size)+(nVar+5)+ii]);
              ii++;
            }
          }
        }
        
        if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
        if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
        
        if (iMarker > nMarker) cout << "ERROR" <<  endl;
        if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
        
        for (iVar = 0; iVar < nPrimVar; iVar++)
          SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
        
        SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
        
      }
      
      /*--- Delete memory for recv the point stuff ---*/
      
      delete [] Buffer_Receive_PrimVar;
      
#endif
      
    }
    
  }
  
  /*--- Wait for the non-blocking sends to complete. ---*/
  
#ifdef HAVE_MPI
  
  SU2_MPI::Barrier(MPI_COMM_WORLD);
  
#endif
  
  /*--- Free all of the memory used for communicating points and elements ---*/
  
  delete[] Buffer_Send_PrimVar;
  
  /*--- Release all of the temporary memory ---*/
  
  delete [] nPointTotal_s;
  delete [] nPointTotal_r;
  delete [] iPrimVar;
  
  unsigned long GlobalIndex_iPoint, GlobalIndex_jPoint;
  
  su2double *Normal    = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  /*--- For the first pass through the periodic BC, make all nodes on
   the periodic surfaces act as point implicit so that we maintain a
   consistent implicit solution on both sides of the periodic boundary. ---*/
  
  if (val_periodic == 1) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (implicit) {
            for (iVar = 0; iVar < nVar; iVar++) {
              unsigned long total_index = iPoint*nVar+iVar;
              Jacobian.SetPointImplicit(total_index);
            }
          }
        }
      }
    }
  }
  
  /*--- Now perform the residual & Jacobian updates with the recv data. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
          GlobalIndex_jPoint = GetDonorGlobalIndex(iMarker, iVertex);
          
          if ((geometry->node[iPoint]->GetDomain()) &&
              (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            /*--- Update the volume and time step using the donor info. ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
              PrimVar_j[iVar] = GetDonorPrimVar(iMarker, iVertex, iVar);
            }
            
            su2double Vol = geometry->node[iPoint]->GetVolume();
            geometry->node[iPoint]->SetVolume(Vol + PrimVar_j[nVar+3]);
            
            su2double dt = node[iPoint]->GetDelta_Time();
            if (PrimVar_j[nVar+4] < dt)
              node[iPoint]->SetDelta_Time(PrimVar_j[nVar+4]);
            
            /*--- Access the residual from the donor. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              Residual[iVar] = GetDonorPrimVar(iMarker, iVertex, iVar);
            }
            
            /*--- Access the Jacobian from the donor. ---*/
            
            if (implicit) {
              int ii = 0;
              for (iVar = 0; iVar < nVar; iVar++) {
                for (jVar = 0; jVar < nVar; jVar++) {
                  Jacobian_i[iVar][jVar] = GetDonorJacobian(iMarker, iVertex, iVar, jVar);
                  ii++;
                }
              }
            }
            
            /*--- Add contributions to total residual and Jacobian. ---*/
            
            LinSysRes.AddBlock(iPoint, Residual);
            if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
            
          }
        }
      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CIncEulerSolver::BC_Custom(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) { }

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
  
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement    = config->GetGrid_Movement();
  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool energy           = config->GetEnergy_Equation();
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
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
      
      V_time_nM1 = node[iPoint]->GetSolution_time_n1();
      V_time_n   = node[iPoint]->GetSolution_time_n();
      V_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- Access the density and Cp at this node (constant for now). ---*/
      
      Density     = node[iPoint]->GetDensity();
      Cp          = node[iPoint]->GetSpecificHeatCp();
      
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
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }
      
      if (!energy) Residual[nDim+1] = 0.0;
      
      //cout << " iPoint: " << iPoint << "  " << Volume_nP1<< "  " << Residual[0]<< "  " << Residual[1]<< "  " << Residual[2]<< "  " << Residual[3]<<endl;
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      if (implicit) {
        
        unsigned short iDim, jDim;
        
        su2double  BetaInc2, Density, dRhodT, Temperature, oneOverCp, Cp;
        su2double  Velocity[3] = {0.0,0.0,0.0};
        
        /*--- Access the primitive variables at this node. ---*/
        
        Density     = node[iPoint]->GetDensity();
        BetaInc2    = node[iPoint]->GetBetaInc2();
        Cp          = node[iPoint]->GetSpecificHeatCp();
        oneOverCp   = 1.0/Cp;
        Temperature = node[iPoint]->GetTemperature();
        
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
        
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
        
          /*--- For implicit calculations, we multiply the preconditioner
           by the cell volume over the time step and add to the Jac diagonal. ---*/
          
          Jacobian_i[0][0] = 1.0/BetaInc2;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_i[iDim+1][0] = Velocity[iDim]/BetaInc2;
          
          if (energy) Jacobian_i[nDim+1][0] = Cp*Temperature/BetaInc2;
          else        Jacobian_i[nDim+1][0] = 0.0;
          
          for (jDim = 0; jDim < nDim; jDim++) {
            Jacobian_i[0][jDim+1] = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              if (iDim == jDim) Jacobian_i[iDim+1][jDim+1] = Density;
              else Jacobian_i[iDim+1][jDim+1] = 0.0;
            }
            Jacobian_i[nDim+1][jDim+1] = 0.0;
          }
          
          Jacobian_i[0][nDim+1] = dRhodT;
          for (iDim = 0; iDim < nDim; iDim++)
            Jacobian_i[iDim+1][nDim+1] = Velocity[iDim]*dRhodT;
          
          if (energy) Jacobian_i[nDim+1][nDim+1] = Cp*(dRhodT*Temperature + Density);
          else        Jacobian_i[nDim+1][nDim+1] = 1.0;
          
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
              Jacobian_i[iVar][jVar] *= Volume_nP1 / TimeStep;
            if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
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
      
      V_time_n = node[iPoint]->GetSolution_time_n();
      
      /*--- Access the density and Cp at this node (constant for now). ---*/
      
      Density     = node[iPoint]->GetDensity();
      Cp          = node[iPoint]->GetSpecificHeatCp();
      
      /*--- Compute the conservative variable vector for all time levels. ---*/
      
      U_time_n[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_n[iDim+1] = Density*V_time_n[iDim+1];
      }
      U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];
      
      for (iVar = 1; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      V_time_n = node[jPoint]->GetSolution_time_n();
      
      U_time_n[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_n[iDim+1] = Density*V_time_n[iDim+1];
      }
      U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];
      
      for (iVar = 1; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
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
        
        V_time_n = node[iPoint]->GetSolution_time_n();
        
        /*--- Access the density and Cp at this node (constant for now). ---*/
        
        Density     = node[iPoint]->GetDensity();
        Cp          = node[iPoint]->GetSpecificHeatCp();
        
        U_time_n[0] = Density;
        for (iDim = 0; iDim < nDim; iDim++) {
          U_time_n[iDim+1] = Density*V_time_n[iDim+1];
        }
        U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];
        
        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
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
      
      V_time_nM1 = node[iPoint]->GetSolution_time_n1();
      V_time_n   = node[iPoint]->GetSolution_time_n();
      V_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- Access the density and Cp at this node (constant for now). ---*/
      
      Density     = node[iPoint]->GetDensity();
      Cp          = node[iPoint]->GetSpecificHeatCp();
      
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
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 1; iVar < nVar; iVar++) {
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_i[iDim+1][iDim+1] = Density*Jacobian_i[iDim+1][iDim+1];
        Jacobian_i[nDim+1][nDim+1] = Density*Cp*Jacobian_i[nDim+1][nDim+1];
        
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CIncEulerSolver::ComputeResidual_BGS(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
      SetRes_BGS(iVar,0.0);
      SetRes_Max_BGS(iVar,0.0,0);
  }

  /*--- Set the residuals ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
          residual = node[iPoint]->GetSolution(iVar) - node[iPoint]->Get_BGSSolution_k(iVar);
          AddRes_BGS(iVar,residual*residual);
          AddRes_Max_BGS(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_BGS(geometry, config);

}


void CIncEulerSolver::UpdateSolution_BGS(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  /*--- To nPoint: The solution must be communicated beforehand ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++){

    node[iPoint]->Set_BGSSolution_k();

  }

}

void CIncEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {
  
  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool grid_movement  = config->GetGrid_Movement();
  bool static_fsi = ((config->GetUnsteady_Simulation() == STEADY) &&
                     (config->GetFSI_Simulation()));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool steady_restart = config->GetSteadyRestart();
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  string UnstExt, text_line;
  ifstream restart_file;
  
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string restart_filename = config->GetSolution_FlowFileName();

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;
  
  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();
  
  /*--- Adjust the number of solution variables in the restart. We always
   carry a space in nVar for the energy equation in the solver, but we only
   write it to the restart if it is active. Therefore, we must reduce nVar
   here if energy is inactive so that the restart is read correctly. ---*/
  
  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();
  
  unsigned short nVar_Restart = nVar;
  if ((!energy) && (!weakly_coupled_heat)) nVar_Restart--;
  Solution[nVar-1] = GetTemperature_Inf();
  
  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
  restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/
  
  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

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
      node[iPoint_Local]->SetSolution(Solution);
      iPoint_Global_Local++;

      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/

      if (grid_movement && val_update_geo) {

        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/

        if (turb_model == SA || turb_model == SA_NEG) {
          index++;
        } else if (turb_model == SST) {
          index+=2;
        }

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/

        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar_Restart;
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
  
  solver[MESH_0][FLOW_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

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
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  }
  
  /*--- Update the geometry for flows on dynamic meshes ---*/
  
  if (grid_movement && val_update_geo) {
    
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
  
  /*--- Update the geometry for flows on static FSI problems with moving meshes ---*/
  
  if (static_fsi && val_update_geo) {
    
    /*--- Communicate the new coordinates and grid velocities at the halos ---*/
    
    geometry[MESH_0]->Set_MPI_Coord(config);
    
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
    node[iPoint]->SetSolution(0, Pressure_Inf);
    for (iDim = 0; iDim < nDim; iDim++){
      node[iPoint]->SetSolution(iDim+1, Velocity_Inf[iDim]);
    }
    node[iPoint]->SetSolution(nDim+1, Temperature_Inf);
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
  
  DonorPrimVar = NULL; DonorGlobalIndex = NULL; DonorJacobian = NULL;
  
}

CIncNSSolver::CIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CIncEulerSolver() {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  string filename_ = config->GetSolution_FlowFileName();
  bool fsi     = config->GetFSI_Simulation();

  unsigned short direct_diff = config->GetDirectDiff();

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/

    Read_SU2_Restart_Metadata(geometry, config, false, filename_);
    
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
  CSkinFriction = NULL;    Cauchy_Serie = NULL; HF_Visc = NULL;
  
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
  
  /*--- Allocate the node variables ---*/
  node = new CVariable*[nPoint];
  
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

  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/
  
  DonorPrimVar = new su2double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      if (rans) {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar*2];
        for (iVar = 0; iVar < nPrimVar*2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar*2];
        for (iVar = 0; iVar < nPrimVar*2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the Jacobian from the donor for the periodic
   boundary condition (point implicit). ---*/
  
  DonorJacobian = new su2double*** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorJacobian[iMarker] = new su2double** [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorJacobian[iMarker][iVertex] = new su2double* [nPrimVarGrad*nDim];
      for (iVar = 0; iVar < nPrimVarGrad*nDim; iVar++) {
        DonorJacobian[iMarker][iVertex][iVar] = new su2double [nPrimVarGrad*nDim];
        for (unsigned short jVar = 0; jVar < nPrimVarGrad*nDim; jVar++) {
          DonorJacobian[iMarker][iVertex][iVar][jVar] = 0.0;
        }
      }
    }
  }
  
  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/
  
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
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
  AoA_Prev       = 0.0;
  Total_CL_Prev  = 0.0;  Total_CD_Prev      = 0.0;
  Total_CMx_Prev = 0.0;  Total_CMy_Prev     = 0.0;  Total_CMz_Prev       = 0.0;

  /*--- Coefficients for fixed lift mode. ---*/
  
  AoA_Prev = 0.0;
  Total_CL_Prev = 0.0; Total_CD_Prev = 0.0;
  Total_CMx_Prev = 0.0; Total_CMy_Prev = 0.0; Total_CMz_Prev = 0.0;

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

  /*--- Initialize the cauchy critera array for fixed CL mode ---*/

  if (config->GetFixed_CL_Mode())
    Cauchy_Serie = new su2double [config->GetCauchy_Elems()+1];

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CIncNSVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, config);

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (fsi){
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

  /*--- Store volume for periodic if necessary (MG too). ---*/
  
  if (config->GetnMarker_Periodic() != 0) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetAuxVar(geometry->node[iPoint]->GetVolume());
  }
  
  /*--- Perform the MPI communication of the solution ---*/

  Set_MPI_Solution(geometry, config);

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

  if (Cauchy_Serie != NULL) delete [] Cauchy_Serie;
  
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
  
  unsigned long ExtIter     = config->GetExtIter();
  bool cont_adjoint         = config->GetContinuous_Adjoint();
  bool disc_adjoint         = config->GetDiscrete_Adjoint();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center               = ((config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED));
  bool center_jst           = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow         = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_turb         = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_adjflow      = (cont_adjoint && (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool fixed_cl             = config->GetFixed_CL_Mode();
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;

  /*--- Store the original volume for periodic cells on the boundaries,
   since this will be increased as we complete the CVs during our
   updates. ---*/
  
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Reset to the original volume. ---*/
        geometry->node[iPoint]->SetVolume(node[iPoint]->GetAuxVar());
        
      }
    }
  }
  
  /*--- Update the angle of attack at the far-field for fixed CL calculations (only direct problem). ---*/
  
  if ((fixed_cl) && (!disc_adjoint) && (!cont_adjoint)) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Artificial dissipation ---*/
  
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
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
  
  /*--- Update the beta value based on the maximum velocity / viscosity. ---*/

  SetBeta_Parameter(geometry, solver_container, config, iMesh);

  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  
  StrainMag_Max = 0.0; Omega_Max = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    solver_container[FLOW_SOL]->node[iPoint]->SetVorticity();
    solver_container[FLOW_SOL]->node[iPoint]->SetStrainMag();
    
    StrainMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    Vorticity = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity();
    Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);
    
    StrainMag_Max = max(StrainMag_Max, StrainMag);
    Omega_Max = max(Omega_Max, Omega);
    
  }
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (implicit && !disc_adjoint) Jacobian.SetValZero();

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

unsigned long CIncNSSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0, DES_LengthScale = 0.0;
  unsigned short turb_model = config->GetKind_Turb_Model();
  bool physical = true;
  
  bool tkeNeeded = (turb_model == SST);
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Retrieve the value of the kinetic energy (if needed) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
      
      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
        DES_LengthScale = solver_container[TURB_SOL]->node[iPoint]->GetDES_LengthScale();
      }
    }
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Incompressible flow, primitive variables --- */

    physical = node[iPoint]->SetPrimVar(eddy_visc, turb_ke, FluidModel);
    
    /*--- Record any non-physical points. ---*/

    if (!physical) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }

    /*--- Set the DES length scale ---*/
    
    node[iPoint]->SetDES_LengthScale(DES_LengthScale);    
    
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
  bool grid_movement = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool energy = config->GetEnergy_Equation();

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
    
    Mean_ProjVel    = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_BetaInc2   = 0.5 * (node[iPoint]->GetBetaInc2()      + node[jPoint]->GetBetaInc2());
    Mean_Density    = 0.5 * (node[iPoint]->GetDensity()       + node[jPoint]->GetDensity());
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
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
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
    /*--- Viscous contribution ---*/
    
    Mean_LaminarVisc          = 0.5*(node[iPoint]->GetLaminarViscosity()    + node[jPoint]->GetLaminarViscosity());
    Mean_EddyVisc             = 0.5*(node[iPoint]->GetEddyViscosity()       + node[jPoint]->GetEddyViscosity());
    Mean_Density              = 0.5*(node[iPoint]->GetDensity()             + node[jPoint]->GetDensity());
    Mean_Thermal_Conductivity = 0.5*(node[iPoint]->GetThermalConductivity() + node[jPoint]->GetThermalConductivity());
    Mean_Cv                   = 0.5*(node[iPoint]->GetSpecificHeatCv()      + node[jPoint]->GetSpecificHeatCv());

    Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
    Lambda_2 = 0.0;
    if (energy) Lambda_2 = (1.0/Mean_Cv)*Mean_Thermal_Conductivity;
    Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);
    
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
      
      Mean_ProjVel    = node[iPoint]->GetProjVel(Normal);
      Mean_BetaInc2   = node[iPoint]->GetBetaInc2();
      Mean_Density    = node[iPoint]->GetDensity();
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

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

      Mean_LaminarVisc          = node[iPoint]->GetLaminarViscosity();
      Mean_EddyVisc             = node[iPoint]->GetEddyViscosity();
      Mean_Density              = node[iPoint]->GetDensity();
      Mean_Thermal_Conductivity = node[iPoint]->GetThermalConductivity();
      Mean_Cv                   = node[iPoint]->GetSpecificHeatCv();

      Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
      Lambda_2 = 0.0;
      if (energy) Lambda_2 = (1.0/Mean_Cv)*Mean_Thermal_Conductivity;
      Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
      
    }
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
    
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(),
                           node[jPoint]->GetPrimitive());
    numerics->SetSecondary(node[iPoint]->GetSecondary(),
                           node[jPoint]->GetSecondary());
    
    /*--- Gradient and limiters ---*/
    
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                 node[jPoint]->GetGradient_Primitive());
    
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

  su2double Alpha     = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta      = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea   = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
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

    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {

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
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
          }
          Grad_Temp[iDim] = node[iPoint]->GetGradient_Primitive(nDim+1, iDim);
        }

        Viscosity = node[iPoint]->GetLaminarViscosity();
        Density = node[iPoint]->GetDensity();

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
        
        thermal_conductivity       = node[iPoint]->GetThermalConductivity();
        HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature;
        HF_Visc[iMarker]          += HeatFlux[iMarker][iVertex]*Area;
        MaxHF_Visc[iMarker]       += pow(HeatFlux[iMarker][iVertex], MaxNorm);

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
  bool grid_movement = config->GetGrid_Movement();
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
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
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
  bool grid_movement = config->GetGrid_Movement();
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
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
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
        
        dTdn = -(node[Point_Normal]->GetTemperature() - Twall)/dist_ij;
        
        /*--- Get thermal conductivity ---*/

        thermal_conductivity = node[iPoint]->GetThermalConductivity();

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
  unsigned long iVertex, iPoint, Point_Normal, total_index;

  su2double *GridVel;
  su2double *Normal, *Coord_i, *Coord_j, Area, dist_ij;
  su2double Tconjugate, dTdn;
  su2double thermal_conductivity;
  su2double Temperature_Ref = config->GetTemperature_Ref();

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
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

      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }

      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/

      node[iPoint]->SetVelocity_Old(Vector);

      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();

      if (energy) {

        Tconjugate = GetConjugateHeatVariable(val_marker, iVertex, 0)/Temperature_Ref;

//        node[iPoint]->SetSolution_Old(nDim+1, Tconjugate);
//        node[iPoint]->SetEnergy_ResTruncError_Zero();

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

        dTdn = -(node[Point_Normal]->GetTemperature() - Tconjugate)/dist_ij;

        /*--- Get thermal conductivity ---*/

        thermal_conductivity = node[iPoint]->GetThermalConductivity();

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
//        if(energy) {
//          total_index = iPoint*nVar+nDim+1;
//          Jacobian.DeleteValsRowi(total_index);
//        }
      }

    }
  }
}
