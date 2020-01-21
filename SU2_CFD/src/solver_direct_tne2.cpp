/*!
 * \file solver_direct_tne2.cpp
 * \brief Main subrotuines for solving flows in thermochemical nonequilibrium.
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
#include "../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../include/variables/CTNE2EulerVariable.hpp"
#include "../include/variables/CTNE2NSVariable.hpp"
#include <math.h>

CTNE2EulerSolver::CTNE2EulerSolver(void) : CSolver() {

  /*--- Array initialization ---*/
  CD_Inv   = NULL; CL_Inv   = NULL; CSF_Inv  = NULL;  CEff_Inv = NULL;
  CMx_Inv  = NULL; CMy_Inv  = NULL; CMz_Inv  = NULL;
  CFx_Inv  = NULL; CFy_Inv  = NULL; CFz_Inv  = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt   = NULL; CL_Mnt   = NULL; CSF_Mnt  = NULL;  CEff_Mnt = NULL;
  CMx_Mnt  = NULL; CMy_Mnt  = NULL; CMz_Mnt  = NULL;
  CFx_Mnt  = NULL; CFy_Mnt  = NULL; CFz_Mnt  = NULL;
  CoPx_Mnt = NULL; CoPy_Mnt = NULL; CoPz_Mnt = NULL;

  CPressure      = NULL; CPressureTarget = NULL; YPlus = NULL;
  HeatFluxTarget = NULL; HeatFlux        = NULL;
  ForceInviscid  = NULL; MomentInviscid  = NULL;
  ForceMomentum  = NULL; MomentMomentum  = NULL;

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

  /*--- Numerical methods array initialization ---*/
  iPoint_UndLapl = NULL; jPoint_UndLapl = NULL;
  Primitive      = NULL; CharacPrimVar  = NULL;
  Primitive_i    = NULL; Primitive_j    = NULL;
  LowMach_Precontioner = NULL;

  DonorPrimVar   = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT   = NULL;

  Smatrix     = NULL; Cvector     = NULL;
  Secondary   = NULL;
  Secondary_i = NULL; Secondary_j = NULL;

  nodes = nullptr;
  node_infty = nullptr;
}

CTNE2EulerSolver::CTNE2EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

  unsigned long iPoint, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, iSpecies, nLineLets;
  unsigned short nZone = geometry->GetnZone();
  unsigned short iZone = config->GetiZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == TNE2_RANS )|| (config->GetKind_Solver() == DISC_ADJ_TNE2_RANS));
  unsigned short direct_diff = config->GetDirectDiff();
  int Unst_RestartIter;
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;

  bool adjoint = config->GetDiscrete_Adjoint();
  bool multizone = config->GetMultizone_Problem();

  string filename_ = "flow";

  su2double *Mvec_Inf;
  su2double Alpha, Beta;
  bool check_infty, nonPhys;

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the AoA
  before non-dimensionalizing ---*/
  if (!(!restart || (iMesh != MESH_0) || nZone >1 )) {

    /*--- Modify file name for a dual-time unsteady restart ---*/
    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/
    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
    }

    /*--- Read and store the restart metadata ---*/
    Read_SU2_Restart_Metadata(geometry, config, false, filename_);

  }

  /*--- Array initialization ---*/

  /*--- Basic array initialization ---*/
  CD_Inv   = NULL; CL_Inv   = NULL; CSF_Inv  = NULL; CEff_Inv = NULL;
  CMx_Inv  = NULL; CMy_Inv  = NULL; CMz_Inv  = NULL;
  CFx_Inv  = NULL; CFy_Inv  = NULL; CFz_Inv  = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt   = NULL; CL_Mnt   = NULL; CSF_Mnt  = NULL; CEff_Mnt= NULL;
  CMx_Mnt  = NULL; CMy_Mnt  = NULL; CMz_Mnt  = NULL;
  CFx_Mnt  = NULL; CFy_Mnt  = NULL; CFz_Mnt  = NULL;
  CoPx_Mnt = NULL; CoPy_Mnt = NULL; CoPz_Mnt = NULL;

  CPressure     = NULL; CPressureTarget = NULL; YPlus = NULL;
  HeatFlux      = NULL; HeatFluxTarget  = NULL;
  ForceInviscid = NULL; MomentInviscid  = NULL;
  ForceMomentum = NULL; MomentMomentum  = NULL;

  /*--- Surface based array initialization ---*/
  Surface_CL_Inv  = NULL; Surface_CD_Inv  = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt  = NULL; Surface_CD_Mnt  = NULL; Surface_CSF_Mnt = NULL; Surface_CEff_Mnt= NULL;
  Surface_CFx_Mnt = NULL; Surface_CFy_Mnt = NULL; Surface_CFz_Mnt = NULL;
  Surface_CMx_Mnt = NULL; Surface_CMy_Mnt = NULL; Surface_CMz_Mnt = NULL;

  Surface_CL  = NULL; Surface_CD  = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Supersonic simulation array initialization ---*/
  CEquivArea_Inv   = NULL;
  CNearFieldOF_Inv = NULL;

  /*--- Numerical methods array initialization ---*/
  iPoint_UndLapl = NULL; jPoint_UndLapl   = NULL;
  Primitive      = NULL; CharacPrimVar    = NULL;
  Primitive_i    = NULL; Primitive_j      = NULL;
  DonorPrimVar   = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT   = NULL;

  LowMach_Precontioner = NULL;
  Smatrix     = NULL; Cvector     = NULL;
  Secondary   = NULL;
  Secondary_i = NULL; Secondary_j = NULL;

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometric constants in the solver structure ---*/
  nSpecies     = config->GetnSpecies();
  nMarker      = config->GetnMarker_All();
  nDim         = geometry->GetnDim();

  /*--- Set size of the conserved and primitive vectors ---*/
  //     U: [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  //     V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradV: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  nVar         = nSpecies + nDim + 2;
  nPrimVar     = nSpecies + nDim + 8;
  nPrimVarGrad = nSpecies + nDim + 8;

  // THIS ARE WTRONG
  nSecondaryVar     = nPrimVarGrad;
  nSecondaryVarGrad = nPrimVarGrad; //These may be used for AD support?

  /*--- Initialize nVarGrad for deallocation ---*/
  nVarGrad     = nPrimVarGrad;
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Store the number of vertices on each marker for deallocation ---*/
  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /*--- Perform the non-dimensionalization for the flow equations using the
    specified reference values. ---*/
  SetNondimensionalization(geometry, config, iMesh);

  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]     = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]   = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]   = 0.0;
  Res_Conv     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]     = 0.0;
  Res_Visc     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]     = 0.0;
  Res_Sour     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]     = 0.0;

  /*--- Define some structure for locating max residuals ---*/
  Point_Max       = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++){
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

  /*--- Define some auxiliary vectors related to the Source term evolution ---*/
  Source      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Source[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  if (config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }

  /*--- Allocate arrays for conserved variable limits ---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    lowerlimit[iSpecies] = 0.0;
    upperlimit[iSpecies] = 1E16;
  }
  for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
    lowerlimit[iVar] = -1E16;
    upperlimit[iVar] = 1E16;
  }
  for (iVar = nSpecies+nDim; iVar < nSpecies+nDim+2; iVar++) {
    lowerlimit[iVar] = 0.0;
    upperlimit[iVar] = 1E16;
  }

  /*--- Initialize the solution & residual CVectors ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Create the structure for storing extra information ---*/
  if (config->GetExtraOutput()) {
    nOutputVariables = nVar;
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
  }

  /*--- Allocate Jacobians for implicit time-stepping ---*/
  if (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT) {
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Jacobians and vector  structures for implicit computations ---*/
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
  }
  else {
    if (rank == MASTER_NODE)  cout<< "Explicit Scheme. No Jacobian structure (Euler). MG level: " << iMesh <<"."<<endl;
  }

  /*--- Allocate arrays for gradient computation by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];

    /*--- c vector := transpose(WA)*(Wb) ---*/
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
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar+2];
        for (iVar = 0; iVar < nPrimVar + 2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
        for (iVar = 0; iVar < nPrimVar ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
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

  /*--- Allocate force & coefficient arrays on boundaries ---*/
  CPressure = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      CPressure[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Non dimensional coefficients ---*/
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

  /*--- Rotocraft coefficients ---*/
  CT_Inv           = new su2double[nMarker];
  CQ_Inv           = new su2double[nMarker];
  CMerit_Inv       = new su2double[nMarker];

  CT_Mnt           = new su2double[nMarker];
  CQ_Mnt           = new su2double[nMarker];
  CMerit_Mnt       = new su2double[nMarker];

  /*--- Supersonic coefficients ---*/
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];

  /*--- Initialize total coefficients ---*/
  Total_CD        = 0.0;    Total_CL           = 0.0;    Total_CSF            = 0.0;
  Total_CMx       = 0.0;    Total_CMy          = 0.0;    Total_CMz            = 0.0;
  Total_CoPx      = 0.0;    Total_CoPy         = 0.0;    Total_CoPz           = 0.0;
  Total_CEff      = 0.0;    Total_CEquivArea   = 0.0;    Total_CNearFieldOF   = 0.0;
  Total_CFx       = 0.0;    Total_CFy          = 0.0;    Total_CFz            = 0.0;
  Total_CT        = 0.0;    Total_CQ           = 0.0;    Total_CMerit         = 0.0;
  Total_MaxHeat   = 0.0;    Total_Heat         = 0.0;    Total_ComboObj       = 0.0;
  Total_CpDiff    = 0.0;    Total_HeatFluxDiff = 0.0;    Total_Custom_ObjFunc = 0.0;
  Total_NetThrust = 0.0;
  Total_Power     = 0.0;    AoA_Prev           = 0.0;
  Total_CL_Prev   = 0.0;    Total_CD_Prev      = 0.0;
  Total_CMx_Prev  = 0.0;    Total_CMy_Prev     = 0.0;    Total_CMz_Prev       = 0.0;
  Total_AeroCD    = 0.0;    Total_SolidCD      = 0.0;
  Total_IDR       = 0.0;    Total_IDC          = 0.0;

  /*--- Read farfield conditions from the config file ---*/
  Density_Inf        = config->GetDensity_FreeStreamND();
  Pressure_Inf       = config->GetPressure_FreeStreamND();
  Velocity_Inf       = config->GetVelocity_FreeStreamND();
  Temperature_Inf    = config->GetTemperature_FreeStreamND();
  Mach_Inf           = config->GetMach();
  Temperature_ve_Inf = config->GetTemperature_ve_FreeStreamND();
  MassFrac_Inf       = config->GetMassFrac_FreeStream();

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

  /*--- Create a CVariable that stores the free-stream values ---*/
  node_infty = new CTNE2EulerVariable(Pressure_Inf, MassFrac_Inf, Mvec_Inf, Temperature_Inf,
                                      Temperature_ve_Inf, 1, nDim, nVar,
                                      nPrimVar, nPrimVarGrad, config);
  check_infty = node_infty->SetPrimVar_Compressible(0,config);

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  nodes = new CTNE2EulerVariable(Pressure_Inf, MassFrac_Inf, Mvec_Inf, Temperature_Inf,
                                 Temperature_ve_Inf, nPoint, nDim, nVar,
                                 nPrimVar, nPrimVarGrad, config);
  SetBaseClassPointerToNodes();

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nonPhys = nodes->SetPrimVar_Compressible(iPoint, config);

    if (nonPhys) {
      bool ionization;
      unsigned short iEl, nHeavy, nEl, *nElStates;
      su2double RuSI, Ru, T, Tve, rhoCvtr, sqvel, rhoE, rhoEve, num, denom, conc;
      su2double rho, rhos, Ef, Ev, Ee, soundspeed;
      su2double *xi, *Ms, *thetav, **thetae, **g, *Tref, *hf;

      /*--- Determine the number of heavy species ---*/
      ionization = config->GetIonization();
      if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
      else            { nHeavy = nSpecies;   nEl = 0; }

      /*--- Load variables from the config class --*/
      xi        = config->GetRotationModes();      // Rotational modes of energy storage
      Ms        = config->GetMolar_Mass();         // Species molar mass
      thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
      thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
      g         = config->GetElDegeneracy();       // Degeneracy of electron states
      nElStates = config->GetnElStates();          // Number of electron states
      Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
      hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]

      /*--- Rename & initialize for convenience ---*/
      RuSI    = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(mol*K)]
      Ru      = 1000.0*RuSI;                    // Universal gas constant [J/(kmol*K)]
      Tve     = Temperature_ve_Inf;             // Vibrational temperature [K]
      T       = Temperature_Inf;                // Translational-rotational temperature [K]
      sqvel   = 0.0;                            // Velocity^2 [m2/s2]
      rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
      rhoEve  = 0.0;                            // Mixture vib-el energy per mass [J/kg]
      denom   = 0.0;
      conc    = 0.0;
      rhoCvtr = 0.0;

      /*--- Calculate mixture density from supplied primitive quantities ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
        denom += MassFrac_Inf[iSpecies] * (Ru/Ms[iSpecies]) * T;
      for (iSpecies = 0; iSpecies < nEl; iSpecies++)
        denom += MassFrac_Inf[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Tve;
      rho = Pressure_Inf / denom;

      /*--- Calculate sound speed and extract velocities ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
        conc += MassFrac_Inf[iSpecies]*rho/Ms[iSpecies];
        rhoCvtr += rho*MassFrac_Inf[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
      }
      soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * Pressure_Inf/rho);
      for (iDim = 0; iDim < nDim; iDim++){
        sqvel += Mvec_Inf[iDim]*soundspeed * Mvec_Inf[iDim]*soundspeed;
      }
      /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
        // Species density
        rhos = MassFrac_Inf[iSpecies]*rho;

        // Species formation energy
        Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

        // Species vibrational energy
        if (thetav[iSpecies] != 0.0)
          Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
        else
          Ev = 0.0;

        // Species electronic energy
        num = 0.0;
        denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
          denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
        }
        Ee = Ru/Ms[iSpecies] * (num/denom);

        // Mixture total energy
        rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                        + Ev + Ee + Ef + 0.5*sqvel);

        // Mixture vibrational-electronic energy
        rhoEve += rhos * (Ev + Ee);
      }
      for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
        // Species formation energy
        Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];

        // Electron t-r mode contributes to mixture vib-el energy
        rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
      }

      /*--- Initialize Solution & Solution_Old vectors ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Solution[iSpecies]      = rho*MassFrac_Inf[iSpecies];
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        Solution[nSpecies+iDim] = rho*Mvec_Inf[iDim]*soundspeed;
      }
      Solution[nSpecies+nDim]     = rhoE;
      Solution[nSpecies+nDim+1]   = rhoEve;

      nodes->SetSolution(iPoint,Solution);
      nodes->SetSolution_Old(iPoint,Solution);

      counter_local++;
    }
  }

  /*--- Warning message about non-physical points ---*/
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }

  /*--- Define solver parameters needed for execution of destructor ---*/
  if (config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;


  /*--- MPI solution ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Deallocate arrays ---*/
  delete [] Mvec_Inf;
}

CTNE2EulerSolver::~CTNE2EulerSolver(void) {
  unsigned short iVar, iMarker;

  /*--- Array deallocation ---*/
  if (CD_Inv != NULL)           delete [] CD_Inv;
  if (CL_Inv != NULL)           delete [] CL_Inv;
  if (CSF_Inv != NULL)          delete [] CSF_Inv;
  if (CMx_Inv != NULL)          delete [] CMx_Inv;
  if (CMy_Inv != NULL)          delete [] CMy_Inv;
  if (CMz_Inv != NULL)          delete [] CMz_Inv;
  if (CFx_Inv != NULL)          delete [] CFx_Inv;
  if (CFy_Inv != NULL)          delete [] CFy_Inv;
  if (CFz_Inv != NULL)          delete [] CFz_Inv;
  if (Surface_CL_Inv != NULL)   delete [] Surface_CL_Inv;
  if (Surface_CD_Inv != NULL)   delete [] Surface_CD_Inv;
  if (Surface_CSF_Inv != NULL)  delete [] Surface_CSF_Inv;
  if (Surface_CEff_Inv != NULL) delete [] Surface_CEff_Inv;
  if (Surface_CFx_Inv != NULL)  delete [] Surface_CFx_Inv;
  if (Surface_CFy_Inv != NULL)  delete [] Surface_CFy_Inv;
  if (Surface_CFz_Inv != NULL)  delete [] Surface_CFz_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;

  if (lowerlimit != NULL)  		delete [] lowerlimit;
  if (upperlimit != NULL) 		delete [] upperlimit;
  if (Surface_CL_Mnt != NULL)   delete [] Surface_CL_Mnt;
  if (Surface_CD_Mnt != NULL)   delete [] Surface_CD_Mnt;
  if (Surface_CSF_Mnt != NULL)  delete [] Surface_CSF_Mnt;
  if (Surface_CEff_Mnt != NULL) delete [] Surface_CEff_Mnt;
  if (Surface_CFx_Mnt != NULL)  delete [] Surface_CFx_Mnt;
  if (Surface_CFy_Mnt != NULL)  delete [] Surface_CFy_Mnt;
  if (Surface_CFz_Mnt != NULL)  delete [] Surface_CFz_Mnt;
  if (Surface_CMx_Mnt != NULL)  delete [] Surface_CMx_Mnt;
  if (Surface_CMy_Mnt != NULL)  delete [] Surface_CMy_Mnt;
  if (Surface_CMz_Mnt != NULL)  delete [] Surface_CMz_Mnt;

  if (Surface_CL != NULL)       delete [] Surface_CL;
  if (Surface_CD != NULL)       delete [] Surface_CD;
  if (Surface_CSF != NULL)      delete [] Surface_CSF;
  if (Surface_CEff != NULL)     delete [] Surface_CEff;
  if (Surface_CFx != NULL)      delete [] Surface_CFx;
  if (Surface_CFy != NULL)      delete [] Surface_CFy;
  if (Surface_CFz != NULL)      delete [] Surface_CFz;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)         delete [] CEff_Inv;
  if (CMerit_Inv != NULL)       delete [] CMerit_Inv;
  if (CT_Inv != NULL)           delete [] CT_Inv;
  if (CQ_Inv != NULL)           delete [] CQ_Inv;
  if (CEquivArea_Inv != NULL)   delete [] CEquivArea_Inv;
  if (CNearFieldOF_Inv != NULL) delete [] CNearFieldOF_Inv;

  if (Primitive != NULL)        delete [] Primitive;
  if (Primitive_i != NULL)      delete [] Primitive_i;
  if (Primitive_j != NULL)      delete [] Primitive_j;

  if (Secondary != NULL)        delete [] Secondary;
  if (Secondary_i != NULL)      delete [] Secondary_i;
  if (Secondary_j != NULL)      delete [] Secondary_j;

  if (LowMach_Precontioner != NULL) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete [] LowMach_Precontioner[iVar];
    delete [] LowMach_Precontioner;
  }

  if (CPressure != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CPressure[iMarker];
    delete [] CPressure;
  }

  if (HeatFlux != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete HeatFlux[iMarker];
    }
    delete [] HeatFlux;
  }
  if (nVertex != NULL) delete [] nVertex;
}

void CTNE2EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  unsigned long iPoint;
  unsigned short iMesh;
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == TNE2_RANS) ||
               (config->GetKind_Solver() == DISC_ADJ_TNE2_RANS));
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));


  /*--- Make sure that the solution is well initialized for unsteady
   calculations with dual time-stepping (load additional restarts for 2nd-order). ---*/

  if (dual_time && (TimeIter == 0 || (restart && (long)TimeIter == (long)config->GetRestart_Iter()))) {

    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][TNE2_SOL]->GetNodes()->Set_Solution_time_n();
        solver_container[iMesh][TNE2_SOL]->GetNodes()->Set_Solution_time_n1();
        if (rans) {
          solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
          solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1();
        }
      }
    }

    if ((restart && (long)TimeIter == config->GetRestart_Iter()) &&
        (config->GetTime_Marching() == DT_STEPPING_2ND)) {

      /*--- Load an additional restart file for a 2nd-order restart ---*/
      solver_container[MESH_0][TNE2_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetRestart_Iter()-1), true);

      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetRestart_Iter()-1), false);

      /*--- Push back this new solution to time level N. ---*/

      for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][TNE2_SOL]->GetNodes()->Set_Solution_time_n();
          if (rans) {
            solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
          }
        }
      }
    }
  }
}

void CTNE2EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solution_container,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep,
                                     unsigned short RunTime_EqSystem, bool Output) {
  unsigned long iPoint, ErrorCounter = 0;

  unsigned long InnerIter = config->GetInnerIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool muscl            = config->GetMUSCL_TNE2();
  bool limiter          = ((config->GetKind_SlopeLimit_TNE2() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool center           = config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED;
  bool center_jst       = center && (config->GetKind_Centered_TNE2() == JST);
  bool van_albada       = config->GetKind_SlopeLimit_TNE2() == VAN_ALBADA_EDGE;
  bool nonPhys;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the non-physical points vector ---*/
    
    nodes->SetNon_Physical(iPoint,false);

    /*--- Primitive variables [rho1,...,rhoNs,T,Tve,u,v,w,P,rho,h,c] ---*/
    nonPhys = nodes->SetPrimVar_Compressible(iPoint, config);
    if (nonPhys) {
      nodes->SetNon_Physical(iPoint,true);
      ErrorCounter++;
    }

    /*--- Initialize the convective residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);

  }

  /*--- Allowing for Primitive Variables to be passed ---*/
  // if(ErrorCounter > 0) {
  //   InitiateComms(geometry, config, SOLUTION);
  //   CompleteComms(geometry, config, SOLUTION);
  // }
  // InitiateComms(geometry, config, PRIMITIVE);
  // CompleteComms(geometry, config, PRIMITIVE);

  /*--- Upwind second order reconstruction ---*/
  if ((muscl && !center) && (iMesh == MESH_0) && !Output) {

    /*--- Calculate the gradients ---*/
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
      SetSolution_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
      SetSolution_Gradient_LS(geometry, config);
    }

    /*--- Limiter computation ---*/
    if ((limiter) && (iMesh == MESH_0) && !Output && !van_albada) {
      SetPrimitive_Limiter(geometry, config);
    }

  }

  /*--- Artificial dissipation ---*/
  if (center && !Output){
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Initialize the jacobian matrices ---*/
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

void CTNE2EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solution_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time,
      Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool time_steping = config->GetTime_Marching() == TIME_STEPPING;
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));

  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetMax_Lambda_Inv(iPoint, 0.0);

  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (nodes->GetProjVel(iPoint, Normal) + nodes->GetProjVel(jPoint,Normal));
    Mean_SoundSpeed = 0.5 * (nodes->GetSoundSpeed(iPoint) + nodes->GetSoundSpeed(jPoint)) * Area;

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
    if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(iPoint,Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(jPoint,Lambda);

  }

  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- Point identification, Normal vector and area ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

        /*--- Mean Values ---*/
        Mean_ProjVel = nodes->GetProjVel(iPoint,Normal);
        Mean_SoundSpeed = nodes->GetSoundSpeed(iPoint) * Area;

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
          nodes->AddMax_Lambda_Inv(iPoint,Lambda);
        }

      }
  }

  /*--- Each element uses their own speed, steady state simulation ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();

    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / nodes->GetMax_Lambda_Inv(iPoint);
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
  if (time_steping) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Sets the regular CFL equal to the unsteady CFL ---*/
      config->SetCFL(iMesh,config->GetUnst_CFL());

      /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
             it computes the time step based on the unsteady CFL ---*/
      if (config->GetCFL(iMesh) == 0.0) {
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
        nodes->SetDelta_Time(iPoint, Local_Delta_Time);
      }
    }

}

void CTNE2EulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {

  su2double *Normal, Area, Mean_SoundSpeed, Mean_ProjVel, Lambda;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetLambda(iPoint, 0.0);

  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
    Mean_SoundSpeed = 0.5 * (nodes->GetSoundSpeed(iPoint) + nodes->GetSoundSpeed(jPoint)) * Area;

    /*--- Inviscid contribution ---*/
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) nodes->AddLambda(iPoint,Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddLambda(jPoint,Lambda);

  }

  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- Point identification, Normal vector and area ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

        /*--- Mean Values ---*/
        Mean_ProjVel = nodes->GetProjVel(iPoint,Normal);
        Mean_SoundSpeed = nodes->GetSoundSpeed(iPoint) * Area;

        /*--- Inviscid contribution ---*/
        Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
        if (geometry->node[iPoint]->GetDomain()) {
          nodes->AddLambda(iPoint,Lambda);
        }
      }
  }

  /*--- Call the MPI routine ---*/
  InitiateComms(geometry, config, MAX_EIGENVALUE);
  CompleteComms(geometry, config, MAX_EIGENVALUE);

}

void CTNE2EulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, jVar;
  bool err;

  /*--- Set booleans based on config settings ---*/
  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);

  //Unused at the moment
  //bool centered = ((config->GetKind_Centered_TNE2() == JST) && (iMesh == MESH_0));

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  numerics->SetPIndex      ( nodes->GetPIndex()       );
  numerics->SetTIndex      ( nodes->GetTIndex()       );
  numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  numerics->SetHIndex      ( nodes->GetHIndex()       );
  numerics->SetAIndex      ( nodes->GetAIndex()       );
  numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge, set normal vectors, and number of neighbors ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(),
                          geometry->node[jPoint]->GetnNeighbor());

    /*--- Pass conservative & primitive variables w/o reconstruction to CNumerics ---*/
    numerics->SetConservative(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));
    numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(jPoint));

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(nodes->GetdPdU(iPoint), nodes->GetdPdU(jPoint));
    numerics->SetdTdU(nodes->GetdTdU(iPoint), nodes->GetdTdU(jPoint));
    numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));

    /*--- Set the largest convective eigenvalue ---*/
    numerics->SetLambda(nodes->GetLambda(iPoint), nodes->GetLambda(jPoint));

    /*--- Compute residuals, and Jacobians ---*/
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if (Res_Conv[iVar] != Res_Conv[iVar])
        err = true;
    if (implicit)
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          if ((Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) ||
              (Jacobian_j[iVar][jVar] != Jacobian_j[iVar][jVar])   )
            err = true;

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(jPoint, Res_Conv);
      if (implicit) {
        Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
        Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
        Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
        Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
      }
    }
  }
}

void CTNE2EulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solution_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh) {

  unsigned long iEdge, iPoint, jPoint, counter_local=0, counter_global=0;
  unsigned long InnerIter = config->GetInnerIter();
  unsigned short RHO_INDEX, RHOS_INDEX, P_INDEX, TVE_INDEX;
  unsigned short iDim, iVar, jVar, iSpecies;

  su2double *U_i, *U_j, *V_i, *V_j;
  su2double **GradU_i, **GradU_j;
  su2double **GradV_i, **GradV_j;
  su2double ProjGradU_i, ProjGradU_j;
  su2double ProjGradV_i, ProjGradV_j;
  su2double *Limiter_i, *Limiter_j;
  su2double *Conserved_i, *Conserved_j, *Primitive_i, *Primitive_j;
  su2double *dPdU_i, *dPdU_j, *dTdU_i, *dTdU_j, *dTvedU_i, *dTvedU_j;
  su2double *Eve_i, *Eve_j, *Cvve_i, *Cvve_j;
  su2double lim_i, lim_j;
  su2double Non_Physical = 1.0;


  /*--- Set booleans based on config settings ---*/
  bool implicit   	= (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool muscl    	  = (config->GetMUSCL_TNE2() && (iMesh == MESH_0));
  bool disc_adjoint = config->GetDiscrete_Adjoint();
  bool limiter      = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool chk_err_i, chk_err_j, err;

  /*--- Allocate arrays ---*/
  Primitive_i = new su2double[nPrimVar];
  Primitive_j = new su2double[nPrimVar];
  Conserved_i = new su2double[nVar];
  Conserved_j = new su2double[nVar];
  dPdU_i      = new su2double[nVar];
  dPdU_j      = new su2double[nVar];
  dTdU_i      = new su2double[nVar];
  dTdU_j      = new su2double[nVar];
  dTvedU_i    = new su2double[nVar];
  dTvedU_j    = new su2double[nVar];
  Eve_i       = new su2double[nSpecies];
  Eve_j       = new su2double[nSpecies];
  Cvve_i      = new su2double[nSpecies];
  Cvve_j      = new su2double[nSpecies];

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  numerics->SetPIndex      ( nodes->GetPIndex()       );
  numerics->SetTIndex      ( nodes->GetTIndex()       );
  numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  numerics->SetHIndex      ( nodes->GetHIndex()       );
  numerics->SetAIndex      ( nodes->GetAIndex()       );
  numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  RHO_INDEX  = nodes->GetRhoIndex();
  RHOS_INDEX = nodes->GetRhosIndex();
  P_INDEX    = nodes->GetPIndex();
  TVE_INDEX  = nodes->GetTveIndex();

  /*--- Loop over edges and calculate convective fluxes ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Retrieve node numbers and pass edge normal to CNumerics ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Get conserved & primitive variables from CVariable ---*/
    U_i = nodes->GetSolution(iPoint);   U_j = nodes->GetSolution(jPoint);
    V_i = nodes->GetPrimitive(iPoint);  V_j = nodes->GetPrimitive(jPoint);

    /*--- High order reconstruction using MUSCL strategy ---*/
    if (muscl) {

      /*--- Assign i-j and j-i to projection vectors ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) -
                              geometry->node[iPoint]->GetCoord(iDim)   );
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) -
                              geometry->node[jPoint]->GetCoord(iDim)   );
      }

      /*---+++ Conserved variable reconstruction & limiting +++---*/

      /*--- Retrieve gradient information & limiter ---*/
      GradU_i = nodes->GetGradient(iPoint);
      GradU_j = nodes->GetGradient(jPoint);

      if (limiter) {
        Limiter_i = nodes->GetLimiter_Primitive(iPoint);
        Limiter_j = nodes->GetLimiter_Primitive(jPoint);
      }

      /*--- Reconstruct conserved variables at the edge interface ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        ProjGradU_i = 0.0; ProjGradU_j = 0.0;
        Non_Physical = nodes->GetNon_Physical(iPoint)*nodes->GetNon_Physical(jPoint);
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjGradU_i += Vector_i[iDim]*GradU_i[iVar][iDim];
          ProjGradU_j += Vector_j[iDim]*GradU_j[iVar][iDim];
        }
        if (limiter) {
          Conserved_i[iVar] = U_i[iVar] + Limiter_i[iVar]*ProjGradU_i;
          Conserved_j[iVar] = U_j[iVar] + Limiter_j[iVar]*ProjGradU_j;
        }
        else {
          Conserved_i[iVar] = U_i[iVar] + ProjGradU_i;
          Conserved_j[iVar] = U_j[iVar] + ProjGradU_j;
        }
      }

      chk_err_i = nodes->Cons2PrimVar(config, Conserved_i, Primitive_i,
                                      dPdU_i, dTdU_i, dTvedU_i, Eve_i, Cvve_i);
      chk_err_j = nodes->Cons2PrimVar(config, Conserved_j, Primitive_j,
                                      dPdU_j, dTdU_j, dTvedU_j, Eve_j, Cvve_j);

      /*--- Check for physical solutions in the reconstructed values ---*/
      // Note: If non-physical, revert to first order
      if ( chk_err_i || chk_err_j) {
        numerics->SetPrimitive   (V_i, V_j);
        numerics->SetConservative(U_i, U_j);
        numerics->SetdPdU  (nodes->GetdPdU(iPoint),   nodes->GetdPdU(jPoint));
        numerics->SetdTdU  (nodes->GetdTdU(iPoint),   nodes->GetdTdU(jPoint));
        numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));
        numerics->SetEve   (nodes->GetEve(iPoint),    nodes->GetEve(jPoint));
        numerics->SetCvve  (nodes->GetCvve(iPoint),   nodes->GetCvve(jPoint));
        if(chk_err_i && !chk_err_j) cout << "Non-physical reconstruction at " << iPoint << "." << endl;
        if(chk_err_j && !chk_err_i) cout << "Non-physical reconstruction at " << jPoint << "." << endl;
        if(chk_err_i && chk_err_j) cout << "Non-physical reconstruction at (" << iPoint << ", " << jPoint << ")." << endl;
      } else {
        numerics->SetConservative(Conserved_i, Conserved_j);
        numerics->SetPrimitive   (Primitive_i, Primitive_j);
        numerics->SetdPdU  (dPdU_i,   dPdU_j  );
        numerics->SetdTdU  (dTdU_i,   dTdU_j  );
        numerics->SetdTvedU(dTvedU_i, dTvedU_j);
        numerics->SetEve   (Eve_i,    Eve_j   );
        numerics->SetCvve  (Cvve_i,   Cvve_j  );
      }

    } else {

      /*--- Set variables without reconstruction ---*/
      numerics->SetPrimitive   (V_i, V_j);
      numerics->SetConservative(U_i, U_j);
      numerics->SetdPdU  (nodes->GetdPdU(iPoint),   nodes->GetdPdU(jPoint));
      numerics->SetdTdU  (nodes->GetdTdU(iPoint),   nodes->GetdTdU(jPoint));
      numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));
      numerics->SetEve   (nodes->GetEve(iPoint),    nodes->GetEve(jPoint));
      numerics->SetCvve  (nodes->GetCvve(iPoint),   nodes->GetCvve(jPoint));
    }

    /*--- Compute the upwind residual ---*/
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if (Res_Conv[iVar] != Res_Conv[iVar]) err = true;
    if (implicit)
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          if ((Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) ||
              (Jacobian_j[iVar][jVar] != Jacobian_j[iVar][jVar])   )
            err = true;

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(jPoint, Res_Conv);
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
      }
    }
  }

  /*--- Warning message about non-physical reconstructions ---*/
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);
  }


  delete [] Conserved_i;
  delete [] Conserved_j;
  delete [] Primitive_i;
  delete [] Primitive_j;
  delete [] dPdU_i;
  delete [] dPdU_j;
  delete [] dTdU_i;
  delete [] dTdU_j;
  delete [] dTvedU_i;
  delete [] dTvedU_j;
  //delete [] Eve_i;
  //delete [] Eve_j;
  delete [] Cvve_i;
  delete [] Cvve_j;
}

void CTNE2EulerSolver::Source_Residual(CGeometry *geometry, CSolver **solution_container, CNumerics *numerics,
                                       CNumerics *second_solver, CConfig *config, unsigned short iMesh) {

  unsigned short iVar, jVar;
  unsigned long iPoint;
  unsigned long eAxi_local, eChm_local, eVib_local;
  unsigned long eAxi_global, eChm_global, eVib_global;
  int rank = MASTER_NODE;

  /*--- Assign booleans ---*/
  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool err = false;

  /*--- Initialize the error counter ---*/
  eAxi_local = 0;
  eChm_local = 0;
  eVib_local = 0;

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  numerics->SetPIndex      ( nodes->GetPIndex()       );
  numerics->SetTIndex      ( nodes->GetTIndex()       );
  numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  numerics->SetHIndex      ( nodes->GetHIndex()       );
  numerics->SetAIndex      ( nodes->GetAIndex()       );
  numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );


  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  for (iVar = 0; iVar < nVar; iVar++) Source[iVar] = 0.0;

  /*--- loop over interior points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Set conserved & primitive variables  ---*/
    numerics->SetConservative(nodes->GetSolution(iPoint),   nodes->GetSolution(iPoint));
    numerics->SetPrimitive   (nodes->GetPrimitive(iPoint),  nodes->GetPrimitive(iPoint) );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(nodes->GetdPdU(iPoint), nodes->GetdPdU(iPoint));
    numerics->SetdTdU(nodes->GetdTdU(iPoint), nodes->GetdTdU(iPoint));
    numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(iPoint));
    numerics->SetEve(nodes->GetEve(iPoint), nodes->GetEve(iPoint));
    numerics->SetCvve(nodes->GetCvve(iPoint), nodes->GetCvve(iPoint));

    /*--- Set volume of the dual grid cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[iPoint]->GetCoord() );


    /*--- Compute axisymmetric source terms (if needed) ---*/
    if (config->GetAxisymmetric()) {
      numerics->ComputeAxisymmetric(Residual, Source, Jacobian_i, config);

      /*--- Check for errors before applying source to the linear system ---*/
      err = false;
      for (iVar = 0; iVar < nVar; iVar++)
        if (Residual[iVar] != Residual[iVar]) err = true;
      if (implicit)
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) err = true;

      /*--- Apply the update to the linear system ---*/
      if (!err) {
        LinSysRes.AddBlock(iPoint, Residual);
        if (implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        eAxi_local++;
    }

    /*--- Compute the non-equilibrium chemistry ---*/
    numerics->ComputeChemistry(Residual, Source, Jacobian_i, config);

    /*--- Check for errors before applying source to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if (Residual[iVar] != Residual[iVar]) err = true;
    if (implicit)
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) err = true;

    /*--- Apply the chemical sources to the linear system ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, Residual);
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    } else
      eChm_local++;

    /*--- Compute vibrational energy relaxation ---*/
    // NOTE: Jacobians don't account for relaxation time derivatives
    numerics->ComputeVibRelaxation(Residual, Source, Jacobian_i, config);

    /*--- Check for errors before applying source to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if (Residual[iVar] != Residual[iVar]) err = true;
    if (implicit)
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) err = true;

    /*--- Apply the vibrational relaxation terms to the linear system ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, Residual);
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    } else
      eVib_local++;

    /*---Set source term ---*/
    nodes->SetSource(iPoint, Source);
  }

  /*--- Checking for NaN ---*/
  eAxi_global = eAxi_local;
  eChm_global = eChm_local;
  eVib_global = eVib_local;

  //THIS IS NO FUN
  if ((rank != MASTER_NODE) &&
      (
        (eAxi_global != 0) ||
        (eChm_global != 0) ||
        (eVib_global != 0)
        )
      ) {
    cout << "Warning!! Instances of NaN in the following source terms: " << endl;
    cout << "Axisymmetry: " << eAxi_global << endl;
    cout << "Chemical:    " << eChm_global << endl;
    cout << "Vib. Relax:  " << eVib_global << endl;
  }
}

void CTNE2EulerSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double Pressure = 0.0, *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord,
  factor, NFPressOF, RefVel2, RefTemp, RefDensity, RefPressure, Mach2Vel, Mach_Motion,
  Force[3] = {0.0,0.0,0.0};
  string Marker_Tag, Monitoring_Tag;
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0},
  MomentZ_Force[3] = {0.0,0.0,0.0};
  su2double AxiFactor;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Inv, MyAllBound_CL_Inv, MyAllBound_CSF_Inv, MyAllBound_CMx_Inv,
  MyAllBound_CMy_Inv, MyAllBound_CMz_Inv, MyAllBound_CoPx_Inv, MyAllBound_CoPy_Inv,
  MyAllBound_CoPz_Inv, MyAllBound_CFx_Inv, MyAllBound_CFy_Inv, MyAllBound_CFz_Inv,
  MyAllBound_CT_Inv, MyAllBound_CQ_Inv, MyAllBound_CNearFieldOF_Inv, *MySurface_CL_Inv = NULL,
  *MySurface_CD_Inv = NULL, *MySurface_CSF_Inv = NULL, *MySurface_CEff_Inv = NULL,
  *MySurface_CFx_Inv = NULL, *MySurface_CFy_Inv = NULL, *MySurface_CFz_Inv = NULL,
  *MySurface_CMx_Inv = NULL, *MySurface_CMy_Inv = NULL, *MySurface_CMz_Inv = NULL;
#endif

  su2double Alpha     = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta      = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea   = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  su2double *Origin   = NULL;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }
  bool axisymmetric         = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream
   values, which is the standard convention. ---*/

  /*--- Evaluate reference values for non-dimensionalization. ---*/
  RefVel2     = node_infty->GetVelocity2(0);
  RefDensity  = node_infty->GetDensity(0);
  RefPressure = node_infty->GetPressure(0);
  factor      = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*-- Initialization ---*/
  Total_CD           = 0.0; Total_CL   = 0.0; Total_CSF     = 0.0; Total_CEff = 0.0;
  Total_CMx          = 0.0; Total_CMy  = 0.0; Total_CMz     = 0.0;
  Total_CoPx         = 0.0; Total_CoPy = 0.0; Total_CoPz    = 0.0;
  Total_CFx          = 0.0; Total_CFy  = 0.0; Total_CFz     = 0.0;
  Total_CT           = 0.0; Total_CQ   = 0.0; Total_CMerit  = 0.0;
  Total_CNearFieldOF = 0.0; Total_Heat = 0.0; Total_MaxHeat = 0.0;

  AllBound_CD_Inv           = 0.0; AllBound_CL_Inv   = 0.0; AllBound_CSF_Inv    = 0.0;
  AllBound_CMx_Inv          = 0.0; AllBound_CMy_Inv  = 0.0; AllBound_CMz_Inv    = 0.0;
  AllBound_CoPx_Inv         = 0.0; AllBound_CoPy_Inv = 0.0; AllBound_CoPz_Inv   = 0.0;
  AllBound_CFx_Inv          = 0.0; AllBound_CFy_Inv  = 0.0; AllBound_CFz_Inv    = 0.0;
  AllBound_CT_Inv           = 0.0; AllBound_CQ_Inv   = 0.0; AllBound_CMerit_Inv = 0.0;
  AllBound_CNearFieldOF_Inv = 0.0; AllBound_CEff_Inv = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Inv[iMarker_Monitoring]  = 0.0; Surface_CD_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring] = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring] = 0.0; Surface_CFy_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring] = 0.0; Surface_CMx_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring] = 0.0; Surface_CMz_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CL[iMarker_Monitoring]      = 0.0; Surface_CD[iMarker_Monitoring]       = 0.0;
    Surface_CSF[iMarker_Monitoring]     = 0.0; Surface_CEff[iMarker_Monitoring]     = 0.0;
    Surface_CFx[iMarker_Monitoring]     = 0.0; Surface_CFy[iMarker_Monitoring]      = 0.0;
    Surface_CFz[iMarker_Monitoring]     = 0.0; Surface_CMx[iMarker_Monitoring]      = 0.0;
    Surface_CMy[iMarker_Monitoring]     = 0.0; Surface_CMz[iMarker_Monitoring]      = 0.0;
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

    if ((Boundary == EULER_WALL)              || (Boundary == HEAT_FLUX)               ||
        (Boundary == HEAT_FLUX_CATALYTIC)     || (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
        (Boundary == ISOTHERMAL)              || (Boundary == ISOTHERMAL_CATALYTIC)    ||
        (Boundary == ISOTHERMAL_NONCATALYTIC) || (Boundary == NEARFIELD_BOUNDARY)) {

      /*--- Force initialization on each marker ---*/
      CD_Inv[iMarker]   = 0.0;         CL_Inv[iMarker]   = 0.0; CSF_Inv[iMarker]    = 0.0;
      CMx_Inv[iMarker]  = 0.0;         CMy_Inv[iMarker]  = 0.0; CMz_Inv[iMarker]    = 0.0;
      CoPx_Inv[iMarker] = 0.0;         CoPy_Inv[iMarker] = 0.0; CoPz_Inv[iMarker]   = 0.0;
      CFx_Inv[iMarker]  = 0.0;         CFy_Inv[iMarker]  = 0.0; CFz_Inv[iMarker]    = 0.0;
      CT_Inv[iMarker]   = 0.0;         CQ_Inv[iMarker]   = 0.0; CMerit_Inv[iMarker] = 0.0;
      CNearFieldOF_Inv[iMarker] = 0.0; CEff_Inv[iMarker] = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
      MomentInviscid[0] = 0.0;  MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
      MomentX_Force[0]  = 0.0;  MomentX_Force[1]  = 0.0; MomentX_Force[2]  = 0.0;
      MomentY_Force[0]  = 0.0;  MomentY_Force[1]  = 0.0; MomentY_Force[2]  = 0.0;
      MomentZ_Force[0]  = 0.0;  MomentZ_Force[1]  = 0.0; MomentZ_Force[2]  = 0.0;

      NFPressOF = 0.0;

      /*--- Loop over vertices to compute forces ---*/
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Calculate pressure ---*/
        Pressure = nodes->GetPressure(iPoint);
        CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefArea;

        /*--- Note that the pressure coefficient is computed at the
        halo cells (for visualization purposes), but not the forces ---*/
        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();

          /*--- Quadratic objective function for the near field.
          This uses the infinity pressure regardless of Mach number. ---*/
          NFPressOF += 0.5*((Pressure - Pressure_Inf)*
                            (Pressure - Pressure_Inf))*Normal[nDim-1];

          for (iDim = 0; iDim < nDim; iDim++) {
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

          /*--- Axisymmetric simulations ---*/
          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Compute force, note minus sign due to outward orientation of
           the normal vector ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf)*Normal[iDim]*factor*AxiFactor;
            ForceInviscid[iDim] += Force[iDim];
          }

          /*--- Compute moment w.r.t. reference axis ---*/
          if (iDim == 3) {
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


      /*--- Project forces and store non-dimensional coefficients ---*/
      if (Monitoring == YES) {
        if (Boundary != NEARFIELD_BOUNDARY) {
          if (nDim == 2) {
            CD_Inv[iMarker]     =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            CL_Inv[iMarker]     = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
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
            CSF_Inv[iMarker]     = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            CEff_Inv[iMarker]    = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
            CMx_Inv[iMarker]     = MomentInviscid[0];
            CMy_Inv[iMarker]     = MomentInviscid[1];
            CMz_Inv[iMarker]     = MomentInviscid[2];
            CoPx_Inv[iMarker]    = -MomentY_Force[0];
            CoPz_Inv[iMarker]    = MomentY_Force[2];
            CFx_Inv[iMarker]     = ForceInviscid[0];
            CFy_Inv[iMarker]     = ForceInviscid[1];
            CFz_Inv[iMarker]     = ForceInviscid[2];
            CT_Inv[iMarker]      = -CFz_Inv[iMarker];
            CQ_Inv[iMarker]      = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker]  = CT_Inv[iMarker] / (CQ_Inv[iMarker] + EPS);
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
              Surface_CL_Inv[iMarker_Monitoring]    += CL_Inv[iMarker];
              Surface_CD_Inv[iMarker_Monitoring]    += CD_Inv[iMarker];
              Surface_CSF_Inv[iMarker_Monitoring]   += CSF_Inv[iMarker];
              Surface_CEff_Inv[iMarker_Monitoring]   = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
              Surface_CFx_Inv[iMarker_Monitoring]   += CFx_Inv[iMarker];
              Surface_CFy_Inv[iMarker_Monitoring]   += CFy_Inv[iMarker];
              Surface_CFz_Inv[iMarker_Monitoring]   += CFz_Inv[iMarker];
              Surface_CMx_Inv[iMarker_Monitoring]   += CMx_Inv[iMarker];
              Surface_CMy_Inv[iMarker_Monitoring]   += CMy_Inv[iMarker];
              Surface_CMz_Inv[iMarker_Monitoring]   += CMz_Inv[iMarker];
            }
          }
        }

        /*--- At the Nearfield SU2 only cares about the pressure coeffient ---*/
        else {
          CNearFieldOF_Inv[iMarker]  = NFPressOF;
          AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
        }
      }
    }
  }

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/
  MyAllBound_CD_Inv        = AllBound_CD_Inv;    AllBound_CD_Inv = 0.0;
  MyAllBound_CL_Inv        = AllBound_CL_Inv;    AllBound_CL_Inv = 0.0;
  MyAllBound_CSF_Inv       = AllBound_CSF_Inv;   AllBound_CSF_Inv = 0.0;
  AllBound_CEff_Inv        = 0.0;
  MyAllBound_CMx_Inv       = AllBound_CMx_Inv;   AllBound_CMx_Inv = 0.0;
  MyAllBound_CMy_Inv       = AllBound_CMy_Inv;   AllBound_CMy_Inv = 0.0;
  MyAllBound_CMz_Inv       = AllBound_CMz_Inv;   AllBound_CMz_Inv = 0.0;
  MyAllBound_CoPx_Inv      = AllBound_CoPx_Inv;  AllBound_CoPx_Inv = 0.0;
  MyAllBound_CoPy_Inv      = AllBound_CoPy_Inv;  AllBound_CoPy_Inv = 0.0;
  MyAllBound_CoPz_Inv      = AllBound_CoPz_Inv;  AllBound_CoPz_Inv = 0.0;
  MyAllBound_CFx_Inv       = AllBound_CFx_Inv;   AllBound_CFx_Inv = 0.0;
  MyAllBound_CFy_Inv       = AllBound_CFy_Inv;   AllBound_CFy_Inv = 0.0;
  MyAllBound_CFz_Inv       = AllBound_CFz_Inv;   AllBound_CFz_Inv = 0.0;
  MyAllBound_CT_Inv        = AllBound_CT_Inv;    AllBound_CT_Inv = 0.0;
  MyAllBound_CQ_Inv        = AllBound_CQ_Inv;    AllBound_CQ_Inv = 0.0;
  AllBound_CMerit_Inv      = 0.0;
  MyAllBound_CNearFieldOF_Inv = AllBound_CNearFieldOF_Inv; AllBound_CNearFieldOF_Inv = 0.0;

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
    SU2_MPI::Allreduce(&MyAllBound_CNearFieldOF_Inv, &AllBound_CNearFieldOF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  /*--- Add the forces on the surfaces using all the nodes ---*/
  MySurface_CL_Inv   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Inv   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Inv  = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Inv = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Inv  = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Inv  = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Inv  = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Inv  = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Inv  = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Inv  = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CL_Inv[iMarker_Monitoring]   = Surface_CL_Inv[iMarker_Monitoring];
    MySurface_CD_Inv[iMarker_Monitoring]   = Surface_CD_Inv[iMarker_Monitoring];
    MySurface_CSF_Inv[iMarker_Monitoring]  = Surface_CSF_Inv[iMarker_Monitoring];
    MySurface_CEff_Inv[iMarker_Monitoring] = Surface_CEff_Inv[iMarker_Monitoring];
    MySurface_CFx_Inv[iMarker_Monitoring]  = Surface_CFx_Inv[iMarker_Monitoring];
    MySurface_CFy_Inv[iMarker_Monitoring]  = Surface_CFy_Inv[iMarker_Monitoring];
    MySurface_CFz_Inv[iMarker_Monitoring]  = Surface_CFz_Inv[iMarker_Monitoring];
    MySurface_CMx_Inv[iMarker_Monitoring]  = Surface_CMx_Inv[iMarker_Monitoring];
    MySurface_CMy_Inv[iMarker_Monitoring]  = Surface_CMy_Inv[iMarker_Monitoring];
    MySurface_CMz_Inv[iMarker_Monitoring]  = Surface_CMz_Inv[iMarker_Monitoring];

    Surface_CL_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CD_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CEff_Inv[iMarker_Monitoring] = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CFy_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]  = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]  = 0.0;
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
  Total_CD           = AllBound_CD_Inv;
  Total_CL           = AllBound_CL_Inv;
  Total_CSF          = AllBound_CSF_Inv;
  Total_CEff         = Total_CL / (Total_CD + EPS);
  Total_CFx          = AllBound_CFx_Inv;
  Total_CFy          = AllBound_CFy_Inv;
  Total_CFz          = AllBound_CFz_Inv;
  Total_CMx          = AllBound_CMx_Inv;
  Total_CMy          = AllBound_CMy_Inv;
  Total_CMz          = AllBound_CMz_Inv;
  Total_CoPx         = AllBound_CoPx_Inv;
  Total_CoPy         = AllBound_CoPy_Inv;
  Total_CoPz         = AllBound_CoPz_Inv;
  Total_CT           = AllBound_CT_Inv;
  Total_CQ           = AllBound_CQ_Inv;
  Total_CMerit       = Total_CT / (Total_CQ + EPS);
  Total_CNearFieldOF = AllBound_CNearFieldOF_Inv;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]   = Surface_CL_Inv[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]   = Surface_CD_Inv[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring]  = Surface_CSF_Inv[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring] = Surface_CL_Inv[iMarker_Monitoring] / (Surface_CD_Inv[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]  = Surface_CFx_Inv[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]  = Surface_CFy_Inv[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]  = Surface_CFz_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]  = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]  = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]  = Surface_CMz_Inv[iMarker_Monitoring];
  }
}

void CTNE2EulerSolver::Momentum_Forces(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  su2double *Normal = NULL, MomentDist[3] = {0.0,0.0,0.0}, *Coord, Area,
      factor, RefVel2, RefTemp, RefPressure, RefDensity,  Mach2Vel, Mach_Motion,
      Force[3] = {0.0,0.0,0.0}, Velocity[3], MassFlow, Density;
  string Marker_Tag, Monitoring_Tag;
  su2double MomentX_Force[3] = {0.0,0.0,0.0}, MomentY_Force[3] = {0.0,0.0,0.0},
  MomentZ_Force[3] = {0.0,0.0,0.0};
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

  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double RefArea     = config->GetRefArea();
  su2double RefLength  = config->GetRefLength();
  su2double *Origin = NULL;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }
  bool axisymmetric          = config->GetAxisymmetric();

  /*--- Evaluate reference values for non-dimensionalization.
   For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  RefVel2     = node_infty->GetVelocity2(0);
  RefDensity  = node_infty->GetDensity(0);
  RefPressure = node_infty->GetPressure(0);
  factor      = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*-- Variables initialization ---*/
  AllBound_CD_Mnt = 0.0;        AllBound_CL_Mnt = 0.0; AllBound_CSF_Mnt = 0.0;
  AllBound_CMx_Mnt = 0.0;          AllBound_CMy_Mnt = 0.0;   AllBound_CMz_Mnt = 0.0;
  AllBound_CoPx_Mnt = 0.0;          AllBound_CoPy_Mnt = 0.0;   AllBound_CoPz_Mnt = 0.0;
  AllBound_CFx_Mnt = 0.0;          AllBound_CFy_Mnt = 0.0;   AllBound_CFz_Mnt = 0.0;
  AllBound_CT_Mnt = 0.0;           AllBound_CQ_Mnt = 0.0;    AllBound_CMerit_Mnt = 0.0;
  AllBound_CEff_Mnt = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Mnt[iMarker_Monitoring]  = 0.0; Surface_CD_Mnt[iMarker_Monitoring]   = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring] = 0.0; Surface_CEff_Mnt[iMarker_Monitoring] = 0.0;
    Surface_CFx_Mnt[iMarker_Monitoring] = 0.0; Surface_CFy_Mnt[iMarker_Monitoring]  = 0.0;
    Surface_CFz_Mnt[iMarker_Monitoring] = 0.0;
    Surface_CMx_Mnt[iMarker_Monitoring] = 0.0; Surface_CMy_Mnt[iMarker_Monitoring]  = 0.0; Surface_CMz_Mnt[iMarker_Monitoring] = 0.0;
  }

  /*--- Loop over the Inlet -Outlet Markers  ---*/

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

      CD_Mnt[iMarker]   = 0.0;  CL_Mnt[iMarker]   = 0.0;  CSF_Mnt[iMarker]    = 0.0;
      CFx_Mnt[iMarker]  = 0.0;  CFy_Mnt[iMarker]  = 0.0;  CFz_Mnt[iMarker]    = 0.0;
      CMx_Mnt[iMarker]  = 0.0;  CMy_Mnt[iMarker]  = 0.0;  CMz_Mnt[iMarker]    = 0.0;
      CoPx_Mnt[iMarker] = 0.0;  CoPy_Mnt[iMarker] = 0.0;  CoPz_Mnt[iMarker]   = 0.0;
      CT_Mnt[iMarker]   = 0.0;  CQ_Mnt[iMarker]   = 0.0;  CMerit_Mnt[iMarker] = 0.0;
      CEff_Mnt[iMarker] = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) ForceMomentum[iDim] = 0.0;
      MomentMomentum[0] = 0.0; MomentMomentum[1] = 0.0; MomentMomentum[2] = 0.0;
      MomentX_Force[0] = 0.0;  MomentX_Force[1] = 0.0;  MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0;  MomentY_Force[1] = 0.0;  MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0;  MomentZ_Force[1] = 0.0;  MomentZ_Force[2] = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/

        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {

          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          Density   = nodes->GetDensity(iPoint);

          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/

          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

          MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim]  = nodes->GetVelocity(iPoint,iDim);
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
          CD_Mnt[iMarker]     =  ForceMomentum[0]*cos(Alpha) + ForceMomentum[1]*sin(Alpha);
          CL_Mnt[iMarker]     = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[1]*cos(Alpha);
          CEff_Mnt[iMarker]   = CL_Mnt[iMarker] / (CD_Mnt[iMarker]+EPS);
          CFx_Mnt[iMarker]    = ForceMomentum[0];
          CFy_Mnt[iMarker]    = ForceMomentum[1];
          CMz_Mnt[iMarker]    = MomentMomentum[2];
          CoPx_Mnt[iMarker]   = MomentZ_Force[1];
          CoPy_Mnt[iMarker]   = -MomentZ_Force[0];
          CT_Mnt[iMarker]     = -CFx_Mnt[iMarker];
          CQ_Mnt[iMarker]     = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker] = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }
        if (nDim == 3) {
          CD_Mnt[iMarker]         =  ForceMomentum[0]*cos(Alpha)*cos(Beta) + ForceMomentum[1]*sin(Beta) + ForceMomentum[2]*sin(Alpha)*cos(Beta);
          CL_Mnt[iMarker]         = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[2]*cos(Alpha);
          CSF_Mnt[iMarker]        = -ForceMomentum[0]*sin(Beta)*cos(Alpha) + ForceMomentum[1]*cos(Beta) - ForceMomentum[2]*sin(Beta)*sin(Alpha);
          CEff_Mnt[iMarker]       = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
          CFx_Mnt[iMarker]        = ForceMomentum[0];
          CFy_Mnt[iMarker]        = ForceMomentum[1];
          CFz_Mnt[iMarker]        = ForceMomentum[2];
          CMx_Mnt[iMarker]        = MomentMomentum[0];
          CMy_Mnt[iMarker]        = MomentMomentum[1];
          CMz_Mnt[iMarker]        = MomentMomentum[2];
          CoPx_Mnt[iMarker]       = -MomentY_Force[0];
          CoPz_Mnt[iMarker]       =  MomentY_Force[2];
          CT_Mnt[iMarker]         = -CFz_Mnt[iMarker];
          CQ_Mnt[iMarker]         = -CMz_Mnt[iMarker];
          CMerit_Mnt[iMarker]     = CT_Mnt[iMarker] / (CQ_Mnt[iMarker] + EPS);
        }

        AllBound_CD_Mnt           += CD_Mnt[iMarker];
        AllBound_CL_Mnt           += CL_Mnt[iMarker];
        AllBound_CSF_Mnt          += CSF_Mnt[iMarker];
        AllBound_CEff_Mnt          = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
        AllBound_CFx_Mnt          += CFx_Mnt[iMarker];
        AllBound_CFy_Mnt          += CFy_Mnt[iMarker];
        AllBound_CFz_Mnt          += CFz_Mnt[iMarker];
        AllBound_CMx_Mnt          += CMx_Mnt[iMarker];
        AllBound_CMy_Mnt          += CMy_Mnt[iMarker];
        AllBound_CMz_Mnt          += CMz_Mnt[iMarker];
        AllBound_CoPy_Mnt         += CoPy_Mnt[iMarker];
        AllBound_CoPz_Mnt         += CoPz_Mnt[iMarker];
        AllBound_CT_Mnt           += CT_Mnt[iMarker];
        AllBound_CQ_Mnt           += CQ_Mnt[iMarker];
        AllBound_CMerit_Mnt       += AllBound_CT_Mnt / (AllBound_CQ_Mnt + EPS);

        /*--- Compute the coefficients per surface ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Mnt[iMarker_Monitoring]      += CL_Mnt[iMarker];
            Surface_CD_Mnt[iMarker_Monitoring]      += CD_Mnt[iMarker];
            Surface_CSF_Mnt[iMarker_Monitoring]     += CSF_Mnt[iMarker];
            Surface_CEff_Mnt[iMarker_Monitoring]     = CL_Mnt[iMarker] / (CD_Mnt[iMarker] + EPS);
            Surface_CFx_Mnt[iMarker_Monitoring]     += CFx_Mnt[iMarker];
            Surface_CFy_Mnt[iMarker_Monitoring]     += CFy_Mnt[iMarker];
            Surface_CFz_Mnt[iMarker_Monitoring]     += CFz_Mnt[iMarker];
            Surface_CMx_Mnt[iMarker_Monitoring]     += CMx_Mnt[iMarker];
            Surface_CMy_Mnt[iMarker_Monitoring]     += CMy_Mnt[iMarker];
            Surface_CMz_Mnt[iMarker_Monitoring]     += CMz_Mnt[iMarker];
          }
        }

      }


    }
  }

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  MyAllBound_CD_Mnt           = AllBound_CD_Mnt;           AllBound_CD_Mnt = 0.0;
  MyAllBound_CL_Mnt           = AllBound_CL_Mnt;           AllBound_CL_Mnt = 0.0;
  MyAllBound_CSF_Mnt          = AllBound_CSF_Mnt;          AllBound_CSF_Mnt = 0.0;
  MyAllBound_CFx_Mnt          = AllBound_CFx_Mnt;          AllBound_CFx_Mnt = 0.0;
  MyAllBound_CFy_Mnt          = AllBound_CFy_Mnt;          AllBound_CFy_Mnt = 0.0;
  MyAllBound_CFz_Mnt          = AllBound_CFz_Mnt;          AllBound_CFz_Mnt = 0.0;
  MyAllBound_CMx_Mnt          = AllBound_CMx_Mnt;          AllBound_CMx_Mnt = 0.0;
  MyAllBound_CMy_Mnt          = AllBound_CMy_Mnt;          AllBound_CMy_Mnt = 0.0;
  MyAllBound_CMz_Mnt          = AllBound_CMz_Mnt;          AllBound_CMz_Mnt = 0.0;
  MyAllBound_CoPx_Mnt         = AllBound_CoPx_Mnt;         AllBound_CoPx_Mnt = 0.0;
  MyAllBound_CoPy_Mnt         = AllBound_CoPy_Mnt;         AllBound_CoPy_Mnt = 0.0;
  MyAllBound_CoPz_Mnt         = AllBound_CoPz_Mnt;         AllBound_CoPz_Mnt = 0.0;
  MyAllBound_CT_Mnt           = AllBound_CT_Mnt;           AllBound_CT_Mnt = 0.0;
  MyAllBound_CQ_Mnt           = AllBound_CQ_Mnt;           AllBound_CQ_Mnt = 0.0;

  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(&MyAllBound_CD_Mnt, &AllBound_CD_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CL_Mnt, &AllBound_CL_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CSF_Mnt, &AllBound_CSF_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AllBound_CEff_Mnt = AllBound_CL_Mnt / (AllBound_CD_Mnt + EPS);
    SU2_MPI::Allreduce(&MyAllBound_CFx_Mnt, &AllBound_CFx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFy_Mnt, &AllBound_CFy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CFz_Mnt, &AllBound_CFz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMx_Mnt, &AllBound_CMx_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMy_Mnt, &AllBound_CMy_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyAllBound_CMz_Mnt, &AllBound_CMz_Mnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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

    Surface_CL_Mnt[iMarker_Monitoring]         = 0.0;
    Surface_CD_Mnt[iMarker_Monitoring]         = 0.0;
    Surface_CSF_Mnt[iMarker_Monitoring]        = 0.0;
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
  Total_CFx           += AllBound_CFx_Mnt;
  Total_CFy           += AllBound_CFy_Mnt;
  Total_CFz           += AllBound_CFz_Mnt;
  Total_CMx           += AllBound_CMx_Mnt;
  Total_CMy           += AllBound_CMy_Mnt;
  Total_CMz           += AllBound_CMz_Mnt;
  Total_CoPx          += AllBound_CoPx_Mnt;
  Total_CoPy          += AllBound_CoPy_Mnt;
  Total_CoPz          += AllBound_CoPz_Mnt;
  Total_CT            += AllBound_CT_Mnt;
  Total_CQ            += AllBound_CQ_Mnt;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]         += Surface_CL_Mnt[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]         += Surface_CD_Mnt[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring]        += Surface_CSF_Mnt[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]       += Surface_CL_Mnt[iMarker_Monitoring] / (Surface_CD_Mnt[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]        += Surface_CFx_Mnt[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]        += Surface_CFy_Mnt[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]        += Surface_CFz_Mnt[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]        += Surface_CMx_Mnt[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]        += Surface_CMy_Mnt[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]        += Surface_CMz_Mnt[iMarker_Monitoring];
  }

}

void CTNE2EulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
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
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        nodes->AddSolution(iPoint, iVar, -Res*Delta);
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

}

void CTNE2EulerSolver::ExplicitRK_Iteration(CGeometry *geometry,CSolver **solver_container, CConfig *config, unsigned short iRKStep) {

  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;

  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry-> node[iPoint]->GetVolume();
    Delta = nodes->GetDelta_Time(iPoint) / Vol;

    Res_TruncError = nodes->GetResTruncError(iPoint);
    Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        nodes->AddSolution(iPoint,iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry-> node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
    }

  }

  /*--- MPI solution ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);


}

void CTNE2EulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solution_container, CConfig *config) {
  unsigned short iVar;
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
    Vol = geometry-> node[iPoint]->GetVolume();

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    if (nodes->GetDelta_Time(iPoint) != 0.0) {
      Delta = Vol / nodes->GetDelta_Time(iPoint);
      Jacobian.AddVal2Diag(iPoint, Delta);
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
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry-> node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
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

  /*--- MPI solution ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}

void CTNE2EulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iSpecies, iVar, iMarker, RHOS_INDEX, RHO_INDEX;
  su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
      Partial_Gradient, Partial_Res, *Normal;
  su2double rho_i,rho_j;

  /*--- Initialize arrays ---*/

  /*--- Primitive variables: [Y1, ..., YNs, T, Tve, u, v, w]^T ---*/
  PrimVar_Vertex = new su2double [nPrimVarGrad];
  PrimVar_i = new su2double [nPrimVarGrad];
  PrimVar_j = new su2double [nPrimVarGrad];

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Set Gradient_Primitive to zero ---*/
  nodes->SetGradient_PrimitiveZero();

  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Pull primitives from CVariable ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
      PrimVar_j[iVar] = nodes->GetPrimitive(jPoint,iVar);
    }

    Normal = geometry-> edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          nodes->AddGradient_Primitive(iPoint, iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          nodes->SubtractGradient_Primitive(jPoint, iVar, iDim, Partial_Res);
      }
    }
  }

  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {

        /*--- Get primitives from CVariable ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          PrimVar_Vertex[iVar] = nodes->GetPrimitive(iPoint, iVar);

        // /*--- Modify species density to mass concentration ---*/
        // rho_i = nodes->GetPrimVar(RHO_INDEX);
        // for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        //   PrimVar_Vertex[RHOS_INDEX+iSpecies] = PrimVar_Vertex[RHOS_INDEX+iSpecies]/rho_i;

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
            nodes->SubtractGradient_Primitive(iPoint, iVar, iDim, Partial_Res);
          }
      }
    }
  }

  /*--- Update gradient value ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = nodes->GetGradient_Primitive(iPoint, iVar,iDim) / (geometry->node[iPoint]->GetVolume());
        nodes->SetGradient_Primitive(iPoint, iVar, iDim, Partial_Gradient);

      }
    }
  }

  delete [] PrimVar_Vertex;
  delete [] PrimVar_i;
  delete [] PrimVar_j;

  InitiateComms(geometry, config, PRIMITIVE_GRADIENT);
  CompleteComms(geometry, config, PRIMITIVE_GRADIENT);

}

void CTNE2EulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {

  unsigned short iVar, iDim, jDim, iNeigh, RHOS_INDEX, RHO_INDEX;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
      r23_b, r33, rho_i, rho_j, weight, product, detR2, z11, z12, z13, z22, z23, z33;
  bool singular;

  //Unused at the momemnt
  //unsigned short iSpecies;
  //unsigned long

  /*--- Initialize arrays, Primitive variables:
   [Y1, ..., YNs, T, Tve, u, v, w, P]^T ---*/
  PrimVar_i      = new su2double [nPrimVarGrad];
  PrimVar_j      = new su2double [nPrimVarGrad];

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Loop over points of the grid ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Set the value of singulare ---*/
    singular = false;

    /*--- Get coordinates ---*/
    Coord_i = geometry->node[iPoint]->GetCoord();

    /*--- Get primitives from CVariable ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      PrimVar_i[iVar] = nodes->GetPrimVar(iVar);

    /*--- Inizialization of variables ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;

    r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();

      for (iVar = 0; iVar < nPrimVarGrad; iVar++)
        PrimVar_j[iVar] = nodes->GetPrimitive(jPoint, iVar);

      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

      /*--- Sumations for entries of upper triangular matrix R ---*/
      if (weight != 0.0){
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
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim]) *
                (PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
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
        for (jDim = 0; jDim < nDim; jDim++)
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        nodes->SetGradient_Primitive(iPoint, iVar, iDim, product);
      }
    }
  }

  delete [] PrimVar_i;
  delete [] PrimVar_j;

  InitiateComms(geometry, config, PRIMITIVE_GRADIENT);
  CompleteComms(geometry, config, PRIMITIVE_GRADIENT);
}

void CTNE2EulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry,
                                                CConfig *config,
                                                unsigned long val_Point) {

  unsigned short iSpecies, iVar, iDim, jDim, iNeigh, RHOS_INDEX, RHO_INDEX;
  unsigned long iPoint, jPoint;
  su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
      r23_b, r33, rho_i, rho_j, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular=false;

  /*--- Initialize arrays Primitive variables:
   [Y1, ..., YNs, T, Tve, u, v, w]^T ---*/
  PrimVar_i = new su2double [nPrimVarGrad];
  PrimVar_j = new su2double [nPrimVarGrad];

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  iPoint = val_Point;
  Coord_i = geometry->node[iPoint]->GetCoord();

  /*--- Get primitives from CVariable ---*/
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    PrimVar_i[iVar] = nodes->GetPrimVar(iVar);

  /*--- Inizialization of variables ---*/
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Cvector[iVar][iDim] = 0.0;

  r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
  r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

  for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
    jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
    Coord_j = geometry->node[jPoint]->GetCoord();

    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      PrimVar_j[iVar] = nodes->GetPrimitive(jPoint,iVar);

    weight = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

    /*--- Sumations for entries of upper triangular matrix R ---*/
    if (weight != 0.0) {
      r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
      r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
      r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
      if (nDim == 3) {
        r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
        r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
        r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
        r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
      }

      /*--- Entries of c:= transpose(A)*b ---*/
      for (iVar = 0; iVar < nPrimVarGrad; iVar++)
        for (iDim = 0; iDim < nDim; iDim++)
          Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/(weight);
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
      Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
      Smatrix[0][1] = -r11*r12/(detR2);
      Smatrix[1][0] = Smatrix[0][1];
      Smatrix[1][1] = r11*r11/(detR2);
    }
    else {
      z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
      z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
      Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
      Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
      Smatrix[0][2] = (z13*z33)/(detR2);
      Smatrix[1][0] = Smatrix[0][1];
      Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
      Smatrix[1][2] = (z23*z33)/(detR2);
      Smatrix[2][0] = Smatrix[0][2];
      Smatrix[2][1] = Smatrix[1][2];
      Smatrix[2][2] = (z33*z33)/(detR2);
    }
  }

  /*--- Computation of the gradient: S*c ---*/
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      product = 0.0;
      for (jDim = 0; jDim < nDim; jDim++)
        product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
      nodes->SetGradient_Primitive(iPoint, iVar, iDim, product);
    }
  }

  AD::SetPreaccOut(nodes->GetGradient_Primitive(iPoint), nPrimVarGrad, nDim);
  AD::EndPreacc();

  delete [] PrimVar_i;
  delete [] PrimVar_j;

  /*--- Communicate the gradient values via MPI. ---*/
  InitiateComms(geometry, config, PRIMITIVE_GRADIENT);
  CompleteComms(geometry, config, PRIMITIVE_GRADIENT);

}

void CTNE2EulerSolver::SetPrimitive_Gradient(CConfig *config) {
  unsigned long iPoint;
  unsigned short iVar;
  su2double *U, *V, **GradU, **GradV;

  /*--- Loop over all points ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    U = nodes->GetSolution(iPoint);
    V = nodes->GetPrimitive(iPoint);
    GradU = nodes->GetGradient(iPoint);
    GradV = nodes->GetGradient_Primitive(iPoint);

    nodes->GradCons2GradPrimVar(config, U, V, GradU, GradV);
  }
}

void CTNE2EulerSolver::SetGradient_L2Proj2(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iFlux, iSpecies;
  unsigned short nVarMetr = nVar, nFluxMetr = 2;  //--- TODO: adjust size of grad vector later for goal vs. feature
  su2double density, *densitys, velocity[2], pressure, enthalpy, rhoeve;
  su2double vnx[3], vny[3];
  su2double graTri[2], graTriVisc[2], graTriSource[2];
  su2double Crd[3][2], Sens[3][nVarMetr][nFluxMetr], SensVisc[3][nVarMetr][nFluxMetr], SensSource[3][nVarMetr];
  bool dummy_bool;
  bool viscous = config->GetViscous(), source = config->GetAdap_Source();

  densitys = new su2double[nSpecies];

  //--- note: currently only implemented for Tri

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<3; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (unsigned short iDim = 0; iDim<2; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }
      //--- store sensors (goal-oriented)
      density = nodes->GetDensity(kNode);
      for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        densitys[iSpecies] = nodes->GetSolution(kNode, iSpecies);
      }
      velocity[0] = nodes->GetVelocity(kNode, 0);
      velocity[1] = nodes->GetVelocity(kNode, 1);
      pressure    = nodes->GetPressure(kNode);
      enthalpy    = nodes->GetEnthalpy(kNode);
      rhoeve      = nodes->GetSolution(kNode, nVar-1);

      for(iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Sens[iNode][iSpecies][0] = densitys[iSpecies]*velocity[0];
        Sens[iNode][iSpecies][1] = densitys[iSpecies]*velocity[1];
      }
      
      Sens[iNode][nSpecies][0]   = density*velocity[0]*velocity[0]+pressure;
      Sens[iNode][nSpecies+1][0] = density*velocity[0]*velocity[1];
      Sens[iNode][nSpecies+2][0] = density*velocity[0]*enthalpy;
      Sens[iNode][nSpecies+3][0] = rhoeve*velocity[0];

      Sens[iNode][nSpecies][1]   = density*velocity[1]*velocity[0];
      Sens[iNode][nSpecies+1][1] = density*velocity[1]*velocity[1]+pressure;
      Sens[iNode][nSpecies+2][1] = density*velocity[1]*enthalpy;
      Sens[iNode][nSpecies+3][1] = rhoeve*velocity[1];

      if(viscous) {

        //--- add viscous flux contribution

        const su2double laminar_viscosity = nodes->GetLaminarViscosity(kNode);
        const su2double eddy_viscosity    = nodes->GetEddyViscosity(kNode);
        const su2double total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        su2double Grad_Vel[2][2], tau[2][2], delta[2][2] = {{1.0, 0.0},{0.0,1.0}};
        su2double Grad_Temp[3], Heat_Flux[3];
        
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(kNode, iDim+1, jDim);
          }
          Grad_Temp[iDim] = nodes->GetGradient_Primitive(kNode, 0, iDim);
        }
        
        /*--- Divergence of the velocity ---*/
        
        su2double div_vel = 0.0; for (unsigned short iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          for (unsigned short jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim]+Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
        }

        /*--- Compute the heat flux vector ---*/
        const su2double Gas_Constant = config->GetGas_ConstantND();
        const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
        const su2double prandtl_lam = config->GetPrandtl_Lam();
        const su2double prandtl_turb = config->GetPrandtl_Turb();
        const su2double thermal_conductivity = Cp * (laminar_viscosity/prandtl_lam + eddy_viscosity/prandtl_turb);
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Heat_Flux[iDim] = thermal_conductivity*Grad_Temp[iDim];
        }

        SensVisc[iNode][0][0] = 0.;
        SensVisc[iNode][1][0] = tau[0][0];
        SensVisc[iNode][2][0] = tau[0][1];
        SensVisc[iNode][3][0] = velocity[0]*tau[0][0] + velocity[1]*tau[0][1] + Heat_Flux[0];

        SensVisc[iNode][0][1] = 0.;
        SensVisc[iNode][1][1] = tau[1][0];
        SensVisc[iNode][2][1] = tau[1][1];
        SensVisc[iNode][3][1] = velocity[0]*tau[1][0] + velocity[1]*tau[1][1] + Heat_Flux[1];

      }

      if(source) {
        for(iVar = 0; iVar < nVar; iVar++) {
          SensSource[iNode][iVar] = nodes->GetSource(kNode)[iVar];
        }
      }
    }

    //--- inward edge's normals : edg[0]=P1P2, edg[1]=P2P0, edg[2]=P0P1
    vnx[0] = Crd[1][1]-Crd[2][1];
    vny[0] = Crd[2][0]-Crd[1][0];

    vnx[1] = Crd[2][1]-Crd[0][1];
    vny[1] = Crd[0][0]-Crd[2][0];

    vnx[2] = Crd[0][1]-Crd[1][1];
    vny[2] = Crd[1][0]-Crd[0][0];

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 3; ++iNode) {
      su2double CrdAvg[2] = {0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 3; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
        }
      }
      CrdAvg[0] /= 2.;
      CrdAvg[1] /= 2.;
      su2double u[2] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
      }
    }

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- gradient at the element ( graTri = 2*|T|*gradT ) 
        graTri[0] = Sens[0][iVar][iFlux]*vnx[0] + Sens[1][iVar][iFlux]*vnx[1] + Sens[2][iVar][iFlux]*vnx[2];
        graTri[1] = Sens[0][iVar][iFlux]*vny[0] + Sens[1][iVar][iFlux]*vny[1] + Sens[2][iVar][iFlux]*vny[2];

        graTriVisc[0] = SensVisc[0][iVar][iFlux]*vnx[0] + SensVisc[1][iVar][iFlux]*vnx[1] + SensVisc[2][iVar][iFlux]*vnx[2];
        graTriVisc[1] = SensVisc[0][iVar][iFlux]*vny[0] + SensVisc[1][iVar][iFlux]*vny[1] + SensVisc[2][iVar][iFlux]*vny[2];
    
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
        for (unsigned short iNode=0; iNode<3; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Area = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Area*6.);
          nodes->AddAnisoGrad(kNode, i+0, graTri[0] * rap);
          nodes->AddAnisoGrad(kNode, i+1, graTri[1] * rap);

          if(viscous){
            nodes->AddAnisoViscGrad(kNode, i+0, graTriVisc[0] * rap);
            nodes->AddAnisoViscGrad(kNode, i+1, graTriVisc[1] * rap);
          }
        }
      }

      if(source) {
        graTriSource[0] = SensSource[0][iVar]*vnx[0] + SensSource[1][iVar]*vnx[1] + SensSource[2][iVar]*vnx[2];
        graTriSource[1] = SensSource[0][iVar]*vny[0] + SensSource[1][iVar]*vny[1] + SensSource[2][iVar]*vny[2];
        //--- assembling
        const unsigned short i = iVar*nDim;
        for (unsigned short iNode=0; iNode<3; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Area = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Area*6.);
          nodes->AddAnisoSourceGrad(kNode, i+0, graTriSource[0] * rap);
          nodes->AddAnisoSourceGrad(kNode, i+1, graTriSource[1] * rap);
        }
      }
    }
  }

  /*--- Communicate the gradient values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_GRADIENT);
  CompleteComms(geometry, config, ANISO_GRADIENT);
  if(viscous) {
    InitiateComms(geometry, config, ANISO_GRADIENT_VISC);
    CompleteComms(geometry, config, ANISO_GRADIENT_VISC);
  }
  if(source) {
    InitiateComms(geometry, config, ANISO_GRADIENT_SOURCE);
    CompleteComms(geometry, config, ANISO_GRADIENT_SOURCE);
  }
}

void CTNE2EulerSolver::SetHessian_L2Proj2(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint(), nPointDomain = geometry->GetnPointDomain(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iFlux;
  unsigned short nVarMetr = nVar, nFluxMetr = 2;  //--- TODO: adjust size of grad vector later for goal vs. feature
  unsigned short nMetr = 3;
  su2double vnx[3], vny[3];
  su2double hesTri[3], hesTriVisc[3], hesTriSource[3];
  su2double Crd[3][2], Grad[3][2][nVarMetr][nFluxMetr], GradVisc[3][2][nVarMetr][nFluxMetr], GradSource[3][2][nVarMetr];
  bool viscous = config->GetViscous(), source = config->GetAdap_Source();

  su2double **A      = new su2double*[nDim],
            **EigVec = new su2double*[nDim], 
            *EigVal  = new su2double[nDim];

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    A[iDim]      = new su2double[nDim];
    EigVec[iDim] = new su2double[nDim];
  }
  
  //--- note: currently only implemented for Tri

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<3; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (unsigned short iDim = 0; iDim<2; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }
      //--- store gradient
      for(iVar = 0; iVar < nVarMetr; iVar++){
        for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
          const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
          Grad[iNode][0][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+0);
          Grad[iNode][1][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+1);

          if(viscous){
            GradVisc[iNode][0][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+0);
            GradVisc[iNode][1][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+1);
          }
        }

        if(source){
          const unsigned short i = iVar*nDim;
          GradSource[iNode][0][iVar] = nodes->GetAnisoSourceGrad(kNode, i+0);
          GradSource[iNode][1][iVar] = nodes->GetAnisoSourceGrad(kNode, i+1);
        }
      }
    }

    //--- inward edge's normals : edg[0]=P1P2, edg[1]=P2P0, edg[2]=P0P1
    vnx[0] = Crd[1][1]-Crd[2][1];
    vny[0] = Crd[2][0]-Crd[1][0];

    vnx[1] = Crd[2][1]-Crd[0][1];
    vny[1] = Crd[0][0]-Crd[2][0];

    vnx[2] = Crd[0][1]-Crd[1][1];
    vny[2] = Crd[1][0]-Crd[0][0];

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 3; ++iNode) {
      su2double CrdAvg[2] = {0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 3; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
        }
      }
      CrdAvg[0] /= 2.;
      CrdAvg[1] /= 2.;
      su2double u[2] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
      }
    }

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- hessian at the element ( hesTri = 2*|T|*hessienT ) 
        hesTri[0] =         Grad[0][0][iVar][iFlux]*vnx[0] 
                          + Grad[1][0][iVar][iFlux]*vnx[1] 
                          + Grad[2][0][iVar][iFlux]*vnx[2];
        hesTri[1] = 0.5 * ( Grad[0][0][iVar][iFlux]*vny[0] 
                          + Grad[1][0][iVar][iFlux]*vny[1] 
                          + Grad[2][0][iVar][iFlux]*vny[2]
                          + Grad[0][1][iVar][iFlux]*vnx[0] 
                          + Grad[1][1][iVar][iFlux]*vnx[1] 
                          + Grad[2][1][iVar][iFlux]*vnx[2] );
        hesTri[2] =         Grad[0][1][iVar][iFlux]*vny[0] 
                          + Grad[1][1][iVar][iFlux]*vny[1] 
                          + Grad[2][1][iVar][iFlux]*vny[2];

        if(viscous){
          hesTriVisc[0] =         GradVisc[0][0][iVar][iFlux]*vnx[0] 
                                + GradVisc[1][0][iVar][iFlux]*vnx[1] 
                                + GradVisc[2][0][iVar][iFlux]*vnx[2];
          hesTriVisc[1] = 0.5 * ( GradVisc[0][0][iVar][iFlux]*vny[0] 
                                + GradVisc[1][0][iVar][iFlux]*vny[1] 
                                + GradVisc[2][0][iVar][iFlux]*vny[2]
                                + GradVisc[0][1][iVar][iFlux]*vnx[0] 
                                + GradVisc[1][1][iVar][iFlux]*vnx[1] 
                                + GradVisc[2][1][iVar][iFlux]*vnx[2] );
          hesTriVisc[2] =         GradVisc[0][1][iVar][iFlux]*vny[0] 
                                + GradVisc[1][1][iVar][iFlux]*vny[1] 
                                + GradVisc[2][1][iVar][iFlux]*vny[2];
        }
        
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;
        for (unsigned short iNode=0; iNode<3; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Area = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Area*6.);
          nodes->AddAnisoHess(kNode, i+0, hesTri[0] * rap);
          nodes->AddAnisoHess(kNode, i+1, hesTri[1] * rap);
          nodes->AddAnisoHess(kNode, i+2, hesTri[2] * rap);

          if(viscous) {
            nodes->AddAnisoViscHess(kNode, i+0, hesTriVisc[0] * rap);
            nodes->AddAnisoViscHess(kNode, i+1, hesTriVisc[1] * rap);
            nodes->AddAnisoViscHess(kNode, i+2, hesTriVisc[2] * rap);
          }
        }
      }

      if(source){
        hesTriSource[0] =         GradSource[0][0][iVar]*vnx[0] 
                                + GradSource[1][0][iVar]*vnx[1] 
                                + GradSource[2][0][iVar]*vnx[2];
        hesTriSource[1] = 0.5 * ( GradSource[0][0][iVar]*vny[0] 
                                + GradSource[1][0][iVar]*vny[1] 
                                + GradSource[2][0][iVar]*vny[2]
                                + GradSource[0][1][iVar]*vnx[0] 
                                + GradSource[1][1][iVar]*vnx[1] 
                                + GradSource[2][1][iVar]*vnx[2] );
        hesTriSource[2] =         GradSource[0][1][iVar]*vny[0] 
                                + GradSource[1][1][iVar]*vny[1] 
                                + GradSource[2][1][iVar]*vny[2];
      
        //--- assembling
        const unsigned short i = iVar*nMetr;
        for (unsigned short iNode=0; iNode<3; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Area = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Area*6.);
          nodes->AddAnisoSourceHess(kNode, i+0, hesTriSource[0] * rap);
          nodes->AddAnisoSourceHess(kNode, i+1, hesTriSource[1] * rap);
          nodes->AddAnisoSourceHess(kNode, i+2, hesTriSource[2] * rap);
        }
      }
    }
  }

  /*--- Communicate the Hessian values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_HESSIAN);
  CompleteComms(geometry, config, ANISO_HESSIAN);
  if(viscous) {
    InitiateComms(geometry, config, ANISO_HESSIAN_VISC);
    CompleteComms(geometry, config, ANISO_HESSIAN_VISC);
  }
  if(source) {
    InitiateComms(geometry, config, ANISO_HESSIAN_SOURCE);
    CompleteComms(geometry, config, ANISO_HESSIAN_SOURCE);
  }

  CorrectBoundAnisoHess(geometry, config);

  //--- Make positive definite matrix
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    for(iVar = 0; iVar < nVarMetr; iVar++){
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

        const su2double a = nodes->GetAnisoHess(iPoint, i+0);
        const su2double b = nodes->GetAnisoHess(iPoint, i+1);
        const su2double c = nodes->GetAnisoHess(iPoint, i+2);
        
        A[0][0] = a; A[0][1] = b;
        A[1][0] = b; A[1][1] = c;

        CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

        for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

        CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

        nodes->SetAnisoHess(iPoint, i+0, A[0][0]);
        nodes->SetAnisoHess(iPoint, i+1, A[0][1]);
        nodes->SetAnisoHess(iPoint, i+2, A[1][1]);
      }
    }
  }

  if(viscous) {
    for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      //--- viscous mass flux is 0, so start with momentum
      for(iVar = 1; iVar < nVarMetr; iVar++){
        for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
          const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

          const su2double a = nodes->GetAnisoViscHess(iPoint, i+0);
          const su2double b = nodes->GetAnisoViscHess(iPoint, i+1);
          const su2double c = nodes->GetAnisoViscHess(iPoint, i+2);

          A[0][0] = a; A[0][1] = b;
          A[1][0] = b; A[1][1] = c;

          CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

          for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

          CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

          nodes->SetAnisoViscHess(iPoint, i+0, A[0][0]);
          nodes->SetAnisoViscHess(iPoint, i+1, A[0][1]);
          nodes->SetAnisoViscHess(iPoint, i+2, A[1][1]);
        }
      }
    }
  }

  if(source) {
    for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      //--- viscous mass flux is 0, so start with momentum
      for(iVar = 1; iVar < nVarMetr; iVar++){
        const unsigned short i = iVar*nMetr;

        const su2double a = nodes->GetAnisoSourceHess(iPoint, i+0);
        const su2double b = nodes->GetAnisoSourceHess(iPoint, i+1);
        const su2double c = nodes->GetAnisoSourceHess(iPoint, i+2);

        A[0][0] = a; A[0][1] = b;
        A[1][0] = b; A[1][1] = c;

        CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

        for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

        CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

        nodes->SetAnisoSourceHess(iPoint, i+0, A[0][0]);
        nodes->SetAnisoSourceHess(iPoint, i+1, A[0][1]);
        nodes->SetAnisoSourceHess(iPoint, i+2, A[1][1]);
      }
    }
  }
}

void CTNE2EulerSolver::SetGradient_L2Proj3(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iDim, iFlux;
  unsigned short nVarMetr = nVar, nFluxMetr = 3;
  su2double density, velocity[3], pressure, enthalpy;
  su2double vnx[4], vny[4], vnz[4];
  su2double graTet[3], graTetVisc[3];
  su2double Crd[4][3], Sens[4][nVarMetr][nFluxMetr], SensVisc[4][nVarMetr][nFluxMetr];
  bool dummy_bool;
  bool viscous = config->GetViscous(), source = config->GetAdap_Source();

  //--- note: currently only implemented for Tet

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<4; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (iDim = 0; iDim<3; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }
      //--- store sensors (goal-oriented)
      density     = nodes->GetDensity(kNode);
      velocity[0] = nodes->GetVelocity(kNode, 0);
      velocity[1] = nodes->GetVelocity(kNode, 1);
      velocity[2] = nodes->GetVelocity(kNode, 2);
      pressure    = nodes->GetPressure(kNode);
      enthalpy    = nodes->GetEnthalpy(kNode);

      Sens[iNode][0][0] = density*velocity[0];
      Sens[iNode][1][0] = Sens[iNode][0][0]*velocity[0]+pressure;
      Sens[iNode][2][0] = Sens[iNode][0][0]*velocity[1];
      Sens[iNode][3][0] = Sens[iNode][0][0]*velocity[2];
      Sens[iNode][4][0] = Sens[iNode][0][0]*enthalpy;

      Sens[iNode][0][1] = density*velocity[1];
      Sens[iNode][1][1] = Sens[iNode][0][1]*velocity[0];
      Sens[iNode][2][1] = Sens[iNode][0][1]*velocity[1]+pressure;
      Sens[iNode][3][1] = Sens[iNode][0][1]*velocity[2];
      Sens[iNode][4][1] = Sens[iNode][0][1]*enthalpy;

      Sens[iNode][0][2] = density*velocity[2];
      Sens[iNode][1][2] = Sens[iNode][0][2]*velocity[0];
      Sens[iNode][2][2] = Sens[iNode][0][2]*velocity[1];
      Sens[iNode][3][2] = Sens[iNode][0][2]*velocity[2]+pressure;
      Sens[iNode][4][2] = Sens[iNode][0][2]*enthalpy;

      if(viscous) {

        //--- add viscous flux contribution

        const su2double laminar_viscosity = nodes->GetLaminarViscosity(kNode);
        const su2double eddy_viscosity    = nodes->GetEddyViscosity(kNode);
        const su2double total_viscosity   = laminar_viscosity + eddy_viscosity;
        
        su2double Grad_Vel[3][3], tau[3][3], delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
        su2double Grad_Temp[3], Heat_Flux[3];
        
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
            Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(kNode, iDim+1, jDim);
          }
          Grad_Temp[iDim] = nodes->GetGradient_Primitive(kNode, 0, iDim);
        }
        
        /*--- Divergence of the velocity ---*/
        
        su2double div_vel = 0.0; for (unsigned short iDim = 0 ; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          for (unsigned short jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim]+Grad_Vel[iDim][jDim] ) - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
          }
        }

        /*--- Compute the heat flux vector ---*/
        const su2double Gas_Constant = config->GetGas_ConstantND();
        const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
        const su2double prandtl_lam = config->GetPrandtl_Lam();
        const su2double prandtl_turb = config->GetPrandtl_Turb();
        const su2double thermal_conductivity = Cp * (laminar_viscosity/prandtl_lam + eddy_viscosity/prandtl_turb);
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Heat_Flux[iDim] = thermal_conductivity*Grad_Temp[iDim];
        }

        SensVisc[iNode][0][0] = 0.;
        SensVisc[iNode][1][0] = tau[0][0];
        SensVisc[iNode][2][0] = tau[0][1];
        SensVisc[iNode][3][0] = tau[0][2];
        SensVisc[iNode][4][0] = velocity[0]*tau[0][0] + velocity[1]*tau[0][1] + velocity[2]*tau[0][2] + Heat_Flux[0];

        SensVisc[iNode][0][1] = 0.;
        SensVisc[iNode][1][1] = tau[1][0];
        SensVisc[iNode][2][1] = tau[1][1];
        SensVisc[iNode][3][1] = tau[1][2];
        SensVisc[iNode][4][1] = velocity[0]*tau[1][0] + velocity[1]*tau[1][1] + velocity[2]*tau[1][2] + Heat_Flux[1];

        SensVisc[iNode][0][2] = 0.;
        SensVisc[iNode][1][2] = tau[2][0];
        SensVisc[iNode][2][2] = tau[2][1];
        SensVisc[iNode][3][2] = tau[2][2];
        SensVisc[iNode][4][2] = velocity[0]*tau[2][0] + velocity[1]*tau[2][1] + velocity[2]*tau[2][2] + Heat_Flux[2];
      }
    }

    //--- inward face's normals : fac[0]=P1P2P3, fac[1]=P2P3P0, fac[2]=P3P0P1, fac[3]=P0P1P2
    vnx[0] = (Crd[2][1]-Crd[1][1])*(Crd[3][2]-Crd[1][2]) - (Crd[2][2]-Crd[1][2])*(Crd[3][1]-Crd[1][1]);
    vny[0] = (Crd[2][2]-Crd[1][2])*(Crd[3][0]-Crd[1][0]) - (Crd[2][0]-Crd[1][0])*(Crd[3][2]-Crd[1][2]);
    vnz[0] = (Crd[2][0]-Crd[1][0])*(Crd[3][1]-Crd[1][1]) - (Crd[2][1]-Crd[1][1])*(Crd[3][0]-Crd[1][0]);

    vnx[1] = (Crd[3][1]-Crd[2][1])*(Crd[0][2]-Crd[2][2]) - (Crd[3][2]-Crd[2][2])*(Crd[0][1]-Crd[2][1]);
    vny[1] = (Crd[3][2]-Crd[2][2])*(Crd[0][0]-Crd[2][0]) - (Crd[3][0]-Crd[2][0])*(Crd[0][2]-Crd[2][2]);
    vnz[1] = (Crd[3][0]-Crd[2][0])*(Crd[0][1]-Crd[2][1]) - (Crd[3][1]-Crd[2][1])*(Crd[0][0]-Crd[2][0]);

    vnx[2] = (Crd[0][1]-Crd[3][1])*(Crd[1][2]-Crd[3][2]) - (Crd[0][2]-Crd[3][2])*(Crd[1][1]-Crd[3][1]);
    vny[2] = (Crd[0][2]-Crd[3][2])*(Crd[1][0]-Crd[3][0]) - (Crd[0][0]-Crd[3][0])*(Crd[1][2]-Crd[3][2]);
    vnz[2] = (Crd[0][0]-Crd[3][0])*(Crd[1][1]-Crd[3][1]) - (Crd[0][1]-Crd[3][1])*(Crd[1][0]-Crd[3][0]);

    vnx[3] = (Crd[1][1]-Crd[0][1])*(Crd[2][2]-Crd[0][2]) - (Crd[1][2]-Crd[0][2])*(Crd[2][1]-Crd[0][1]);
    vny[3] = (Crd[1][2]-Crd[0][2])*(Crd[2][0]-Crd[0][0]) - (Crd[1][0]-Crd[0][0])*(Crd[2][2]-Crd[0][2]);
    vnz[3] = (Crd[1][0]-Crd[0][0])*(Crd[2][1]-Crd[0][1]) - (Crd[1][1]-Crd[0][1])*(Crd[2][0]-Crd[0][0]);

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 4; ++iNode) {
      su2double CrdAvg[3] = {0.0, 0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 4; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
          CrdAvg[2] += Crd[jNode][2];
        }
      }
      CrdAvg[0] /= 3.;
      CrdAvg[1] /= 3.;
      CrdAvg[2] /= 3.;
      su2double u[3] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1],
                        CrdAvg[2]-Crd[iNode][2]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1] + vnz[iNode]*u[2]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
        vnz[iNode] *= -1.0;
      }
    }

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- gradient at the element ( graTet = 6*|T|*gradT ) 
        graTet[0] = Sens[0][iVar][iFlux]*vnx[0] + Sens[1][iVar][iFlux]*vnx[1] + Sens[2][iVar][iFlux]*vnx[2] + Sens[3][iVar][iFlux]*vnx[3];
        graTet[1] = Sens[0][iVar][iFlux]*vny[0] + Sens[1][iVar][iFlux]*vny[1] + Sens[2][iVar][iFlux]*vny[2] + Sens[3][iVar][iFlux]*vny[3];
        graTet[2] = Sens[0][iVar][iFlux]*vnz[0] + Sens[1][iVar][iFlux]*vnz[1] + Sens[2][iVar][iFlux]*vnz[2] + Sens[3][iVar][iFlux]*vnz[3];

        if(viscous) {
          graTetVisc[0] = SensVisc[0][iVar][iFlux]*vnx[0] + SensVisc[1][iVar][iFlux]*vnx[1] + SensVisc[2][iVar][iFlux]*vnx[2] + SensVisc[3][iVar][iFlux]*vnx[3];
          graTetVisc[1] = SensVisc[0][iVar][iFlux]*vny[0] + SensVisc[1][iVar][iFlux]*vny[1] + SensVisc[2][iVar][iFlux]*vny[2] + SensVisc[3][iVar][iFlux]*vny[3];
          graTetVisc[2] = SensVisc[0][iVar][iFlux]*vnz[0] + SensVisc[1][iVar][iFlux]*vnz[1] + SensVisc[2][iVar][iFlux]*vnz[2] + SensVisc[3][iVar][iFlux]*vnz[3];
        }
    
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
        for (unsigned short iNode=0; iNode<4; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Vol = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Vol*24.);
          nodes->AddAnisoGrad(kNode, i+0, graTet[0] * rap);
          nodes->AddAnisoGrad(kNode, i+1, graTet[1] * rap);
          nodes->AddAnisoGrad(kNode, i+2, graTet[2] * rap);

          if(viscous) {
            nodes->AddAnisoViscGrad(kNode, i+0, graTetVisc[0] * rap);
            nodes->AddAnisoViscGrad(kNode, i+1, graTetVisc[1] * rap);
            nodes->AddAnisoViscGrad(kNode, i+2, graTetVisc[2] * rap);
          }
        }
      }
    }
  }

  /*--- Communicate the gradient values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_GRADIENT);
  CompleteComms(geometry, config, ANISO_GRADIENT);
  if(viscous) {
    InitiateComms(geometry, config, ANISO_GRADIENT_VISC);
    CompleteComms(geometry, config, ANISO_GRADIENT_VISC);
  }
  if(source) {
    InitiateComms(geometry, config, ANISO_GRADIENT_SOURCE);
    CompleteComms(geometry, config, ANISO_GRADIENT_SOURCE);
  }
}

void CTNE2EulerSolver::SetHessian_L2Proj3(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint(), nPointDomain = geometry->GetnPointDomain(), iElem, nElem = geometry->GetnElem();
  unsigned short iVar, iFlux;
  unsigned short nVarMetr = nVar, nFluxMetr = 3;  //--- TODO: adjust size of grad vector later for goal vs. feature
  unsigned short nMetr = 6;
  su2double vnx[4], vny[4], vnz[4];
  su2double hesTet[6], hesTetVisc[6];
  su2double Crd[4][3], Grad[4][3][nVarMetr][nFluxMetr], GradVisc[4][3][nVarMetr][nFluxMetr];
  bool viscous = config->GetViscous(), source = config->GetAdap_Source();

  //--- note: currently only implemented for Tri

  for (iElem=0; iElem<nElem; ++iElem) {
    for (unsigned short iNode=0; iNode<4; ++iNode) {
      const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
      //--- store coordinates
      for (unsigned short iDim = 0; iDim<3; ++iDim) {
        Crd[iNode][iDim] = geometry->node[kNode]->GetCoord(iDim);
      }
      //--- store gradient
      for(iVar = 0; iVar < nVarMetr; iVar++){
        for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
          const unsigned short i = iFlux*nVarMetr*nDim + iVar*nDim;
          Grad[iNode][0][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+0);
          Grad[iNode][1][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+1);
          Grad[iNode][2][iVar][iFlux] = nodes->GetAnisoGrad(kNode, i+2);

          if(viscous) {
            GradVisc[iNode][0][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+0);
            GradVisc[iNode][1][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+1);
            GradVisc[iNode][2][iVar][iFlux] = nodes->GetAnisoViscGrad(kNode, i+2);
          }
        }
      }
    }

    //--- inward face's normals : fac[0]=P1P2P3, fac[1]=P2P3P0, fac[2]=P3P0P1, fac[3]=P0P1P2
    vnx[0] = (Crd[2][1]-Crd[1][1])*(Crd[3][2]-Crd[1][2]) - (Crd[2][2]-Crd[1][2])*(Crd[3][1]-Crd[1][1]);
    vny[0] = (Crd[2][2]-Crd[1][2])*(Crd[3][0]-Crd[1][0]) - (Crd[2][0]-Crd[1][0])*(Crd[3][2]-Crd[1][2]);
    vnz[0] = (Crd[2][0]-Crd[1][0])*(Crd[3][1]-Crd[1][1]) - (Crd[2][1]-Crd[1][1])*(Crd[3][0]-Crd[1][0]);

    vnx[1] = (Crd[3][1]-Crd[2][1])*(Crd[0][2]-Crd[2][2]) - (Crd[3][2]-Crd[2][2])*(Crd[0][1]-Crd[2][1]);
    vny[1] = (Crd[3][2]-Crd[2][2])*(Crd[0][0]-Crd[2][0]) - (Crd[3][0]-Crd[2][0])*(Crd[0][2]-Crd[2][2]);
    vnz[1] = (Crd[3][0]-Crd[2][0])*(Crd[0][1]-Crd[2][1]) - (Crd[3][1]-Crd[2][1])*(Crd[0][0]-Crd[2][0]);

    vnx[2] = (Crd[0][1]-Crd[3][1])*(Crd[1][2]-Crd[3][2]) - (Crd[0][2]-Crd[3][2])*(Crd[1][1]-Crd[3][1]);
    vny[2] = (Crd[0][2]-Crd[3][2])*(Crd[1][0]-Crd[3][0]) - (Crd[0][0]-Crd[3][0])*(Crd[1][2]-Crd[3][2]);
    vnz[2] = (Crd[0][0]-Crd[3][0])*(Crd[1][1]-Crd[3][1]) - (Crd[0][1]-Crd[3][1])*(Crd[1][0]-Crd[3][0]);

    vnx[3] = (Crd[1][1]-Crd[0][1])*(Crd[2][2]-Crd[0][2]) - (Crd[1][2]-Crd[0][2])*(Crd[2][1]-Crd[0][1]);
    vny[3] = (Crd[1][2]-Crd[0][2])*(Crd[2][0]-Crd[0][0]) - (Crd[1][0]-Crd[0][0])*(Crd[2][2]-Crd[0][2]);
    vnz[3] = (Crd[1][0]-Crd[0][0])*(Crd[2][1]-Crd[0][1]) - (Crd[1][1]-Crd[0][1])*(Crd[2][0]-Crd[0][0]);

    //--- check if inward normal
    for(unsigned short iNode = 0; iNode < 4; ++iNode) {
      su2double CrdAvg[3] = {0.0, 0.0, 0.0};
      for(unsigned short jNode = 0; jNode < 4; ++jNode) {
        if(iNode != jNode) {
          CrdAvg[0] += Crd[jNode][0];
          CrdAvg[1] += Crd[jNode][1];
          CrdAvg[2] += Crd[jNode][2];
        }
      }
      CrdAvg[0] /= 3.;
      CrdAvg[1] /= 3.;
      CrdAvg[2] /= 3.;
      su2double u[3] = {CrdAvg[0]-Crd[iNode][0],
                        CrdAvg[1]-Crd[iNode][1],
                        CrdAvg[2]-Crd[iNode][2]};
      if((vnx[iNode]*u[0] + vny[iNode]*u[1] + vnz[iNode]*u[2]) > 0.) {
        vnx[iNode] *= -1.0;
        vny[iNode] *= -1.0;
        vnz[iNode] *= -1.0;
      }
    }

    //--- loop over conservative variables
    for(iVar = 0; iVar < nVarMetr; iVar++){

      //--- loop over directions
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){

        //--- hessian at the element ( hesTet = 6*|T|*hessienT ) 
        hesTet[0] =         Grad[0][0][iVar][iFlux]*vnx[0] 
                          + Grad[1][0][iVar][iFlux]*vnx[1] 
                          + Grad[2][0][iVar][iFlux]*vnx[2]
                          + Grad[3][0][iVar][iFlux]*vnx[3];

        hesTet[1] = 0.5 * ( Grad[0][0][iVar][iFlux]*vny[0] 
                          + Grad[1][0][iVar][iFlux]*vny[1] 
                          + Grad[2][0][iVar][iFlux]*vny[2]
                          + Grad[3][0][iVar][iFlux]*vny[3]
                          + Grad[0][1][iVar][iFlux]*vnx[0] 
                          + Grad[1][1][iVar][iFlux]*vnx[1] 
                          + Grad[2][1][iVar][iFlux]*vnx[2]
                          + Grad[3][1][iVar][iFlux]*vnx[3] );

        hesTet[2] = 0.5 * ( Grad[0][0][iVar][iFlux]*vnz[0] 
                          + Grad[1][0][iVar][iFlux]*vnz[1] 
                          + Grad[2][0][iVar][iFlux]*vnz[2]
                          + Grad[3][0][iVar][iFlux]*vnz[3]
                          + Grad[0][2][iVar][iFlux]*vnx[0] 
                          + Grad[1][2][iVar][iFlux]*vnx[1] 
                          + Grad[2][2][iVar][iFlux]*vnx[2]
                          + Grad[3][2][iVar][iFlux]*vnx[3] );

        hesTet[3] =         Grad[0][1][iVar][iFlux]*vny[0] 
                          + Grad[1][1][iVar][iFlux]*vny[1] 
                          + Grad[2][1][iVar][iFlux]*vny[2]
                          + Grad[3][1][iVar][iFlux]*vny[3];

        hesTet[4] = 0.5 * ( Grad[0][1][iVar][iFlux]*vnz[0] 
                          + Grad[1][1][iVar][iFlux]*vnz[1] 
                          + Grad[2][1][iVar][iFlux]*vnz[2]
                          + Grad[3][1][iVar][iFlux]*vnz[3]
                          + Grad[0][2][iVar][iFlux]*vny[0] 
                          + Grad[1][2][iVar][iFlux]*vny[1] 
                          + Grad[2][2][iVar][iFlux]*vny[2]
                          + Grad[3][2][iVar][iFlux]*vny[3] );

        hesTet[5] =         Grad[0][2][iVar][iFlux]*vnz[0] 
                          + Grad[1][2][iVar][iFlux]*vnz[1] 
                          + Grad[2][2][iVar][iFlux]*vnz[2]
                          + Grad[3][2][iVar][iFlux]*vnz[3];

        if(viscous) {
          hesTetVisc[0] =         GradVisc[0][0][iVar][iFlux]*vnx[0] 
                                + GradVisc[1][0][iVar][iFlux]*vnx[1] 
                                + GradVisc[2][0][iVar][iFlux]*vnx[2]
                                + GradVisc[3][0][iVar][iFlux]*vnx[3];

          hesTetVisc[1] = 0.5 * ( GradVisc[0][0][iVar][iFlux]*vny[0] 
                                + GradVisc[1][0][iVar][iFlux]*vny[1] 
                                + GradVisc[2][0][iVar][iFlux]*vny[2]
                                + GradVisc[3][0][iVar][iFlux]*vny[3]
                                + GradVisc[0][1][iVar][iFlux]*vnx[0] 
                                + GradVisc[1][1][iVar][iFlux]*vnx[1] 
                                + GradVisc[2][1][iVar][iFlux]*vnx[2]
                                + GradVisc[3][1][iVar][iFlux]*vnx[3] );

          hesTetVisc[2] = 0.5 * ( GradVisc[0][0][iVar][iFlux]*vnz[0] 
                                + GradVisc[1][0][iVar][iFlux]*vnz[1] 
                                + GradVisc[2][0][iVar][iFlux]*vnz[2]
                                + GradVisc[3][0][iVar][iFlux]*vnz[3]
                                + GradVisc[0][2][iVar][iFlux]*vnx[0] 
                                + GradVisc[1][2][iVar][iFlux]*vnx[1] 
                                + GradVisc[2][2][iVar][iFlux]*vnx[2]
                                + GradVisc[3][2][iVar][iFlux]*vnx[3] );

          hesTetVisc[3] =         GradVisc[0][1][iVar][iFlux]*vny[0] 
                                + GradVisc[1][1][iVar][iFlux]*vny[1] 
                                + GradVisc[2][1][iVar][iFlux]*vny[2]
                                + GradVisc[3][1][iVar][iFlux]*vny[3];

          hesTetVisc[4] = 0.5 * ( GradVisc[0][1][iVar][iFlux]*vnz[0] 
                                + GradVisc[1][1][iVar][iFlux]*vnz[1] 
                                + GradVisc[2][1][iVar][iFlux]*vnz[2]
                                + GradVisc[3][1][iVar][iFlux]*vnz[3]
                                + GradVisc[0][2][iVar][iFlux]*vny[0] 
                                + GradVisc[1][2][iVar][iFlux]*vny[1] 
                                + GradVisc[2][2][iVar][iFlux]*vny[2]
                                + GradVisc[3][2][iVar][iFlux]*vny[3] );

          hesTetVisc[5] =         GradVisc[0][2][iVar][iFlux]*vnz[0] 
                                + GradVisc[1][2][iVar][iFlux]*vnz[1] 
                                + GradVisc[2][2][iVar][iFlux]*vnz[2]
                                + GradVisc[3][2][iVar][iFlux]*vnz[3];
        }
        
        //--- assembling
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;
        for (unsigned short iNode=0; iNode<4; ++iNode) {
          const unsigned long kNode = geometry->elem[iElem]->GetNode(iNode);
          const su2double Vol = geometry->node[kNode]->GetVolume();
          const su2double rap = 1./(Vol*24.);
          nodes->AddAnisoHess(kNode, i+0, hesTet[0] * rap);
          nodes->AddAnisoHess(kNode, i+1, hesTet[1] * rap);
          nodes->AddAnisoHess(kNode, i+2, hesTet[2] * rap);
          nodes->AddAnisoHess(kNode, i+3, hesTet[3] * rap);
          nodes->AddAnisoHess(kNode, i+4, hesTet[4] * rap);
          nodes->AddAnisoHess(kNode, i+5, hesTet[5] * rap);

          if(viscous) {
            nodes->AddAnisoViscHess(kNode, i+0, hesTetVisc[0] * rap);
            nodes->AddAnisoViscHess(kNode, i+1, hesTetVisc[1] * rap);
            nodes->AddAnisoViscHess(kNode, i+2, hesTetVisc[2] * rap);
            nodes->AddAnisoViscHess(kNode, i+3, hesTetVisc[3] * rap);
            nodes->AddAnisoViscHess(kNode, i+4, hesTetVisc[4] * rap);
            nodes->AddAnisoViscHess(kNode, i+5, hesTetVisc[5] * rap);
          }
        }
      }
    }
  }

  /*--- Communicate the Hessian values via MPI. ---*/
  
  InitiateComms(geometry, config, ANISO_HESSIAN);
  CompleteComms(geometry, config, ANISO_HESSIAN);
  if(viscous) {
    InitiateComms(geometry, config, ANISO_HESSIAN_VISC);
    CompleteComms(geometry, config, ANISO_HESSIAN_VISC);
  }
  if(source) {
    InitiateComms(geometry, config, ANISO_HESSIAN_SOURCE);
    CompleteComms(geometry, config, ANISO_HESSIAN_SOURCE);
  }

  CorrectBoundAnisoHess(geometry, config);

  //--- Make positive definite matrix
  su2double **A      = new su2double*[nDim],
            **EigVec = new su2double*[nDim], 
            *EigVal  = new su2double[nDim];

  for(unsigned short iDim = 0; iDim < nDim; ++iDim){
    A[iDim]      = new su2double[nDim];
    EigVec[iDim] = new su2double[nDim];
  }

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    for(iVar = 0; iVar < nVarMetr; iVar++){
      for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
        const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

        const su2double a = nodes->GetAnisoHess(iPoint, i+0);
        const su2double b = nodes->GetAnisoHess(iPoint, i+1);
        const su2double c = nodes->GetAnisoHess(iPoint, i+2);
        const su2double d = nodes->GetAnisoHess(iPoint, i+3);
        const su2double e = nodes->GetAnisoHess(iPoint, i+4);
        const su2double f = nodes->GetAnisoHess(iPoint, i+5);

        A[0][0] = a; A[0][1] = b; A[0][2] = c;
        A[1][0] = b; A[1][1] = d; A[1][2] = e;
        A[2][0] = c; A[2][1] = e; A[2][2] = f;

        CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

        for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

        CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

        nodes->SetAnisoHess(iPoint, i+0, A[0][0]);
        nodes->SetAnisoHess(iPoint, i+1, A[0][1]);
        nodes->SetAnisoHess(iPoint, i+2, A[0][2]);
        nodes->SetAnisoHess(iPoint, i+3, A[1][1]);
        nodes->SetAnisoHess(iPoint, i+4, A[1][2]);
        nodes->SetAnisoHess(iPoint, i+5, A[2][2]);
      }
    }
  }

  if(viscous) {
    for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      //--- viscous mass flux is 0, so start with momentum
      for(iVar = 1; iVar < nVarMetr; iVar++){
        for(iFlux = 0; iFlux < nFluxMetr; iFlux++){
          const unsigned short i = iFlux*nVarMetr*nMetr + iVar*nMetr;

          const su2double a = nodes->GetAnisoViscHess(iPoint, i+0);
          const su2double b = nodes->GetAnisoViscHess(iPoint, i+1);
          const su2double c = nodes->GetAnisoViscHess(iPoint, i+2);
          const su2double d = nodes->GetAnisoViscHess(iPoint, i+3);
          const su2double e = nodes->GetAnisoViscHess(iPoint, i+4);
          const su2double f = nodes->GetAnisoViscHess(iPoint, i+5);

          A[0][0] = a; A[0][1] = b; A[0][2] = c;
          A[1][0] = b; A[1][1] = d; A[1][2] = e;
          A[2][0] = c; A[2][1] = e; A[2][2] = f;

          CNumerics::EigenDecomposition(A, EigVec, EigVal, nDim);

          for(unsigned short iDim = 0; iDim < nDim; ++iDim) EigVal[iDim] = abs(EigVal[iDim]);

          CNumerics::EigenRecomposition(A, EigVec, EigVal, nDim);

          nodes->SetAnisoViscHess(iPoint, i+0, A[0][0]);
          nodes->SetAnisoViscHess(iPoint, i+1, A[0][1]);
          nodes->SetAnisoViscHess(iPoint, i+2, A[0][2]);
          nodes->SetAnisoViscHess(iPoint, i+3, A[1][1]);
          nodes->SetAnisoViscHess(iPoint, i+4, A[1][2]);
          nodes->SetAnisoViscHess(iPoint, i+5, A[2][2]);
        }
      }
    }
  }
}

void CTNE2EulerSolver::SetPrimitive_Limiter(CGeometry *geometry,
                                            CConfig *config) {

  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double dave, LimK, eps1, eps2, dm, dp, du, y, limiter;
  su2double *Primitive, *Primitive_i, *Primitive_j, *LocalMinPrimitive, *LocalMaxPrimitive,
      *GlobalMinPrimitive, *GlobalMaxPrimitive;
  su2double *Coord_i, *Coord_j;
  su2double **Gradient_i, **Gradient_j;


  dave = config->GetRefElemLength();
  LimK = config->GetVenkat_LimiterCoeff();

  if (config->GetKind_SlopeLimit_TNE2() == NO_LIMITER) {

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        nodes->SetLimiter_Primitive(iPoint, iVar,1.0);
      }
    }
  }
  else {
    /*--- Initialize solution max, solution min and limiter in entire domain ---*/
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        nodes->SetSolution_Max(iPoint, iVar, -EPS);
        nodes->SetSolution_Min(iPoint, iVar, EPS);
        nodes->SetLimiter_Primitive(iPoint, iVar, 2.0);
      }
    }

    /*--- Establish bounts for Spekreijse monotonicity by finding max/min values
    of neighbor variables ---*/
    for (iEdge = 0; iEdge < geometry-> GetnEdge(); iEdge++) {

      /*--- Point identification, Normal vector and area ---*/
      iPoint = geometry-> edge[iEdge]->GetNode(0);
      jPoint = geometry-> edge[iEdge]->GetNode(1);

      /*--- Get primitive variables ---*/
      Primitive_i = nodes->GetPrimitive(iPoint);
      Primitive_j = nodes->GetPrimitive(jPoint);

      /*--- Compute the max and min values for nodes i & j ---*/
      for (iVar = 0; iVar < nVar; iVar++){
        du = (Primitive_j[iVar]-Primitive_i[iVar]);
        nodes->SetSolution_Min(iPoint, iVar, min(nodes->GetSolution_Min(iPoint, iVar), du));
        nodes->SetSolution_Max(iPoint, iVar, max(nodes->GetSolution_Max(jPoint, iVar), du));
        nodes->SetSolution_Min(jPoint, iVar, min(nodes->GetSolution_Min(jPoint, iVar), -du));
        nodes->SetSolution_Max(jPoint, iVar, max(nodes->GetSolution_Max(jPoint, iVar), -du));
      }
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

        AD::SetPreaccIn(nodes->GetSolution_Max(iPoint, iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(iPoint, iVar));
        AD::SetPreaccIn(nodes->GetSolution_Max(jPoint, iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(jPoint, iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = nodes->GetSolution_Max(iPoint, iVar);
          else dp = nodes->GetSolution_Min(iPoint, iVar);
          limiter = dp/dm;
        }

        if (limiter < nodes->GetLimiter_Primitive(iPoint, iVar)) {
          nodes->SetLimiter_Primitive(iPoint,  iVar, limiter);
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

        if (limiter < nodes->GetLimiter_Primitive(jPoint, iVar)) {
          nodes->SetLimiter_Primitive(jPoint, iVar, limiter);
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

    /*--- Allocate memory for the max and min primitive value --*/
    LocalMinPrimitive = new su2double [nPrimVarGrad]; GlobalMinPrimitive = new su2double [nPrimVarGrad];
    LocalMaxPrimitive = new su2double [nPrimVarGrad]; GlobalMaxPrimitive = new su2double [nPrimVarGrad];

    /*--- Compute the max value and min value of the solution ---*/
    Primitive = nodes->GetPrimitive(iPoint);
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
      AD::SetPreaccIn(eps2);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {

        if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
          eps1 = LimK * (GlobalMaxPrimitive[iVar] - GlobalMinPrimitive[iVar]);
          eps2 = eps1*eps1;
        }
        else {
          eps1 = LimK*dave;
          eps2 = eps1*eps1*eps1;
        }

        AD::SetPreaccIn(nodes->GetSolution_Max(iPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(iPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Max(jPoint,iVar));
        AD::SetPreaccIn(nodes->GetSolution_Min(jPoint,iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        if ( dm > 0.0 ) dp = nodes->GetSolution_Max(iPoint,iVar);
        else dp = nodes->GetSolution_Min(iPoint,iVar);

        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if (limiter < nodes->GetLimiter_Primitive(iPoint, iVar)) {
          nodes->SetLimiter_Primitive(iPoint, iVar, limiter);
          AD::SetPreaccOut(nodes->GetLimiter_Primitive(iPoint)[iVar]);
        }

        /*-- Repeat for point j on the edge ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if ( dm > 0.0 ) dp = nodes->GetSolution_Max(jPoint,iVar);
        else dp = nodes->GetSolution_Min(jPoint,iVar);

        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if (limiter < nodes->GetLimiter_Primitive(jPoint,iVar)) {
          nodes->SetLimiter_Primitive(jPoint,iVar, limiter);
          AD::SetPreaccOut(nodes->GetLimiter_Primitive(jPoint)[iVar]);
        }

      }

      AD::EndPreacc();

    }

    delete [] LocalMinPrimitive; delete [] GlobalMinPrimitive;
    delete [] LocalMaxPrimitive; delete [] GlobalMaxPrimitive;

  }

  /*--- Limiter MPI ---*/
  InitiateComms(geometry, config, PRIMITIVE_LIMITER);
  CompleteComms(geometry, config, PRIMITIVE_LIMITER);

}

void CTNE2EulerSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {

  unsigned short iSpecies, iEl;
  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
      Energy_FreeStream      = 0.0, ModVel_FreeStreamND = 0.0, Velocity_Reynolds = 0.0,
      Omega_FreeStream       = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
      Density_FreeStream     = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0;

  su2double Length_Ref       = 0.0, Density_Ref   = 0.0, Pressure_Ref     = 0.0, Velocity_Ref = 0.0,
      Temperature_Ref  = 0.0, Time_Ref      = 0.0, Omega_Ref        = 0.0, Force_Ref    = 0.0,
      Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref   = 0.0;

  su2double Pressure_FreeStreamND    = 0.0, Density_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
      Temperature_FreeStreamND = 0.0, Temperature_ve_FreeStreamND = 0.0, Gas_Constant         = 0.0, Gas_ConstantND      = 0.0,
      Viscosity_FreeStreamND   = 0.0, Tke_FreeStreamND     = 0.0,
      Total_UnstTimeND         = 0.0, Delta_UnstTimeND     = 0.0,
      Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0};

  su2double Mass    = 0.0, soundspeed = 0.0, GasConstant_Inf = 0.0, Froude = 0.0,
      rhoCvtr = 0.0, rhoE       = 0.0, sqvel           = 0.0,
      denom   = 0.0, num        = 0.0, conc            = 0.0;

  unsigned short iDim,nEl,nHeavy,*nElStates;

  /*--- Local variables ---*/
  su2double Alpha         = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta          = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach          = config->GetMach();
  su2double Reynolds      = config->GetReynolds();

  bool unsteady           = (config->GetTime_Marching() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);
  bool tkeNeeded          = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);
  bool ionization         = config->GetIonization();

  su2double *Ms, *xi, *Tref, *hf, *thetav, **thetae, **g, rhos, Ef, Ee, Ev;
  su2double RuSI          = UNIVERSAL_GAS_CONSTANT;
  su2double Ru            = 1000.0*RuSI;
  su2double T             = config->GetTemperature_FreeStream();
  su2double Tve           = config->GetTemperature_ve_FreeStream();

  unsigned short nSpecies = config->GetnSpecies();
  Tref                    = config->GetRefTemperature();
  Ms                      = config->GetMolar_Mass();
  xi                      = config->GetRotationModes();
  nElStates               = config->GetnElStates();
  hf                      = config->GetEnthalpy_Formation();
  thetav                  = config->GetCharVibTemp();
  thetae                  = config->GetCharElTemp();
  g                       = config->GetElDegeneracy();
  ionization = config->GetIonization();

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Compressible non dimensionalization ---*/

  /*--- Compute Gas Constant ---*/
  // This needs work for Ionization and such
  MassFrac_Inf  = config->GetMassFrac_FreeStream();
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    Mass += MassFrac_Inf[iSpecies] * Ms[iSpecies];
  GasConstant_Inf = Ru / Mass; config->SetGas_Constant(GasConstant_Inf);

  /*--- Compute the Free Stream Pressure, Temperatrue, and Density ---*/
  Pressure_FreeStream     = config->GetPressure_FreeStream();
  Temperature_FreeStream  = T;

  /*--- Compute the density using mixtures ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    denom += MassFrac_Inf[iSpecies] * (Ru/Ms[iSpecies]) * T;
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    denom += MassFrac_Inf[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Tve;
  Density_FreeStream = Pressure_FreeStream / denom;

  /*--- Calculate sound speed and extract velocities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    conc += MassFrac_Inf[iSpecies]*Density_FreeStream/Ms[iSpecies];
    rhoCvtr += Density_FreeStream*MassFrac_Inf[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
  }
  soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * Pressure_FreeStream/Density_FreeStream);

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
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  sqvel = ModVel_FreeStream;
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    // Species density
    rhos = MassFrac_Inf[iSpecies]*Density_FreeStream;

    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

    // Species vibrational energy
    if (thetav[iSpecies] != 0.0)
      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
    }
    Ee = Ru/Ms[iSpecies] * (num/denom);

    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*sqvel);
  }

  /*--- Viscous initialization ---*/
  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
          To accomplish this, simply set the non-dimensional coefficients to the
          dimensional ones. This will be overruled later.---*/
    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());
    config->SetMu_ConstantND(config->GetMu_Constant());

    /*--- Reynolds based initialization ---*/
    if (reynolds_init) {

      /*--- First, check if there is mesh motion. If yes, use the Mach
         number relative to the body to initialize the flow. ---*/

      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;

      /*--- For viscous flows, pressure will be computed from a density
            that is found from the Reynolds number. The viscosity is computed
            from the dimensional version of Sutherland's law or the constant
            viscosity, depending on the input option.---*/

      // THIS NEEDS TO BE REVISITED
      Viscosity_FreeStream = 1.853E-5*(pow(Temperature_FreeStream/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_FreeStream+110.3));
      Density_FreeStream   = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds* config->GetLength_Reynolds());
      Pressure_FreeStream  = Density_FreeStream*GasConstant_Inf*Temperature_FreeStream;
      Energy_FreeStream    = Pressure_FreeStream/(Density_FreeStream*Gamma_Minus_One)+0.5*ModVel_FreeStream *ModVel_FreeStream;

      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      config->SetPressure_FreeStream(Pressure_FreeStream);
    }

    /*--- Thermodynamics quantities based initialization ---*/
    else {
      Viscosity_FreeStream = 1.853E-5*(pow(Temperature_FreeStream/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_FreeStream+110.3));
      Density_FreeStream   = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds* config->GetLength_Reynolds());
      Pressure_FreeStream  = Density_FreeStream*GasConstant_Inf*Temperature_FreeStream;
      Energy_FreeStream    = Pressure_FreeStream/(Density_FreeStream*Gamma_Minus_One)+0.5*ModVel_FreeStream *ModVel_FreeStream ;
    }

    /*--- Turbulence kinetic energy ---*/
    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }
  else {

    /*--- For inviscid flow, energy is calculated from the specified
       FreeStream quantities using the proper gas law. ---*/
    Energy_FreeStream    = rhoE/Density_FreeStream;
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
    Density_Ref       = Density_FreeStream;                  // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;              // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = config->GetDensity_Ref()*Velocity_Ref*Velocity_Ref*Length_Ref*Length_Ref; config->SetForce_Ref(Force_Ref);
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

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  Temperature_ve_FreeStreamND = Tve/config->GetTemperature_Ref(); config->SetTemperature_ve_FreeStreamND(Temperature_ve_FreeStreamND);
  Gas_ConstantND           = config->GetGas_Constant()/Gas_Constant_Ref;          config->SetGas_ConstantND(Gas_ConstantND);

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

  Energy_FreeStreamND = Pressure_FreeStreamND/(Density_FreeStreamND*Gamma_Minus_One)+0.5*ModVel_FreeStreamND *ModVel_FreeStreamND;

  if (viscous) {

    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

    /*--- Sutherland's model ---*/

    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());

    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);

  }

  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is the master node and first domain ---*/

  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {

    cout.precision(6);

    if (viscous) {
      cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
      cout << "based on the free-stream temperature and a density computed" << endl;
      cout << "from the Reynolds number." << endl;
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "temperature and pressure using the ideal gas law." << endl;
    }

    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;

    cout <<"-- Input conditions:"<< endl;

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
        cout << "Laminar Viscosity: " << config->GetMu_Constant();
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
      switch (config->GetKind_ConductivityModel()) {

      case CONSTANT_PRANDTL:
        cout << "Conductivity Model: CONSTANT_PRANDTL  "<< endl;
        cout << "Prandtl: " << config->GetPrandtl_Lam()<< endl;
        break;

      case CONSTANT_CONDUCTIVITY:
        cout << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< endl;
        cout << "Molecular Conductivity: " << config->GetKt_Constant()<< " W/m^2.K." << endl;
        cout << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< endl;
        break;

      }
    }

    cout << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+Mach*Mach*0.5*(Gamma-1.0), Gamma/(Gamma-1.0) );
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Free-stream temperature: " << config->GetTemperature_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;

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

    cout << "Magnitude: "   << config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m/s (" << config->GetModVel_FreeStream()*1.94384 << " KTS)." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s (" << config->GetModVel_FreeStream()*0.592484 << " KTS)." << endl;

    cout << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;

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

    cout << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;

    cout << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;

    cout << "Reference temperature: " << config->GetTemperature_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " K." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " R." << endl;

    cout << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " kg/m^3." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " slug/ft^3." << endl;

    cout << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m/s." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft/s." << endl;

    cout << "Reference energy per unit mass: " << config->GetEnergy_Ref();
    if (config->GetSystemMeasurements() == SI) cout << " m^2/s^2." << endl;
    else if (config->GetSystemMeasurements() == US) cout << " ft^2/s^2." << endl;

    if (viscous) {
      cout << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " N.s/m^2." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.s/ft^2." << endl;
      cout << "Reference conductivity: " << config->GetConductivity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << " W/m^2.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf/ft.s.R." << endl;
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

    cout << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
    cout << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << endl;

    cout << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << endl;

    cout << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << endl;

    if (nDim == 2) {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      cout << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      cout << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    cout << "Magnitude: "    << config->GetModVel_FreeStreamND() << endl;

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

void CTNE2EulerSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) {
  unsigned short iDim, jDim, iVar, jVar;
  su2double Beta, local_Mach, Beta2, rho, enthalpy, soundspeed, sq_vel;
  su2double *U_i = NULL;
  su2double Beta_max = config->GetmaxTurkelBeta();
  su2double Mach_infty2, Mach_lim2, aux, parameter;


  cout << "This dont work" << endl;
  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(nodes->GetVelocity2(iPoint))/nodes->GetSoundSpeed(iPoint);

  /*--- Weiss and Smith Preconditionng ---*/
  Mach_infty2 = pow(config->GetMach(),2.0);
  Mach_lim2   = pow(0.00001,2.0);
  aux         = max(pow(local_Mach,2.0),Mach_lim2);
  parameter   = min(1.0, max(aux,Beta_max*Mach_infty2));

  U_i = nodes->GetSolution(iPoint);

  rho = U_i[0];
  enthalpy = nodes->GetEnthalpy(iPoint);
  soundspeed = nodes->GetSoundSpeed(iPoint);
  sq_vel = nodes->GetVelocity2(iPoint);

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

void CTNE2EulerSolver::Evaluate_ObjFunc(CConfig *config) {

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
  case EQUIVALENT_AREA:
    Total_ComboObj+=Weight_ObjFunc*Total_CEquivArea;
    break;
  case NEARFIELD_PRESSURE:
    Total_ComboObj+=Weight_ObjFunc*Total_CNearFieldOF;
    break;
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
  case SURFACE_MACH:
    Total_ComboObj+=Weight_ObjFunc*config->GetSurface_Mach(0);
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
  case TOTAL_AVG_TEMPERATURE:
    Total_ComboObj+=Weight_ObjFunc*config->GetSurface_Temperature(0);
    break;
  case CUSTOM_OBJFUNC:
    Total_ComboObj+=Weight_ObjFunc*Total_Custom_ObjFunc;
    break;
  default:
    break;
  }

}

void CTNE2EulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                     CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim, jDim, iSpecies, iVar, jVar, kVar;
  unsigned long iPoint, iVertex;

  su2double *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3], *NormalArea,
      ProjGridVel = 0.0, turb_ke;
  su2double Density_b, StaticEnergy_b, Enthalpy_b, *Velocity_b, Kappa_b, Chi_b, Energy_b, VelMagnitude2_b, Pressure_b;
  su2double Density_i, *Velocity_i, ProjVelocity_i = 0.0, Energy_i, VelMagnitude2_i;
  su2double **Jacobian_b, **DubDu;
  su2double rhoCvtr, rhoCvve, rho_el, RuSI, Ru;
  su2double rho, cs, P, rhoE, rhoEve, conc, Beta;
  su2double *u, *Ms, *dPdU;
  su2double *V_domain,*V_reflected;
  su2double *U_domain,*U_reflected;
  su2double *dPdU_reflected, *dTdU_reflected, *dTvedU_reflected;

  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));

  /*--- Allocate arrays ---*/
  Normal     = new su2double[nDim];
  NormalArea = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_i = new su2double[nDim];
  Jacobian_b = new su2double*[nVar];
  DubDu      = new su2double*[nVar];
  u          = new su2double[nDim];

  U_reflected = new su2double[nVar];
  V_reflected = new su2double[nPrimVar];

  dPdU_reflected   = new su2double[nVar];
  dTdU_reflected   = new su2double[nVar];
  dTvedU_reflected = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_b[iVar] = new su2double[nVar];
    DubDu[iVar] = new su2double[nVar];
  }

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  conv_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  conv_numerics->SetPIndex      ( nodes->GetPIndex()       );
  conv_numerics->SetTIndex      ( nodes->GetTIndex()       );
  conv_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  conv_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  conv_numerics->SetHIndex      ( nodes->GetHIndex()       );
  conv_numerics->SetAIndex      ( nodes->GetAIndex()       );
  conv_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  conv_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  visc_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  visc_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  visc_numerics->SetPIndex      ( nodes->GetPIndex()       );
  visc_numerics->SetTIndex      ( nodes->GetTIndex()       );
  visc_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  visc_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  visc_numerics->SetHIndex      ( nodes->GetHIndex()       );
  visc_numerics->SetAIndex      ( nodes->GetAIndex()       );
  visc_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  visc_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  /*--- Load parameters from the config class ---*/
  Ms = config->GetMolar_Mass();

  /*--- Rename for convenience ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negative for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      /*--- Calculate parameters from the geometry ---*/
      Area   = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++){
        Normal[iDim]     = -Normal[iDim];
        UnitNormal[iDim] = -Normal[iDim]/Area;
      }
      conv_numerics->SetNormal(Normal);

      // /*--- Retrieve the pressure on the vertex ---*/
      // P   = nodes->GetPressure(iPoint);

      /*--- Compute the residual ---*/
      turb_ke = 0.0;
      if (tkeNeeded) turb_ke=solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint, 0);

      /*--- Retrieve solution at the boundary node & allocate reflected state ---*/
      U_domain = nodes->GetSolution(iPoint);
      V_domain = nodes->GetPrimitive(iPoint);
      for(iVar = 0; iVar < nVar; iVar++)
        U_reflected[iVar] = nodes->GetSolution(iPoint, iVar);
      for(iVar = 0; iVar < nPrimVar; iVar++)
        V_reflected[iVar] = nodes->GetPrimitive(iPoint, iVar);

      ProjVelocity_i = nodes->GetProjVel(iPoint,UnitNormal);
      Density_i      = nodes->GetDensity(iPoint);
      for (iDim = 0; iDim < nDim; iDim++) {
        V_reflected[nSpecies+2+iDim] = nodes->GetVelocity(iPoint,iDim) - 2.*ProjVelocity_i*UnitNormal[iDim];
        U_reflected[nSpecies+iDim] = Density_i*V_reflected[nSpecies+2+iDim];
        // nodes->SetSolution(iPoint,nSpecies+iDim, Density_i*Velocity_i[iDim]);
      }
      nodes->CalcdPdU(  V_reflected, nodes->GetEve(iPoint), config, dPdU_reflected  );
      nodes->CalcdTdU(  V_reflected, config, dTdU_reflected  );
      nodes->CalcdTvedU(V_reflected, nodes->GetEve(iPoint), config, dTvedU_reflected);

      /*--- Pass conserved & primitive variables to CNumerics ---*/
      conv_numerics->SetConservative(U_domain, U_reflected);
      conv_numerics->SetPrimitive(V_domain, V_reflected);

      /*--- Pass supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU(nodes->GetdPdU(iPoint),     dPdU_reflected);
      conv_numerics->SetdTdU(nodes->GetdTdU(iPoint),     dTdU_reflected);
      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), dTvedU_reflected);
      conv_numerics->SetEve(nodes->GetEve(iPoint),       nodes->GetEve(iPoint));
      conv_numerics->SetCvve(nodes->GetCvve(iPoint),     nodes->GetCvve(iPoint));

      /*--- Compute the convective residual (and Jacobian) ---*/
      // Note: This uses the specified boundary num. method specified in driver_structure.cpp
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      // /*--- Apply the flow-tangency b.c. to the convective flux ---*/
      // for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      //   Residual[iSpecies] = 0.0;
      // for (iDim = 0; iDim < nDim; iDim++)
      //   Residual[nSpecies+iDim] = P * UnitNormal[iDim] * Area;
      // Residual[nSpecies+nDim]   = 0.0;
      // Residual[nSpecies+nDim+1] = 0.0;

      // /*--- Add the Reynolds stress tensor contribution ---*/
      // if (tkeNeeded) {
      //   for (iDim = 0; iDim < nDim; iDim++)
      //     Residual[nSpecies+iDim+1] += (2.0/3.0)*Density_b*turb_ke*NormalArea[iDim];
      // }

      /*--- Add value to the residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration. ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }

      // /*--- If using implicit time-stepping, calculate b.c. contribution to Jacobian ---*/
      // if (implicit) {

      //   /*--- Initialize Jacobian ---*/
      //   for (iVar = 0; iVar < nVar; iVar++)
      //     for (jVar = 0; jVar < nVar; jVar++)
      //       Jacobian_i[iVar][jVar] = 0.0;

      //   /*--- Calculate state i ---*/
      //   rho     = nodes->GetDensity(iPoint);
      //   rhoCvtr = nodes->GetRhoCv_tr(iPoint);
      //   rhoCvve = nodes->GetRhoCv_ve(iPoint);
      //   rhoE    = nodes->GetSolution(iPoint,nSpecies+nDim);
      //   rhoEve  = nodes->GetSolution(iPoint,nSpecies+nDim+1);
      //   dPdU    = nodes->GetdPdU(iPoint);
      //   for (iDim = 0; iDim < nDim; iDim++)
      //     u[iDim] = nodes->GetVelocity(iPoint,iDim);

      //   /*--- If free electrons are present, retrieve the electron gas density ---*/
      //   if (config->GetIonization()) rho_el = nodes->GetMassFraction(iPoint,nSpecies-1) * rho;
      //   else                         rho_el = 0.0;

      //   conc = 0.0;
      //   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      //     cs    = nodes->GetMassFraction(iPoint,iSpecies);
      //     conc += cs * rho/Ms[iSpecies];

      //     /////// NEW //////
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Jacobian_i[nSpecies+iDim][iSpecies] = dPdU[iSpecies] * UnitNormal[iDim];
      //       Jacobian_i[iSpecies][nSpecies+iDim] = cs * UnitNormal[iDim];
      //     }
      //   }

      //   Beta = Ru*conc/rhoCvtr;

      //   for (iDim = 0; iDim < nDim; iDim++) {
      //     for (jDim = 0; jDim < nDim; jDim++) {
      //       Jacobian_i[nSpecies+iDim][nSpecies+jDim] = u[iDim]*UnitNormal[jDim]
      //           + dPdU[nSpecies+jDim]*UnitNormal[iDim];
      //     }
      //     Jacobian_i[nSpecies+iDim][nSpecies+nDim]   = dPdU[nSpecies+nDim]  *UnitNormal[iDim];
      //     Jacobian_i[nSpecies+iDim][nSpecies+nDim+1] = dPdU[nSpecies+nDim+1]*UnitNormal[iDim];

      //     Jacobian_i[nSpecies+nDim][nSpecies+iDim]   = (rhoE+P)/rho * UnitNormal[iDim];
      //     Jacobian_i[nSpecies+nDim+1][nSpecies+iDim] = rhoEve/rho   * UnitNormal[iDim];
      //   }

      //   /*--- Integrate over the dual-grid area ---*/
      //   for (iVar = 0; iVar < nVar; iVar++)
      //     for (jVar = 0; jVar < nVar; jVar++)
      //       Jacobian_i[iVar][jVar] = Jacobian_i[iVar][jVar] * Area;

      //   /*--- Apply the contribution to the system ---*/
      //   Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);

      // }
    }
  }
  delete [] u;
}

// void CTNE2EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solution_container,
//                                     CNumerics *conv_numerics, CNumerics *visc_numerics,
//                                     CConfig *config, unsigned short val_marker) {
//   unsigned short iDim;
//   unsigned long iVertex, iPoint, Point_Normal;

//   su2double *V_infty, *V_domain;
//   su2double *U_domain,*U_infty;

//   /*--- Set booleans from configuration parameters ---*/
//   bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
//   bool viscous  = config->GetViscous();
//   //bool tkeNeeded = (((config->GetKind_Solver() == RANS ) ||
//   //                   (config->GetKind_Solver() == DISC_ADJ_RANS))
//   //                  && (config->GetKind_Turb_Model() == SST));

//   /*--- Allocate arrays ---*/
//   su2double *Normal = new su2double[nDim];

//   /*--- Pass structure of the primitive variable vector to CNumerics ---*/
//   conv_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
//   conv_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
//   conv_numerics->SetPIndex      ( nodes->GetPIndex()       );
//   conv_numerics->SetTIndex      ( nodes->GetTIndex()       );
//   conv_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
//   conv_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
//   conv_numerics->SetHIndex      ( nodes->GetHIndex()       );
//   conv_numerics->SetAIndex      ( nodes->GetAIndex()       );
//   conv_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
//   conv_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

//   visc_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
//   visc_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
//   visc_numerics->SetPIndex      ( nodes->GetPIndex()       );
//   visc_numerics->SetTIndex      ( nodes->GetTIndex()       );
//   visc_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
//   visc_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
//   visc_numerics->SetHIndex      ( nodes->GetHIndex()       );
//   visc_numerics->SetAIndex      ( nodes->GetAIndex()       );
//   visc_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
//   visc_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

//   /*--- Loop over all the vertices on this boundary (val_marker) ---*/
//   for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//     iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

//     /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//     if (geometry->node[iPoint]->GetDomain()) {

//       /*--- Retrieve index of the closest interior node ---*/
//       Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

//       /*--- Pass boundary node normal to CNumerics ---*/
//       geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
//       for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
//       conv_numerics->SetNormal(Normal);

//       /*--- Retrieve solution at the boundary node & free-stream ---*/
//       U_domain = nodes->GetSolution(iPoint);
//       V_domain = nodes->GetPrimitive(iPoint);
//       U_infty  = node_infty->GetSolution(0);
//       V_infty  = node_infty->GetPrimitive(0);

//       /*--- Pass conserved & primitive variables to CNumerics ---*/
//       conv_numerics->SetConservative(U_domain, U_infty);
//       conv_numerics->SetPrimitive(V_domain, V_infty);

//       /*--- Pass supplementary information to CNumerics ---*/
//       conv_numerics->SetdPdU(nodes->GetdPdU(iPoint),     node_infty->GetdPdU(0));
//       conv_numerics->SetdTdU(nodes->GetdTdU(iPoint),     node_infty->GetdTdU(0));
//       conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));

//       /*--- Compute the convective residual (and Jacobian) ---*/
//       // Note: This uses the specified boundary num. method specified in driver_structure.cpp
//       conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

//       /*--- Apply contribution to the linear system ---*/
//       LinSysRes.AddBlock(iPoint, Residual);
//       if (implicit)
//         Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

//       /*--- Viscous contribution ---*/
//       if (viscous) {
//         visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
//                                 geometry->node[Point_Normal]->GetCoord() );
//         visc_numerics->SetNormal(Normal);

//         /*--- Primitive variables, and gradient ---*/
//         visc_numerics->SetConservative(nodes->GetSolution(iPoint),
//                                        node_infty->GetSolution(0) );
//         visc_numerics->SetConsVarGradient(nodes->GetGradient(iPoint),
//                                           node_infty->GetGradient(0) );
//         visc_numerics->SetPrimitive(nodes->GetPrimitive(iPoint),
//                                     node_infty->GetPrimitive(0) );
//         visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
//                                           node_infty->GetGradient_Primitive(0) );

//         /*--- Pass supplementary information to CNumerics ---*/
//         visc_numerics->SetdPdU(nodes->GetdPdU(iPoint), node_infty->GetdPdU(0));
//         visc_numerics->SetdTdU(nodes->GetdTdU(iPoint), node_infty->GetdTdU(0));
//         visc_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));

//         /*--- Species diffusion coefficients ---*/
//         visc_numerics->SetDiffusionCoeff(nodes->GetDiffusionCoeff(iPoint),
//                                          node_infty->GetDiffusionCoeff(0) );

//         /*--- Laminar viscosity ---*/
//         visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint),
//                                            node_infty->GetLaminarViscosity(0) );

//         /*--- Thermal conductivity ---*/
//         visc_numerics->SetThermalConductivity(nodes->GetThermalConductivity(iPoint),
//                                               node_infty->GetThermalConductivity(0));

//         /*--- Vib-el. thermal conductivity ---*/
//         visc_numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint),
//                                                  node_infty->GetThermalConductivity_ve(0) );

//         /*--- Compute and update residual ---*/
//         visc_numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
//         LinSysRes.SubtractBlock(iPoint, Res_Visc);
//         if (implicit) {
//           Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//         }
//       }
//     }
//   }

//   /*--- Free locally allocated memory ---*/
//   delete [] Normal;
// }

void CTNE2EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solution_container,
                                    CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iSpecies;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double *GridVel;
  su2double Area, UnitNormal[3] = {0.0,0.0,0.0};
  su2double Density, Pressure, Temperature, Temperature_ve, Energy,  Velocity[3] = {0.0,0.0,0.0};
  su2double Density_Bound, Pressure_Bound, Vel_Bound[3] = {0.0,0.0,0.0};
  su2double Density_Infty, Pressure_Infty, Vel_Infty[3] = {0.0,0.0,0.0};
  su2double SoundSpeed, Entropy, Velocity2, Vn;
  su2double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound, Mach_Bound;
  su2double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  su2double RiemannPlus, RiemannMinus;

  su2double Alpha, Beta, *Mvec, *Mvec_Inf;
  su2double *Ys;

  su2double *V_infty, *V_domain;
  su2double *U_domain,*U_infty;

  su2double Gas_Constant = config->GetGas_ConstantND();

  /*--- Set booleans from configuration parameters ---*/
  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool viscous  = config->GetViscous();
  bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);

  /*--- Allocate arrays ---*/
  su2double *Normal = new su2double[nDim];

  U_infty = new su2double[nVar];
  V_infty = new su2double[nPrimVar];
  Ys = new su2double[nSpecies];

  /*--- Vectorize free stream Mach number based on AoA & AoS ---*/
  Mvec     = new su2double[nDim];
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

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  conv_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  conv_numerics->SetPIndex      ( nodes->GetPIndex()       );
  conv_numerics->SetTIndex      ( nodes->GetTIndex()       );
  conv_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  conv_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  conv_numerics->SetHIndex      ( nodes->GetHIndex()       );
  conv_numerics->SetAIndex      ( nodes->GetAIndex()       );
  conv_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  conv_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  visc_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  visc_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  visc_numerics->SetPIndex      ( nodes->GetPIndex()       );
  visc_numerics->SetTIndex      ( nodes->GetTIndex()       );
  visc_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  visc_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  visc_numerics->SetHIndex      ( nodes->GetHIndex()       );
  visc_numerics->SetAIndex      ( nodes->GetAIndex()       );
  visc_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  visc_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Retrieve index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Pass boundary node normal to CNumerics ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at the boundary node & free-stream ---*/
      U_domain = nodes->GetSolution(iPoint);
      V_domain = nodes->GetPrimitive(iPoint);

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

      Density_Bound = V_domain[nSpecies+nDim+3];
      Vel2_Bound = 0.0; Vn_Bound = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Bound[iDim] = V_domain[nSpecies+2+nDim];
        Vel2_Bound     += Vel_Bound[iDim]*Vel_Bound[iDim];
        Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
      }
      Pressure_Bound   = nodes->GetPressure(iPoint);
      SoundSpeed_Bound = nodes->GetSoundSpeed(iPoint);
      Entropy_Bound    = pow(Density_Bound, Gamma)/Pressure_Bound;

      /*--- Store the primitive variable state for the freestream. Project
         the freestream velocity vector into the local normal direction,
         i.e. compute v_infty.n. ---*/

      Density_Infty = node_infty->GetDensity(0);
      Vel2_Infty = 0.0; Vn_Infty = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Infty[iDim] = node_infty->GetVelocity(0, iDim);
        Vel2_Infty     += Vel_Infty[iDim]*Vel_Infty[iDim];
        Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
      }
      Pressure_Infty   = node_infty->GetPressure(0);
      SoundSpeed_Infty = node_infty->GetSoundSpeed(0);
      Entropy_Infty    = pow(Density_Infty, Gamma)/Pressure_Infty;

      /*--- Adjust the normal freestream velocity for grid movement ---*/

      Qn_Infty = Vn_Infty;
      if (dynamic_grid) {
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
        Temperature_ve = V_domain[nSpecies+1];
        for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
          Ys[iSpecies] = V_domain[iSpecies]/V_domain[nSpecies+nDim+3];
        }
      } else  {
        /*--- Inflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Infty[iDim] + (Vn-Vn_Infty)*UnitNormal[iDim];
        Entropy = Entropy_Infty;
        Temperature_ve = Temperature_ve_Inf;
        for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
          Ys[iSpecies] = MassFrac_Inf[iSpecies];
        }
      }

      /*--- Recompute the primitive variables. ---*/

      Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      Temperature = Pressure/(Gas_Constant*Density);
      // TODO: Fix these later!
      for (iDim = 0; iDim < nDim; iDim++) {
        Mvec[iDim] = Velocity[iDim]/SoundSpeed;
      }

      // if (Qn_Infty > 0){
      //   // Subsonic outflow
      //   // TODO: check these (mainly Tve)
      //   if (Qn_Infty <= SoundSpeed_Infty) {
      //     Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      //     Velocity2 = 0.0;
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Velocity2 += Velocity[iDim]*Velocity[iDim];
      //     }
      //     Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      //     Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      //     Temperature = Pressure/(Gas_Constant*Density);
      //     Temperature_ve = V_domain[nSpecies+1];
      //     for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
      //       Ys[iSpecies] = MassFrac_Inf[iSpecies];
      //     }
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Mvec[iDim] = Velocity[iDim]/SoundSpeed;
      //     }
      //   }
      //   // Supersonic outflow
      //   // All values set to interior
      //   else {
      //     Density = V_domain[nSpecies+nDim+3];
      //     Velocity2 = 0.0;
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Velocity2 += V_domain[nSpecies+2+iDim]*V_domain[nSpecies+2+iDim];
      //     }
      //     Pressure = V_domain[nSpecies+nDim+2];
      //     SoundSpeed = V_domain[nSpecies+nDim+5];
      //     Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      //     Temperature = V_domain[nSpecies];
      //     Temperature_ve = V_domain[nSpecies+1];
      //     for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
      //       Ys[iSpecies] = V_domain[iSpecies]/Density;
      //     }
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Mvec[iDim] = V_domain[nSpecies+2+iDim]/SoundSpeed;
      //     }
      //   }
      // }
      // else {
      //   // Supersonic inflow
      //   // All values set to freestream
      //   if (Qn_Infty < -SoundSpeed_Infty) {
      //     Density = node_infty->GetDensity(0);
      //     Velocity2 = 0.0;
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Velocity2 += Vel_Infty[iDim]*Vel_Infty[iDim];
      //     }
      //     Pressure = node_infty->GetPressure(0);
      //     SoundSpeed = node_infty->GetSoundSpeed(0);
      //     Energy = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      //     Temperature = Pressure/(Gas_Constant*Density);
      //     Temperature_ve = node_infty->GetTemperature_ve(0);
      //     for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
      //       Ys[iSpecies] = MassFrac_Inf[iSpecies];
      //     }
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Mvec[iDim] = Mvec_Inf[iDim];
      //     }
      //   }
      //   // Subsonic inflow
      //   // TODO: check these (mainly Tve)
      //   else {
      //     Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      //     Velocity2 = 0.0;
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Velocity2 += Velocity[iDim]*Velocity[iDim];
      //     }
      //     Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      //     Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      //     Temperature = Density*SoundSpeed*SoundSpeed/Gamma;
      //     Temperature_ve = Temperature_ve_Inf;
      //     for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
      //       Ys[iSpecies] = V_domain[iSpecies]/Density;
      //     }
      //     for (iDim = 0; iDim < nDim; iDim++) {
      //       Mvec[iDim] = Velocity[iDim]/SoundSpeed;
      //     }
      //   }
      // }
      if (tkeNeeded) Energy += GetTke_Inf();

      /*--- Store new primitive state for computing the flux. ---*/

      CTNE2EulerVariable *node_bc = new CTNE2EulerVariable(Pressure, Ys, Mvec, Temperature,
                                                           Temperature_ve, 1, nDim, nVar,
                                                           nPrimVar, nPrimVarGrad, config);

      const bool check_bc = node_bc->SetPrimVar_Compressible(0,config);

      /*--- Pass conserved & primitive variables to CNumerics ---*/
      conv_numerics->SetConservative(U_domain, node_bc->GetSolution(0));
      conv_numerics->SetPrimitive(V_domain, node_bc->GetPrimitive(0));

      /*--- Pass supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU(nodes->GetdPdU(iPoint),     node_bc->GetdPdU(0));
      conv_numerics->SetdTdU(nodes->GetdTdU(iPoint),     node_bc->GetdTdU(0));
      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_bc->GetdTvedU(0));
      conv_numerics->SetEve(nodes->GetEve(iPoint),       node_bc->GetEve(0));
      conv_numerics->SetCvve(nodes->GetCvve(iPoint),     node_bc->GetCvve(0));

      /*--- Compute the convective residual (and Jacobian) ---*/
      // Note: This uses the specified boundary num. method specified in driver_structure.cpp
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Apply contribution to the linear system ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      if (viscous) {
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord() );
        visc_numerics->SetNormal(Normal);

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetConservative(nodes->GetSolution(iPoint),
                                       node_infty->GetSolution(0) );
        visc_numerics->SetConsVarGradient(nodes->GetGradient(iPoint),
                                          node_infty->GetGradient(0) );
        visc_numerics->SetPrimitive(nodes->GetPrimitive(iPoint),
                                    node_infty->GetPrimitive(0) );
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          node_infty->GetGradient_Primitive(0) );

        /*--- Pass supplementary information to CNumerics ---*/
        visc_numerics->SetdPdU(nodes->GetdPdU(iPoint), node_infty->GetdPdU(0));
        visc_numerics->SetdTdU(nodes->GetdTdU(iPoint), node_infty->GetdTdU(0));
        visc_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));

        /*--- Species diffusion coefficients ---*/
        visc_numerics->SetDiffusionCoeff(nodes->GetDiffusionCoeff(iPoint),
                                         node_infty->GetDiffusionCoeff(0) );

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint),
                                           node_infty->GetLaminarViscosity(0) );

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(nodes->GetThermalConductivity(iPoint),
                                              node_infty->GetThermalConductivity(0));

        /*--- Vib-el. thermal conductivity ---*/
        visc_numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint),
                                                 node_infty->GetThermalConductivity_ve(0) );

        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Res_Visc);
        if (implicit) {
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;
}

void CTNE2EulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solution_container,
                                CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  cout << "This dont work" << endl;
  unsigned short iVar, iDim, iSpecies, RHO_INDEX, nSpecies;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
      Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
      alpha, aa, bb, cc, dd, Area, UnitaryNormal[3];

  bool implicit             = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  su2double Two_Gamma_M1    = 2.0/Gamma_Minus_One;
  su2double Gas_Constant    = config->GetGas_ConstantND();
  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded            = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                               (config->GetKind_Turb_Model() == SST));
  bool gravity              = (config->GetGravityForce());
  bool viscous              = config->GetViscous();

  su2double *U_domain = new su2double[nVar];      su2double *U_inlet = new su2double[nVar];
  su2double *V_domain = new su2double[nPrimVar];  su2double *V_inlet = new su2double[nPrimVar];
  su2double *Normal   = new su2double[nDim];
  su2double UnitNormal[3];

  nSpecies = config->GetnSpecies();
  su2double Spec_Density[nSpecies]; //QT THROWING AN ERROR WHEN USING nSPECIES

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
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
      case TOTAL_CONDITIONS:

        /*--- Retrieve the specified total conditions for this inlet. ---*/
        P_Total  = config->GetInlet_Ptotal(Marker_Tag);
        T_Total  = config->GetInlet_Ttotal(Marker_Tag);
        Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_Total /= config->GetPressure_Ref();
        T_Total /= config->GetTemperature_Ref();

        /*--- Store primitives and set some variables for clarity. ---*/
        Density = V_domain[RHO_INDEX];
        Velocity2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = U_domain[nSpecies+iDim]/Density;
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy      = U_domain[nVar-2]/Density;
        Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
        H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
        SoundSpeed2 = Gamma*Pressure/Density;

        /*--- Compute the acoustic Riemann invariant that is extrapolated
           from the domain interior. ---*/
        Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
        for (iDim = 0; iDim < nDim; iDim++)
          Riemann += Velocity[iDim]*UnitaryNormal[iDim];

        /*--- Total speed of sound ---*/
        SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

        /*--- Dot product of normal and flow direction. This should
           be negative due to outward facing boundary normal convention. ---*/
        alpha = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];

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
          Velocity[iDim] = nodes->GetVelocity(iPoint, iDim);
        Pressure    = nodes->GetPressure(iPoint);
        SoundSpeed2 = Gamma*Pressure/U_domain[0];

        /*--- Compute the acoustic Riemann invariant that is extrapolated
           from the domain interior. ---*/
        Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
        for (iDim = 0; iDim < nDim; iDim++)
          Riemann += Velocity[iDim]*UnitaryNormal[iDim];

        /*--- Speed of sound squared for fictitious inlet state ---*/
        SoundSpeed2 = Riemann;
        for (iDim = 0; iDim < nDim; iDim++)
          SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitaryNormal[iDim];

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
      }

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetConservative(U_domain, U_inlet);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

//      /*--- Viscous contribution ---*/
//      if (viscous) {

//        /*--- Set the normal vector and the coordinates ---*/
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

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
}

void CTNE2EulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solution_container,
                                 CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iVar, iDim, iSpecies,iEl, nHeavy, nEl, *nElStates;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Pressure, P_Exit, Velocity[3], Temperature, Tve,
      Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
      Area, UnitaryNormal[3];
  su2double *Ms, *xi, *thetav, **thetae, *Tref, **g, *hf, RuSI, Ru;
  su2double Ef, Ev, Ee;
  su2double thoTve, exptv, thsqr, Cvvs, Cves;
  su2double num, num2, num3, denom;

  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool implicit           = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool grid_movement      = config->GetGrid_Movement();
  bool viscous            = config->GetViscous();
  bool gravity            = config->GetGravityForce();
  bool ionization         = config->GetIonization();

  su2double *U_domain = new su2double[nVar];      su2double *U_outlet = new su2double[nVar];
  su2double *V_domain = new su2double[nPrimVar];  su2double *V_outlet = new su2double[nPrimVar];
  su2double *Normal   = new su2double[nDim];
  su2double *Ys       = new su2double[nSpecies];

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  conv_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  conv_numerics->SetPIndex      ( nodes->GetPIndex()       );
  conv_numerics->SetTIndex      ( nodes->GetTIndex()       );
  conv_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  conv_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  conv_numerics->SetHIndex      ( nodes->GetHIndex()       );
  conv_numerics->SetAIndex      ( nodes->GetAIndex()       );
  conv_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  conv_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  visc_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  visc_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  visc_numerics->SetPIndex      ( nodes->GetPIndex()       );
  visc_numerics->SetTIndex      ( nodes->GetTIndex()       );
  visc_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  visc_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  visc_numerics->SetHIndex      ( nodes->GetHIndex()       );
  visc_numerics->SetAIndex      ( nodes->GetAIndex()       );
  visc_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  visc_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  unsigned short T_INDEX       = nodes->GetTIndex();
  unsigned short TVE_INDEX     = nodes->GetTveIndex();
  unsigned short VEL_INDEX     = nodes->GetVelIndex();
  unsigned short PRESS_INDEX   = nodes->GetPIndex();
  unsigned short RHO_INDEX     = nodes->GetRhoIndex();
  unsigned short H_INDEX       = nodes->GetHIndex();
  unsigned short A_INDEX       = nodes->GetAIndex();
  unsigned short RHOCVTR_INDEX = nodes->GetRhoCvtrIndex();
  unsigned short RHOCVVE_INDEX = nodes->GetRhoCvveIndex();

  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));

  Ms   = config->GetMolar_Mass();
  xi   = config->GetRotationModes();
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  nElStates = config->GetnElStates();

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else { nHeavy = nSpecies; nEl = 0; }

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
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
        UnitaryNormal[iDim] = Normal[iDim]/Area;

      /*--- Current solution at this boundary node ---*/
      for (iVar = 0; iVar < nVar; iVar++)     U_domain[iVar] = nodes->GetSolution(iPoint, iVar);
      for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = nodes->GetPrimitive(iPoint, iVar);

      /*--- Initialize solution at outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++)     U_outlet[iVar] = 0.0;
      for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = 0.0;

      /*--- Build the fictitious intlet state based on characteristics ---*/

      /*--- Retrieve the specified back pressure for this outlet. ---*/
      if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;
      else P_Exit = config->GetOutlet_Pressure(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
       of boundary update depends on this. ---*/
      Density = V_domain[RHO_INDEX];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[VEL_INDEX+iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitaryNormal[iDim];
      }
      Energy      = U_domain[nVar-2]/Density;
      Temperature = V_domain[T_INDEX];
      Tve         = V_domain[TVE_INDEX];
      Pressure    = V_domain[PRESS_INDEX];
      SoundSpeed  = V_domain[A_INDEX];
      Mach_Exit   = sqrt(Velocity2)/SoundSpeed;

      /*--- Compute Species Concentrations ---*/
      //Using partial pressures, maybe not
      for (iSpecies =0; iSpecies<nSpecies;iSpecies++){
        Ys[iSpecies] = V_domain[iSpecies]/Density;
      }

      /*--- Recompute boundary state depending Mach number ---*/
      if (Mach_Exit >= 1.0) {

        /*--- Supersonic exit flow: there are no incoming characteristics,
         so no boundary condition is necessary. Set outlet state to current
         state so that upwinding handles the direction of propagation. ---*/
        for (iVar = 0; iVar < nVar; iVar++) U_outlet[iVar] = U_domain[iVar];
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
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitaryNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }

        /*--- Primitive variables, using the derived quantities ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
          V_outlet[iSpecies]= Ys[iSpecies]*Density;
        }

        V_outlet[T_INDEX]     = V_domain[T_INDEX];
        V_outlet[TVE_INDEX]   = V_domain[TVE_INDEX];

        for (iDim = 0; iDim < nDim; iDim++){
          V_outlet[VEL_INDEX+iDim] = Velocity[iDim];
        }

        V_outlet[PRESS_INDEX] = Pressure;
        V_outlet[RHO_INDEX]   = Density;
        V_outlet[A_INDEX]     = SoundSpeed;

        for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
          V_outlet[RHOCVTR_INDEX ] += Density* Ys[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
        }

        for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {

          /*--- Vibrational energy ---*/
          if (thetav[iSpecies] != 0.0) {

            /*--- Rename for convenience ---*/
            thoTve = thetav[iSpecies]/Tve;
            exptv = exp(thetav[iSpecies]/Tve);
            thsqr = thetav[iSpecies]*thetav[iSpecies];

            /*--- Calculate species vibrational specific heats ---*/
            Cvvs  = Ru/Ms[iSpecies] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));

            /*--- Add contribution ---*/
            V_outlet[RHOCVVE_INDEX] += V_outlet[iSpecies] * Cvvs;
          }

          /*--- Electronic energy ---*/
          if (nElStates[iSpecies] != 0) {
            num = 0.0; num2 = 0.0;
            denom = g[iSpecies][0] * exp(-thetae[iSpecies][0]/Tve);
            num3  = g[iSpecies][0] * (thetae[iSpecies][0]/(Tve*Tve))*exp(-thetae[iSpecies][0]/Tve);
            for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
              thoTve = thetae[iSpecies][iEl]/Tve;
              exptv = exp(-thetae[iSpecies][iEl]/Tve);

              num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exptv;
              denom += g[iSpecies][iEl] * exptv;
              num2  += g[iSpecies][iEl] * (thoTve*thoTve) * exptv;
              num3  += g[iSpecies][iEl] * thoTve/Tve * exptv;
            }
            Cves = Ru/Ms[iSpecies] * (num2/denom - num*num3/(denom*denom));
            V_outlet[RHOCVVE_INDEX] += V_outlet[iSpecies] * Cves;
          }

        }

        for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
          cout << "THIS MAY BE WRONG" << endl;
          Cves = 3.0/2.0 * Ru/Ms[nSpecies-1];
          V_outlet[RHOCVVE_INDEX] += V_outlet[nSpecies-1] * Cves;
        }

        /*--- Conservative variables, using the derived quantities ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++){
          U_outlet[iSpecies] = V_outlet[iSpecies];
        }

        for (iDim = 0; iDim < nDim; iDim++)
          U_outlet[nSpecies+iDim] = Velocity[iDim]*Density;

        /*--- Calculate energy (RRHO and Electronic) ---*/
        for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {

          // Species formation energy
          Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

          // Species vibrational energy
          if (thetav[iSpecies] != 0.0)
            Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
          else
            Ev = 0.0;

          // Species electronic energy
          num = 0.0;
          denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
          for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
            num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
            denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
          }
          Ee = Ru/Ms[iSpecies] * (num/denom);

          // Mixture total Energy
          U_outlet[nVar-2] += U_outlet[iSpecies] * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies]*
                  (V_outlet[T_INDEX]-Tref[iSpecies]) + Ev + Ee + Ef +
                  0.5*Velocity2);

          // Mixture vibrational-electronic energy
          U_outlet[nVar-1] += U_outlet[iSpecies] * (Ev + Ee);


        }

        for (iSpecies = 0; iSpecies < nEl; iSpecies++) {

          // Species formation energy
          Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];

          // Electron t-r mode contributes to mixture vib-el energy
          U_outlet[nVar-1] += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
        }

      }

      /*--- Setting Last remaining variables ---*/
      V_outlet[H_INDEX]= (U_outlet[nVar-2]+Pressure)/Density;

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetConservative(U_domain, U_outlet);
      conv_numerics->SetPrimitive(V_domain,V_outlet);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Passing supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU(nodes->GetdPdU(iPoint), node_infty->GetdPdU(0));
      conv_numerics->SetdTdU(nodes->GetdTdU(iPoint), node_infty->GetdTdU(0));
      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
//      if (viscous) {

//        /*--- Set the normal vector and the coordinates ---*/
//        visc_numerics->SetNormal(Normal);
//        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

//        /*--- Primitive variables, and gradient ---*/
//        visc_numerics->SetPrimitive(V_domain, V_outlet);
//        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), nodes->GetGradient_Primitive());

//        /*--- Conservative variables, and gradient ---*/
//        visc_numerics->SetConservative(U_domain, U_outlet);
//        visc_numerics->SetConsVarGradient(nodes->GetGradient(), node_infty->GetGradient() );


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
//          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//      }
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] U_domain;
  delete [] U_outlet;
  delete [] V_domain;
  delete [] V_outlet;
  delete [] Normal;
  delete [] Ys;

}

void CTNE2EulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solution_container,
                                           CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, iSpecies, iEl;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Density, Pressure, Temperature, Temperature_ve, Energy, *Velocity, *Mass_Frac, Velocity2, soundspeed;
  su2double Gas_Constant = config->GetGas_ConstantND();

  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous              = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double *Ms, *xi, *Tref, *hf, *thetav, **thetae, **g, rhos, Ef, Ee, Ev;
  unsigned short nEl, nHeavy, *nElStates;
  Tref                    = config->GetRefTemperature();
  Ms                      = config->GetMolar_Mass();
  xi                      = config->GetRotationModes();
  nElStates               = config->GetnElStates();
  hf                      = config->GetEnthalpy_Formation();
  thetav                  = config->GetCharVibTemp();
  thetae                  = config->GetCharElTemp();
  g                       = config->GetElDegeneracy();
  cout << "This doesnt work" << endl;
  su2double denom = 0.0, conc   = 0.0, rhoCvtr = 0.0, num = 0.0,
      rhoE  = 0.0, rhoEve = 0.0;
  su2double RuSI  = UNIVERSAL_GAS_CONSTANT;
  su2double Ru = 1000.0*RuSI;

  su2double *U_inlet = new su2double[nVar];     su2double *U_domain = new su2double[nVar];
  su2double *V_inlet = new su2double[nPrimVar]; su2double *V_domain = new su2double[nPrimVar];
  su2double *Normal = new su2double[nDim];

  bool ionization = config->GetIonization();

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else { nHeavy = nSpecies; nEl = 0; }

  /*--- Supersonic inlet flow: there are no outgoing characteristics,
   so all flow variables can be imposed at the inlet.
   First, retrieve the specified values for the primitive variables. ---*/
  //ASSUME TVE = T for the time being
  Mass_Frac      = config->GetInlet_MassFrac(Marker_Tag);
  Temperature    = config->GetInlet_Temperature(Marker_Tag);
  Pressure       = config->GetInlet_Pressure(Marker_Tag);
  Velocity       = config->GetInlet_Velocity(Marker_Tag);
  Temperature_ve = Temperature;

  /*--- Compute Density and Species Densities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    denom += Mass_Frac[iSpecies] * (Ru/Ms[iSpecies]) * Temperature;
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    denom += Mass_Frac[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Temperature_ve;
  Density = Pressure / denom;

  /*--- Compute Soundspeed and Velocity squared ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    conc += Mass_Frac[iSpecies]*Density/Ms[iSpecies];
    rhoCvtr += Density*Mass_Frac[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
  }
  soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * Pressure/Density);

  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity2 += Velocity[iDim]*Velocity[iDim];

  /*--- Non-dim. the inputs if necessary. ---*/
  // Need to update this portion
  //Temperature = Temperature/config->GetTemperature_Ref();
  //Pressure    = Pressure/config->GetPressure_Ref();
  //Density     = Density/config->GetDensity_Ref();
  //for (iDim = 0; iDim < nDim; iDim++)
  //  Velocity[iDim] = Velocity[iDim]/config->GetVelocity_Ref();

  /*--- Compute energy (RRHO) from supplied primitive quanitites ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {

    // Species density
    rhos = Mass_Frac[iSpecies]*Density;

    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

    // Species vibrational energy
    if (thetav[iSpecies] != 0.0)
      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Temperature_ve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Temperature_ve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Temperature_ve);
      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Temperature_ve);
    }
    Ee = Ru/Ms[iSpecies] * (num/denom);

    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (Temperature-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*Velocity2);

    // Mixture vibrational-electronic energy
    rhoEve += rhos * (Ev + Ee);
  }

  /*--- Setting Conservative Variables ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    U_inlet[iSpecies] = Mass_Frac[iSpecies]*Density;
  for (iDim = 0; iDim < nDim; iDim++)
    U_inlet[nSpecies+iDim] = Density*Velocity[iDim];
  U_inlet[nVar-2] = rhoE;
  U_inlet[nVar-1] = rhoEve;

  /*--- Setting Primitive Vaariables ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    V_inlet[iSpecies] = Mass_Frac[iSpecies]*Density;
  V_inlet[nSpecies] = Temperature;
  V_inlet[nSpecies+1] = Temperature_ve;
  for (iDim = 0; iDim < nDim; iDim++)
    V_inlet[nSpecies+2+iDim] = Velocity[iDim];
  V_inlet[nSpecies+2+nDim] = Pressure;
  V_inlet[nSpecies+3+nDim] = Density;
  V_inlet[nSpecies+4+nDim] = rhoE+Pressure/Density;
  V_inlet[nSpecies+5+nDim] = soundspeed;
  V_inlet[nSpecies+6+nDim] = rhoCvtr;

  //This requires Newtown Raphson.....So this is not currently operational (See Deathstar)
  //V_inlet[nSpecies+7+nDim] = rhoCvve;

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Current solution at this boundary node ---*/
      for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = nodes->GetSolution(iPoint,iVar);
      for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      su2double Area = 0.0;su2double UnitaryNormal[3];
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitaryNormal[iDim] = Normal[iDim]/Area;

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);
      conv_numerics->SetConservative(U_domain, U_inlet);
      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Pass supplementary info to CNumerics ---*/
      conv_numerics->SetdPdU(nodes->GetdPdU(iPoint), node_infty->GetdPdU(0));
      conv_numerics->SetdTdU(nodes->GetdTdU(iPoint), node_infty->GetdTdU(0));
      conv_numerics->SetdTvedU(nodes->GetdTvedU(iPoint), node_infty->GetdTvedU(0));

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

        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint), nodes->GetLaminarViscosity(iPoint));

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
  delete [] U_domain;
  delete [] U_inlet;
  delete [] V_domain;
  delete [] V_inlet;
  delete [] Normal;

}

void CTNE2EulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solution_container,
                                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_outlet, *V_domain;
  su2double *U_outlet, *U_domain;

  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double *Normal = new su2double[nDim];

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  conv_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  conv_numerics->SetPIndex      ( nodes->GetPIndex()       );
  conv_numerics->SetTIndex      ( nodes->GetTIndex()       );
  conv_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  conv_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  conv_numerics->SetHIndex      ( nodes->GetHIndex()       );
  conv_numerics->SetAIndex      ( nodes->GetAIndex()       );
  conv_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  conv_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  /*--- Supersonic outlet flow: there are no ingoing characteristics,
   so all flow variables can should be interpolated from the domain. ---*/

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Current solution at this boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);
      U_domain = nodes->GetSolution(iPoint);

      /*--- Allocate the value at the outlet ---*/
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
      conv_numerics->SetdPdU(   nodes->GetdPdU(iPoint),   nodes->GetdPdU(iPoint));
      conv_numerics->SetdTdU(   nodes->GetdTdU(iPoint),   nodes->GetdTdU(iPoint));
      conv_numerics->SetdTvedU( nodes->GetdTvedU(iPoint), nodes->GetdTvedU(iPoint));

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}

void CTNE2EulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Call the Euler wall routine ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTNE2EulerSolver::SetResidual_DualTime(CGeometry *geometry,
                                            CSolver **solution_container,
                                            CConfig *config,
                                            unsigned short iRKStep,
                                            unsigned short iMesh,
                                            unsigned short RunTime_EqSystem) {
  unsigned short iVar, jVar;
  unsigned long iPoint;
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_n, Volume_nP1, TimeStep;

  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool Grid_Movement = config->GetGrid_Movement();

  /*--- loop over points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Solution at time n-1, n and n+1 ---*/
    U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
    U_time_n   = nodes->GetSolution_time_n(iPoint);
    U_time_nP1 = nodes->GetSolution(iPoint);

    /*--- Volume at time n-1 and n ---*/
    if (Grid_Movement) {
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_n = geometry->node[iPoint]->GetVolume_n();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
    }
    else {
      Volume_nM1 = geometry->node[iPoint]->GetVolume();
      Volume_n = geometry->node[iPoint]->GetVolume();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
    }

    /*--- Time Step ---*/
    TimeStep = config->GetDelta_UnstTimeND();

    /*--- Compute Residual ---*/
    for(iVar = 0; iVar < nVar; iVar++) {
      if (config->GetTime_Marching() == DT_STEPPING_1ST)
        Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
      if (config->GetTime_Marching() == DT_STEPPING_2ND)
        Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
                           +  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
    }

    /*--- Add Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual);

    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;

        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
      }
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }

}

//void CTNE2EulerSolver::GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone) {

//  /*--- Restart the solution from file information ---*/
//  string restart_filename = config->GetSolution_FlowFileName();
//  unsigned long iPoint, index, nFlowIter, adjIter, flowIter;
//  char buffer[50];
//  string UnstExt, text_line;
//  ifstream restart_file;
//  bool grid_movement = config->GetGrid_Movement();
//  unsigned short nZone = geometry->GetnZone();
//  unsigned short iVar;

//  /*--- Multi-zone restart files. ---*/
//  if (nZone > 1 && !(config->GetTime_Marching() == DT_STEPPING_1ST)) {
//    restart_filename.erase(restart_filename.end()-4, restart_filename.end());
//    sprintf (buffer, "_%d.dat", int(val_iZone));
//    UnstExt = string(buffer);
//    restart_filename.append(UnstExt);
//  }

//  /*--- For the unsteady adjoint, we integrate backwards through
//   physical time, so load in the direct solution files in reverse. ---*/
//  if (config->GetTime_Marching() == DT_STEPPING_1ST) {
//    flowIter = val_iZone;
//    restart_filename.erase(restart_filename.end()-4, restart_filename.end());
//    if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.dat", int(val_iZone));
//    if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.dat", int(val_iZone));
//    if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.dat", int(val_iZone));
//    if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.dat", int(val_iZone));
//    if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.dat", int(val_iZone));
//    UnstExt = string(buffer);
//    restart_filename.append(UnstExt);
//  } else if (config->GetTime_Marching() && config->GetWrt_Unsteady()) {
//    nFlowIter = config->GetnExtIter() - 1;
//    adjIter   = config->GetExtIter();
//    flowIter  = nFlowIter - adjIter;
//    restart_filename.erase (restart_filename.end()-4, restart_filename.end());
//    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
//    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
//    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
//    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
//    if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
//    UnstExt = string(buffer);
//    restart_filename.append(UnstExt);
//  } else {
//    flowIter = config->GetExtIter();
//    restart_filename.erase (restart_filename.end()-4, restart_filename.end());
//    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
//    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
//    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
//    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
//    if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
//    UnstExt = string(buffer);
//    restart_filename.append(UnstExt);
//  }

//  /*--- Open the flow solution from the restart file ---*/
//  if (rank == MASTER_NODE && val_iZone == ZONE_0)
//    cout << "Reading in the direct flow solution from iteration " << flowIter << "." << endl;
//  restart_file.open(restart_filename.data(), ios::in);

//  /*--- In case there is no file ---*/
//  if (restart_file.fail()) {
//    cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
//    exit(1);
//  }

//  /*--- In case this is a parallel simulation, we need to perform the
//   Global2Local index transformation first. ---*/
//  long *Global2Local = NULL;
//  Global2Local = new long[geometry->GetGlobal_nPointDomain()];
//  /*--- First, set all indices to a negative value by default ---*/
//  for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
//    Global2Local[iPoint] = -1;
//  }

//  /*--- Now fill array with the transform values only for local points ---*/
//  for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
//    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
//  }

//  /*--- Read all lines in the restart file ---*/
//  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

//  /*--- The first line is the header ---*/
//  getline (restart_file, text_line);

//  while (getline (restart_file,text_line)) {
//    istringstream point_line(text_line);

//    /*--- Retrieve local index. If this node from the restart file lives
//     on a different processor, the value of iPoint_Local will be -1, as
//     initialized above. Otherwise, the local index for this node on the
//     current processor will be returned and used to instantiate the vars. ---*/
//    iPoint_Local = Global2Local[iPoint_Global];
//    if (iPoint_Local >= 0) {

//      /*--- First value is the point index, then the conservative variables ---*/
//      point_line >> index;

//      for (iVar = 0; iVar < nVar; iVar++)
//        point_line >> Solution[iVar];

//      node[iPoint_Local]->SetSolution(iPoint,Solution);

//      /*--- If necessary, read in the grid velocities for the unsteady adjoint ---*/
//      if (config->GetTime_Marching() && config->GetWrt_Unsteady() && grid_movement) {
//        su2double Volume, GridVel[3];
//        if (nDim == 2) point_line >> Volume >> GridVel[0] >> GridVel[1];
//        if (nDim == 3) point_line >> Volume >> GridVel[0] >> GridVel[1] >> GridVel[2];
//        if (iPoint_Local >= 0)
//          for (unsigned short iDim = 0; iDim < nDim; iDim++)
//            geometry->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
//      }

//    }
//    iPoint_Global++;
//  }


//  /*--- Close the restart file ---*/
//  restart_file.close();

//  /*--- Free memory needed for the transformation ---*/
//  delete [] Global2Local;

//}

void CTNE2EulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {


  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool static_fsi = ((config->GetTime_Marching() == STEADY) &&
                     (config->GetFSI_Simulation()));
  bool steady_restart = config->GetSteadyRestart();
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  bool turbulent     = (config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS);

  string UnstExt, text_line;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string  restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

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
    if (turb_model == SST) turbVars = 2;
    else turbVars = 1;
  }

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
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);
      iPoint_Global_Local++;

      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/
      if (grid_movement && val_update_geo) {

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/
        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar + turbVars;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[MESH_0]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
      }

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
  solver[MESH_0][TNE2_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][TNE2_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][TNE2_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_TNE2_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][TNE2_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][TNE2_SOL]->GetNodes()->SetSolution(iPoint,Solution);
    }
    solver[MESH_0][TNE2_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
    solver[MESH_0][TNE2_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);
    solver[iMesh][TNE2_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_TNE2_SYS, false);
  }

  /*--- Update the geometry for flows on dynamic meshes ---*/
  if (grid_movement && val_update_geo) {

    /*--- Communicate the new coordinates and grid velocities at the halos ---*/

    geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
    geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);

    geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, GRID_VELOCITY);
    geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, GRID_VELOCITY);

    /*--- Recompute the edges and dual mesh control volumes in the
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
      geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine], config);
      geometry[iMesh]->SetMaxLength(config);
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
  if (dual_time && grid_movement)
    Restart_OldGeometry(geometry[MESH_0], config);

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

void CTNE2EulerSolver::SetVolume_Output(CConfig *config, CGeometry *geometry, su2double **data_container, unsigned short nOutput_Vars) {

#ifdef DEBUG_TDE

  unsigned short iVar;
  unsigned long iPoint;

  /*--- Add up total number of output variables to be written. ---*/
  nOutput_Vars = nVar;

  for (iVar = 0; iVar < config->GetnOutput_Vars_Vol(); iVar++ ) {

    switch(config->GetOutput_Vars_Vol(iVar)) {
    case PRESSURE:
      nOutput_Vars++;
      break;
    case MACH:
      nOutput_Vars++;
      break;
    }

  }

  // NEEDS TO BE MAX NUMBER OF POINTS ON ANY PARTITION ?
  data_container = new su2double*[nOutput_Vars];
  for (iVar = 0; iVar < nOutput_Vars; iVar++ ) {
    data_container[iVar] = new su2double[nPointDomain];
  }

  for (iVar = 0; iVar < config->GetnOutput_Vars_Vol(); iVar++ ) {

    switch(config->GetOutput_Vars_Vol(iVar)) {
    case PRESSURE:
      nOutput_Vars++;
      break;
    case MACH:
      nOutput_Vars++;
      break;
    }

  }
#endif
}

void CTNE2EulerSolver::ResetNodeInfty(su2double pressure_inf, su2double *massfrac_inf, su2double *mvec_inf, su2double temperature_inf,
                                        su2double temperature_ve_inf, CConfig *config){
  su2double check_infty;
  if (node_infty != nullptr) delete node_infty;

  node_infty = new CTNE2EulerVariable(pressure_inf, massfrac_inf, mvec_inf, temperature_inf,
                                      temperature_ve_inf, 1, nDim, nVar,
                                      nPrimVar, nPrimVarGrad, config);

  check_infty = node_infty->SetPrimVar_Compressible(0,config);
}

CTNE2NSSolver::CTNE2NSSolver(void) : CTNE2EulerSolver() {

  /*--- Array initialization ---*/
  CD_Visc   = NULL; CL_Visc   = NULL; CSF_Visc  = NULL; CEff_Visc = NULL;
  CFx_Visc  = NULL; CFy_Visc  = NULL; CFz_Visc  = NULL;
  CMx_Visc  = NULL; CMy_Visc  = NULL; CMz_Visc  = NULL;
  CoPx_Visc = NULL; CoPy_Visc = NULL; CoPz_Visc = NULL;

  ForceViscous = NULL; MomentViscous = NULL; CSkinFriction = NULL;

  Buffet_Sensor = NULL; Buffet_Metric = NULL;

  /*--- Surface based array initialization ---*/
  Surface_CL_Visc  = NULL; Surface_CD_Visc    = NULL; Surface_CSF_Visc      = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL; Surface_CFy_Visc   = NULL; Surface_CFz_Visc      = NULL;
  Surface_CMx_Visc = NULL; Surface_CMy_Visc   = NULL; Surface_CMz_Visc      = NULL;
  Surface_HF_Visc  = NULL; Surface_MaxHF_Visc = NULL; Surface_Buffet_Metric = NULL;

  /*--- Rotorcraft simulation array initialization ---*/
  CMerit_Visc = NULL; CT_Visc = NULL; CQ_Visc = NULL;

  /*--- Inlet Variables ---*/
  Inlet_Ttotal = NULL;
  Inlet_Ptotal = NULL;
  Inlet_FlowDir = NULL;

}

CTNE2NSSolver::CTNE2NSSolver(CGeometry *geometry, CConfig *config,
                             unsigned short iMesh) : CTNE2EulerSolver() {

  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iSpecies, iMarker, nLineLets;
  su2double Density, Velocity2, Pressure, Temperature, StaticEnergy;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart    = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;

  bool low_mach_prec = config->Low_Mach_Preconditioning();

  bool adjoint = (config->GetDiscrete_Adjoint());
  string filename_ = "flow";

  unsigned short direct_diff = config->GetDirectDiff();
  bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_TNE2_RANS));
  bool multizone = config->GetMultizone_Problem();

  bool check_infty, check;
  su2double *Mvec_Inf, Alpha, Beta, dull_val;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
     before computing all the non-dimesional quantities. ---*/
  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Modify file name for a dual-time unsteady restart ---*/
    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/
    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
    }

    filename_ = config->GetFilename(filename_, ".meta", Unst_RestartIter);

    /*--- Read and store the restart metadata. ---*/
    Read_SU2_Restart_Metadata(geometry, config, false, filename_);

  }

  /*--- Array initialization ---*/
  CD_Visc   = NULL; CL_Visc   = NULL; CSF_Visc  = NULL; CEff_Visc = NULL;
  CFx_Visc  = NULL; CFy_Visc  = NULL; CFz_Visc  = NULL;
  CMx_Visc  = NULL; CMy_Visc  = NULL; CMz_Visc  = NULL;
  CoPx_Visc = NULL; CoPy_Visc = NULL; CoPz_Visc = NULL;

  Buffet_Sensor = NULL; Buffet_Metric = NULL;

  Surface_CL_Visc  = NULL; Surface_CD_Visc    = NULL; Surface_CSF_Visc = NULL; Surface_CEff_Visc = NULL;
  Surface_CFx_Visc = NULL; Surface_CFy_Visc   = NULL; Surface_CFz_Visc = NULL;
  Surface_CMx_Visc = NULL; Surface_CMy_Visc   = NULL; Surface_CMz_Visc = NULL;
  Surface_HF_Visc  = NULL; Surface_MaxHF_Visc = NULL;

  Surface_Buffet_Metric = NULL;

  CMerit_Visc      = NULL; CT_Visc      = NULL; CQ_Visc       = NULL;
  MaxHF_Visc       = NULL; ForceViscous = NULL; MomentViscous = NULL;
  CSkinFriction    = NULL; Cauchy_Serie = NULL; HF_Visc       = NULL;
  HeatConjugateVar = NULL;

  /*--- Initialize quantities for the average process for internal flow ---*/
  AverageVelocity 		    = NULL;
  AverageTurboVelocity 	    = NULL;
  OldAverageTurboVelocity   = NULL;
  ExtAverageTurboVelocity   = NULL;
  AverageFlux 			    = NULL;
  SpanTotalFlux 		    = NULL;
  AveragePressure  			= NULL;
  OldAveragePressure        = NULL;
  RadialEquilibriumPressure = NULL;
  ExtAveragePressure  		= NULL;
  AverageDensity   			= NULL;
  OldAverageDensity         = NULL;
  ExtAverageDensity   		= NULL;
  AverageNu                 = NULL;
  AverageKine               = NULL;
  AverageOmega              = NULL;
  ExtAverageNu              = NULL;
  ExtAverageKine            = NULL;
  ExtAverageOmega           = NULL;

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure ---*/
  nSpecies     = config->GetnSpecies();
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim         = geometry->GetnDim();

  /*--- Set the size of the primitive and conserve vectors ---*/
  //     U: [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  //     V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradV: [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  nVar         = nSpecies+nDim+2;
  nPrimVar     = nSpecies+nDim+8;
  nPrimVarGrad = nSpecies+nDim+8;

  /*--- Initialize nVarGrad for deallocation ---*/
  nVarGrad     = nPrimVarGrad;

  /*--- Store the number of vertices on each marker for deallocation later ---*/
  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /*--- Perform the non-dimensionalization for the flow equations using the
    specified reference values. ---*/
  SetNondimensionalization(geometry, config, iMesh);

  /*--- Define auxiliary vectors to store residual-related quantities ---*/
  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]     = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]   = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]   = 0.0;
  Res_Conv     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]     = 0.0;
  Res_Visc     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]     = 0.0;
  Res_Sour     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]     = 0.0;

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

  /*--- Allocate arrays for conserved variable limits ---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    lowerlimit[iSpecies] = 0.0;
    upperlimit[iSpecies] = 1E16;
  }
  for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
    lowerlimit[iVar] = -1E16;
    upperlimit[iVar] = 1E16;
  }
  for (iVar = nSpecies+nDim; iVar < nSpecies+nDim+2; iVar++) {
    lowerlimit[iVar] = 1E-4;
    upperlimit[iVar] = 1E16;
  }

  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
  if (config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }

  /*--- Initialize the solution & residual CVectors ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Create the structure for storing extra information ---*/
  if (config->GetExtraOutput()) {
    nOutputVariables = nVar;
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
  }

  /*--- Allocate Jacobians for implicit time-stepping ---*/
  if (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT) {
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the global Jacobian ---*/
    if (rank == MASTER_NODE) cout << "Initialize jacobian structure (TNE2 Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No jacobian structure (TNE2 Navier-Stokes). MG level: "
           << iMesh <<"." << endl;
  }

  /*--- Computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];

    /*--- c vector := transpose(WA)*(Wb) ---*/
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
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar+2];
        for (iVar = 0; iVar < nPrimVar + 2 ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
      }
      else {
        DonorPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
        for (iVar = 0; iVar < nPrimVar ; iVar++) {
          DonorPrimVar[iMarker][iVertex][iVar] = 0.0;
        }
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

  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      DonorGlobalIndex[iMarker][iVertex] = 0;
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

  /*--- Buffet sensor in all the markers ---*/
  if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){
    Buffet_Sensor = new su2double*[nMarker];
    for(iMarker = 0; iMarker < nMarker; iMarker++) {
      Buffet_Sensor[iMarker] = new su2double[geometry->nVertex[iMarker]];
    }
  }

  /*--- Non dimensional coefficients ---*/
  ForceInviscid  = new su2double[3];
  MomentInviscid = new su2double[3];
  CD_Inv         = new su2double[nMarker];
  CL_Inv         = new su2double[nMarker];
  CSF_Inv        = new su2double[nMarker];
  CEff_Inv       = new su2double[nMarker];
  CFx_Inv        = new su2double[nMarker];
  CFy_Inv        = new su2double[nMarker];
  CFz_Inv        = new su2double[nMarker];
  CMx_Inv        = new su2double[nMarker];
  CMy_Inv        = new su2double[nMarker];
  CMz_Inv        = new su2double[nMarker];
  CoPx_Inv       = new su2double[nMarker];
  CoPy_Inv       = new su2double[nMarker];
  CoPz_Inv       = new su2double[nMarker];

  ForceMomentum  = new su2double[3];
  MomentMomentum = new su2double[3];
  CD_Mnt         = new su2double[nMarker];
  CL_Mnt         = new su2double[nMarker];
  CSF_Mnt        = new su2double[nMarker];
  CEff_Mnt       = new su2double[nMarker];
  CFx_Mnt        = new su2double[nMarker];
  CFy_Mnt        = new su2double[nMarker];
  CFz_Mnt        = new su2double[nMarker];
  CMx_Mnt        = new su2double[nMarker];
  CMy_Mnt        = new su2double[nMarker];
  CMz_Mnt        = new su2double[nMarker];
  CoPx_Mnt       = new su2double[nMarker];
  CoPy_Mnt       = new su2double[nMarker];
  CoPz_Mnt       = new su2double[nMarker];

  ForceViscous   = new su2double[3];
  MomentViscous  = new su2double[3];
  CD_Visc        = new su2double[nMarker];
  CL_Visc        = new su2double[nMarker];
  CSF_Visc       = new su2double[nMarker];
  CEff_Visc      = new su2double[nMarker];
  CFx_Visc       = new su2double[nMarker];
  CFy_Visc       = new su2double[nMarker];
  CFz_Visc       = new su2double[nMarker];
  CMx_Visc       = new su2double[nMarker];
  CMy_Visc       = new su2double[nMarker];
  CMz_Visc       = new su2double[nMarker];
  CoPx_Visc      = new su2double[nMarker];
  CoPy_Visc      = new su2double[nMarker];
  CoPz_Visc      = new su2double[nMarker];

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

  Surface_CL       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD       = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff     = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy      = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz      = new su2double[config->GetnMarker_Monitoring()];

  Surface_CL_Visc    = new su2double[config->GetnMarker_Monitoring()];
  Surface_CD_Visc    = new su2double[config->GetnMarker_Monitoring()];
  Surface_CSF_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CEff_Visc  = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFx_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFy_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CFz_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc   = new su2double[config->GetnMarker_Monitoring()];
  Surface_HF_Visc    = new su2double[config->GetnMarker_Monitoring()];
  Surface_MaxHF_Visc = new su2double[config->GetnMarker_Monitoring()];

  if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){
    Buffet_Metric          = new su2double[nMarker];
    Surface_Buffet_Metric = new su2double[config->GetnMarker_Monitoring()];
  }

  /*--- Rotational coefficients ---*/
  CMerit_Inv = new su2double[nMarker];
  CT_Inv     = new su2double[nMarker];
  CQ_Inv     = new su2double[nMarker];

  CMerit_Mnt = new su2double[nMarker];
  CT_Mnt     = new su2double[nMarker];
  CQ_Mnt     = new su2double[nMarker];

  CMerit_Visc = new su2double[nMarker];
  CT_Visc     = new su2double[nMarker];
  CQ_Visc     = new su2double[nMarker];

  /*--- Heat based coefficients ---*/
  HF_Visc    = new su2double[nMarker];
  MaxHF_Visc = new su2double[nMarker];

  /*--- Supersonic coefficients ---*/
  CEquivArea_Inv   = new su2double[nMarker];
  CNearFieldOF_Inv = new su2double[nMarker];

  /*--- Init total coefficients ---*/
  Total_CD        = 0.0; Total_CL           = 0.0; Total_CSF          = 0.0;
  Total_CMx       = 0.0; Total_CMy          = 0.0; Total_CMz          = 0.0;
  Total_CoPx      = 0.0; Total_CoPy         = 0.0; Total_CoPz          = 0.0;
  Total_CEff      = 0.0; Total_CEquivArea   = 0.0; Total_CNearFieldOF = 0.0;
  Total_CFx       = 0.0; Total_CFy          = 0.0; Total_CFz          = 0.0;
  Total_CT        = 0.0; Total_CQ           = 0.0; Total_CMerit       = 0.0;
  Total_MaxHeat   = 0.0; Total_Heat         = 0.0; Total_ComboObj     = 0.0;
  Total_CpDiff    = 0.0; Total_HeatFluxDiff = 0.0;
  Total_NetThrust = 0.0; Total_CL_Prev      = 0.0;
  Total_Power     = 0.0; AoA_Prev           = 0.0; Total_CD_Prev      = 0.0;
  Total_CMx_Prev  = 0.0; Total_CMy_Prev     = 0.0; Total_CMz_Prev     = 0.0;
  Total_AeroCD    = 0.0; Total_SolidCD      = 0.0; Total_IDR          = 0.0;
  Total_IDC       = 0.0;
  Total_Custom_ObjFunc = 0.0;

  /*--- Read farfield conditions from config ---*/
  Density_Inf        = config->GetDensity_FreeStreamND();
  Pressure_Inf       = config->GetPressure_FreeStreamND();
  Temperature_Inf    = config->GetTemperature_FreeStreamND();
  Temperature_ve_Inf = config->GetTemperature_ve_FreeStreamND();
  MassFrac_Inf       = config->GetMassFrac_FreeStream();
  Mach_Inf           = config->GetMach();
  Viscosity_Inf      = config->GetViscosity_FreeStreamND();
  Prandtl_Lam        = config->GetPrandtl_Lam();
  Prandtl_Turb       = config->GetPrandtl_Turb();

  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  switch(direct_diff) {
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


  /*--- Create a CVariable that stores the free-stream values ---*/
  node_infty = new CTNE2NSVariable(Pressure_Inf, MassFrac_Inf,
                                   Mvec_Inf, Temperature_Inf,
                                   Temperature_ve_Inf, 1, nDim, nVar,
                                   nPrimVar, nPrimVarGrad, config);
  check_infty = node_infty->SetPrimVar_Compressible(0, config);

  Velocity_Inf = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_Inf[iDim] = node_infty->GetVelocity(iPoint, iDim);

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  nodes = new CTNE2NSVariable(Pressure_Inf, MassFrac_Inf,
                                       Mvec_Inf, Temperature_Inf,
                                       Temperature_ve_Inf, nPoint, nDim, nVar,
                                       nPrimVar, nPrimVarGrad, config);

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    check = nodes->SetPrimVar_Compressible(iPoint, config);

    if (check) {
      bool ionization;
      unsigned short iEl, nHeavy, nEl, *nElStates;
      su2double RuSI, Ru, T, Tve, rhoCvtr, sqvel, rhoE, rhoEve, num, denom, conc;
      su2double rho, rhos, Ef, Ev, Ee, soundspeed;
      su2double *xi, *Ms, *thetav, **thetae, **g, *Tref, *hf;

      /*--- Determine the number of heavy species ---*/
      ionization = config->GetIonization();
      if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
      else            { nHeavy = nSpecies;   nEl = 0; }

      /*--- Load variables from the config class --*/
      xi        = config->GetRotationModes();      // Rotational modes of energy storage
      Ms        = config->GetMolar_Mass();         // Species molar mass
      thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
      thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
      g         = config->GetElDegeneracy();       // Degeneracy of electron states
      nElStates = config->GetnElStates();          // Number of electron states
      Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
      hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]

      /*--- Rename & initialize for convenience ---*/
      RuSI    = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(mol*K)]
      Ru      = 1000.0*RuSI;                    // Universal gas constant [J/(kmol*K)]
      Tve     = Temperature_ve_Inf;             // Vibrational temperature [K]
      T       = Temperature_Inf;                // Translational-rotational temperature [K]
      sqvel   = 0.0;                            // Velocity^2 [m2/s2]
      rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
      rhoEve  = 0.0;                            // Mixture vib-el energy per mass [J/kg]
      denom   = 0.0;
      conc    = 0.0;
      rhoCvtr = 0.0;

      /*--- Calculate mixture density from supplied primitive quantities ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
        denom += MassFrac_Inf[iSpecies] * (Ru/Ms[iSpecies]) * T;
      for (iSpecies = 0; iSpecies < nEl; iSpecies++)
        denom += MassFrac_Inf[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Tve;
      rho = Pressure_Inf / denom;

      /*--- Calculate sound speed and extract velocities ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
        conc += MassFrac_Inf[iSpecies]*rho/Ms[iSpecies];
        rhoCvtr += rho*MassFrac_Inf[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
      }
      soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * Pressure_Inf/rho);
      for (iDim = 0; iDim < nDim; iDim++){
        sqvel += Mvec_Inf[iDim]*soundspeed * Mvec_Inf[iDim]*soundspeed;
      }
      /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
      for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
        // Species density
        rhos = MassFrac_Inf[iSpecies]*rho;

        // Species formation energy
        Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

        // Species vibrational energy
        if (thetav[iSpecies] != 0.0)
          Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
        else
          Ev = 0.0;

        // Species electronic energy
        num = 0.0;
        denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
          denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
        }
        Ee = Ru/Ms[iSpecies] * (num/denom);

        // Mixture total energy
        rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                        + Ev + Ee + Ef + 0.5*sqvel);

        // Mixture vibrational-electronic energy
        rhoEve += rhos * (Ev + Ee);
      }
      for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
        // Species formation energy
        Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];

        // Electron t-r mode contributes to mixture vib-el energy
        rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
      }

      /*--- Initialize Solution & Solution_Old vectors ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Solution[iSpecies]      = rho*MassFrac_Inf[iSpecies];
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        Solution[nSpecies+iDim] = rho*Mvec_Inf[iDim]*soundspeed;
      }
      Solution[nSpecies+nDim]     = rhoE;
      Solution[nSpecies+nDim+1]   = rhoEve;

      nodes->SetSolution(iPoint,Solution);
      nodes->SetSolution_Old(iPoint,Solution);

      counter_local++;
    }
  }

  /*--- Warning message about non-physical points ---*/
  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }

  /*--- Define some structures for locating max residuals ---*/
  Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar]  = 0;
  Point_Max_Coord_BGS = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord_BGS[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
  }


  /*--- Define solver parameters needed for execution of destructor ---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Perform the MPI communication of the solution ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Deallocate arrays ---*/
  delete [] Mvec_Inf;
}

CTNE2NSSolver::~CTNE2NSSolver(void) {
  unsigned short iMarker, iDim;

  unsigned long iVertex;

  if (CD_Visc != NULL)          delete [] CD_Visc;
  if (CL_Visc != NULL)          delete [] CL_Visc;
  if (CSF_Visc != NULL)         delete [] CSF_Visc;
  if (CFx_Visc != NULL)         delete [] CFx_Visc;
  if (CFy_Visc != NULL)         delete [] CFy_Visc;
  if (CFz_Visc != NULL)         delete [] CFz_Visc;
  if (CMx_Visc != NULL)         delete [] CMx_Visc;
  if (CMy_Visc != NULL)         delete [] CMy_Visc;
  if (CMz_Visc != NULL)         delete [] CMz_Visc;
  if (CoPx_Visc != NULL)        delete [] CoPx_Visc;
  if (CoPy_Visc != NULL)        delete [] CoPy_Visc;
  if (CoPz_Visc != NULL)        delete [] CoPz_Visc;
  if (CEff_Visc != NULL)        delete [] CEff_Visc;
  if (CMerit_Visc != NULL)      delete [] CMerit_Visc;
  if (Buffet_Metric != NULL)    delete [] Buffet_Metric;
  if (CT_Visc != NULL)          delete [] CT_Visc;
  if (CQ_Visc != NULL)          delete [] CQ_Visc;
  if (HF_Visc != NULL)          delete [] HF_Visc;
  if (MaxHF_Visc != NULL)       delete [] MaxHF_Visc;
  if (ForceViscous != NULL)     delete [] ForceViscous;
  if (MomentViscous != NULL)    delete [] MomentViscous;

  if (Surface_CL_Visc != NULL)      delete [] Surface_CL_Visc;
  if (Surface_CD_Visc != NULL)      delete [] Surface_CD_Visc;
  if (Surface_CSF_Visc != NULL)     delete [] Surface_CSF_Visc;
  if (Surface_CEff_Visc != NULL)    delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != NULL)     delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != NULL)     delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != NULL)     delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != NULL)     delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)     delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)     delete [] Surface_CMz_Visc;
  if (Surface_HF_Visc != NULL)      delete [] Surface_HF_Visc;
  if (Surface_MaxHF_Visc != NULL)   delete [] Surface_MaxHF_Visc;
  if (Surface_Buffet_Metric != NULL) delete [] Surface_Buffet_Metric;

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

void CTNE2NSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh, unsigned short iRKStep,
                                  unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;

  unsigned long InnerIter = config->GetInnerIter();
  bool disc_adjoint       = config->GetDiscrete_Adjoint();
  bool implicit           = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool center             = (config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED);
  bool center_jst         = (center && config->GetKind_Centered_TNE2() == JST);
  bool limiter_flow       = ((config->GetKind_SlopeLimit_TNE2() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()));
  bool limiter_turb       = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()));
  bool van_albada         = (config->GetKind_SlopeLimit_TNE2() == VAN_ALBADA_EDGE);
  bool nonPhys;

  su2double *errU, *errV;
  su2double StrainMag = 0.0, Omega = 0.0, *Vorticity;

  errU = new su2double[nVar];
  errV = new su2double[nPrimVar];

  /*--- Set the primitive variables ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    /*--- Set the primitive variables incompressible (dens, vx, vy, vz, beta)
          and compressible (temp, vx, vy, vz, press, dens, enthal, sos)---*/
    nonPhys = nodes->SetPrimVar_Compressible(iPoint, config);
    if (nonPhys) {
      ErrorCounter++;
    }

    /*--- Initialize the convective, source and viscous residual vector ---*/
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
  }

  /*--- Allowing for Primitive Variables to be passed ---*/
  InitiateComms(geometry, config, PRIMITIVE);
  CompleteComms(geometry, config, PRIMITIVE);

  /*--- Artificial dissipation ---*/
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Gradient computation ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
    SetSolution_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
    SetSolution_Gradient_LS(geometry, config);

  }

  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/
  if ((iMesh == MESH_0) && (limiter_flow || limiter_turb)
      && !Output && !van_albada) { SetPrimitive_Limiter(geometry, config); }


  // THIS ISNT SET UP YET
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  //StrainMag_Max = 0.0; Omega_Max = 0.0;
  //for (iPoint = 0; iPoint < nPoint; iPoint++) {

  //  solver_container[TNE2_SOL]->nodes->SetVorticity();
  //  solver_container[TNE2_SOL]->nodes->SetStrainMag();

  //  StrainMag = solver_container[TNE2_SOL]->nodes->GetStrainMag();
  //  Vorticity = solver_container[TNE2_SOL]->nodes->GetVorticity();
  //  Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);

  //  StrainMag_Max = max(StrainMag_Max, StrainMag);
  //  Omega_Max = max(Omega_Max, Omega);
  //}

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
      solver_container[TNE2_SOL]->SetStrainMag_Max(StrainMag_Max);
      solver_container[TNE2_SOL]->SetOmega_Max(Omega_Max);
    }
  }
}

void CTNE2NSSolver::SetTime_Step(CGeometry *geometry,
                                 CSolver **solution_container,
                                 CConfig *config,
                                 unsigned short iMesh,
                                 unsigned long Iteration) {

  unsigned short iDim, iMarker, iSpecies;
  unsigned short VEL_INDEX, RHO_INDEX, RHOS_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  su2double *Normal, Area, Vol;
  su2double Mean_SoundSpeed, Mean_ProjVel;
  su2double Lambda, Local_Delta_Time, Local_Delta_Time_Visc, Global_Delta_Time;
  su2double Mean_LaminarVisc, Mean_ThermalCond, Mean_ThermalCond_ve, Mean_Density, Mean_Tve;
  su2double cp, cv, cvve, RuSI, Ru, *xi, *Ms, Na, Mmix, Rmix;
  su2double Lambda_1, Lambda_2, K_v, Global_Delta_UnstTimeND;
  su2double *V_i, *V_j, *X;
  su2double UnitNormal[3];
  su2double tmp;

  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));

  /*--- Initialize parameters ---*/
  Global_Delta_Time = 1E6;
  Min_Delta_Time    = 1.E6;
  Max_Delta_Time    = 0.0;
  K_v    = 0.5;
  iPoint = 0;
  jPoint = 0;
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  Na   = AVOGAD_CONSTANT;

  A_INDEX       = nodes->GetAIndex();
  VEL_INDEX     = nodes->GetVelIndex();
  RHO_INDEX     = nodes->GetRhoIndex();
  RHOS_INDEX    = nodes->GetRhosIndex();
  RHOCVTR_INDEX = nodes->GetRhoCvtrIndex();
  RHOCVVE_INDEX = nodes->GetRhoCvveIndex();

  X = new su2double[nSpecies];

  /*--- Get from config ---*/
  xi = config->GetRotationModes();
  Ms = config->GetMolar_Mass();

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    nodes->SetMax_Lambda_Inv(iPoint,0.0);
    nodes->SetMax_Lambda_Visc(iPoint,0.0);
  }

  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Calculate geometric quantities ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    Normal = geometry->edge[iEdge]->GetNormal();
    Area   = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);
    for (iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    /*--- Acquire the primitive variable information at each node ---*/
    V_i = nodes->GetPrimitive(iPoint);
    V_j = nodes->GetPrimitive(jPoint);

    /*--- Calculate the required mean values ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_ProjVel      = 0.5*( V_i[VEL_INDEX+iDim]
          +V_j[VEL_INDEX+iDim] )*UnitNormal[iDim];
    Mean_SoundSpeed     = 0.5*(V_i[A_INDEX]   + V_j[A_INDEX]);
    Mean_Density        = 0.5*(V_i[RHO_INDEX] + V_j[RHO_INDEX]);

    /*--- Calculate the maximum spectral radius from convection ---*/
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
    if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(iPoint,Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(jPoint,Lambda);

    /*--- Calculate mean viscous quantities ---*/
    Mean_LaminarVisc    = 0.5*(nodes->GetLaminarViscosity(iPoint) +
                               nodes->GetLaminarViscosity(jPoint)  );
    Mean_ThermalCond    = 0.5*(nodes->GetThermalConductivity(iPoint) +
                               nodes->GetThermalConductivity(jPoint)  );
    Mean_ThermalCond_ve = 0.5*(nodes->GetThermalConductivity_ve(iPoint) +
                               nodes->GetThermalConductivity_ve(jPoint)  );
    Mean_Density        = 0.5*(nodes->GetDensity(iPoint) +
                               nodes->GetDensity(jPoint)  );
    cv = 0.5*(nodes->GetRhoCv_tr(iPoint) + nodes->GetRhoCv_ve(iPoint) +
              nodes->GetRhoCv_tr(jPoint) + nodes->GetRhoCv_ve(jPoint)  )/ Mean_Density;

    /*--- Determine the viscous spectral radius and apply it to the control volume ---*/
    Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
    Lambda_2 = (Mean_ThermalCond+Mean_ThermalCond_ve)/cv;
    Lambda   = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;

    if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Visc(iPoint, Lambda);
    if (geometry->node[jPoint]->GetDomain()) nodes->AddMax_Lambda_Visc(jPoint, Lambda);
  }

  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area   = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Acquire the primitive variable information at each node ---*/
      V_i = nodes->GetPrimitive(iPoint);

      /*--- Calculate the required mean values ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Mean_ProjVel      = V_i[VEL_INDEX+iDim]*UnitNormal[iDim];
      Mean_SoundSpeed     = V_i[A_INDEX];
      Mean_Density        = V_i[RHO_INDEX];
      Mean_ThermalCond    = nodes->GetThermalConductivity(iPoint);
      Mean_ThermalCond_ve = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Calculate the maximum spectral radius from convection ---*/
      Lambda = (fabs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
      if (geometry->node[iPoint]->GetDomain())
        nodes->AddMax_Lambda_Inv(iPoint,Lambda);

      /*--- Calculate viscous mean quantities ---*/
      Mean_LaminarVisc    = nodes->GetLaminarViscosity(iPoint);
      Mean_ThermalCond    = nodes->GetThermalConductivity(iPoint);
      Mean_ThermalCond_ve = nodes->GetThermalConductivity_ve(iPoint);
      Mean_Density        = nodes->GetDensity(iPoint);
      cv = (nodes->GetRhoCv_tr(iPoint) +
            nodes->GetRhoCv_ve(iPoint)  ) / Mean_Density;

      Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
      Lambda_2 = (Mean_ThermalCond+Mean_ThermalCond_ve)/cv;
      Lambda   = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;

      if (geometry->node[iPoint]->GetDomain())
        nodes->AddMax_Lambda_Visc(iPoint,Lambda);
    }
  }

  /*--- Each element uses their own speed ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();

    /*--- Calculate local inv. and visc. dTs, take the minimum of the two ---*/
    Local_Delta_Time      = config->GetCFL(iMesh)*Vol / nodes->GetMax_Lambda_Inv(iPoint);
    Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ nodes->GetMax_Lambda_Visc(iPoint);
    Local_Delta_Time      = min(Local_Delta_Time, Local_Delta_Time_Visc);
    Global_Delta_Time     = min(Global_Delta_Time, Local_Delta_Time);

    /*--- Store minimum and maximum dt's within the grid for printing ---*/
    Min_Delta_Time        = min(Min_Delta_Time, Local_Delta_Time);
    Max_Delta_Time        = max(Max_Delta_Time, Local_Delta_Time);

    /*--- Set the time step ---*/
    nodes->SetDelta_Time(iPoint,Local_Delta_Time);
  }

  /*--- Communicate minimum and maximum time steps ---*/
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

  /*--- Check if there is any element with only one neighbor...
   a CV that is inside another CV ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (geometry->node[iPoint]->GetnPoint() == 1)
      nodes->SetDelta_Time(iPoint,Min_Delta_Time);
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
      nodes->SetDelta_Time(iPoint, Global_Delta_Time);
  }

  /*--- Recompute the unsteady time step for the dual time stratey
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
        /*--- Check if there is any element with only one neighbor...
         a CV that is inside another CV ---*/
        if (geometry->node[iPoint]->GetnPoint() == 1) Local_Delta_Time = 0.0;
        nodes->SetDelta_Time(iPoint, Local_Delta_Time);
      }
    }
  delete [] X;
}

void CTNE2NSSolver::Viscous_Residual(CGeometry *geometry,
                                     CSolver **solution_container,
                                     CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep) {
  bool implicit, err;
  unsigned short iVar, jVar;
  unsigned long iPoint, jPoint, iEdge;

  /*--- Determine time integration scheme ---*/
  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  numerics->SetPIndex      ( nodes->GetPIndex()       );
  numerics->SetTIndex      ( nodes->GetTIndex()       );
  numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  numerics->SetHIndex      ( nodes->GetHIndex()       );
  numerics->SetAIndex      ( nodes->GetAIndex()       );
  numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );


  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord() );
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Primitive variables, and gradient ---*/
    numerics->SetConservative   (nodes->GetSolution(iPoint),
                                 nodes->GetSolution(jPoint) );
    numerics->SetConsVarGradient(nodes->GetGradient(iPoint),
                                 nodes->GetGradient(jPoint) );
    numerics->SetPrimitive      (nodes->GetPrimitive(iPoint),
                                 nodes->GetPrimitive(jPoint) );
    numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                 nodes->GetGradient_Primitive(jPoint) );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU  (nodes->GetdPdU(iPoint),   nodes->GetdPdU(jPoint));
    numerics->SetdTdU  (nodes->GetdTdU(iPoint),   nodes->GetdTdU(jPoint));
    numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));
    numerics->SetEve   (nodes->GetEve(iPoint),    nodes->GetEve(jPoint));
    numerics->SetCvve  (nodes->GetCvve(iPoint),   nodes->GetCvve(jPoint));

    /*--- Species diffusion coefficients ---*/
    numerics->SetDiffusionCoeff(nodes->GetDiffusionCoeff(iPoint),
                                nodes->GetDiffusionCoeff(jPoint) );

    /*--- Laminar viscosity ---*/
    numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint),
                                  nodes->GetLaminarViscosity(jPoint) );

    /*--- Thermal conductivity ---*/
    numerics->SetThermalConductivity(nodes->GetThermalConductivity(iPoint),
                                     nodes->GetThermalConductivity(jPoint));

    /*--- Vib-el. thermal conductivity ---*/
    numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint),
                                        nodes->GetThermalConductivity_ve(jPoint) );

    /*--- Compute and update residual ---*/
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if (Res_Visc[iVar] != Res_Visc[iVar]) err = true;
    if (implicit)
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          if ((Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) ||
              (Jacobian_j[iVar][jVar] != Jacobian_j[iVar][jVar])   )
            err = true;

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      LinSysRes.AddBlock(jPoint, Res_Visc);
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
      }
    }
  } //iEdge
}

void CTNE2NSSolver::Source_Residual(CGeometry *geometry,
                                    CSolver **solution_container,
                                    CNumerics *numerics,
                                    CNumerics *second_solver, CConfig *config,
                                    unsigned short iMesh) {
  bool implicit, err;
  unsigned short iMarker, iVar, jVar;
  unsigned long iPoint, iVertex;
  unsigned long eAxi_local,  eChm_local,  eVib_local;
  unsigned long eAxi_global, eChm_global, eVib_global;
  int rank = MASTER_NODE;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  err = false;

  /*--- Initialize error counters ---*/
  eAxi_local = 0;
  eChm_local = 0;
  eVib_local = 0;

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  numerics->SetPIndex      ( nodes->GetPIndex()       );
  numerics->SetTIndex      ( nodes->GetTIndex()       );
  numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  numerics->SetHIndex      ( nodes->GetHIndex()       );
  numerics->SetAIndex      ( nodes->GetAIndex()       );
  numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  /*--- loop over interior points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    //    if (geometry->node[iPoint]->GetSolidBoundary()) {

    /*--- Set conserved & primitive variables  ---*/
    numerics->SetConservative(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));
    numerics->SetPrimitive   (nodes->GetPrimitive(iPoint),  nodes->GetPrimitive(iPoint) );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU  (nodes->GetdPdU(iPoint),   nodes->GetdPdU(iPoint)  );
    numerics->SetdTdU  (nodes->GetdTdU(iPoint),   nodes->GetdTdU(iPoint)  );
    numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(iPoint));
    numerics->SetEve   (nodes->GetEve(iPoint),    nodes->GetEve(iPoint)   );
    numerics->SetCvve  (nodes->GetCvve(iPoint),   nodes->GetCvve(iPoint)  );

    /*--- Set volume of the dual grid cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[iPoint]->GetCoord() );

    /*--- Store the value of the source term residuals (only for visualization and debugging) ---*/
    if (config->GetExtraOutput()) {
      for (iVar = 0; iVar < nVar; iVar++) {
        OutputVariables[iPoint* (unsigned long) nOutputVariables + iVar] = 0.0;
      }
    }

    /*--- Initialize source vector ---*/
    for (iVar = 0; iVar< nVar; iVar++) Source[iVar] = 0.0;

    /*--- Compute axisymmetric source terms (if needed) ---*/
    if (config->GetAxisymmetric()) {
      numerics->ComputeAxisymmetric(Residual, Source, Jacobian_i, config);

      /*--- Check for errors in the axisymmetric source ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        if (err)
          break;
        if (Residual[iVar] != Residual[iVar])
          err = true;
        if (implicit)
          for (jVar = 0; jVar < nVar; jVar++)
            if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) {
              err = true;
              break;
            }
      }

      /*--- Apply the update to the linear system ---*/
      if (!err) {
        LinSysRes.AddBlock(iPoint, Residual);
        if (implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        eAxi_local++;
    }

    /*--- Compute the non-equilibrium chemistry ---*/
    numerics->ComputeChemistry(Residual, Source, Jacobian_i, config);

    /*--- Check for errors in the Chemical source term ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++) {
      if (err)
        break;
      if (Residual[iVar] != Residual[iVar])
        err = true;
      if (implicit)
        for (jVar = 0; jVar < nVar; jVar++)
          if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) {
            err = true;
            break;
          }
    }

    /*--- Apply the chemical sources to the linear system ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, Residual);
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    } else
      eChm_local++;

    /*--- Store the value of the source term residuals (only for visualization and debugging) ---*/
    if (config->GetExtraOutput()) {
      for (iVar = 0; iVar < nVar; iVar++) {
        OutputVariables[iPoint* (unsigned long) nOutputVariables + iVar] += Residual[iVar];
      }
    }

    /*--- Compute vibrational energy relaxation ---*/
    // NOTE: Jacobians don't account for relaxation time derivatives
    numerics->ComputeVibRelaxation(Residual, Source, Jacobian_i, config);

    /*--- Check for errors in the relaxation source term ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++) {
      if (err)
        break;
      if (Residual[iVar] != Residual[iVar])
        err = true;
      if (implicit)
        for (jVar = 0; jVar < nVar; jVar++)
          if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) {
            err = true;
            break;
          }
    }

    /*--- Apply the vibrational relaxation terms to the linear system ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, Residual);
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    } else
      eVib_local++;

    /*--- Store the value of the source term residuals (only for visualization and debugging) ---*/
    if (config->GetExtraOutput()) {
      for (iVar = 0; iVar < nVar; iVar++) {
        OutputVariables[iPoint* (unsigned long) nOutputVariables + iVar] += Residual[iVar];
      }
    }
  }

  /*-#ifndef NO_MPI
  rank = MPI::COMM_WORLD.Get_rank();
  MPI::COMM_WORLD.Reduce(&eAxi_local, &eAxi_global, 1, MPI::UNSIGNED_LONG,
                         MPI::SUM, MASTER_NODE                             );
  MPI::COMM_WORLD.Reduce(&eChm_local, &eChm_global, 1, MPI::UNSIGNED_LONG,
                         MPI::SUM, MASTER_NODE                             );
  MPI::COMM_WORLD.Reduce(&eVib_local, &eVib_global, 1, MPI::UNSIGNED_LONG,
                         MPI::SUM, MASTER_NODE                             );
  #else
  eAxi_global = eAxi_local;
  eChm_global = eChm_local;
  eVib_global = eVib_local;
  #endif

  if ((rank == MASTER_NODE) &&
      (
       (eAxi_global != 0) ||
       (eChm_global != 0) ||
       (eVib_global != 0)
       )
      ) {
    cout << "Warning!! Instances of NaN in the following source terms: " << endl;
    cout << "Axisymmetry: " << eAxi_global << endl;
    cout << "Chemical:    " << eChm_global << endl;
    cout << "Vib. Relax:  " << eVib_global << endl;
  }
  -*/

  /*--- Loop over boundaries ---*/
  //  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
  //		switch (config->GetMarker_All_Boundary(iMarker)) {
  //      case EULER_WALL: case SYMMETRY_PLANE: case FAR_FIELD:
  //      case HEAT_FLUX: case ISOTHERMAL:
  //
  //        /*--- Loop over all of the vertices on this boundary marker ---*/
  //        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
  //          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  //
  //          /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
  //          if (geometry->node[iPoint]->GetDomain()) {
  //            /*--- Set conserved & primitive variables  ---*/
  //            numerics->SetConservative(node[iPoint]->GetSolution(),
  //                                      node[iPoint]->GetSolution());
  //            numerics->SetPrimitive   (node[iPoint]->GetPrimVar(),
  //                                      node[iPoint]->GetPrimVar() );
  //            numerics->SetdPdU        (node[iPoint]->GetdPdU(),
  //                                      node[iPoint]->GetdPdU());
  //            numerics->SetdTdU        (node[iPoint]->GetdTdU(),
  //                                      node[iPoint]->GetdTdU());
  //            numerics->SetdTvedU      (node[iPoint]->GetdTvedU(),
  //                                      node[iPoint]->GetdTvedU());
  //
  //            /*--- Set volume of the dual grid cell ---*/
  //            numerics->SetVolume(geometry->node[iPoint]->GetVolume());
  //
  //            /*--- Compute the non-equilibrium chemistry ---*/
  //            numerics->ComputeChemistry(Residual, Jacobian_i, config);
  //            LinSysRes.SubtractBlock(iPoint, Residual);
  //            if (implicit)
  //              Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
  //
  //            /*--- Error checking ---*/
  //            for (iVar = 0; iVar < nVar; iVar++)
  //              if (Residual[iVar] != Residual[iVar])
  //                cout << "NaN in Chemistry Residual" << endl;
  //            if (implicit) {
  //              for (iVar = 0; iVar < nVar; iVar++) {
  //                for (jVar = 0; jVar < nVar; jVar++) {
  //                  if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar])
  //                    cout << "NaN in Chemistry Jacobian i" << endl;
  //                }
  //              }
  //            }
  //            /*--- Compute vibrational energy relaxation ---*/
  //              // NOTE: Jacobians don't account for relaxation time derivatives
  //            numerics->ComputeVibRelaxation(Residual, Jacobian_i, config);
  //            LinSysRes.SubtractBlock(iPoint, Residual);
  //            if (implicit)
  //              Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
  //
  //            /*--- Error checking ---*/
  //            for (iVar = 0; iVar < nVar; iVar++)
  //              if (Residual[iVar] != Residual[iVar])
  //                cout << "NaN in vibrational Residual" << endl;
  //            if (implicit) {
  //              for (iVar = 0; iVar < nVar; iVar++) {
  //                for (jVar = 0; jVar < nVar; jVar++) {
  //                  if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar])
  //                    cout << "NaN in vibrational Jacobian i" << endl;
  //                }
  //              }
  //            }
  //          }
  //        }
  //        break;
  //
  //      case HEAT_FLUX_NONCATALYTIC: case HEAT_FLUX_CATALYTIC:
  //
  //        /*--- Loop over all of the vertices on this boundary marker ---*/
  //        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
  //          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  //
  //          /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
  //          if (geometry->node[iPoint]->GetDomain()) {
  //            /*--- Set conserved & primitive variables  ---*/
  //            numerics->SetConservative(node[iPoint]->GetSolution(),
  //                                      node[iPoint]->GetSolution());
  //            numerics->SetPrimitive   (node[iPoint]->GetPrimVar(),
  //                                      node[iPoint]->GetPrimVar() );
  //            numerics->SetdPdU        (node[iPoint]->GetdPdU(),
  //                                      node[iPoint]->GetdPdU());
  //            numerics->SetdTdU        (node[iPoint]->GetdTdU(),
  //                                      node[iPoint]->GetdTdU());
  //            numerics->SetdTvedU      (node[iPoint]->GetdTvedU(),
  //                                      node[iPoint]->GetdTvedU());
  //
  //            /*--- Set volume of the dual grid cell ---*/
  //            numerics->SetVolume(geometry->node[iPoint]->GetVolume());
  //
  //            /*--- Compute vibrational energy relaxation ---*/
  //              // NOTE: Jacobians don't account for relaxation time derivatives
  //            numerics->ComputeVibRelaxation(Residual, Jacobian_i, config);
  //            LinSysRes.SubtractBlock(iPoint, Residual);
  //            if (implicit)
  //              Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
  //
  //            /*--- Error checking ---*/
  //            for (iVar = 0; iVar < nVar; iVar++)
  //              if (Residual[iVar] != Residual[iVar])
  //                cout << "NaN in vibrational Residual" << endl;
  //            if (implicit) {
  //              for (iVar = 0; iVar < nVar; iVar++) {
  //                for (jVar = 0; jVar < nVar; jVar++) {
  //                  if (Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar])
  //                    cout << "NaN in vibrational Jacobian i" << endl;
  //                }
  //              }
  //            }
  //          }
  //        }
  //        break;
  //
  //      case ISOTHERMAL_NONCATALYTIC: case ISOTHERMAL_CATALYTIC:
  //          // NO SOURCE TERMS TO BE SET HERE!
  //        break;
  //    }
  //  }
  //  }
}

void CTNE2NSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {

  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  unsigned short VEL_INDEX, T_INDEX, TVE_INDEX;
  unsigned long iVertex, iPoint, iPointNormal;
  su2double **Grad_PrimVar;
  su2double Delta, Viscosity, ThermalCond, ThermalCond_ve;
  su2double TauNormal;
  su2double FrictionVel;
  su2double *Normal, *Coord, *Coord_Normal, Area;
  su2double Force[3];
  su2double MomentDist[3];
  su2double RefDensity, Density;
  su2double div_vel, RefVel2;
  su2double dTn, dTven, pnorm;
  su2double Alpha, Beta, RefLength, RefArea, *Origin;
  su2double factor;
  su2double MaxNorm = 8.0;

  su2double WallShearStress, WallDistMod, WallDist[3];

  su2double Vel[3], Velocity_Inf[3];

  su2double UnitNormal[3];
  su2double TauElem[3];
  su2double TauTangent[3];
  su2double Tau[3][3];
  su2double MomentX_Force[3], MomentY_Force[3], MomentZ_Force[3];
  string Marker_Tag, Monitoring_Tag;

#ifdef HAVE_MPI
  su2double MyAllBound_CD_Visc, MyAllBound_CL_Visc, MyAllBound_CSF_Visc, MyAllBound_CMx_Visc, MyAllBound_CMy_Visc, MyAllBound_CMz_Visc, MyAllBound_CoPx_Visc, MyAllBound_CoPy_Visc, MyAllBound_CoPz_Visc, MyAllBound_CFx_Visc, MyAllBound_CFy_Visc, MyAllBound_CFz_Visc, MyAllBound_CT_Visc, MyAllBound_CQ_Visc, MyAllBound_HF_Visc, MyAllBound_MaxHF_Visc, *MySurface_CL_Visc = NULL, *MySurface_CD_Visc = NULL, *MySurface_CSF_Visc = NULL, *MySurface_CEff_Visc = NULL, *MySurface_CFx_Visc = NULL, *MySurface_CFy_Visc = NULL, *MySurface_CFz_Visc = NULL, *MySurface_CMx_Visc = NULL, *MySurface_CMy_Visc = NULL, *MySurface_CMz_Visc = NULL, *MySurface_HF_Visc = NULL, *MySurface_MaxHF_Visc;
#endif

  /*--- Retrieve index information from CVariable ---*/
  VEL_INDEX = nodes->GetVelIndex();
  T_INDEX   = nodes->GetTIndex();
  TVE_INDEX = nodes->GetTveIndex();

  /*--- Retrieve data from CConfig ---*/
  pnorm = config->GetPnormHeat();

  /*--- Calculate angle of attack & sideslip ---*/
  Alpha = config->GetAoA()*PI_NUMBER/180.0;
  Beta  = config->GetAoS()*PI_NUMBER/180.0;

  /*--- Determine reference geometrical parameters ---*/
  RefArea    = config->GetRefArea();
  RefLength  = config->GetRefLength();
  Origin     = config->GetRefOriginMoment(0);

  /*--- Get reference values from the freestream node. ---*/
  RefVel2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_Inf[iDim] = node_infty->GetVelocity(1,iDim);
    RefVel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  RefDensity  = node_infty->GetDensity(1);

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*-- Initialization --*/
  AllBound_CMx_Visc  = 0.0; AllBound_CMy_Visc   = 0.0; AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc  = 0.0; AllBound_CFy_Visc   = 0.0; AllBound_CFz_Visc = 0.0;
  AllBound_CD_Visc   = 0.0; AllBound_CL_Visc    = 0.0;
  AllBound_HF_Visc   = 0.0; AllBound_MaxHF_Visc = 0.0;
  AllBound_CEff_Visc = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Visc[iMarker_Monitoring]  = 0.0; Surface_CD_Visc[iMarker_Monitoring]    = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring]  = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring] = 0.0; Surface_CFy_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring] = 0.0; Surface_CMx_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring] = 0.0; Surface_CMz_Visc[iMarker_Monitoring]   = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]  = 0.0; Surface_MaxHF_Visc[iMarker_Monitoring] = 0.0;
  }

  /*--- Loop over the Navier-Stokes markers ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- Identify boundary information ---*/
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

    /*--- Solve for the forces ---*/
    if ((Boundary == HEAT_FLUX              ) ||
        (Boundary == HEAT_FLUX_CATALYTIC    ) ||
        (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == ISOTHERMAL_CATALYTIC   ) ||
        (Boundary == ISOTHERMAL_NONCATALYTIC)) {

      /*--- Forces initialization at each Marker ---*/
      CD_Visc[iMarker]   = 0.0; CL_Visc[iMarker]    = 0.0; CSF_Visc[iMarker]    = 0.0;
      CFx_Visc[iMarker]  = 0.0; CFy_Visc[iMarker]   = 0.0; CFz_Visc[iMarker]    = 0.0;
      CMx_Visc[iMarker]  = 0.0; CMy_Visc[iMarker]   = 0.0; CMz_Visc[iMarker]    = 0.0;
      CoPx_Visc[iMarker] = 0.0; CoPy_Visc[iMarker]  = 0.0; CoPz_Visc[iMarker]   = 0.0;
      CT_Visc[iMarker]   = 0.0; CQ_Visc[iMarker]    = 0.0; CMerit_Visc[iMarker] = 0.0;
      HF_Visc[iMarker]   = 0.0; MaxHF_Visc[iMarker] = 0.0; CEff_Visc[iMarker]   = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) {
        ForceViscous[iDim]  = 0.0; MomentViscous[iDim] = 0.0;
        MomentX_Force[iDim] = 0.0; MomentY_Force[iDim] = 0.0; MomentZ_Force[iDim] = 0.0;
      }

      /*--- Loop over the boundary points ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Acquire & calculate geometric parameters ---*/
        iPoint       = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
        Coord        = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPointNormal]->GetCoord();
        Normal       = geometry->vertex[iMarker][iVertex]->GetNormal();

        Area         = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);

        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
          MomentDist[iDim] = Coord[iDim] - Origin[iDim];
        }

        /*--- Get vertex flow parameters ---*/
        Grad_PrimVar   = nodes->GetGradient_Primitive(iPoint);
        Viscosity      = nodes->GetLaminarViscosity(iPoint);
        ThermalCond    = nodes->GetThermalConductivity(iPoint);
        ThermalCond_ve = nodes->GetThermalConductivity_ve(iPoint);
        Density        = nodes->GetDensity(iPoint);

        /*--- Calculate the viscous stress tensor ---*/
        div_vel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          div_vel += Grad_PrimVar[VEL_INDEX+iDim][iDim];

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Delta = 0.0; if (iDim == jDim) Delta = 1.0;
            Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[VEL_INDEX+jDim][iDim] +
                Grad_PrimVar[VEL_INDEX+iDim][jDim]  )
                - TWO3*Viscosity*div_vel*Delta;
          }
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++)
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
        }

        /*--- Compute wall shear stress (using the stress tensor) ---*/
        TauNormal = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          TauNormal += TauElem[iDim] * UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++)
          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];

        WallShearStress = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          WallShearStress += TauTangent[iDim]*TauTangent[iDim];
        WallShearStress = sqrt(WallShearStress);

        for (iDim = 0; iDim < nDim; iDim++)
          WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          WallDistMod += WallDist[iDim]*WallDist[iDim];
        WallDistMod = sqrt(WallDistMod);

        /*--- Compute wall skin friction coefficient, and heat flux on the wall ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          CSkinFriction[iMarker][iDim][iVertex] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);

        /*--- Compute y+ and non-dimensional velocity ---*/
        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);

        /*--- Compute heat flux on the wall ---*/
        dTn = 0.0; dTven = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          dTn   += Grad_PrimVar[T_INDEX][iDim]*UnitNormal[iDim];
          dTven += Grad_PrimVar[TVE_INDEX][iDim]*UnitNormal[iDim];
        }

        HeatFlux[iMarker][iVertex] = ThermalCond*dTn + ThermalCond_ve*dTven;
        HF_Visc[iMarker] += HeatFlux[iMarker][iVertex]*Area;
        MaxHF_Visc[iMarker] += pow(HeatFlux[iMarker][iVertex], pnorm)*Area;

        /*--- Compute viscous forces, and moment using the stress tensor ---*/
        if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {

          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim]*Area*factor;
            ForceViscous[iDim] += Force[iDim];
          }

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

      /*--- Transform ForceInviscid into CLift and CDrag ---*/
      if  (Monitoring == YES) {

        if (nDim == 2) {
          CD_Visc[iMarker]     =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CL_Visc[iMarker]     = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker]   = CL_Visc[iMarker] / (CD_Visc[iMarker]+EPS);
          CFx_Visc[iMarker]    = ForceViscous[0];
          CFy_Visc[iMarker]    = ForceViscous[1];
          CMz_Visc[iMarker]    = MomentViscous[2];
          CoPx_Visc[iMarker]   = MomentZ_Force[1];
          CoPy_Visc[iMarker]   = -MomentZ_Force[0];
          CT_Visc[iMarker]     = -CFx_Visc[iMarker];
          CQ_Visc[iMarker]     = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker] = CT_Visc[iMarker] / (CQ_Visc[iMarker]+EPS);
          MaxHF_Visc[iMarker]  = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }

        if (nDim == 3) {
          CD_Visc[iMarker]     =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          CL_Visc[iMarker]     = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSF_Visc[iMarker]    = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker]   = CL_Visc[iMarker]/(CD_Visc[iMarker] + EPS);
          CFx_Visc[iMarker]    = ForceViscous[0];
          CFy_Visc[iMarker]    = ForceViscous[1];
          CFz_Visc[iMarker]    = ForceViscous[2];
          CMx_Visc[iMarker]    = MomentViscous[0];
          CMy_Visc[iMarker]    = MomentViscous[1];
          CMz_Visc[iMarker]    = MomentViscous[2];
          CoPx_Visc[iMarker]   =  -MomentY_Force[0];
          CoPz_Visc[iMarker]   = MomentY_Force[2];
          CT_Visc[iMarker]     = -CFz_Visc[iMarker];
          CQ_Visc[iMarker]     = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker] = CT_Visc[iMarker] / (CQ_Visc[iMarker] + EPS);
          MaxHF_Visc[iMarker]  = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
        }

        AllBound_CD_Visc    += CD_Visc[iMarker];
        AllBound_CL_Visc    += CL_Visc[iMarker];
        AllBound_CSF_Visc   += CSF_Visc[iMarker];
        AllBound_CFx_Visc   += CFx_Visc[iMarker];
        AllBound_CFy_Visc   += CFy_Visc[iMarker];
        AllBound_CFz_Visc   += CFz_Visc[iMarker];
        AllBound_CMx_Visc   += CMx_Visc[iMarker];
        AllBound_CMy_Visc   += CMy_Visc[iMarker];
        AllBound_CMz_Visc   += CMz_Visc[iMarker];
        AllBound_CoPx_Visc  += CoPx_Visc[iMarker];
        AllBound_CoPy_Visc  += CoPy_Visc[iMarker];
        AllBound_CoPz_Visc  += CoPz_Visc[iMarker];
        AllBound_CT_Visc    += CT_Visc[iMarker];
        AllBound_CQ_Visc    += CQ_Visc[iMarker];
        AllBound_HF_Visc    += HF_Visc[iMarker];
        AllBound_MaxHF_Visc += pow(MaxHF_Visc[iMarker], MaxNorm);

        /*--- Compute the coefficients per surface ---*/
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CL_Visc[iMarker_Monitoring]    += CL_Visc[iMarker];
            Surface_CD_Visc[iMarker_Monitoring]    += CD_Visc[iMarker];
            Surface_CSF_Visc[iMarker_Monitoring]   += CSF_Visc[iMarker];
            Surface_CEff_Visc[iMarker_Monitoring]  += CEff_Visc[iMarker];
            Surface_CFx_Visc[iMarker_Monitoring]   += CFx_Visc[iMarker];
            Surface_CFy_Visc[iMarker_Monitoring]   += CFy_Visc[iMarker];
            Surface_CFz_Visc[iMarker_Monitoring]   += CFz_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]   += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]   += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]   += CMz_Visc[iMarker];
            Surface_HF_Visc[iMarker_Monitoring]    += HF_Visc[iMarker];
            Surface_MaxHF_Visc[iMarker_Monitoring] += pow(MaxHF_Visc[iMarker],MaxNorm);
          }
        }
      }
    }
  }

  /*--- Update some global coeffients ---*/
  AllBound_CEff_Visc   = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  AllBound_MaxHF_Visc  = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/
  MyAllBound_CD_Visc    = AllBound_CD_Visc;   AllBound_CD_Visc = 0.0;
  MyAllBound_CL_Visc    = AllBound_CL_Visc;   AllBound_CL_Visc = 0.0;
  MyAllBound_CSF_Visc   = AllBound_CSF_Visc;  AllBound_CSF_Visc = 0.0;
  AllBound_CEff_Visc    = 0.0;
  MyAllBound_CMx_Visc   = AllBound_CMx_Visc;  AllBound_CMx_Visc = 0.0;
  MyAllBound_CMy_Visc   = AllBound_CMy_Visc;  AllBound_CMy_Visc = 0.0;
  MyAllBound_CMz_Visc   = AllBound_CMz_Visc;  AllBound_CMz_Visc = 0.0;
  MyAllBound_CoPx_Visc  = AllBound_CoPx_Visc; AllBound_CoPx_Visc = 0.0;
  MyAllBound_CoPy_Visc  = AllBound_CoPy_Visc; AllBound_CoPy_Visc = 0.0;
  MyAllBound_CoPz_Visc  = AllBound_CoPz_Visc; AllBound_CoPz_Visc = 0.0;
  MyAllBound_CFx_Visc   = AllBound_CFx_Visc;  AllBound_CFx_Visc = 0.0;
  MyAllBound_CFy_Visc   = AllBound_CFy_Visc;  AllBound_CFy_Visc = 0.0;
  MyAllBound_CFz_Visc   = AllBound_CFz_Visc;  AllBound_CFz_Visc = 0.0;
  MyAllBound_CT_Visc    = AllBound_CT_Visc;   AllBound_CT_Visc = 0.0;
  MyAllBound_CQ_Visc    = AllBound_CQ_Visc;   AllBound_CQ_Visc = 0.0;
  AllBound_CMerit_Visc  = 0.0;
  MyAllBound_HF_Visc    = AllBound_HF_Visc;   AllBound_HF_Visc = 0.0;
  MyAllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, MaxNorm);
  AllBound_MaxHF_Visc = 0.0;

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

  MySurface_CL_Visc    = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CD_Visc    = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CSF_Visc   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CEff_Visc  = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFx_Visc   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFy_Visc   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CFz_Visc   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Visc   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Visc   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Visc   = new su2double[config->GetnMarker_Monitoring()];
  MySurface_HF_Visc    = new su2double[config->GetnMarker_Monitoring()];
  MySurface_MaxHF_Visc = new su2double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {

    MySurface_CL_Visc[iMarker_Monitoring]    = Surface_CL_Visc[iMarker_Monitoring];
    MySurface_CD_Visc[iMarker_Monitoring]    = Surface_CD_Visc[iMarker_Monitoring];
    MySurface_CSF_Visc[iMarker_Monitoring]   = Surface_CSF_Visc[iMarker_Monitoring];
    MySurface_CEff_Visc[iMarker_Monitoring]  = Surface_CEff_Visc[iMarker_Monitoring];
    MySurface_CFx_Visc[iMarker_Monitoring]   = Surface_CFx_Visc[iMarker_Monitoring];
    MySurface_CFy_Visc[iMarker_Monitoring]   = Surface_CFy_Visc[iMarker_Monitoring];
    MySurface_CFz_Visc[iMarker_Monitoring]   = Surface_CFz_Visc[iMarker_Monitoring];
    MySurface_CMx_Visc[iMarker_Monitoring]   = Surface_CMx_Visc[iMarker_Monitoring];
    MySurface_CMy_Visc[iMarker_Monitoring]   = Surface_CMy_Visc[iMarker_Monitoring];
    MySurface_CMz_Visc[iMarker_Monitoring]   = Surface_CMz_Visc[iMarker_Monitoring];
    MySurface_HF_Visc[iMarker_Monitoring]    = Surface_HF_Visc[iMarker_Monitoring];
    MySurface_MaxHF_Visc[iMarker_Monitoring] = Surface_MaxHF_Visc[iMarker_Monitoring];

    Surface_CL_Visc[iMarker_Monitoring]    = 0.0;
    Surface_CD_Visc[iMarker_Monitoring]    = 0.0;
    Surface_CSF_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CEff_Visc[iMarker_Monitoring]  = 0.0;
    Surface_CFx_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CFy_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CFz_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]   = 0.0;
    Surface_HF_Visc[iMarker_Monitoring]    = 0.0;
    Surface_MaxHF_Visc[iMarker_Monitoring] = 0.0;
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

  delete [] MySurface_CL_Visc;   delete [] MySurface_CD_Visc;  delete [] MySurface_CSF_Visc;
  delete [] MySurface_CEff_Visc; delete [] MySurface_CFx_Visc; delete [] MySurface_CFy_Visc;
  delete [] MySurface_CFz_Visc;  delete [] MySurface_CMx_Visc; delete [] MySurface_CMy_Visc;
  delete [] MySurface_CMz_Visc;  delete [] MySurface_HF_Visc;  delete [] MySurface_MaxHF_Visc;

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value)---*/

  Total_CD      += AllBound_CD_Visc;
  Total_CL      += AllBound_CL_Visc;
  Total_CSF     += AllBound_CSF_Visc;
  Total_CEff     = Total_CL / (Total_CD + EPS);
  Total_CFx     += AllBound_CFx_Visc;
  Total_CFy     += AllBound_CFy_Visc;
  Total_CFz     += AllBound_CFz_Visc;
  Total_CMx     += AllBound_CMx_Visc;
  Total_CMy     += AllBound_CMy_Visc;
  Total_CMz     += AllBound_CMz_Visc;
  Total_CoPx    += AllBound_CoPx_Visc;
  Total_CoPy    += AllBound_CoPy_Visc;
  Total_CoPz    += AllBound_CoPz_Visc;
  Total_CT      += AllBound_CT_Visc;
  Total_CQ      += AllBound_CQ_Visc;
  Total_CMerit   = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  Total_Heat     = AllBound_HF_Visc;
  Total_MaxHeat  = AllBound_MaxHF_Visc;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL[iMarker_Monitoring]   += Surface_CL_Visc[iMarker_Monitoring];
    Surface_CD[iMarker_Monitoring]   += Surface_CD_Visc[iMarker_Monitoring];
    Surface_CSF[iMarker_Monitoring]  += Surface_CSF_Visc[iMarker_Monitoring];
    Surface_CEff[iMarker_Monitoring]  = Surface_CL[iMarker_Monitoring] / (Surface_CD[iMarker_Monitoring] + EPS);
    Surface_CFx[iMarker_Monitoring]  += Surface_CFx_Visc[iMarker_Monitoring];
    Surface_CFy[iMarker_Monitoring]  += Surface_CFy_Visc[iMarker_Monitoring];
    Surface_CFz[iMarker_Monitoring]  += Surface_CFz_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]  += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]  += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]  += Surface_CMz_Visc[iMarker_Monitoring];
  }
}

void CTNE2NSSolver::BC_Sym_Plane(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker) {

  /*--- Call the Euler wall routine ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, visc_numerics,
                config, val_marker);

}

void CTNE2NSSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                     CSolver **solution_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *sour_numerics,
                                     CConfig *config,
                                     unsigned short val_marker) {

  /*--- Local variables ---*/
  bool implicit;
  unsigned short iDim, iVar;
  unsigned short T_INDEX, TVE_INDEX;
  unsigned long iVertex, iPoint, total_index;
  su2double Wall_HeatFlux, dTdn, dTvedn, ktr, kve, pcontrol;
  su2double *Normal, Area;
  su2double **GradV;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);

  /*--- Set "Proportional control" coefficient ---*/
  pcontrol = 1.0;

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config ---*/
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /*--- Get the locations of the primitive variables ---*/
  T_INDEX    = nodes->GetTIndex();
  TVE_INDEX  = nodes->GetTveIndex();

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
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
        Res_Visc[iVar] = 0.0;
      }

      /*--- Set the residual on the boundary with the specified heat flux ---*/
      // Note: Contributions from qtr and qve are used for proportional control
      //       to drive the solution toward the specified heatflux more quickly.
      GradV  = nodes->GetGradient_Primitive(iPoint);
      dTdn   = 0.0;
      dTvedn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
        dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
      }
      ktr = nodes->GetThermalConductivity(iPoint);
      //      kve = nodes->GetThermalConductivity_ve();
      //			Res_Visc[nSpecies+nDim]   += pcontrol*(ktr*dTdn+kve*dTvedn) +
      //                                   Wall_HeatFlux*Area;
      //      Res_Visc[nSpecies+nDim+1] += pcontrol*(kve*dTvedn) +
      //                                   Wall_HeatFlux*Area;
      //
      //			/*--- Apply viscous residual to the linear system ---*/
      //      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Apply the no-slip condition in a strong way ---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      nodes->SetVelocity_Old(iPoint,Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }
      if (implicit) {
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
    }
  }
}

void CTNE2NSSolver::BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry,
                                                 CSolver **solution_container,
                                                 CNumerics *conv_numerics,
                                                 CNumerics *sour_numerics,
                                                 CConfig *config,
                                                 unsigned short val_marker) {

  /*--- Call standard "HeatFlux" wall to apply no-slip & energy b.c.'s ---*/
  BC_HeatFlux_Wall(geometry, solution_container, conv_numerics,
                   sour_numerics, config, val_marker);

  //	/*--- Local variables ---*/
  //  bool implicit;
  //	unsigned short iDim, iSpecies, iVar;
  //  unsigned short RHOS_INDEX, RHO_INDEX, T_INDEX, TVE_INDEX;
  //	unsigned long iVertex, iPoint;
  //	double pcontrol;
  //  su2double rho, Ys, eves, hs;
  //	double *Normal, Area;
  //  su2double *Ds, *V, *dYdn, SdYdn;
  //  su2double **GradV, **GradY;
  //
  //  /*--- Assign booleans ---*/
  //	implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  //
  //  /*--- Set "Proportional control" coefficient ---*/
  //  pcontrol = 0.6;
  //
  //  /*--- Get the locations of the primitive variables ---*/
  //  RHOS_INDEX = nodes->GetRhosIndex();
  //  RHO_INDEX  = nodes->GetRhoIndex();
  //  T_INDEX    = nodes->GetTIndex();
  //  TVE_INDEX  = nodes->GetTveIndex();
  //
  //  /*--- Allocate arrays ---*/
  //  dYdn = new su2double[nSpecies];
  //  GradY = new su2double*[nSpecies];
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    GradY[iSpecies] = new su2double[nDim];
  //
  //	/*--- Loop over all of the vertices on this boundary marker ---*/
  //	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
  //		if (geometry->nodes->GetDomain()) {
  //
  //			/*--- Compute dual-grid area and boundary normal ---*/
  //			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
  //			Area = 0.0;
  //			for (iDim = 0; iDim < nDim; iDim++)
  //				Area += Normal[iDim]*Normal[iDim];
  //			Area = sqrt (Area);
  //
  //			/*--- Initialize the convective & viscous residuals to zero ---*/
  //			for (iVar = 0; iVar < nVar; iVar++)
  //				Res_Visc[iVar] = 0.0;
  //
  //      /*--- Get temperature gradient information ---*/
  //      V = nodes->GetPrimVar();
  //      GradV  = nodes->GetGradient_Primitive();
  //
  //      /*--- Rename for convenience ---*/
  //      rho = V[RHO_INDEX];
  //      Ds  = nodes->GetDiffusionCoeff();
  //
  //      /*--- Calculate normal derivative of mass fraction ---*/
  //      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //        Ys = V[RHOS_INDEX+iSpecies]/rho;
  //        dYdn[iSpecies] = 0.0;
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          dYdn[iSpecies] += 1.0/rho * (GradV[RHOS_INDEX+iSpecies][iDim] -
  //                                       Ys*GradV[RHO_INDEX][iDim])*Normal[iDim];
  //      }
  //
  //      /*--- Calculate supplementary quantities ---*/
  //      SdYdn = 0.0;
  //      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //        SdYdn += rho*Ds[iSpecies]*dYdn[iSpecies];
  //
  //      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //        Ys   = V[RHOS_INDEX+iSpecies]/rho;
  //        eves = nodes->CalcEve(config, V[TVE_INDEX], iSpecies);
  //        hs   = nodes->CalcHs(config, V[T_INDEX], eves, iSpecies);
  //        Res_Visc[iSpecies] = rho*Ds[iSpecies]*dYdn[iSpecies] - Ys*SdYdn;
  //        Res_Visc[nSpecies+nDim]   += Res_Visc[iSpecies]*hs;
  //        Res_Visc[nSpecies+nDim+1] += Res_Visc[iSpecies]*eves;
  //      }
  //
  //			/*--- Viscous contribution to the residual at the wall ---*/
  //      LinSysRes.SubtractBlock(iPoint, Res_Visc);
  //		}
  //	}
  //
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    delete [] GradY[iSpecies];
  //  delete [] GradY;
  //  delete [] dYdn;
}

void CTNE2NSSolver::BC_HeatFluxCatalytic_Wall(CGeometry *geometry,
                                              CSolver **solution_container,
                                              CNumerics *conv_numerics,
                                              CNumerics *sour_numerics,
                                              CConfig *config,
                                              unsigned short val_marker) {

  /*--- Local variables ---*/
  bool implicit, catalytic;
  unsigned short iDim, iSpecies, iVar;
  unsigned short T_INDEX, TVE_INDEX, RHOS_INDEX, RHO_INDEX;
  unsigned long iVertex, iPoint, total_index;
  su2double Wall_HeatFlux, dTdn, dTvedn, ktr, kve, pcontrol;
  su2double rho, Ys, eves, hs;
  su2double *Normal, Area;
  su2double *Ds, *V, *dYdn, SdYdn;
  su2double **GradV, **GradY;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  catalytic = false;

  /*--- Set "Proportional control" coefficient ---*/
  pcontrol = 0.6;

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config ---*/
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /*--- Get the locations of the primitive variables ---*/
  T_INDEX    = nodes->GetTIndex();
  TVE_INDEX  = nodes->GetTveIndex();
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Allocate arrays ---*/
  dYdn = new su2double[nSpecies];
  GradY = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    GradY[iSpecies] = new su2double[nDim];

  //  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  //  sour_numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  //  sour_numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  //  sour_numerics->SetPIndex      ( nodes->GetPIndex()       );
  //  sour_numerics->SetTIndex      ( nodes->GetTIndex()       );
  //  sour_numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  //  sour_numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  //  sour_numerics->SetHIndex      ( nodes->GetHIndex()       );
  //  sour_numerics->SetAIndex      ( nodes->GetAIndex()       );
  //  sour_numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  //  sour_numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
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
        Res_Visc[iVar] = 0.0;
        Res_Sour[iVar] = 0.0;
      }

      /*--- Assign wall velocity to "Vector" array ---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;

      /*--- Set the residual, truncation error, and velocity value ---*/
      nodes->SetVelocity_Old(iPoint,Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Get temperature gradient information ---*/
      V = nodes->GetPrimVar();
      GradV  = nodes->GetGradient_Primitive(iPoint);
      dTdn   = 0.0;
      dTvedn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
        dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
      }

      if (catalytic) {
        cout << "NEED TO IMPLEMENT CATALYTIC BOUNDARIES IN HEATFLUX!!!" << endl;
        exit(1);
      }
      else {

        /*--- Rename for convenience ---*/
        rho = V[RHO_INDEX];
        Ds  = nodes->GetDiffusionCoeff(iPoint);

        /*--- Calculate normal derivative of mass fraction ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          Ys = V[RHOS_INDEX+iSpecies]/rho;
          dYdn[iSpecies] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            dYdn[iSpecies] += 1.0/rho * (GradV[RHOS_INDEX+iSpecies][iDim] -
                Ys*GradV[RHO_INDEX][iDim])*Normal[iDim];
        }

        /*--- Calculate supplementary quantities ---*/
        SdYdn = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          SdYdn += rho*Ds[iSpecies]*dYdn[iSpecies];

        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          Ys   = V[RHOS_INDEX+iSpecies]/rho;
          eves = nodes->CalcEve(config, V[TVE_INDEX], iSpecies);
          hs   = nodes->CalcHs(config, V[T_INDEX], eves, iSpecies);
          //          Res_Visc[iSpecies] = -rho*Ds[iSpecies]*dYdn[iSpecies] + Ys*SdYdn;
          //          Res_Visc[nSpecies+nDim]   += Res_Visc[iSpecies]*hs;
          //          Res_Visc[nSpecies+nDim+1] += Res_Visc[iSpecies]*eves;
        }
      }

      /*--- Get node thermal conductivity ---*/
      ktr = nodes->GetThermalConductivity(iPoint);
      kve = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Set the residual on the boundary with the specified heat flux ---*/
      // Note: Contributions from qtr and qve are used for proportional control
      //       to drive the solution toward the specified heatflux more quickly.
      Res_Visc[nSpecies+nDim]   += pcontrol*(ktr*dTdn+kve*dTvedn) +
          Wall_HeatFlux*Area;
      Res_Visc[nSpecies+nDim+1] += pcontrol*(kve*dTvedn) +
          Wall_HeatFlux*Area;

      /*--- Viscous contribution to the residual at the wall ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      //      /*--- Apply the non-catalytic wall boundary ---*/
      //      // Note: We are re-calculating the chemistry residual and adding it to
      //      //       the linear system to eliminate the contribution from the solution
      //      //       (convention is to subtract sources)
      //      sour_numerics->SetConservative(nodes->GetSolution(),
      //                                     nodes->GetSolution() );
      //      sour_numerics->SetPrimitive   (nodes->GetPrimVar() ,
      //                                     nodes->GetPrimVar()  );
      //      sour_numerics->SetdPdU        (nodes->GetdPdU()    ,
      //                                     nodes->GetdPdU()     );
      //      sour_numerics->SetdTdU        (nodes->GetdTdU()    ,
      //                                     nodes->GetdTdU()     );
      //      sour_numerics->SetdTvedU      (nodes->GetdTvedU()  ,
      //                                     nodes->GetdTvedU()   );
      //      sour_numerics->SetVolume      (geometry->nodes->GetVolume());
      //      sour_numerics->ComputeChemistry(Res_Sour, Jacobian_i, config);
      //      LinSysRes.AddBlock(iPoint, Res_Sour);
      //      if (implicit)
      //        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal)/
       Note that we need to add a contribution for moving walls to the Jacobian. ---*/
      if (implicit) {
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }

    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] GradY[iSpecies];
  delete [] GradY;
  delete [] dYdn;
}

void CTNE2NSSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                       CSolver **solution_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *sour_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) {

  unsigned short iDim, iVar, jVar;
  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr, kve;
  su2double Ti, Tvei, Tj, Tvej, *dTdU, *dTvedU;
  su2double Twall, dij, theta;
  su2double Area, *Normal, UnitNormal[3];
  su2double *Coord_i, *Coord_j;
  su2double C;

  bool implicit   = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool ionization = config->GetIonization();

  if (ionization) {
    cout << "BC_ISOTHERMAL: NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION" << endl;
    exit(1);
  }

  /*--- Define 'proportional control' constant ---*/
  C = 5;

  /*--- Identify the boundary ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature ---*/
  Twall = config->GetIsothermal_Temperature(Marker_Tag);

  /*--- Extract necessary indices ---*/
  RHOS_INDEX    = nodes->GetRhosIndex();
  T_INDEX       = nodes->GetTIndex();
  TVE_INDEX     = nodes->GetTveIndex();
  RHOCVTR_INDEX = nodes->GetRhoCvtrIndex();
  RHOCVVE_INDEX = nodes->GetRhoCvveIndex();

  /*--- Loop over boundary points to calculate energy flux ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      Coord_i = geometry->node[iPoint]->GetCoord();
      Coord_j = geometry->node[jPoint]->GetCoord();

      dij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dij += (Coord_j[iDim] - Coord_i[iDim])*(Coord_j[iDim] - Coord_i[iDim]);
      dij = sqrt(dij);

      /*--- Calculate geometrical parameters ---*/
      theta = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        theta += UnitNormal[iDim]*UnitNormal[iDim];
      }

      /*--- Initialize viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar ++)
        Res_Visc[iVar] = 0.0;

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there is grid motion (v = u_wall)---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      nodes->SetVelocity_Old(iPoint,Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Calculate the gradient of temperature ---*/

      Ti   = nodes->GetTemperature(iPoint);
      Tj   = nodes->GetTemperature(jPoint);
      Tvei = nodes->GetTemperature_ve(iPoint);
      Tvej = nodes->GetTemperature_ve(jPoint);

      /*--- Rename variables for convenience ---*/
      ktr     = nodes->GetThermalConductivity(iPoint);
      kve     = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Apply to the linear system ---*/
      Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                   (ktr*(Twall-Ti) + kve*(Twall-Tvei))*C)*Area/dij;
      Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Twall-Tvei) *C)*Area/dij;

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;

        dTdU   = nodes->GetdTdU(iPoint);
        dTvedU = nodes->GetdTvedU(iPoint);
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_i[nSpecies+nDim][iVar]   = -(ktr*theta/dij*dTdU[iVar] +
                                                kve*theta/dij*dTvedU[iVar])*Area;
          Jacobian_i[nSpecies+nDim+1][iVar] = - kve*theta/dij*dTvedU[iVar]*Area;
        }
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      } // implicit
    }
  }
}

void CTNE2NSSolver::BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
                                                   CSolver **solution_container,
                                                   CNumerics *conv_numerics,
                                                   CNumerics *sour_numerics,
                                                   CConfig *config,
                                                   unsigned short val_marker) {

  /*--- Call standard isothermal BC to apply no-slip and energy b.c.'s ---*/
  BC_Isothermal_Wall(geometry, solution_container, conv_numerics,
                     sour_numerics, config, val_marker);

}

void CTNE2NSSolver::BC_IsothermalCatalytic_Wall(CGeometry *geometry,
                                                CSolver **solution_container,
                                                CNumerics *conv_numerics,
                                                CNumerics *sour_numerics,
                                                CConfig *config,
                                                unsigned short val_marker) {

  /*--- Call standard isothermal BC to apply no-slip and energy b.c.'s ---*/
  BC_Isothermal_Wall(geometry, solution_container, conv_numerics,
                     sour_numerics, config, val_marker);

  ///////////// FINITE DIFFERENCE METHOD ///////////////
  /*--- Local variables ---*/
  bool implicit;
  unsigned short iDim, iSpecies, jSpecies, iVar, jVar, kVar;
  unsigned short RHOS_INDEX, RHO_INDEX, T_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double pcontrol;
  su2double rho, *eves, *hs, RuSI, Ru, *Ms, *xi;
  su2double *dTdU, *dTvedU, *Cvtr, *Cvve;
  su2double *Normal, Area, dij, UnitNormal[3];
  su2double *Di, *Dj, *Vi, *Vj, *Yj, *Yst, *dYdn, SdYdn;
  su2double **GradY;
  su2double **dVdU;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);

  /*--- Set "Proportional control" coefficient ---*/
  pcontrol = 0.6;

  /*--- Get species mass fractions at the wall ---*/
  Yst = config->GetWall_Catalycity();

  /*--- Get the locations of the primitive variables ---*/
  RHOS_INDEX    = nodes->GetRhosIndex();
  RHO_INDEX     = nodes->GetRhoIndex();
  T_INDEX       = nodes ->GetTIndex();

  /*--- Allocate arrays ---*/
  hs    = new su2double[nSpecies];
  Yj    = new su2double[nSpecies];
  dYdn  = new su2double[nSpecies];
  GradY = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    GradY[iSpecies] = new su2double[nDim];
  dVdU = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    dVdU[iVar] = new su2double[nVar];
  Cvtr = new su2double[nSpecies];

  /*--- Get universal information ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  Ms   = config->GetMolar_Mass();
  xi   = config->GetRotationModes();

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      dij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dij += (geometry->node[jPoint]->GetCoord(iDim) -
                geometry->node[iPoint]->GetCoord(iDim))
            * (geometry->node[jPoint]->GetCoord(iDim) -
               geometry->node[iPoint]->GetCoord(iDim));
      }
      dij = sqrt(dij);


      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;


      /*--- Initialize the viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Res_Visc[iVar] = 0.0;

      /*--- Get primitive information ---*/
      Vi = nodes->GetPrimitive(iPoint);
      Vj = nodes->GetPrimitive(jPoint);
      Di = nodes->GetDiffusionCoeff(iPoint);
      Dj = nodes->GetDiffusionCoeff(jPoint);
      eves = nodes->GetEve(iPoint);
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        hs[iSpecies] = nodes->CalcHs(config, Vi[T_INDEX],
                                            eves[iSpecies], iSpecies);
        Yj[iSpecies] = Vj[RHOS_INDEX+iSpecies]/Vj[RHO_INDEX];
      }
      rho    = Vi[RHO_INDEX];
      dTdU   = nodes->GetdTdU(iPoint);
      dTvedU = nodes->GetdTvedU(iPoint);


      /*--- Calculate normal derivative of mass fraction ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dYdn[iSpecies] = (Yst[iSpecies]-Yj[iSpecies])/dij;

      /*--- Calculate supplementary quantities ---*/
      SdYdn = 0.0;
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        SdYdn += rho*Di[iSpecies]*dYdn[iSpecies];

      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Res_Visc[iSpecies]         = -(-rho*Di[iSpecies]*dYdn[iSpecies]
                                       +Yst[iSpecies]*SdYdn            )*Area;
        Res_Visc[nSpecies+nDim]   += (Res_Visc[iSpecies]*hs[iSpecies]  )*Area;
        Res_Visc[nSpecies+nDim+1] += (Res_Visc[iSpecies]*eves[iSpecies])*Area;
      }

      /*--- Viscous contribution to the residual at the wall ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if (implicit) {

        /*--- Initialize the transformation matrix ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++) {
            dVdU[iVar][jVar] = 0.0;
            Jacobian_j[iVar][jVar] = 0.0;
            Jacobian_i[iVar][jVar] = 0.0;
          }

        /*--- Populate transformation matrix ---*/
        // dYsdrk
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
            dVdU[iSpecies][jSpecies] += -1.0/rho*Yst[iSpecies];
          dVdU[iSpecies][iSpecies] += 1.0/rho;
        }
        for (iVar = 0; iVar < nVar; iVar++) {
          dVdU[nSpecies+nDim][iVar]   = dTdU[iVar];
          dVdU[nSpecies+nDim+1][iVar] = dTvedU[iVar];
        }

        /*--- Calculate supplementary quantities ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          Cvtr[iSpecies] = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
        }
        Cvve = nodes->GetCvve(iPoint);

        /*--- Take the primitive var. Jacobian & store in Jac. jj ---*/
        // Species mass fraction
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
            Jacobian_j[iSpecies][jSpecies] += -Yst[iSpecies]*rho*Di[jSpecies]/dij;
          Jacobian_j[iSpecies][iSpecies] += rho*Di[iSpecies]/dij - SdYdn;
        }

        // Temperature
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
            Jacobian_j[nSpecies+nDim][iSpecies] += Jacobian_j[jSpecies][iSpecies]*hs[iSpecies];
          }
          Jacobian_j[nSpecies+nDim][nSpecies+nDim] += Res_Visc[iSpecies]/Area*(Ru/Ms[iSpecies] +
                                                                               Cvtr[iSpecies]  );
          Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += Res_Visc[iSpecies]/Area*Cvve[iSpecies];
        }

        // Vib.-El. Temperature
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
            Jacobian_j[nSpecies+nDim+1][iSpecies] += Jacobian_j[jSpecies][iSpecies]*eves[iSpecies];
          Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += Res_Visc[iSpecies]/Area*Cvve[iSpecies];
        }

        /*--- Multiply by the transformation matrix and store in Jac. ii ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            for (kVar = 0; kVar < nVar; kVar++)
              Jacobian_i[iVar][jVar] += Jacobian_j[iVar][kVar]*dVdU[kVar][jVar]*Area;

        /*--- Apply to the linear system ---*/
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] GradY[iSpecies];
  delete [] GradY;
  delete [] dYdn;
  delete [] hs;
  delete [] Yj;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] dVdU[iVar];
  delete [] dVdU;
  delete [] Cvtr;
}
