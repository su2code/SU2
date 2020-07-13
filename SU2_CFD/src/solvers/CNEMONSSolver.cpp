/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
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

#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/solvers/CNEMONSSolver.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"

CNEMONSSolver::CNEMONSSolver(void) : CNEMOEulerSolver() {

  /*--- Array initialization ---*/
  CD_Visc   = nullptr; CL_Visc   = nullptr; CSF_Visc  = nullptr; CEff_Visc = nullptr;
  CFx_Visc  = nullptr; CFy_Visc  = nullptr; CFz_Visc  = nullptr;
  CMx_Visc  = nullptr; CMy_Visc  = nullptr; CMz_Visc  = nullptr;
  CoPx_Visc = nullptr; CoPy_Visc = nullptr; CoPz_Visc = nullptr;

  ForceViscous = nullptr; MomentViscous = nullptr; CSkinFriction = nullptr;

  Buffet_Sensor = nullptr; Buffet_Metric = nullptr;

  /*--- Surface based array initialization ---*/
  Surface_CL_Visc  = nullptr; Surface_CD_Visc    = nullptr; Surface_CSF_Visc      = nullptr; Surface_CEff_Visc = nullptr;
  Surface_CFx_Visc = nullptr; Surface_CFy_Visc   = nullptr; Surface_CFz_Visc      = nullptr;
  Surface_CMx_Visc = nullptr; Surface_CMy_Visc   = nullptr; Surface_CMz_Visc      = nullptr;
  Surface_HF_Visc  = nullptr; Surface_MaxHF_Visc = nullptr; Surface_Buffet_Metric = nullptr;

  /*--- Rotorcraft simulation array initialization ---*/
  CMerit_Visc = nullptr; CT_Visc = nullptr; CQ_Visc = nullptr;

  /*--- Inlet Variables ---*/
  Inlet_Ttotal = nullptr;
  Inlet_Ptotal = nullptr;
  Inlet_FlowDir = nullptr;

}

CNEMONSSolver::CNEMONSSolver(CGeometry *geometry, CConfig *config,
                             unsigned short iMesh) : CNEMOEulerSolver() {

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
  bool rans = false; //((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_NEMO_RANS));
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
  CD_Visc   = nullptr; CL_Visc   = nullptr; CSF_Visc  = nullptr; CEff_Visc = nullptr;
  CFx_Visc  = nullptr; CFy_Visc  = nullptr; CFz_Visc  = nullptr;
  CMx_Visc  = nullptr; CMy_Visc  = nullptr; CMz_Visc  = nullptr;
  CoPx_Visc = nullptr; CoPy_Visc = nullptr; CoPz_Visc = nullptr;

  Buffet_Sensor = nullptr; Buffet_Metric = nullptr;

  Surface_CL_Visc  = nullptr; Surface_CD_Visc    = nullptr; Surface_CSF_Visc = nullptr; Surface_CEff_Visc = nullptr;
  Surface_CFx_Visc = nullptr; Surface_CFy_Visc   = nullptr; Surface_CFz_Visc = nullptr;
  Surface_CMx_Visc = nullptr; Surface_CMy_Visc   = nullptr; Surface_CMz_Visc = nullptr;
  Surface_HF_Visc  = nullptr; Surface_MaxHF_Visc = nullptr;

  Surface_Buffet_Metric = nullptr;

  CMerit_Visc      = nullptr; CT_Visc      = nullptr; CQ_Visc       = nullptr;
  MaxHF_Visc       = nullptr; ForceViscous = nullptr; MomentViscous = nullptr;
  CSkinFriction    = nullptr; Cauchy_Serie = nullptr; HF_Visc       = nullptr;
  HeatConjugateVar = nullptr;

  /*--- Initialize quantities for the average process for internal flow ---*/
  AverageVelocity 		    = nullptr;
  AverageTurboVelocity 	    = nullptr;
  OldAverageTurboVelocity   = nullptr;
  ExtAverageTurboVelocity   = nullptr;
  AverageFlux 			    = nullptr;
  SpanTotalFlux 		    = nullptr;
  AveragePressure  			= nullptr;
  OldAveragePressure        = nullptr;
  RadialEquilibriumPressure = nullptr;
  ExtAveragePressure  		= nullptr;
  AverageDensity   			= nullptr;
  OldAverageDensity         = nullptr;
  ExtAverageDensity   		= nullptr;
  AverageNu                 = nullptr;
  AverageKine               = nullptr;
  AverageOmega              = nullptr;
  ExtAverageNu              = nullptr;
  ExtAverageKine            = nullptr;
  ExtAverageOmega           = nullptr;

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
  SetNondimensionalization(config, iMesh);

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

  /*--- Define some auxiliary vectors related to the Source term evolution ---*/
  Source      = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Source[iVar] = 0.0;

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
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
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
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the global Jacobian ---*/
    if (rank == MASTER_NODE) cout << "Initialize jacobian structure (NEMO Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No jacobian structure (NEMO Navier-Stokes). MG level: "
           << iMesh <<"." << endl;
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
  Pressure_Inf       = config->GetPressure_FreeStream();
  Temperature_Inf    = config->GetTemperature_FreeStream();
  Temperature_ve_Inf = config->GetTemperature_ve_FreeStream();
  MassFrac_Inf       = config->GetMassFrac_FreeStream();
  Mach_Inf           = config->GetMach();
  Velocity_Inf       = config->GetVelocity_FreeStreamND();
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

  node_infty = new CNEMONSVariable(Pressure_Inf, MassFrac_Inf,
                                   Mvec_Inf, Temperature_Inf,
                                   Temperature_ve_Inf, 1, nDim, nVar,
                                   nPrimVar, nPrimVarGrad, config);


  check_infty = node_infty->SetPrimVar_Compressible(0, config); 

 
  /*--- Initialize the solution to the far-field state everywhere. ---*/
  nodes = new CNEMONSVariable(Pressure_Inf, MassFrac_Inf,
                                       Mvec_Inf, Temperature_Inf,
                                       Temperature_ve_Inf, nPoint, nDim, nVar,
                                       nPrimVar, nPrimVarGrad, config); 

  SetBaseClassPointerToNodes();

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    check = nodes->SetPrimVar_Compressible(iPoint, config); 

    if (check) {

      bool ionization;
      unsigned short iEl, nHeavy, nEl;
      const unsigned short *nElStates;
      su2double RuSI, Ru, T, Tve, rhoCvtr, sqvel, rhoE, rhoEve, num, denom, conc;
      su2double rho, rhos, Ef, Ev, Ee, soundspeed;
      const su2double *xi, *Ms, *thetav, *Tref, *hf;

      /*--- Determine the number of heavy species ---*/
      ionization = config->GetIonization();
      if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
      else            { nHeavy = nSpecies;   nEl = 0; }

      /*--- Load variables from the config class --*/
      xi        = config->GetRotationModes();      // Rotational modes of energy storage
      Ms        = config->GetMolar_Mass();         // Species molar mass
      thetav    = config->GetCharVibTemp();        // Species characteristic vib. temperature [K]
      const auto& thetae    = config->GetCharElTemp();         // Characteristic electron temperature [K]
      const auto& g         = config->GetElDegeneracy();       // Degeneracy of electron states
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

  /* Auxiliary vector for storing primitives for gradient computation in viscous flow */
  /* V = [Y1, ... , Yn, T, Tve, ... ] */
  primitives_aux = new su2double[nPrimVar];

  /*--- Perform the MPI communication of the solution ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Deallocate arrays ---*/
  delete [] Mvec_Inf;

}

CNEMONSSolver::~CNEMONSSolver(void) {
  unsigned short iMarker, iDim;

  unsigned long iVertex;

  if (CD_Visc != nullptr)          delete [] CD_Visc;
  if (CL_Visc != nullptr)          delete [] CL_Visc;
  if (CSF_Visc != nullptr)         delete [] CSF_Visc;
  if (CFx_Visc != nullptr)         delete [] CFx_Visc;
  if (CFy_Visc != nullptr)         delete [] CFy_Visc;
  if (CFz_Visc != nullptr)         delete [] CFz_Visc;
  if (CMx_Visc != nullptr)         delete [] CMx_Visc;
  if (CMy_Visc != nullptr)         delete [] CMy_Visc;
  if (CMz_Visc != nullptr)         delete [] CMz_Visc;
  if (CoPx_Visc != nullptr)        delete [] CoPx_Visc;
  if (CoPy_Visc != nullptr)        delete [] CoPy_Visc;
  if (CoPz_Visc != nullptr)        delete [] CoPz_Visc;
  if (CEff_Visc != nullptr)        delete [] CEff_Visc;
  if (CMerit_Visc != nullptr)      delete [] CMerit_Visc;
  if (Buffet_Metric != nullptr)    delete [] Buffet_Metric;
  if (CT_Visc != nullptr)          delete [] CT_Visc;
  if (CQ_Visc != nullptr)          delete [] CQ_Visc;
  if (HF_Visc != nullptr)          delete [] HF_Visc;
  if (MaxHF_Visc != nullptr)       delete [] MaxHF_Visc;
  if (ForceViscous != nullptr)     delete [] ForceViscous;
  if (MomentViscous != nullptr)    delete [] MomentViscous;

  if (Surface_CL_Visc != nullptr)      delete [] Surface_CL_Visc;
  if (Surface_CD_Visc != nullptr)      delete [] Surface_CD_Visc;
  if (Surface_CSF_Visc != nullptr)     delete [] Surface_CSF_Visc;
  if (Surface_CEff_Visc != nullptr)    delete [] Surface_CEff_Visc;
  if (Surface_CFx_Visc != nullptr)     delete [] Surface_CFx_Visc;
  if (Surface_CFy_Visc != nullptr)     delete [] Surface_CFy_Visc;
  if (Surface_CFz_Visc != nullptr)     delete [] Surface_CFz_Visc;
  if (Surface_CMx_Visc != nullptr)     delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != nullptr)     delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != nullptr)     delete [] Surface_CMz_Visc;
  if (Surface_HF_Visc != nullptr)      delete [] Surface_HF_Visc;
  if (Surface_MaxHF_Visc != nullptr)   delete [] Surface_MaxHF_Visc;
  if (Surface_Buffet_Metric != nullptr) delete [] Surface_Buffet_Metric;

  if (CSkinFriction != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        delete [] CSkinFriction[iMarker][iDim];
      }
      delete [] CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }

  if (HeatConjugateVar != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        delete [] HeatConjugateVar[iMarker][iVertex];
      }
      delete [] HeatConjugateVar[iMarker];
    }
    delete [] HeatConjugateVar;
  }

  if (Buffet_Sensor != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      delete [] Buffet_Sensor[iMarker];
    }
    delete [] Buffet_Sensor;
  }

  if (primitives_aux != nullptr) delete [] primitives_aux;

}

void CNEMONSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh, unsigned short iRKStep,
                                  unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;

  unsigned long InnerIter = config->GetInnerIter();
  bool muscl              = config->GetMUSCL_Flow();
  bool disc_adjoint       = config->GetDiscrete_Adjoint();
  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center             = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  bool center_jst         = (center && config->GetKind_Centered_Flow() == JST);
  bool limiter            = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()));
  bool limiter_turb       = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()));
  bool van_albada         = (config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE);
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
  //InitiateComms(geometry, config, PRIMITIVE);
  //CompleteComms(geometry, config, PRIMITIVE);


  /*--- Artificial dissipation ---*/
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Upwind second order reconstruction ---*/
  if ((muscl && !center) && (iMesh == MESH_0) && !Output) {

    /*--- Calculate the gradients ---*/
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetSolution_Gradient_GG(geometry, config, true);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetSolution_Gradient_LS(geometry, config, true);
    }

    /*--- Limiter computation ---*/
    if ((limiter) && (iMesh == MESH_0) && !Output && !van_albada) {
      SetSolution_Limiter(geometry, config);
    }

  }

 
  /*--- Primitive gradient computation ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/
  //if ((iMesh == MESH_0) && (limiter_flow || limiter_turb)
  //    && !Output && !van_albada) { SetSolution_Limiter(geometry, config); }

  // THIS ISNT SET UP YET
  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  //StrainMag_Max = 0.0; Omega_Max = 0.0;
  //for (iPoint = 0; iPoint < nPoint; iPoint++) {

  //  solver_container[NEMO_SOL]->nodes->SetVorticity();
  //  solver_container[NEMO_SOL]->nodes->SetStrainMag();

  //  StrainMag = solver_container[NEMO_SOL]->nodes->GetStrainMag();
  //  Vorticity = solver_container[NEMO_SOL]->nodes->GetVorticity();
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
    //su2double MyStrainMag_Max = StrainMag_Max; StrainMag_Max = 0.0;

    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    //SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    if (iMesh == MESH_0) {
      config->SetNonphysical_Points(ErrorCounter);
      //solver_container[NEMO_SOL]->SetStrainMag_Max(StrainMag_Max);
      //solver_container[NEMO_SOL]->SetOmega_Max(Omega_Max);
    }
  }

}

void CNEMONSSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction) {

  unsigned long iPoint, iVar;
  unsigned short iSpecies, RHO_INDEX, RHOS_INDEX;
  
  auto& gradient = nodes->GetGradient_Primitive();

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Modify species density to mass concentration ---*/
  for ( iPoint = 0; iPoint < nPoint; iPoint++){
    for( iVar = 0; iVar < nPrimVar; iVar++) {
      primitives_aux[iVar] = nodes->GetPrimitive(iPoint, iVar);
    }
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      primitives_aux[RHOS_INDEX+iSpecies] = primitives_aux[RHOS_INDEX+iSpecies]/primitives_aux[RHO_INDEX];
    for( iVar = 0; iVar < nPrimVar; iVar++)
      nodes->SetPrimitive_Aux(iPoint, iVar, primitives_aux[iVar] );
  }

  const auto& primitives = nodes->GetPrimitive_Aux();

  computeGradientsGreenGauss(this, PRIMITIVE_GRADIENT, PERIODIC_PRIM_GG, *geometry,
                             *config, primitives, 0, nPrimVarGrad, gradient);
}

void CNEMONSSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction) {

  /*--- Set a flag for unweighted or weighted least-squares. ---*/
  bool weighted;

  if (reconstruction)
    weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
  else
    weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);

  
  unsigned long iPoint, iVar;
  unsigned short iSpecies, RHO_INDEX, RHOS_INDEX;

  auto& rmatrix = nodes->GetRmatrix();
  auto& gradient = nodes->GetGradient_Primitive();


  PERIODIC_QUANTITIES kindPeriodicComm = weighted? PERIODIC_PRIM_LS : PERIODIC_PRIM_ULS;

  /*--- Get indices of species & mixture density ---*/
  //RHOS_INDEX = nodes->GetRhosIndex();
  //RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Modify species density to mass concentration ---*/
  //for ( iPoint = 0; iPoint < nPoint; iPoint++){
  //  for( iVar = 0; iVar < nPrimVar; iVar++) {
  //    primitives_aux[iVar] = nodes->GetPrimitive(iPoint, iVar);
  //  }
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    primitives_aux[RHOS_INDEX+iSpecies] = primitives_aux[RHOS_INDEX+iSpecies]/primitives_aux[RHO_INDEX];
  //  for( iVar = 0; iVar < nPrimVar; iVar++)
  //    nodes->SetPrimitive_Aux(iPoint, iVar, primitives_aux[iVar] );
  //}

  //const auto& primitives = nodes->GetPrimitive_Aux();
  const auto& primitives = nodes->GetPrimitive();
  
  computeGradientsLeastSquares(this, PRIMITIVE_GRADIENT, kindPeriodicComm, *geometry, *config,
                               weighted, primitives, 0, nPrimVarGrad, gradient, rmatrix);
}

void CNEMONSSolver::Viscous_Residual(CGeometry *geometry,
                                     CSolver **solution_container,
                                     CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep) {

  bool implicit, err;
  unsigned short iVar, jVar;
  unsigned long iPoint, jPoint, iEdge;

  /*--- Determine time integration scheme ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  CNumerics* numerics = numerics_container[VISC_TERM];

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  //numerics->SetRhosIndex   ( nodes->GetRhosIndex()    );
  //numerics->SetRhoIndex    ( nodes->GetRhoIndex()     );
  //numerics->SetPIndex      ( nodes->GetPIndex()       );
  //numerics->SetTIndex      ( nodes->GetTIndex()       );
  //numerics->SetTveIndex    ( nodes->GetTveIndex()     );
  //numerics->SetVelIndex    ( nodes->GetVelIndex()     );
  //numerics->SetHIndex      ( nodes->GetHIndex()       );
  //numerics->SetAIndex      ( nodes->GetAIndex()       );
  //numerics->SetRhoCvtrIndex( nodes->GetRhoCvtrIndex() );
  //numerics->SetRhoCvveIndex( nodes->GetRhoCvveIndex() );


  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edges->GetNode(iEdge, 0);
    jPoint = geometry->edges->GetNode(iEdge, 1);
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                       geometry->nodes->GetCoord(jPoint) );
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

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

void CNEMONSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {

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
  su2double MyAllBound_CD_Visc, MyAllBound_CL_Visc, MyAllBound_CSF_Visc, MyAllBound_CMx_Visc, MyAllBound_CMy_Visc, MyAllBound_CMz_Visc, MyAllBound_CoPx_Visc, MyAllBound_CoPy_Visc, MyAllBound_CoPz_Visc, MyAllBound_CFx_Visc, MyAllBound_CFy_Visc, MyAllBound_CFz_Visc, MyAllBound_CT_Visc, MyAllBound_CQ_Visc, MyAllBound_HF_Visc, MyAllBound_MaxHF_Visc, *MySurface_CL_Visc = nullptr, *MySurface_CD_Visc = nullptr, *MySurface_CSF_Visc = nullptr, *MySurface_CEff_Visc = nullptr, *MySurface_CFx_Visc = nullptr, *MySurface_CFy_Visc = nullptr, *MySurface_CFz_Visc = nullptr, *MySurface_CMx_Visc = nullptr, *MySurface_CMy_Visc = nullptr, *MySurface_CMz_Visc = nullptr, *MySurface_HF_Visc = nullptr, *MySurface_MaxHF_Visc;
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
    Velocity_Inf[iDim] = node_infty->GetVelocity(0,iDim);
    RefVel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  RefDensity  = node_infty->GetDensity(0);

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
        (Boundary == ISOTHERMAL_NONCATALYTIC) ||
        (Boundary == SMOLUCHOWSKI_MAXWELL)) {

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
        Coord        = geometry->nodes->GetCoord(iPoint);
        Coord_Normal = geometry->nodes->GetCoord(iPointNormal);
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
        if ((geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES)) {

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

void CNEMONSSolver::BC_Sym_Plane(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker) {

  /*--- Call the Euler wall routine ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, visc_numerics,
                config, val_marker);

}

void CNEMONSSolver::BC_HeatFlux_Wall(CGeometry *geometry,
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
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

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
    if (geometry->nodes->GetDomain(iPoint)) {

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

void CNEMONSSolver::BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry,
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
  //	implicit = (config->GetKind_TimeIntScheme_NEMO() == EULER_IMPLICIT);
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

void CNEMONSSolver::BC_HeatFluxCatalytic_Wall(CGeometry *geometry,
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
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
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
    if (geometry->nodes->GetDomain(iPoint)) {

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
      V = nodes->GetPrimitive(iPoint);
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

void CNEMONSSolver::BC_Isothermal_Wall(CGeometry *geometry,
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

  bool implicit   = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
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


    if (geometry->nodes->GetDomain(iPoint)) {


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
      Coord_i = geometry->nodes->GetCoord(iPoint);
      Coord_j = geometry->nodes->GetCoord(jPoint);

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

void CNEMONSSolver::BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
                                                   CSolver **solution_container,
                                                   CNumerics *conv_numerics,
                                                   CNumerics *sour_numerics,
                                                   CConfig *config,
                                                   unsigned short val_marker) {

  /*--- Call standard isothermal BC to apply no-slip and energy b.c.'s ---*/
  BC_Isothermal_Wall(geometry, solution_container, conv_numerics,
                     sour_numerics, config, val_marker);

}

void CNEMONSSolver::BC_IsothermalCatalytic_Wall(CGeometry *geometry,
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
  su2double rho, *eves, *hs, RuSI, Ru;
  su2double *dTdU, *dTvedU, *Cvtr, *Cvve;
  su2double *Normal, Area, dij, UnitNormal[3];
  su2double *Di, *Dj, *Vi, *Vj, *Yj, *dYdn, SdYdn;
  su2double **GradY;
  su2double **dVdU;
  const su2double *Yst, *Ms, *xi;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

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
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      dij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dij += (geometry->nodes->GetCoord(jPoint, iDim) -
                geometry->nodes->GetCoord(iPoint, iDim))
            * (geometry->nodes->GetCoord(jPoint, iDim) -
               geometry->nodes->GetCoord(iPoint, iDim));
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

void CNEMONSSolver::BC_Smoluchowski_Maxwell(CGeometry *geometry,
                                CSolver **solution_container,
                                CNumerics *conv_numerics,
                                CNumerics *visc_numerics,
                                CConfig *config,
                                unsigned short val_marker) {


  unsigned short iDim, jDim, iVar, jVar, iSpecies;
  unsigned short T_INDEX, TVE_INDEX, VEL_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr, kve;
  su2double Ti, Tvei, Tj, Tvej, *dTdU, *dTvedU;
  su2double Twall, Tslip, dij, theta;
  su2double Pi;
  su2double Area, *Normal, UnitNormal[3];
  su2double *Coord_i, *Coord_j;
  su2double C;

  su2double TMAC, TAC;
  su2double Viscosity, Lambda;
  su2double Density, GasConstant, Ru;

  su2double **Grad_PrimVar;
  su2double Vector_Tangent_dT[3], Vector_Tangent_dTve[3], Vector_Tangent_HF[3];
  su2double dTn, dTven;

  su2double TauElem[3], TauTangent[3];
  su2double Tau[3][3];
  su2double TauNormal;
  su2double div_vel=0, Delta;
  
  bool implicit   = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool ionization = config->GetIonization();

  if (ionization) {
    cout << "BC_SMOLUCHOWSKI_MAXWELL: NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION" << endl;
    exit(1);
  }

  /*--- Define 'proportional control' constant ---*/
  C = 5;

  /*--- Identify the boundary ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature and accomodation coefficients---*/
  Twall = config->GetIsothermal_Temperature(Marker_Tag);
  TMAC  = 1.0;
  TAC   = 1.0;

  /*--- Extract necessary indices ---*/
  T_INDEX       = nodes->GetTIndex();
  VEL_INDEX     = nodes->GetVelIndex();
  TVE_INDEX     = nodes->GetTveIndex();

  /*--- Loop over boundary points to calculate energy flux ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->nodes->GetDomain(iPoint)) {

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
      Coord_i = geometry->nodes->GetCoord(iPoint);
      Coord_j = geometry->nodes->GetCoord(jPoint);

      dij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dij += (Coord_j[iDim] - Coord_i[iDim])*(Coord_j[iDim] - Coord_i[iDim]);
      dij = sqrt(dij);

      /*--- Calculate geometrical parameters ---*/
      /*theta = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        theta += UnitNormal[iDim]*UnitNormal[iDim];
      }*/

      /*--- Calculate Pressure ---*/
      Pi   = nodes->GetPressure(iPoint);

      /*--- Calculate the gradient of temperature ---*/
      Ti   = nodes->GetTemperature(iPoint);
      Tj   = nodes->GetTemperature(jPoint);
      Tvei = nodes->GetTemperature_ve(iPoint);
      Tvej = nodes->GetTemperature_ve(jPoint);

      /*--- Rename variables for convenience ---*/
      ktr  = nodes->GetThermalConductivity(iPoint);
      kve  = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Retrieve Flow Data ---*/
      Viscosity = nodes->GetLaminarViscosity(iPoint);
      Density = nodes->GetDensity(iPoint);

      /*--- Calculate specific gas constant --- */
      GasConstant=0;
      for(iSpecies=0;iSpecies<nSpecies;iSpecies++)
        GasConstant+=UNIVERSAL_GAS_CONSTANT*1000.0/config->GetMolar_Mass(iSpecies)*nodes->GetMassFraction(iPoint,iSpecies);
      
      /*--- Calculate temperature gradients normal to surface---*/ //Doubt about minus sign
      dTn   = - (Ti-Tj)/dij;
      dTven = - (Tvei-Tvej)/dij;

      /*--- Calculate molecular mean free path ---*/
      Lambda = Viscosity/Density*sqrt(PI_NUMBER/(2*GasConstant*Ti));

      /*--- Calculate Temperature Slip ---*/
      Tslip = ((2-TAC)/TAC)*2*Gamma/(Gamma+1)/Prandtl_Lam*Lambda*dTn+Twall;

      /*--- Retrieve Primitive Gradients ---*/
      Grad_PrimVar = nodes->GetGradient_Primitive(iPoint);

      /*--- Calculate temperature gradients tangent to surface ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_Tangent_dT[iDim]   = Grad_PrimVar[T_INDEX][iDim] - dTn * UnitNormal[iDim];
        Vector_Tangent_dTve[iDim] = Grad_PrimVar[TVE_INDEX][iDim] - dTven * UnitNormal[iDim];
      }

      /*--- Calculate Heatflux tangent to surface ---*/
      for (iDim = 0; iDim < nDim; iDim++) 
        Vector_Tangent_HF[iDim] = ktr*Vector_Tangent_dT[iDim]+kve*Vector_Tangent_dTve[iDim];
      
      /*--- Initialize viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar ++)
        Res_Visc[iVar] = 0.0;

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
      for (iDim = 0; iDim < nDim; iDim++) {
        TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
      }

      /*--- Store the Slip Velocity at the wall */
      for (iDim = 0; iDim < nDim; iDim++)
        Vector[iDim] = - Lambda/Viscosity*(2-TMAC)/TMAC*(TauTangent[iDim])-3/4*(Gamma-1)/Gamma*Prandtl_Lam/Pi*Vector_Tangent_HF[iDim];

      nodes->SetVelocity_Old(iPoint,Vector);

      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Apply to the linear system ---*/
      Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                   (ktr*(Tslip-Ti) + kve*(Tslip-Tvei))*C)*Area/dij;
      Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Tslip-Tvei)*C)*Area/dij;

      LinSysRes.SubtractBlock(iPoint, Res_Visc);
    }
  } 
}