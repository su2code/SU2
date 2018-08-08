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
#include <math.h>

CTNE2EulerSolver::CTNE2EulerSolver(void) : CSolver() {

	/*--- Array initialization ---*/
  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt = NULL; CL_Mnt = NULL; CSF_Mnt = NULL;  CEff_Mnt = NULL;
  CMx_Mnt = NULL; CMy_Mnt = NULL; CMz_Mnt = NULL;
  CFx_Mnt = NULL; CFy_Mnt = NULL; CFz_Mnt = NULL;
  CoPx_Mnt = NULL; CoPy_Mnt = NULL; CoPz_Mnt = NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL; MomentMomentum = NULL;

	/*--- Surface based array initialization ---*/
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt = NULL; Surface_CD_Mnt = NULL; Surface_CSF_Mnt = NULL; Surface_CEff_Mnt = NULL;
  Surface_CFx_Mnt = NULL; Surface_CFy_Mnt = NULL; Surface_CFz_Mnt = NULL;
  Surface_CMx_Mnt = NULL; Surface_CMy_Mnt = NULL; Surface_CMz_Mnt = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

	 /*--- Numerical methods array initialization ---*/
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;

  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;

  Secondary = NULL; Secondary_i = NULL; Secondary_j = NULL;
}

CTNE2EulerSolver::CTNE2EulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {

	unsigned long iPoint, counter_local = 0, counter_global = 0, iVertex;
	unsigned short iVar, iDim, iMarker, iSpecies, nLineLets;
	su2double StaticEnergy, Density, Velocity2, Pressure, Temperature;
	unsigned short nZone = geometry->GetnZone();
	unsigned short iZone = config->GetiZone();
	bool restart   = (config->GetRestart() || config->GetRestart_Flow());
	bool rans = ((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS));
	unsigned short direct_diff = config->GetDirectDiff();
	string filename_ = config->GetSolution_FlowFileName();

	su2double *Mvec_Inf;
  su2double Alpha, Beta;
	bool check_infty, nonPhys;


	/*--- Check for a restart file to evaluate if there is a change in the AoA
		before non-dimensionalizing ---*/
	if (!(!restart || (iMesh != MESH_0) || nZone >1 )) {

		/*--- Multizone problems require number of the zone to be appended ---*/
		if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone);

		/*--- Read and store the restart metadata ---*/
		Read_SU2_Restart_Metadata(geometry, config, false, filename_);

	}

	/*--- Array initialization ---*/

	/*--- Basic array initialization ---*/
  CD_Inv = NULL; CL_Inv = NULL; CSF_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  CoPx_Inv = NULL; CoPy_Inv = NULL; CoPz_Inv = NULL;

  CD_Mnt= NULL; CL_Mnt= NULL; CSF_Mnt= NULL; CEff_Mnt= NULL;
  CMx_Mnt= NULL;   CMy_Mnt= NULL;   CMz_Mnt= NULL;
  CFx_Mnt= NULL;   CFy_Mnt= NULL;   CFz_Mnt= NULL;
  CoPx_Mnt= NULL;   CoPy_Mnt= NULL;   CoPz_Mnt= NULL;

  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  ForceMomentum = NULL;  MomentMomentum = NULL;

  /*--- Surface based array initialization ---*/
  Surface_CL_Inv = NULL; Surface_CD_Inv = NULL; Surface_CSF_Inv = NULL; Surface_CEff_Inv = NULL;
  Surface_CFx_Inv = NULL; Surface_CFy_Inv = NULL; Surface_CFz_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CL_Mnt= NULL; Surface_CD_Mnt= NULL; Surface_CSF_Mnt= NULL; Surface_CEff_Mnt= NULL;
  Surface_CFx_Mnt= NULL;   Surface_CFy_Mnt= NULL;   Surface_CFz_Mnt= NULL;
  Surface_CMx_Mnt= NULL;   Surface_CMy_Mnt= NULL;   Surface_CMz_Mnt = NULL;

  Surface_CL = NULL; Surface_CD = NULL; Surface_CSF = NULL; Surface_CEff = NULL;
  Surface_CFx = NULL; Surface_CFy = NULL; Surface_CFz = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Numerical methods array initialization ---*/
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  DonorPrimVar = NULL; DonorGlobalIndex = NULL;
  ActDisk_DeltaP = NULL; ActDisk_DeltaT = NULL;

  Smatrix = NULL; Cvector = NULL;

  Secondary=NULL; Secondary_i=NULL; Secondary_j=NULL;

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
	//nSecondaryVar     = ????;
	//nSecondaryVarGrad = ????;

	/*--- Initialize nVarGrad for deallocation ---*/
	nVarGrad     = nPrimVarGrad;
	nPoint       = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();

	/*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  SetNondimensionalization(geometry, config, iMesh);

	/*--- Store the number of vertices on each marker for deallocation ---*/
	nVertex = new unsigned long[nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		nVertex[iMarker] = geometry->nVertex[iMarker];

	/*--- Allocate a CVariable array for each node of the mesh ---*/
	node = new CVariable*[nPoint];

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

		if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
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
	ForceInviscid    = new su2double[nDim];
	MomentInviscid   = new su2double[3];
	CD_Inv           = new su2double[nMarker];
	CL_Inv           = new su2double[nMarker];
	CSF_Inv          = new su2double[nMarker];
	CMx_Inv          = new su2double[nMarker];
	CMy_Inv          = new su2double[nMarker];
	CMz_Inv          = new su2double[nMarker];
	CEff_Inv         = new su2double[nMarker];
	CFx_Inv          = new su2double[nMarker];
	CFy_Inv          = new su2double[nMarker];
	CFz_Inv          = new su2double[nMarker];

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
	Temperature_ve_Inf = config->GetTemperature_ve_FreeStream();
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
  node_infty = new CTNE2EulerVariable(Pressure_Inf, MassFrac_Inf,
                                      Mvec_Inf, Temperature_Inf,
                                      Temperature_ve_Inf, nDim, nVar,
                                      nPrimVar, nPrimVarGrad, config);
  check_infty = node_infty->SetPrimVar_Compressible(config);

	/*--- Initialize the solution to the far-field state everywhere. ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CTNE2EulerVariable(Pressure_Inf, MassFrac_Inf, Mvec_Inf, Temperature_Inf,
                                            Temperature_ve_Inf, nDim, nVar, nPrimVar, nPrimVarGrad,
                                            config);

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    nonPhys = node[iPoint]->SetPrimVar_Compressible(config);

    if (nonPhys) {
      bool ionization;
      unsigned short iEl, nHeavy, nEl, *nElStates;
      su2double Ru, T, Tve, rhoCvtr, sqvel, rhoE, rhoEve, num, denom, conc;
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
      Ru      = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(kmol*K)]
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
      for (iDim = 0; iDim < nDim; iDim++)
        sqvel += Mvec_Inf[iDim]*soundspeed * Mvec_Inf[iDim]*soundspeed;

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
        Solution[iSpecies]     = rho*MassFrac_Inf[iSpecies];
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        Solution[nSpecies+iDim]     = rho*Mvec_Inf[iDim]*soundspeed;
      }
      Solution[nSpecies+nDim]       = rhoE;
      Solution[nSpecies+nDim+1]     = rhoEve;

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

  /*--- Define solver parameters needed for execution of destructor ---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;


	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

  /*--- Deallocate arrays ---*/
  delete [] Mvec_Inf;
}

CTNE2EulerSolver::~CTNE2EulerSolver(void) {
	unsigned short iVar, iMarker;

	/*--- Array deallocation ---*/
  if (Velocity_Inf != NULL)     delete [] Velocity_Inf;
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

  if (lowerlimit != NULL)  			delete [] lowerlimit;
  if (upperlimit != NULL) 			delete [] upperlimit;

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
			delete CPressure[iMarker];
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

void CTNE2EulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
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
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;

        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;

        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
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
      Buffer_Send_U    = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sent---*/
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
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_U[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_U[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_U[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  rotMatrix[3][3], *angles, *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;

  su2double *Limiter = new su2double[nVar];

  #ifndef HAVE_MPI
    int send_to, receive_from;
    SU2_MPI::Status status;
  #endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
      (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

  #ifdef HAVE_MPI
    send_to      = config->GetMarker_All_SendRecv(MarkerS) -1;
    receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
  #endif

      nVertexS = geometry->nVertex[MarkerS]; nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;       nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit    = new su2double[nBufferS_Vector];

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
        //1st column
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        //2nd column
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        //3rd column
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex];
          Limiter[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Limiter[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
          Limiter[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
          Limiter[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_Limit[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_Limit[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_Limit[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Primitive(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR, VEL_INDEX;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
   *Buffer_Receive_V = NULL, *Buffer_Send_V = NULL;

  su2double *Primitive;

  #ifndef HAVE_MPI
    int send_to, receive_from;
    SU2_MPI::Status status;
  #endif

  Primitive = new su2double[nPrimVar];
  VEL_INDEX = node_infty->GetVelIndex();

	for (iMarker = 0; iMarker < nMarker; iMarker++) {

		if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

		MarkerS = iMarker;  MarkerR = iMarker+1;

  #ifdef HAVE_MPI
    send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
		receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
  #endif

		nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
		nBufferS_Vector = nVertexS*nPrimVar;    nBufferR_Vector = nVertexR*nPrimVar;

    /*--- Allocate Receive and send buffers  ---*/
    Buffer_Receive_V = new su2double [nBufferR_Vector];
    Buffer_Send_V    = new su2double [nBufferS_Vector];

    /*--- Copy the solution that should be sended ---*/
    for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVar; iVar++)
          Buffer_Send_V[iVar*nVertexS+iVertex] = node[iPoint]->GetPrimVar(iVar);
    }

  #ifdef HAVE_MPI

    /*--- Send/Receive information using Sendrecv ---*/
    SU2_MPI::Sendrecv(Buffer_Send_V, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                      Buffer_Receive_V, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

  #else

      /*--- Receive information without MPI ---*/
    for (iVertex = 0; iVertex < nVertexR; iVertex++) {
      for (iVar = 0; iVar < nVar; iVar++)
        Buffer_Receive_V[iVar*nVertexR+iVertex] = Buffer_Send_V[iVar*nVertexR+iVertex];

    }

  #endif

    /*--- Deallocate send buffer ---*/
    delete [] Buffer_Send_V;

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
      for (iVar = 0; iVar < nPrimVar; iVar++)
        Primitive[iVar] = Buffer_Receive_V[iVar*nVertexR+iVertex];

      /*--- Rotate the velocity components. ---*/
      if (nDim == 2) {
        Primitive[VEL_INDEX]   = rotMatrix[0][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex];
        Primitive[VEL_INDEX+1] = rotMatrix[1][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex];
      }
      else {
        Primitive[VEL_INDEX]   = rotMatrix[0][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_V[(VEL_INDEX+2)*nVertexR+iVertex];
        Primitive[VEL_INDEX+1] = rotMatrix[1][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_V[(VEL_INDEX+2)*nVertexR+iVertex];
        Primitive[VEL_INDEX+2] = rotMatrix[2][0]*Buffer_Receive_V[(VEL_INDEX)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_V[(VEL_INDEX+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_V[(VEL_INDEX+2)*nVertexR+iVertex];
      }

      /*--- Copy transformed conserved variables back into buffer. ---*/
      for (iVar = 0; iVar < nPrimVar; iVar++)
        node[iPoint]->SetPrimVar(iVar, Primitive[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_V;
    }
	}
  delete [] Primitive;
}

void CTNE2EulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;

  #ifndef HAVE_MPI
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
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex];
        }
        else {
          Solution[nSpecies]   = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+1] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
          Solution[nSpecies+2] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[(nSpecies)*nVertexR+iVertex]
                               + rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[(nSpecies+1)*nVertexR+iVertex]
                               + rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[(nSpecies+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_MaxEigenvalue(CGeometry *geometry, CConfig *config) {
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

void CTNE2EulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
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

      /*--- Copy the solution old that should be sent ---*/
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
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Primitive_Gradient(CGeometry *geometry, CConfig *config) {
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
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex]
                              + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
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

void CTNE2EulerSolver::Set_MPI_Primitive_Limiter(CGeometry *geometry,
                                                 CConfig *config) {

  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned short VEL_INDEX;
  unsigned long iVertex, iPoint;
  unsigned long nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  rotMatrix[3][3], *angles, *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;

  /*--- Initialize array to store the limiter ---*/
  su2double *Limiter = new su2double [nPrimVarGrad];

  #ifdef HAVE_MPI
    int send_to, receive_from;
    SU2_MPI::Status status;
  #endif

  /*--- Get the position of the velocity terms in the primitive vector ---*/
  VEL_INDEX = node_infty->GetVelIndex();

  /*--- Loop over all send/receive boundaries ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

  #ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
  #endif

			nVertexS = geometry->nVertex[MarkerS];
      nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad;
      nBufferR_Vector = nVertexR*nPrimVarGrad;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double[nBufferR_Vector];
      Buffer_Send_Limit    = new su2double[nBufferS_Vector];

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
        iPoint          = geometry->vertex[MarkerR][iVertex]->GetNode();
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
        //1st column
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        //2nd column
        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        //3rd column
        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Limiter[VEL_INDEX]   = rotMatrix[0][0]*Buffer_Receive_Limit[(VEL_INDEX  )*nVertexR+iVertex] +
                                 rotMatrix[0][1]*Buffer_Receive_Limit[(VEL_INDEX+1)*nVertexR+iVertex];
          Limiter[VEL_INDEX+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(VEL_INDEX  )*nVertexR+iVertex] +
                                 rotMatrix[1][1]*Buffer_Receive_Limit[(VEL_INDEX+1)*nVertexR+iVertex];
        }
        else {
          Limiter[VEL_INDEX]   = rotMatrix[0][0]*Buffer_Receive_Limit[(VEL_INDEX  )*nVertexR+iVertex] +
                                 rotMatrix[0][1]*Buffer_Receive_Limit[(VEL_INDEX+1)*nVertexR+iVertex] +
                                 rotMatrix[0][2]*Buffer_Receive_Limit[(VEL_INDEX+2)*nVertexR+iVertex];
          Limiter[VEL_INDEX+1] = rotMatrix[1][0]*Buffer_Receive_Limit[(VEL_INDEX  )*nVertexR+iVertex] +
                                 rotMatrix[1][1]*Buffer_Receive_Limit[(VEL_INDEX+1)*nVertexR+iVertex] +
                                 rotMatrix[1][2]*Buffer_Receive_Limit[(VEL_INDEX+2)*nVertexR+iVertex];
          Limiter[VEL_INDEX+2] = rotMatrix[2][0]*Buffer_Receive_Limit[(VEL_INDEX  )*nVertexR+iVertex] +
                                 rotMatrix[2][1]*Buffer_Receive_Limit[(VEL_INDEX+1)*nVertexR+iVertex] +
                                 rotMatrix[2][2]*Buffer_Receive_Limit[(VEL_INDEX+2)*nVertexR+iVertex];
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

void CTNE2EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solution_container,CConfig *config,
                                     unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long ErrorCounter = 0;

	unsigned long ExtIter = config->GetExtIter();
  bool cont_adjoint     = config->GetContinuous_Adjoint();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
	bool implicit         = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool muscl            = (config->GetMUSCL_Flow() || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
	bool limiter          = ((config->GetKind_SlopeLimit_TNE2() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) &&
													!(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool center       		= ((config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED) ||
                      		(cont_adjoint && config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED));
	bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
	bool nonPhys;
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
	bool interface        = (config->GetnMarker_InterfaceBound() != 0);

	/* --- Compute interface MPI --- */
	if (interface) {Set_MPI_Interface(geometry, config); }

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
    if (limiter && (iMesh == MESH_0) && !Output && !van_albada)
				{ SetPrimitive_Limiter(geometry, config); }
  }

	/*--- Artificial dissipation ---*/
  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

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

void CTNE2EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solution_container, CConfig *config,
                                      unsigned short iMesh, unsigned long Iteration) {

	su2double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel=0.0, Lambda, Local_Delta_Time,
	Global_Delta_Time = 1E6, Global_Delta_UnstTimeND;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;

	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement = config->GetGrid_Movement();
  bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;
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
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);

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
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

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

}

void CTNE2EulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {

	su2double *Normal, Area, Mean_SoundSpeed, Mean_ProjVel, Lambda;
	unsigned long iEdge, iVertex, iPoint, jPoint;
	unsigned short iDim, iMarker;

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetLambda(0.0);

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		Normal = geometry->edge[iEdge]->GetNormal();
		Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

		/*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

		/*--- Inviscid contribution ---*/
		Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
		if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
		if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);

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
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

      /*--- Inviscid contribution ---*/
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddLambda(Lambda);
      }
    }
  }

  /*--- Call the MPI routine ---*/
  Set_MPI_MaxEigenvalue(geometry, config);

}

void CTNE2EulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, jVar;
  bool err;

  /*--- Set booleans based on config settings ---*/
	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool centered = ((config->GetKind_Centered_TNE2() == JST) && (iMesh == MESH_0));

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );

	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge, set normal vectors, and number of neighbors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(),
                          geometry->node[jPoint]->GetnNeighbor());

		/*--- Pass conservative & primitive variables w/o reconstruction to CNumerics ---*/
		numerics->SetConservative(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetPrimitive(node[iPoint]->GetPrimVar(), node[jPoint]->GetPrimVar());

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(node[iPoint]->GetdPdU(), node[jPoint]->GetdPdU());
    numerics->SetdTdU(node[iPoint]->GetdTdU(), node[jPoint]->GetdTdU());
    numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node[jPoint]->GetdTvedU());

    /*--- Set the largest convective eigenvalue ---*/
		numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());

		/*--- Compute residuals, and Jacobians ---*/
		numerics->ComputeResidual(Res_Conv, Res_Visc, Jacobian_i, Jacobian_j, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if ((Res_Conv[iVar] != Res_Conv[iVar]) ||
          (Res_Visc[iVar] != Res_Visc[iVar])   )
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
      LinSysRes.AddBlock(iPoint, Res_Visc);
      LinSysRes.SubtractBlock(jPoint, Res_Visc);
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
	unsigned long iEdge, iPoint, jPoint;
	unsigned long ExtIter = config->GetExtIter();
  unsigned short RHO_INDEX, RHOS_INDEX, P_INDEX, TVE_INDEX;
  unsigned short iDim, iSpecies, iVar, jVar;

	su2double *U_i, *U_j, *V_i, *V_j;
  su2double **GradU_i, **GradU_j, ProjGradU_i, ProjGradU_j;
  su2double **GradV_i, **GradV_j, ProjGradV_i, ProjGradV_j;
  su2double *Limiter_i, *Limiter_j;
  su2double *Conserved_i, *Conserved_j, *Primitive_i, *Primitive_j;
  su2double *dPdU_i, *dPdU_j, *dTdU_i, *dTdU_j, *dTvedU_i, *dTvedU_j;
  su2double *Eve_i, *Eve_j, *Cvve_i, *Cvve_j;
  su2double counter_local=0, counter_global=0;
  su2double lim_i, lim_j, lim_ij;

	/*--- Set booleans based on config settings ---*/
	bool implicit 		= (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool muscl    		= (config->GetMUSCL_Flow() && (iMesh == MESH_0));
	bool disc_adjoint = config->GetDiscrete_Adjoint();
  bool limiter      = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) &&
											!(disc_adjoint && config->GetFrozen_Limiter_Disc()));
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
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );

  RHO_INDEX  = node[0]->GetRhoIndex();
  RHOS_INDEX = node[0]->GetRhosIndex();
  P_INDEX    = node[0]->GetPIndex();
  TVE_INDEX  = node[0]->GetTveIndex();

  /*--- Loop over edges and calculate convective fluxes ---*/
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Retrieve node numbers and pass edge normal to CNumerics ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Get conserved & primitive variables from CVariable ---*/
    U_i = node[iPoint]->GetSolution();  U_j = node[jPoint]->GetSolution();
    V_i = node[iPoint]->GetPrimVar();   V_j = node[jPoint]->GetPrimVar();
    //S_i = node[iPoint]->GetSecondary(); S_j = node[jPoint]->GetSecondary();

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
      GradU_i = node[iPoint]->GetGradient();
      GradU_j = node[jPoint]->GetGradient();
      GradV_i = node[iPoint]->GetGradient_Primitive();
      GradV_j = node[jPoint]->GetGradient_Primitive();

			if (limiter) {
        Limiter_i = node[iPoint]->GetLimiter();
        Limiter_j = node[jPoint]->GetLimiter();
      }

      /*--- Reconstruct conserved variables at the edge interface ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        ProjGradU_i = 0.0; ProjGradU_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          ProjGradU_i += Vector_i[iDim]*GradU_i[iVar][iDim];
          ProjGradU_j += Vector_j[iDim]*GradU_j[iVar][iDim];
        }
        if (limiter) {
          Conserved_i[iVar] = U_i[iVar] + lim_i*ProjGradU_i;
          Conserved_j[iVar] = U_j[iVar] + lim_j*ProjGradU_j;
        }
        else {
          Conserved_i[iVar] = U_i[iVar] + ProjGradU_i;
          Conserved_j[iVar] = U_j[iVar] + ProjGradU_j;
        }
      }

      /*--- Calculate corresponding primitive reconstructed variables ---*/
  //      for (iVar = 0; iVar < nPrimVar; iVar++) {
  //        ProjGradV_i = 0.0; ProjGradV_j = 0.0;
  //        for (iDim = 0; iDim < nDim; iDim++) {
  //          ProjGradV_i += Vector_i[iDim]*GradV_i[iVar][iDim];
  //          ProjGradV_j += Vector_j[iDim]*GradV_j[iVar][iDim];
  //        }
  //        Primitive_i[iVar] = V_i[iVar] + lim_ij*ProjGradV_i;
  //        Primitive_j[iVar] = V_j[iVar] + lim_ij*ProjGradV_j;
  //
  //        // Vib.-el. energy
  //        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //          Eve_i[iSpecies] = node[iPoint]->CalcEve(config, Primitive_i[TVE_INDEX], iSpecies);
  //          Eve_j[iSpecies] = node[jPoint]->CalcEve(config, Primitive_j[TVE_INDEX], iSpecies);
  //          Cvve_i[iSpecies] = node[iPoint]->CalcCvve(Primitive_i[TVE_INDEX], config, iSpecies);
  //          Cvve_j[iSpecies] = node[jPoint]->CalcCvve(Primitive_j[TVE_INDEX], config, iSpecies);
  //        }
  //
  //        // Recalculate derivatives of pressure
  //        // NOTE: We need to pass the vib-el. energy, but it is only used when
  //        //       ionized species are present, so, for now, we just load it with
  //        //       the value at i or j, and come back later when ionized is ready.
  //        node[iPoint]->CalcdPdU(Primitive_i, Eve_i, config, dPdU_i);
  //        node[jPoint]->CalcdPdU(Primitive_j, Eve_j, config, dPdU_j);
  //
  //        // Recalculate temperature derivatives
  //        node[iPoint]->CalcdTdU(Primitive_i, config, dTdU_i);
  //        node[jPoint]->CalcdTdU(Primitive_j, config, dTdU_j);
  //
  //        // Recalculate Tve derivatives
  //        // Note: Species vib.-el. energies are required for species density
  //        //       terms.  For now, just pass the values at i and j and hope it works
  //        node[iPoint]->CalcdTvedU(Primitive_i, Eve_i, config, dTvedU_i);
  //        node[jPoint]->CalcdTvedU(Primitive_j, Eve_j, config, dTvedU_j);
  //      }


      chk_err_i = node[iPoint]->Cons2PrimVar(config, Conserved_i, Primitive_i,
                                             dPdU_i, dTdU_i, dTvedU_i, Eve_i, Cvve_i);
      chk_err_j = node[jPoint]->Cons2PrimVar(config, Conserved_j, Primitive_j,
                                             dPdU_j, dTdU_j, dTvedU_j, Eve_j, Cvve_j);


      /*--- Check for physical solutions in the reconstructed values ---*/
      // Note: If non-physical, revert to first order
      if ( chk_err_i || chk_err_j) {
        numerics->SetPrimitive   (V_i, V_j);
        numerics->SetConservative(U_i, U_j);
        numerics->SetdPdU  (node[iPoint]->GetdPdU(),   node[jPoint]->GetdPdU());
        numerics->SetdTdU  (node[iPoint]->GetdTdU(),   node[jPoint]->GetdTdU());
        numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node[jPoint]->GetdTvedU());
        numerics->SetEve   (node[iPoint]->GetEve(),    node[jPoint]->GetEve());
        numerics->SetCvve  (node[iPoint]->GetCvve(),   node[jPoint]->GetCvve());
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
      numerics->SetdPdU  (node[iPoint]->GetdPdU(),   node[jPoint]->GetdPdU());
      numerics->SetdTdU  (node[iPoint]->GetdTdU(),   node[jPoint]->GetdTdU());
      numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node[jPoint]->GetdTvedU());
      numerics->SetEve   (node[iPoint]->GetEve(),    node[jPoint]->GetEve());
      numerics->SetCvve  (node[iPoint]->GetCvve(),   node[jPoint]->GetCvve());
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
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
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
  delete [] Eve_i;
  delete [] Eve_j;
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
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );


	/*--- Initialize the source residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

  /*--- loop over interior points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Set conserved & primitive variables  ---*/
    numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
    numerics->SetPrimitive   (node[iPoint]->GetPrimVar(),  node[iPoint]->GetPrimVar() );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU(node[iPoint]->GetdPdU(), node[iPoint]->GetdPdU());
    numerics->SetdTdU(node[iPoint]->GetdTdU(), node[iPoint]->GetdTdU());
    numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node[iPoint]->GetdTvedU());
    numerics->SetEve(node[iPoint]->GetEve(), node[iPoint]->GetEve());
    numerics->SetCvve(node[iPoint]->GetCvve(), node[iPoint]->GetCvve());

    /*--- Set volume of the dual grid cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[iPoint]->GetCoord() );


    /*--- Compute axisymmetric source terms (if needed) ---*/
    if (config->GetAxisymmetric()) {
      numerics->ComputeAxisymmetric(Residual, Jacobian_i, config);

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
    numerics->ComputeChemistry(Residual, Jacobian_i, config);

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
    numerics->ComputeVibRelaxation(Residual, Jacobian_i, config);

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

  }

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
}

void CTNE2EulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {

	unsigned long iVertex, iPoint;
	unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
	su2double Pressure, *Normal = NULL, *Coord, Area, factor, NFPressOF,
	RefVel2, RefDensity, RefPressure;

  su2double MomentDist[3];
  su2double Force[3];
  su2double UnitNormal[3];
	su2double MomentX_Force[3] = {0.0,0.0,0.0},
						MomentY_Force[3] = {0.0,0.0,0.0},
						MomentZ_Force[3] = {0.0,0.0,0.0};
	su2double AxiFactor;
  string Marker_Tag, Monitoring_Tag;
	bool grid_movement        = config->GetGrid_Movement();
  bool axisymmetric         = config->GetAxisymmetric();

  #ifdef HAVE_MPI
  su2double MyAllBound_CD_Inv, MyAllBound_CL_Inv, MyAllBound_CSF_Inv,
	MyAllBound_CMx_Inv, MyAllBound_CMy_Inv, MyAllBound_CMz_Inv, MyAllBound_CoPx_Inv,
	MyAllBound_CoPy_Inv, MyAllBound_CoPz_Inv, MyAllBound_CFx_Inv, MyAllBound_CFy_Inv,
	MyAllBound_CFz_Inv, MyAllBound_CT_Inv, MyAllBound_CQ_Inv, MyAllBound_CNearFieldOF_Inv,
	*MySurface_CL_Inv = NULL, *MySurface_CD_Inv = NULL, *MySurface_CSF_Inv = NULL,
	*MySurface_CEff_Inv = NULL, *MySurface_CFx_Inv = NULL, *MySurface_CFy_Inv = NULL,
	*MySurface_CFz_Inv = NULL, *MySurface_CMx_Inv = NULL, *MySurface_CMy_Inv = NULL,
	*MySurface_CMz_Inv = NULL;
  #endif

	su2double Alpha           = config->GetAoA()*PI_NUMBER/180.0;
	su2double Beta            = config->GetAoS()*PI_NUMBER/180.0;
	su2double RefArea         = config->GetRefArea();
	su2double RefLength       = config->GetRefLength();
	su2double *Origin         = NULL;
  if (config->GetnMarker_Monitoring() != 0){
    Origin = config->GetRefOriginMoment(0);
  }

	/*--- Evaluate reference values for non-dimensionalization. ---*/
  RefVel2     = node_infty->GetVelocity2();
	RefDensity  = node_infty->GetDensity();
	RefPressure = node_infty->GetPressure();
  factor      = 1.0 / (0.5*RefDensity*RefArea*RefVel2);


	/*-- Initialization ---*/
  Total_CD = 0.0;           Total_CL = 0.0;    Total_CSF = 0.0;     Total_CEff = 0.0;
  Total_CMx = 0.0;          Total_CMy = 0.0;   Total_CMz = 0.0;
  Total_CoPx = 0.0;         Total_CoPy = 0.0;  Total_CoPz = 0.0;
  Total_CFx = 0.0;          Total_CFy = 0.0;   Total_CFz = 0.0;
  Total_CT = 0.0;           Total_CQ = 0.0;    Total_CMerit = 0.0;
  Total_CNearFieldOF = 0.0; Total_Heat = 0.0;  Total_MaxHeat = 0.0;

  AllBound_CD_Inv = 0.0;           AllBound_CL_Inv = 0.0; AllBound_CSF_Inv = 0.0;
  AllBound_CMx_Inv = 0.0;          AllBound_CMy_Inv = 0.0;   AllBound_CMz_Inv = 0.0;
  AllBound_CoPx_Inv = 0.0;         AllBound_CoPy_Inv = 0.0;  AllBound_CoPz_Inv = 0.0;
  AllBound_CFx_Inv = 0.0;          AllBound_CFy_Inv = 0.0;   AllBound_CFz_Inv = 0.0;
  AllBound_CT_Inv = 0.0;           AllBound_CQ_Inv = 0.0;    AllBound_CMerit_Inv = 0.0;
  AllBound_CNearFieldOF_Inv = 0.0; AllBound_CEff_Inv = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CL_Inv[iMarker_Monitoring]      = 0.0; Surface_CD_Inv[iMarker_Monitoring]      = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring] = 0.0; Surface_CEff_Inv[iMarker_Monitoring]       = 0.0;
    Surface_CFx_Inv[iMarker_Monitoring]        = 0.0; Surface_CFy_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CFz_Inv[iMarker_Monitoring]        = 0.0; Surface_CMx_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]        = 0.0; Surface_CMz_Inv[iMarker_Monitoring]        = 0.0;
    Surface_CL[iMarker_Monitoring]          = 0.0; Surface_CD[iMarker_Monitoring]          = 0.0;
    Surface_CSF[iMarker_Monitoring]     = 0.0; Surface_CEff[iMarker_Monitoring]           = 0.0;
    Surface_CFx[iMarker_Monitoring]            = 0.0; Surface_CFy[iMarker_Monitoring]            = 0.0;
    Surface_CFz[iMarker_Monitoring]            = 0.0; Surface_CMx[iMarker_Monitoring]            = 0.0;
    Surface_CMy[iMarker_Monitoring]            = 0.0; Surface_CMz[iMarker_Monitoring]            = 0.0;
  }
	/*--- Loop over the Euler and Navier-Stokes markers ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		Boundary   = config->GetMarker_All_KindBC(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);

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
			MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
      MomentX_Force[0] = 0.0; MomentX_Force[1] = 0.0; MomentX_Force[2] = 0.0;
      MomentY_Force[0] = 0.0; MomentY_Force[1] = 0.0; MomentY_Force[2] = 0.0;
      MomentZ_Force[0] = 0.0; MomentZ_Force[1] = 0.0; MomentZ_Force[2] = 0.0;

			NFPressOF = 0.0;


      /*--- Loop over vertices to compute forces ---*/
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Calculate pressure ---*/
				Pressure = node[iPoint]->GetPressure();
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

          /*--- Geometry terms ---*/
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) {
            UnitNormal[iDim] = Normal[iDim]/Area;
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

					/*--- Axisymmetric simulations ---*/
          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          /*--- Compute force, note minus sign due to outward orientation of
           the normal vector ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf)*Normal[iDim]*factor;
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

          AllBound_CD_Inv           += CD_Inv[iMarker];
          AllBound_CL_Inv           += CL_Inv[iMarker];
          AllBound_CSF_Inv          += CSF_Inv[iMarker];
          AllBound_CEff_Inv          = AllBound_CL_Inv / (AllBound_CD_Inv + EPS);
          AllBound_CMx_Inv          += CMx_Inv[iMarker];
          AllBound_CMy_Inv          += CMy_Inv[iMarker];
          AllBound_CMz_Inv          += CMz_Inv[iMarker];
          AllBound_CoPx_Inv         += CoPx_Inv[iMarker];
          AllBound_CoPy_Inv         += CoPy_Inv[iMarker];
          AllBound_CoPz_Inv         += CoPz_Inv[iMarker];
          AllBound_CFx_Inv          += CFx_Inv[iMarker];
          AllBound_CFy_Inv          += CFy_Inv[iMarker];
          AllBound_CFz_Inv          += CFz_Inv[iMarker];
          AllBound_CT_Inv           += CT_Inv[iMarker];
          AllBound_CQ_Inv           += CQ_Inv[iMarker];
          AllBound_CMerit_Inv        = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);

          /*--- Compute the coefficients per surface ---*/

          for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
            Marker_Tag = config->GetMarker_All_TagBound(iMarker);
            if (Marker_Tag == Monitoring_Tag) {
              Surface_CL_Inv[iMarker_Monitoring]      += CL_Inv[iMarker];
              Surface_CD_Inv[iMarker_Monitoring]      += CD_Inv[iMarker];
              Surface_CSF_Inv[iMarker_Monitoring]     += CSF_Inv[iMarker];
              Surface_CEff_Inv[iMarker_Monitoring]     = CL_Inv[iMarker] / (CD_Inv[iMarker] + EPS);
              Surface_CFx_Inv[iMarker_Monitoring]     += CFx_Inv[iMarker];
              Surface_CFy_Inv[iMarker_Monitoring]     += CFy_Inv[iMarker];
              Surface_CFz_Inv[iMarker_Monitoring]     += CFz_Inv[iMarker];
              Surface_CMx_Inv[iMarker_Monitoring]     += CMx_Inv[iMarker];
              Surface_CMy_Inv[iMarker_Monitoring]     += CMy_Inv[iMarker];
              Surface_CMz_Inv[iMarker_Monitoring]     += CMz_Inv[iMarker];
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
  MyAllBound_CNearFieldOF_Inv = AllBound_CNearFieldOF_Inv; AllBound_CNearFieldOF_Inv = 0.0;

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

    Surface_CL_Inv[iMarker_Monitoring]         = 0.0;
    Surface_CD_Inv[iMarker_Monitoring]         = 0.0;
    Surface_CSF_Inv[iMarker_Monitoring]        = 0.0;
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
  Total_CFx           = AllBound_CFx_Inv;
  Total_CFy           = AllBound_CFy_Inv;
  Total_CFz           = AllBound_CFz_Inv;
  Total_CMx           = AllBound_CMx_Inv;
  Total_CMy           = AllBound_CMy_Inv;
  Total_CMz           = AllBound_CMz_Inv;
  Total_CoPx          = AllBound_CoPx_Inv;
  Total_CoPy          = AllBound_CoPy_Inv;
  Total_CoPz          = AllBound_CoPz_Inv;
  Total_CT            = AllBound_CT_Inv;
  Total_CQ            = AllBound_CQ_Inv;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);
  Total_CNearFieldOF  = AllBound_CNearFieldOF_Inv;

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
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = local_Residual[iVar] + local_Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry-> node[iPoint]->GetGlobalIndex());
      }
    }

  }

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);

}

void CTNE2EulerSolver::ExplicitRK_Iteration(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CConfig *config,
                                            unsigned short iRKStep) {

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
    Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry-> node[iPoint]->GetGlobalIndex());
      }
    }

  }

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);


}

void CTNE2EulerSolver::ImplicitEuler_Iteration(CGeometry *geometry,
                                               CSolver **solution_container,
                                               CConfig *config) {
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
		local_Res_TruncError = node[iPoint]->GetResTruncError();

		/*--- Read the volume ---*/
		Vol = geometry-> node[iPoint]->GetVolume();

		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
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

void CTNE2EulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint, iEdge, iVertex;
	unsigned short iDim, iSpecies, iVar, iMarker, RHOS_INDEX, RHO_INDEX;
	su2double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average, rho_i, rho_j,
	Partial_Gradient, Partial_Res, *Normal;

	/*--- Initialize arrays ---*/

	/*--- Primitive variables: [Y1, ..., YNs, T, Tve, u, v, w]^T ---*/
	PrimVar_Vertex = new su2double [nPrimVarGrad];
	PrimVar_i = new su2double [nPrimVarGrad];
	PrimVar_j = new su2double [nPrimVarGrad];

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = node[0]->GetRhosIndex();
  RHO_INDEX  = node[0]->GetRhoIndex();

	/*--- Set Gradient_Primitive to zero ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

	/*--- Loop interior edges ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Pull primitives from CVariable ---*/
		for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
			PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);
			PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iVar);
		}

		Normal = geometry-> edge[iEdge]->GetNormal();
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

        /*--- Get primitives from CVariable ---*/
				for (iVar = 0; iVar < nPrimVarGrad; iVar++)
					PrimVar_Vertex[iVar] = node[iPoint]->GetPrimVar(iVar);

        /*--- Modify species density to mass concentration ---*/
        rho_i = node[iPoint]->GetPrimVar(RHO_INDEX);
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          PrimVar_Vertex[RHOS_INDEX+iSpecies] = PrimVar_Vertex[RHOS_INDEX+iSpecies]/rho_i;

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
				Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar,iDim) / (geometry->node[iPoint]->GetVolume());
				node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);

			}
		}
	}

	delete [] PrimVar_Vertex;
	delete [] PrimVar_i;
	delete [] PrimVar_j;

	Set_MPI_Primitive_Gradient(geometry, config);

}

void CTNE2EulerSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {

	unsigned short iSpecies, iVar, iDim, jDim, iNeigh, RHOS_INDEX, RHO_INDEX;
	unsigned long iPoint, jPoint;
	su2double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
	r23_b, r33, rho_i, rho_j, weight, product, detR2, z11, z12, z13, z22, z23, z33;
  bool singular;

	/*--- Initialize arrays, Primitive variables:
   [Y1, ..., YNs, T, Tve, u, v, w, P]^T ---*/
	PrimVar_i = new su2double [nPrimVarGrad];
	PrimVar_j = new su2double [nPrimVarGrad];

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = node[0]->GetRhosIndex();
  RHO_INDEX  = node[0]->GetRhoIndex();

	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Set the value of singulare ---*/
		singular = false;

    /*--- Get coordinates ---*/
		Coord_i = geometry->node[iPoint]->GetCoord();

    /*--- Get primitives from CVariable ---*/
		for (iVar = 0; iVar < nPrimVarGrad; iVar++)
			PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);

		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nPrimVarGrad; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				Cvector[iVar][iDim] = 0.0;

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		AD::StartPreacc();
    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    AD::SetPreaccIn(Coord_i, nDim);

		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();

			for (iVar = 0; iVar < nPrimVarGrad; iVar++)
				PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iVar);

			AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

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
				node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
			}
		}

		AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    AD::EndPreacc();
	}

	delete [] PrimVar_i;
	delete [] PrimVar_j;

	Set_MPI_Primitive_Gradient(geometry, config);

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
  RHOS_INDEX = node[0]->GetRhosIndex();
  RHO_INDEX  = node[0]->GetRhoIndex();

  iPoint = val_Point;
  Coord_i = geometry->node[iPoint]->GetCoord();

  /*--- Get primitives from CVariable ---*/
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);


  /*--- Inizialization of variables ---*/
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Cvector[iVar][iDim] = 0.0;

  r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
  r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

  AD::StartPreacc();
  AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
  AD::SetPreaccIn(Coord_i, nDim);

  for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
    jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
    Coord_j = geometry->node[jPoint]->GetCoord();

    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iVar);

		AD::SetPreaccIn(Coord_j, nDim);
    AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

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
        product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
      node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
    }
  }

	AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
  AD::EndPreacc();

	delete [] PrimVar_i;
	delete [] PrimVar_j;

	Set_MPI_Primitive_Gradient(geometry, config);
}

void CTNE2EulerSolver::SetPrimitive_Gradient(CConfig *config) {
  unsigned long iPoint;
  unsigned short iVar;
  su2double *U, *V, **GradU, **GradV;

	/*--- Loop over all points ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    U = node[iPoint]->GetSolution();
    V = node[iPoint]->GetPrimVar();
    GradU = node[iPoint]->GetGradient();
    GradV = node[iPoint]->GetGradient_Primitive();

    node[iPoint]->GradCons2GradPrimVar(config, U, V, GradU, GradV);
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
				node[iPoint]->SetLimiter_Primitive(iVar,1.0);
			}
		}
	}
	else {
		/*--- Initialize solution max, solution min and limiter in entire domain ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
				node[iPoint]->SetSolution_Max(iVar, -EPS);
				node[iPoint]->SetSolution_Min(iVar, EPS);
				node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
			}
		}

		/*--- Establish bounts for Spekreijse monotonicity by finding max/min values
		of neighbor variables ---*/
		for (iEdge = 0; iEdge < geometry-> GetnEdge(); iEdge++) {

			/*--- Point identification, Normal vector and area ---*/
			iPoint = geometry-> edge[iEdge]->GetNode(0);
			jPoint = geometry-> edge[iEdge]->GetNode(1);

			/*--- Get primitive variables ---*/
			Primitive_i = node[iPoint]->GetPrimVar();
			Primitive_j = node[jPoint]->GetPrimVar();

			/*--- Compute the max and min values for nodes i & j ---*/
			for (iVar = 0; iVar < nVar; iVar++){
				du = (Primitive_j[iVar]-Primitive_i[iVar]);
				node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
        node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
        node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
        node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
			}
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

        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
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

        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }

      }

      AD::EndPreacc();

    }

    delete [] LocalMinPrimitive; delete [] GlobalMinPrimitive;
    delete [] LocalMaxPrimitive; delete [] GlobalMaxPrimitive;

  }

  /*--- Limiter MPI ---*/
  Set_MPI_Primitive_Limiter(geometry, config);

}

void CTNE2EulerSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) {
	unsigned short iDim, jDim, iVar, jVar;
	su2double Beta, local_Mach, Beta2, rho, enthalpy, soundspeed, sq_vel;
	su2double *U_i = NULL;
	su2double Beta_max = config->GetmaxTurkelBeta();
  su2double Mach_infty2, Mach_lim2, aux, parameter;

	/*--- Variables to calculate the preconditioner parameter Beta ---*/
	local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();

  /*--- Weiss and Smith Preconditionng ---*/
	Mach_infty2 = pow(config->GetMach(),2.0);
	Mach_lim2   = pow(0.00001,2.0);
	aux         = max(pow(local_Mach,2.0),Mach_lim2);
	parameter   = min(1.0, max(aux,Beta_max*Mach_infty2));

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

void CTNE2EulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, jDim, iSpecies, iVar, jVar, kVar;
	unsigned long iPoint, iVertex;

  su2double *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3], *NormalArea,
  ProjGridVel = 0.0, turb_ke;
  su2double Density_b, StaticEnergy_b, Enthalpy_b, *Velocity_b, Kappa_b, Chi_b, Energy_b, VelMagnitude2_b, Pressure_b;
  su2double Density_i, *Velocity_i, ProjVelocity_i = 0.0, Energy_i, VelMagnitude2_i;
  su2double **Jacobian_b, **DubDu;
	su2double rhoCvtr, rhoCvve, rho_el, Ru;
  su2double rho, cs, P, rhoE, rhoEve, conc, Beta;
  su2double *u, *Ms, *dPdU;

	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));

  /*--- Allocate arrays ---*/
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

  /*--- Load parameters from the config class ---*/
  Ms = config->GetMolar_Mass();

  /*--- Rename for convenience ---*/
  Ru = UNIVERSAL_GAS_CONSTANT;

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
				NormalArea[iDim] = -Normal[iDim];
				UnitNormal[iDim] = -Normal[iDim]/Area;
			}

			/*--- Retrieve the pressure on the vertex ---*/
      P   = node[iPoint]->GetPressure();

			/*--- Compute the boundary state b ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity_b[iDim] = Velocity_i[iDim] - ProjVelocity_i * UnitNormal[iDim];

			VelMagnitude2_b = 0.0;
			for (iDim = 0; iDim <nDim; iDim++)
				VelMagnitude2_b += Velocity_b[iDim] * Velocity_b[iDim];


			/*--- Compute the residual ---*/
			turb_ke = 0.0;
			if (tkeNeeded) turb_ke=solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);

      /*--- Apply the flow-tangency b.c. to the convective flux ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        Residual[iSpecies] = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[nSpecies+iDim] = P * UnitNormal[iDim] * Area;
      Residual[nSpecies+nDim]   = 0.0;
			Residual[nSpecies+nDim+1] = 0.0;

			/*--- Add the Reynolds stress tensor contribution ---*/
      if (tkeNeeded) {
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[nSpecies+iDim+1] += (2.0/3.0)*Density_b*turb_ke*NormalArea[iDim];
      }
			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- If using implicit time-stepping, calculate b.c. contribution to Jacobian ---*/
			if (implicit) {

				/*--- Initialize Jacobian ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;

				/*--- Calculate state i ---*/
        rho     = node[iPoint]->GetDensity();
        rhoCvtr = node[iPoint]->GetRhoCv_tr();
        rhoCvve = node[iPoint]->GetRhoCv_ve();
        rhoE    = node[iPoint]->GetSolution(nSpecies+nDim);
        rhoEve  = node[iPoint]->GetSolution(nSpecies+nDim+1);
        dPdU    = node[iPoint]->GetdPdU();
        for (iDim = 0; iDim < nDim; iDim++)
          u[iDim] = node[iPoint]->GetVelocity(iDim);

        /*--- If free electrons are present, retrieve the electron gas density ---*/
        if (config->GetIonization()) rho_el = node[iPoint]->GetMassFraction(nSpecies-1) * rho;
        else                         rho_el = 0.0;

        conc = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          cs    = node[iPoint]->GetMassFraction(iSpecies);
          conc += cs * rho/Ms[iSpecies];

          /////// NEW //////
          for (iDim = 0; iDim < nDim; iDim++) {
            Jacobian_i[nSpecies+iDim][iSpecies] = dPdU[iSpecies] * UnitNormal[iDim];
            Jacobian_i[iSpecies][nSpecies+iDim] = cs * UnitNormal[iDim];
          }

         Jacobian_i[iSpecies][nSpecies+4] = 0.0;
        }

        Beta = Ru*conc/rhoCvtr;

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            Jacobian_i[nSpecies+iDim][nSpecies+jDim] = u[iDim]*UnitNormal[jDim]
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
        Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);

			}
    }
	}
  delete [] UnitNormal;
  delete [] u;
}

void CTNE2EulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solution_container,
                                    CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
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
  su2double *U_domain,*U_infty;
  su2double Gas_Constant     = config->GetGas_ConstantND();

	/*--- Getting info from config ---*/
	unsigned short VEL_INDEX     = node[0]->GetVelIndex();
	unsigned short RHOCVTR_INDEX = node[0]->GetRhoCvtrIndex() ;
	unsigned short RHO_INDEX     = node[0]->GetRhoIndex();
	unsigned short A_INDEX       = node[0]->GetAIndex();


	/*--- Set booleans from configuration parameters ---*/
  bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool viscous  = config->GetViscous();
  bool tkeNeeded = (((config->GetKind_Solver() == RANS ) ||
                   (config->GetKind_Solver() == DISC_ADJ_RANS))
                   && (config->GetKind_Turb_Model() == SST));

	/*--- Allocate arrays ---*/
	su2double *Normal = new su2double[nDim];

  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  conv_numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  conv_numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  conv_numerics->SetPIndex      ( node[0]->GetPIndex()       );
  conv_numerics->SetTIndex      ( node[0]->GetTIndex()       );
  conv_numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  conv_numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  conv_numerics->SetHIndex      ( node[0]->GetHIndex()       );
  conv_numerics->SetAIndex      ( node[0]->GetAIndex()       );
  conv_numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  conv_numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );

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
      U_domain = node[iPoint]->GetSolution();
      V_domain = node[iPoint]->GetPrimVar();
      U_infty  = node_infty->GetSolution();
      V_infty  = node_infty->GetPrimVar();

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
			for (iDim = 0; iDim < nDim; iDim++){
				Vel_Bound[iDim] = V_domain[VEL_INDEX+iDim];
				Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
			}
			Entropy_Bound    = V_domain[RHOCVTR_INDEX]/V_domain[RHO_INDEX];

			/*--- Store the primitive variable state for the freestream. Project
      	the freestream velocity vector into the local normal direction,
      	i.e. compute v_infty.n. ---*/
      Vn_Infty = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Infty[iDim] = V_infty[VEL_INDEX+iDim];
        Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
      }
			Qn_Infty         = Vn_Infty;
			SoundSpeed_Infty = V_infty[A_INDEX];
			Entropy_Infty    = V_infty[RHOCVTR_INDEX]/V_infty[RHO_INDEX];

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


  // Might need to deal with these Gamma_Minus_Ones business..................
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
        V_infty[nSpecies+iDim+1] = Velocity[iDim];
      V_infty[nSpecies+nDim+1] = Pressure;
      V_infty[nSpecies+nDim+2] = Density;
      V_infty[nSpecies+nDim+3] = Energy + Pressure/Density;

			/*--- Pass conserved & primitive variables to CNumerics ---*/
      conv_numerics->SetConservative(U_domain, U_infty);
      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Pass supplementary information to CNumerics ---*/
      conv_numerics->SetdPdU(node[iPoint]->GetdPdU(), node_infty->GetdPdU());
      conv_numerics->SetdTdU(node[iPoint]->GetdTdU(), node_infty->GetdTdU());
      conv_numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node_infty->GetdTvedU());

			/*--- Compute the convective residual (and Jacobian) ---*/
      // Note: This uses the specified boundary num. method specified in definition_structure.cpp
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Apply contribution to the linear system ---*/
      LinSysRes.AddBlock(iPoint, Residual);
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());

			/*--- Viscous contribution ---*/
			if (viscous) {
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord() );
        visc_numerics->SetNormal(Normal);

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetConservative(node[iPoint]->GetSolution(),
                                       node_infty->GetSolution() );
        visc_numerics->SetConsVarGradient(node[iPoint]->GetGradient(),
                                          node_infty->GetGradient() );
        visc_numerics->SetPrimitive(node[iPoint]->GetPrimVar(),
                                    node_infty->GetPrimVar() );
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                          node_infty->GetGradient_Primitive() );

        /*--- Pass supplementary information to CNumerics ---*/
        visc_numerics->SetdPdU(node[iPoint]->GetdPdU(), node_infty->GetdPdU());
        visc_numerics->SetdTdU(node[iPoint]->GetdTdU(), node_infty->GetdTdU());
        visc_numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node_infty->GetdTvedU());

        /*--- Species diffusion coefficients ---*/
        visc_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(),
                                         node_infty->GetDiffusionCoeff() );

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(),
                                           node_infty->GetLaminarViscosity() );

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(),
                                              node_infty->GetThermalConductivity());

        /*--- Vib-el. thermal conductivity ---*/
        visc_numerics->SetThermalConductivity_ve(node[iPoint]->GetThermalConductivity_ve(),
                                                 node_infty->GetThermalConductivity_ve() );

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
	unsigned short iVar, iDim;
	unsigned long iVertex, iPoint, Point_Normal;
	su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
	Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
	alpha, aa, bb, cc, dd, Area, UnitaryNormal[3];

	bool implicit             = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement        = config->GetGrid_Movement();
	su2double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
	su2double Gas_Constant       = config->GetGas_ConstantND();
	unsigned short Kind_Inlet = config->GetKind_Inlet();
	string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
	bool tkeNeeded = (((config->GetKind_Solver() == RANS )|| (config->GetKind_Solver() == DISC_ADJ_RANS)) &&
                    (config->GetKind_Turb_Model() == SST));
  bool gravity = (config->GetGravityForce());
  bool viscous              = config->GetViscous();

	su2double *U_domain = new su2double[nVar];      su2double *U_inlet = new su2double[nVar];
	su2double *V_domain = new su2double[nPrimVar];  su2double *V_inlet = new su2double[nPrimVar];
	su2double *Normal = new su2double[nDim];
  su2double UnitNormal[3];

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
			for (iVar = 0; iVar < nVar; iVar++)     U_domain[iVar] = node[iPoint]->GetSolution(iVar);
      for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);

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
          Density = V_domain[0];
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = U_domain[iDim+1]/Density;
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Energy      = U_domain[nVar-1]/Density;
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

          /*--- Static pressure using isentropic relation at a point ---*/
          Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);

          /*--- Density at the inlet from the gas law ---*/
          Density = Pressure/(Gas_Constant*Temperature);

          /*--- Using pressure, density, & velocity, compute the energy ---*/
          Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;

          /*--- Conservative variables, using the derived quantities ---*/
          U_inlet[0] = Density;
          for (iDim = 0; iDim < nDim; iDim++)
            U_inlet[iDim+1] = Velocity[iDim]*Density;
          U_inlet[nDim+1] = Energy*Density;

          /*--- Primitive variables, using the derived quantities ---*/
          V_inlet[0] = Temperature;
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Velocity[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;

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
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

			/*--- Roe Turkel preconditioning, set the value of beta ---*/
			if (config->GetKind_Upwind() == TURKEL)
				node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());

			/*--- Viscous contribution ---*/
			if (viscous) {

				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_inlet);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

				/*--- Laminar viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

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

void CTNE2EulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solution_container,
                                   CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iVar, iDim;
	unsigned long iVertex, iPoint, Point_Normal;
	su2double Pressure, P_Exit, Velocity[3],
	Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
	Area, UnitaryNormal[3];

	bool implicit           = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  su2double Gas_Constant     = config->GetGas_ConstantND();
	bool grid_movement      = config->GetGrid_Movement();
	string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
	bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());

	su2double *U_domain = new su2double[nVar];      su2double *U_outlet = new su2double[nVar];
  su2double *V_domain = new su2double[nPrimVar];  su2double *V_outlet = new su2double[nPrimVar];
	su2double *Normal = new su2double[nDim];

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
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
      for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);

			/*--- Build the fictitious intlet state based on characteristics ---*/

      /*--- Retrieve the specified back pressure for this outlet. ---*/
      if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;
      else P_Exit = config->GetOutlet_Pressure(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
       of boundary update depends on this. ---*/
      Density = U_domain[0];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = U_domain[iDim+1]/Density;
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitaryNormal[iDim];
      }
      Energy     = U_domain[nVar-1]/Density;
      Pressure   = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Mach_Exit  = sqrt(Velocity2)/SoundSpeed;

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
         tangential velocity components, are extrapolated. Adapted from an
         original implementation in the Stanford University multi-block
         (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
         der Weide, last modified 09-10-2007. ---*/

        Entropy = Pressure*pow(1.0/Density,Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        Pressure   = P_Exit;
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitaryNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy  = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;

        /*--- Conservative variables, using the derived quantities ---*/
        U_outlet[0] = Density;
        for (iDim = 0; iDim < nDim; iDim++)
          U_outlet[iDim+1] = Velocity[iDim]*Density;
        U_outlet[nDim+1] = Energy*Density;

        /*--- Conservative variables, using the derived quantities ---*/
        V_outlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_outlet[iDim+1] = Velocity[iDim];
        V_outlet[nDim+1] = Pressure;
        V_outlet[nDim+2] = Density;

			}

			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetConservative(U_domain, U_outlet);

			if (grid_movement)
				conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

			/*--- Jacobian contribution for implicit integration ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

			/*--- Roe Turkel preconditioning, set the value of beta ---*/
			if (config->GetKind_Upwind() == TURKEL)
				node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());

			/*--- Viscous contribution ---*/
			if (viscous) {

				/*--- Set the normal vector and the coordinates ---*/
				visc_numerics->SetNormal(Normal);
				visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

				/*--- Primitive variables, and gradient ---*/
				visc_numerics->SetPrimitive(V_domain, V_outlet);
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

				/*--- Laminar viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

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
	delete [] U_outlet;
  delete [] V_domain;
	delete [] V_outlet;
	delete [] Normal;

}

void CTNE2EulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solution_container,
                                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned short iDim, iVar;
	unsigned long iVertex, iPoint, Point_Normal;
	su2double Density, Pressure, Temperature, Energy, *Velocity, Velocity2;
	su2double Gas_Constant = config->GetGas_ConstantND();

	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool grid_movement  = config->GetGrid_Movement();
	bool viscous              = config->GetViscous();
	string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double *U_inlet = new su2double[nVar]; su2double *U_domain = new su2double[nVar];
  su2double *V_inlet = new su2double[nPrimVar]; su2double *V_domain = new su2double[nPrimVar];
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

	/*--- Conservative variables, using the derived quantities ---*/
	U_inlet[0] = Density;
  for (iDim = 0; iDim < nDim; iDim++)
    U_inlet[iDim+1] = Velocity[iDim]*Density;
  U_inlet[nDim+1] = Energy*Density;

  /*--- Primitive variables, using the derived quantities ---*/
	V_inlet[0] = Temperature;
  for (iDim = 0; iDim < nDim; iDim++)
    V_inlet[iDim+1] = Velocity[iDim];
  V_inlet[nDim+1] = Pressure;
  V_inlet[nDim+2] = Density;

	/*--- Loop over all the vertices on this boundary marker ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

			/*--- Current solution at this boundary node ---*/
			for (iVar = 0; iVar < nVar; iVar++) U_domain[iVar] = node[iPoint]->GetSolution(iVar);
			for (iVar = 0; iVar < nPrimVar; iVar++) V_domain[iVar] = node[iPoint]->GetPrimVar(iVar);

			/*--- Normal vector for this vertex (negate for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

			double Area = 0.0;su2double UnitaryNormal[3];
			for (iDim = 0; iDim < nDim; iDim++)
				Area += Normal[iDim]*Normal[iDim];
			Area = sqrt (Area);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = Normal[iDim]/Area;

			/*--- Set various quantities in the solver class ---*/
			conv_numerics->SetNormal(Normal);
			conv_numerics->SetConservative(U_domain, U_inlet);

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
				visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

				/*--- Laminar viscosity ---*/
				visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

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

void CTNE2EulerSolver::BC_Sym_Plane(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *conv_numerics,
                                    CNumerics *visc_numerics, CConfig *config,
                                    unsigned short val_marker) {
  //  bool implicit;
  //  unsigned short iDim;
  //  unsigned long iPoint, iVertex;
  //  su2double *Normal, *UnitNormal, Area;
  //
  //  /*--- Allocate arrays ---*/
  //  UnitNormal = new su2double[3];
  //
  //  /*--- Set booleans based on configuration settings ---*/
  //  implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);

  /*--- Call the Euler wall routine ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);

  //  /*--- Compute the viscous contribution (if any) ---*/
  //	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
  //		if (geometry->node[iPoint]->GetDomain()) {
  //      /*--- Calculate parameters from the geometry ---*/
  //      // Note: The vertex normal points out of the geometry by convention,
  //      //       so to calculate the influence from the boundary condition
  //      //       to the domain, we negate this vector
  //      Area   = 0.0;
  //			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
  //			for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  //			Area = sqrt (Area);
  //			for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
  //
  //    }
  //  }
  //
  //
  //  delete [] UnitNormal;

}

void CTNE2EulerSolver::SetResidual_DualTime(CGeometry *geometry,
                                            CSolver **solution_container,
                                            CConfig *config,
                                            unsigned short iRKStep,
                                            unsigned short iMesh,
                                            unsigned short RunTime_EqSystem) {
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1;
	double Volume_nM1, Volume_n, Volume_nP1, TimeStep;

	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool Grid_Movement = config->GetGrid_Movement();

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n   = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();

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
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
                          +  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
		}

		/*--- Add Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual);

		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar++) {
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;

				if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
					Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
				if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
					Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
			}
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}

}

void CTNE2EulerSolver::GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone) {

  unsigned short iVar;
	int rank = MASTER_NODE;
  #ifndef NO_MPI
  #ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  #else
	rank = MPI::COMM_WORLD.Get_rank();
  #endif
  #endif

	/*--- Restart the solution from file information ---*/
	string restart_filename = config->GetSolution_FlowFileName();
	unsigned long iPoint, index, nFlowIter, adjIter, flowIter;
	char buffer[50];
	string UnstExt, text_line;
	ifstream restart_file;
	bool grid_movement = config->GetGrid_Movement();
	unsigned short nZone = geometry->GetnZone();

	/*--- Multi-zone restart files. ---*/
	if (nZone > 1 && !(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)) {
		restart_filename.erase(restart_filename.end()-4, restart_filename.end());
		sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	}

	/*--- For the unsteady adjoint, we integrate backwards through
   physical time, so load in the direct solution files in reverse. ---*/
	if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST) {
		flowIter = val_iZone;
		restart_filename.erase(restart_filename.end()-4, restart_filename.end());
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.dat", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		nFlowIter = config->GetnExtIter() - 1;
		adjIter   = config->GetExtIter();
		flowIter  = nFlowIter - adjIter;
		restart_filename.erase (restart_filename.end()-4, restart_filename.end());
		if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
		if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
		if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
		if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
		if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	} else {
		flowIter = config->GetExtIter();
		restart_filename.erase (restart_filename.end()-4, restart_filename.end());
		if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
		if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
		if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
		if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
		if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
		UnstExt = string(buffer);
		restart_filename.append(UnstExt);
	}

	/*--- Open the flow solution from the restart file ---*/
	if (rank == MASTER_NODE && val_iZone == ZONE_0)
		cout << "Reading in the direct flow solution from iteration " << flowIter << "." << endl;
	restart_file.open(restart_filename.data(), ios::in);

	/*--- In case there is no file ---*/
	if (restart_file.fail()) {
		cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
		exit(1);
	}

	/*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
	long *Global2Local = NULL;
	Global2Local = new long[geometry->GetGlobal_nPointDomain()];
	/*--- First, set all indices to a negative value by default ---*/
	for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
		Global2Local[iPoint] = -1;
	}

	/*--- Now fill array with the transform values only for local points ---*/
	for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
	}

	/*--- Read all lines in the restart file ---*/
	long iPoint_Local = 0; unsigned long iPoint_Global = 0;

	/*--- The first line is the header ---*/
	getline (restart_file, text_line);

	while (getline (restart_file,text_line)) {
		istringstream point_line(text_line);

		/*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1, as
     initialized above. Otherwise, the local index for this node on the
     current processor will be returned and used to instantiate the vars. ---*/
		iPoint_Local = Global2Local[iPoint_Global];
		if (iPoint_Local >= 0) {

      /*--- First value is the point index, then the conservative variables ---*/
      point_line >> index;

      for (iVar = 0; iVar < nVar; iVar++)
        point_line >> Solution[iVar];

			node[iPoint_Local]->SetSolution(Solution);

			/*--- If necessary, read in the grid velocities for the unsteady adjoint ---*/
			if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady() && grid_movement) {
				double Volume, GridVel[3];
				if (nDim == 2) point_line >> Volume >> GridVel[0] >> GridVel[1];
				if (nDim == 3) point_line >> Volume >> GridVel[0] >> GridVel[1] >> GridVel[2];
				if (iPoint_Local >= 0)
					for (unsigned short iDim = 0; iDim < nDim; iDim++)
						geometry->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
			}

		}
		iPoint_Global++;
	}


	/*--- Close the restart file ---*/
	restart_file.close();

	/*--- Free memory needed for the transformation ---*/
	delete [] Global2Local;

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

CTNE2NSSolver::CTNE2NSSolver(void) : CTNE2EulerSolver() {

	/*--- Array initialization ---*/
	CD_Visc = NULL;
	CL_Visc = NULL;
	CMx_Visc = NULL;
	CMy_Visc = NULL;
	CMz_Visc = NULL;
	CFx_Visc = NULL;
	CFy_Visc = NULL;
	CFz_Visc = NULL;
	CEff_Visc = NULL;

	ForceViscous = NULL;
	MomentViscous = NULL;
	CSkinFriction = NULL;

}

CTNE2NSSolver::CTNE2NSSolver(CGeometry *geometry, CConfig *config,
                             unsigned short iMesh) : CTNE2EulerSolver() {
  bool restart, check_infty, check;
  unsigned short iDim, iMarker, iSpecies, iVar, nZone, nLineLets;
	unsigned long iPoint, index, counter_local, counter_global;
  su2double *Mvec_Inf, Alpha, Beta, dull_val;

  /*--- Get MPI rank ---*/
	int rank = MASTER_NODE;
  #ifndef NO_MPI
  #ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  #else
	rank = MPI::COMM_WORLD.Get_rank();
  #endif
  #endif

	/*--- Array initialization ---*/
	CD_Visc       = NULL;
	CL_Visc       = NULL;
	CMx_Visc         = NULL;
	CMy_Visc         = NULL;
	CMz_Visc         = NULL;
	CFx_Visc         = NULL;
	CFy_Visc         = NULL;
	CFz_Visc         = NULL;
	CEff_Visc        = NULL;
  Heat_Visc        = NULL;
  MaxHeatFlux_Visc = NULL;
	ForceViscous     = NULL;
	MomentViscous    = NULL;
	CSkinFriction    = NULL;
  YPlus            = NULL;

  /*--- Initialize counters ---*/
  counter_local  = 0;
  counter_global = 0;

  /*--- Set booleans from settings in CConfig ---*/
  restart = (config->GetRestart() || config->GetRestart_Flow());

	/*--- Define geometry constants in the solver structure ---*/
  nSpecies     = config->GetnSpecies();
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();
	nDim         = geometry->GetnDim();
  nZone        = geometry->GetnZone();

  /*--- Set the size of the primitive and conserve vectors ---*/
  //     U: [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  //     V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradV: [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
	nVar         = nSpecies+nDim+2;
  nPrimVar     = nSpecies+nDim+8;
  nPrimVarGrad = nSpecies+nDim+8;

	/*--- Allocate an array of CVariable objects for each node in the  mesh ---*/
	node = new CVariable*[nPoint];

	/*--- Define auxiliary vectors to store residual-related quantities ---*/
	Residual      = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS  = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max  = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max  = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
	Residual_i    = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
	Residual_j    = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
	Res_Conv      = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
	Res_Visc      = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
	Res_Sour      = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;

	/*--- Define some auxiliary vectors related to the solution ---*/
	Solution   = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
	Solution_i = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
	Solution_j = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;

	/*--- Define some auxiliary vectors related to the geometry ---*/
	Vector   = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
	Vector_i = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
	Vector_j = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

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
		if (rank == MASTER_NODE)
      cout << "Initialize jacobian structure (TNE2 Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);

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

	/*--- Allocate force & coefficient arrays on boundaries ---*/
	CPressure = new su2double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CPressure[iMarker] = new su2double [geometry->nVertex[iMarker]];

	/*--- Heat tranfer in all the markers ---*/
	HeatFlux = new su2double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		HeatFlux[iMarker] = new su2double [geometry->nVertex[iMarker]];

	/*--- Skin friction in all the markers ---*/
	CSkinFriction = new su2double* [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		CSkinFriction[iMarker] = new su2double [geometry->nVertex[iMarker]];

  /*--- Y+ of all members ---*/
  YPlus = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    YPlus[iMarker] = new su2double[geometry->nVertex[iMarker]];

	/*--- Non dimensional coefficients ---*/
	ForceInviscid  = new su2double[3];
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

	/*--- Initialize total coefficients ---*/
	Total_CD    = 0.0;  Total_CL    = 0.0;  Total_CSF = 0.0;
  Total_CFx   = 0.0;  Total_CFy   = 0.0;  Total_CFz = 0.0;
	Total_CMx   = 0.0;  Total_CMy   = 0.0;  Total_CMz = 0.0;
	Total_CEff  = 0.0;
  Total_Heat     = 0.0;
  Total_MaxHeat  = 0.0;

	ForceViscous  = new su2double[3];
	MomentViscous = new su2double[3];
	CD_Visc       = new su2double[nMarker];
	CL_Visc       = new su2double[nMarker];
	CMx_Visc      = new su2double[nMarker];
	CMy_Visc      = new su2double[nMarker];
	CMz_Visc      = new su2double[nMarker];
	CEff_Visc     = new su2double[nMarker];
	CFx_Visc      = new su2double[nMarker];
	CFy_Visc      = new su2double[nMarker];
	CFz_Visc      = new su2double[nMarker];
  Heat_Visc        = new su2double[nMarker];
  MaxHeatFlux_Visc = new su2double[nMarker];

	/*--- Read farfield conditions from config ---*/
	Pressure_Inf       = config->GetPressure_FreeStream();
  Temperature_Inf    = config->GetTemperature_FreeStream();
  Temperature_ve_Inf = config->GetTemperature_ve_FreeStream();
  MassFrac_Inf       = config->GetMassFrac_FreeStream();
  Mach_Inf           = config->GetMach();

  // Note: May need to investigate these more carefully...
	Viscosity_Inf = config->GetViscosity_FreeStreamND();
	Mach_Inf      = config->GetMach();
	Prandtl_Lam   = config->GetPrandtl_Lam();
	Prandtl_Turb  = config->GetPrandtl_Turb();

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
                                   Temperature_ve_Inf, nDim, nVar,
                                   nPrimVar, nPrimVarGrad, config);
  check_infty = node_infty->SetPrimVar_Compressible(config);

  Velocity_Inf = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_Inf[iDim] = node_infty->GetVelocity(iDim);

	/*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
	if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {

		/*--- Initialize using freestream values ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint] = new CTNE2NSVariable(Pressure_Inf, MassFrac_Inf,
                                         Mvec_Inf, Temperature_Inf,
                                         Temperature_ve_Inf, nDim, nVar,
                                         nPrimVar, nPrimVarGrad, config);
	} else {

		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();
		restart_file.open(filename.data(), ios::in);

		/*--- In case there is no file ---*/
		if (restart_file.fail()) {
			cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
			exit(1);
		}

		/*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
		long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];

		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
			Global2Local[iPoint] = -1;

		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;

		/*--- The first line is the header ---*/
		getline (restart_file, text_line);

		while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);

			/*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
			iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
        point_line >> index;
        for (iDim = 0; iDim < nDim; iDim++)
          point_line >> dull_val;
        for (iVar = 0; iVar < nVar; iVar++) {
          point_line >> Solution[iVar];
        }

        /*--- Call the CVariable constructor with the solution ---*/
				node[iPoint_Local] = new CTNE2NSVariable(Solution, nDim, nVar,
                                                 nPrimVar,nPrimVarGrad, config);
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CTNE2NSVariable(Pressure_Inf, MassFrac_Inf,
                                         Mvec_Inf, Temperature_Inf,
                                         Temperature_ve_Inf, nDim, nVar,
                                         nPrimVar, nPrimVarGrad, config);

			//node[iPoint] = new CTNE2NSVariable(Solution, nDim, nVar, nPrimVar,
        //                                 nPrimVarGrad, config);

		/*--- Close the restart file ---*/
		restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}

  /*--- Check that the initial solution is physical ---*/
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    check = node[iPoint]->SetPrimVar_Compressible(config);

  //    node[iPoint]->SetDensity();
  //    node[iPoint]->SetVelocity2();
  //    check_temp = node[iPoint]->SetTemperature(config);
  //    check_press = node[iPoint]->SetPressure(config);
  //
    if (check) {
      bool ionization;
      unsigned short iEl, nHeavy, nEl, *nElStates;
      su2double Ru, T, Tve, rhoCvtr, sqvel, rhoE, rhoEve, num, denom, conc;
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
      Ru      = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(kmol*K)]
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
      for (iDim = 0; iDim < nDim; iDim++)
        sqvel += Mvec_Inf[iDim]*soundspeed * Mvec_Inf[iDim]*soundspeed;

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
        Solution[iSpecies]     = rho*MassFrac_Inf[iSpecies];
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        Solution[nSpecies+iDim]     = rho*Mvec_Inf[iDim]*soundspeed;
      }
      Solution[nSpecies+nDim]       = rhoE;
      Solution[nSpecies+nDim+1]     = rhoEve;

      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);

      counter_local++;
    }
  }

  #ifndef NO_MPI
  #ifdef WINDOWS
  MPI_Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  #else
  MPI::COMM_WORLD.Reduce(&counter_local, &counter_global, 1, MPI::UNSIGNED_LONG, MPI::SUM, MASTER_NODE);
  #endif
  #else

  counter_global = counter_local;

  #endif


  if ((rank == MASTER_NODE) && (counter_global != 0))
    cout << "Warning. The original solution contains "<< counter_global
         << " points that are not physical." << endl;

	/*--- Define solver parameters needed for execution of destructor ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    least_squares = true;
	else
    least_squares = false;

	/*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);

  /*--- Deallocate arrays ---*/
  delete [] Mvec_Inf;

}

CTNE2NSSolver::~CTNE2NSSolver(void) {
	unsigned short iMarker;

	if (CD_Visc != NULL) delete [] CD_Visc;
	if (CL_Visc != NULL) delete [] CL_Visc;
	if (CMx_Visc != NULL) delete [] CMx_Visc;
	if (CMy_Visc != NULL) delete [] CMy_Visc;
	if (CMz_Visc != NULL) delete [] CMz_Visc;
	if (CFx_Visc != NULL) delete [] CFx_Visc;
	if (CFy_Visc != NULL) delete [] CFy_Visc;
	if (CFz_Visc != NULL) delete [] CFz_Visc;
	if (CEff_Visc != NULL) delete [] CEff_Visc;
  if (Heat_Visc != NULL) delete [] Heat_Visc;
  if (MaxHeatFlux_Visc != NULL) delete [] MaxHeatFlux_Visc;
	if (ForceViscous != NULL) delete [] ForceViscous;
	if (MomentViscous != NULL) delete [] MomentViscous;

  if (YPlus != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] YPlus[iMarker];
    }
    delete [] YPlus;
  }

	if (CSkinFriction != NULL) {
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
			delete CSkinFriction[iMarker];
		}
		delete [] CSkinFriction;
	}

}

void CTNE2NSSolver::Preprocessing(CGeometry *geometry, CSolver **solution_container, CConfig *config,
                                  unsigned short iMesh, unsigned short iRKStep,
                                  unsigned short RunTime_EqSystem, bool Output) {
	unsigned long iPoint, ErrorCounter = 0;
  bool nonPhys;
	unsigned long ExtIter = config->GetExtIter();
  bool adjoint      = config->GetContinuous_Adjoint();
	bool implicit     = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool muscl        = (config->GetMUSCL_Flow());
  bool limiter      = ((config->GetKind_SlopeLimit_TNE2() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()));

  bool center       = ((config->GetKind_ConvNumScheme_TNE2() == SPACE_CENTERED) ||
                       (adjoint && config->GetKind_ConvNumScheme_AdjTNE2() == SPACE_CENTERED));
  int rank;
  su2double *errU, *errV;

  errU = new su2double[nVar];
  errV = new su2double[nPrimVar];

  #ifdef NO_MPI
  rank = MASTER_NODE;
  #else
  #ifdef WINDOWS
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  #else
  rank = MPI::COMM_WORLD.Get_rank();
  #endif
  #endif

	for (iPoint = 0; iPoint < nPoint; iPoint ++) {

		/*--- Set the primitive variables incompressible (dens, vx, vy, vz, beta)
     and compressible (temp, vx, vy, vz, press, dens, enthal, sos)---*/
		nonPhys = node[iPoint]->SetPrimVar_Compressible(config);
    if (nonPhys) {
      ErrorCounter++;
    }

		/*--- Initialize the convective, source and viscous residual vector ---*/
		if (!Output) LinSysRes.SetBlock_Zero(iPoint);

	}

  Set_MPI_Primitive(geometry, config);

	/*--- Compute gradient of the primitive variables ---*/
  switch (config->GetKind_Gradient_Method()) {
    case GREEN_GAUSS:
      SetPrimVar_Gradient_GG(geometry, config);
      SetSolution_Gradient_GG(geometry, config);
      break;
    case WEIGHTED_LEAST_SQUARES:
      SetPrimVar_Gradient_LS(geometry, config);
      SetSolution_Gradient_LS(geometry, config);
      break;
  }

  Set_MPI_Solution_Gradient(geometry, config);
  Set_MPI_Primitive_Gradient(geometry, config);


  if ((muscl) && (iMesh == MESH_0) && limiter) {
    SetSolution_Limiter(geometry, config);
  }

  Set_MPI_Solution_Limiter(geometry, config);

  /*--- Artificial dissipation ---*/
  if (center)
    SetMax_Eigenvalue(geometry, config);

	/*--- Initialize the jacobian matrices ---*/
	if (implicit) Jacobian.SetValZero();

  /*--- Error message ---*/
  #ifndef NO_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
  #ifdef WINDOWS
  MPI_Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  #else
  MPI::COMM_WORLD.Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
  #endif
  #endif
  if ((ErrorCounter != 0) && (rank == MASTER_NODE))
    cout <<"The solution contains "<< ErrorCounter << " non-physical points." << endl;

  delete [] errU;
  delete [] errV;
}

void CTNE2NSSolver::SetTime_Step(CGeometry *geometry,
                                 CSolver **solution_container,
                                 CConfig *config,
                                 unsigned short iMesh,
                                 unsigned long Iteration) {

	unsigned short iDim, iMarker, iSpecies;
  unsigned short VEL_INDEX, RHO_INDEX, RHOS_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;
  unsigned long iEdge, iVertex, iPoint, jPoint;
	double *Normal, Area, Vol;
  su2double Mean_SoundSpeed, Mean_ProjVel;
  su2double Lambda, Local_Delta_Time, Local_Delta_Time_Visc, Global_Delta_Time;
  su2double Mean_LaminarVisc, Mean_ThermalCond, Mean_ThermalCond_ve, Mean_Density, Mean_Tve;
  su2double cp, cv, cvve, Ru, *xi, *Ms, Na, Mmix, Rmix;
  su2double Lambda_1, Lambda_2, K_v, Global_Delta_UnstTimeND;
  su2double *V_i, *V_j, *X;
  su2double UnitNormal[3];
  su2double tmp;

	bool implicit = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
	bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialize parameters ---*/
  Global_Delta_Time = 1E6;
  Min_Delta_Time    = 1.E6;
  Max_Delta_Time    = 0.0;
  K_v    = 0.5;
  iPoint = 0;
  jPoint = 0;
  Ru = UNIVERSAL_GAS_CONSTANT;
  Na = AVOGAD_CONSTANT;

  A_INDEX = node_infty->GetAIndex();
  VEL_INDEX  = node_infty->GetVelIndex();
  RHO_INDEX  = node_infty->GetRhoIndex();
  RHOS_INDEX = node_infty->GetRhosIndex();
  RHOCVTR_INDEX = node_infty->GetRhoCvtrIndex();
  RHOCVVE_INDEX = node_infty->GetRhoCvveIndex();

  X = new su2double[nSpecies];

  /*--- Get from config ---*/
  xi = config->GetRotationModes();
  Ms = config->GetMolar_Mass();

	/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		node[iPoint]->SetMax_Lambda_Inv(0.0);
		node[iPoint]->SetMax_Lambda_Visc(0.0);
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
    V_i = node[iPoint]->GetPrimVar();
    V_j = node[jPoint]->GetPrimVar();

    /*--- Calculate the required mean values ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_ProjVel      = 0.5*( V_i[VEL_INDEX+iDim]
                               +V_j[VEL_INDEX+iDim] )*UnitNormal[iDim];
    Mean_SoundSpeed     = 0.5*(V_i[A_INDEX]   + V_j[A_INDEX]);
    Mean_Density        = 0.5*(V_i[RHO_INDEX] + V_j[RHO_INDEX]);
    Mean_ThermalCond    = 0.5*(node[iPoint]->GetThermalConductivity() +
                               node[jPoint]->GetThermalConductivity()  );
    Mean_ThermalCond_ve = 0.5*(node[iPoint]->GetThermalConductivity_ve() +
                               node[jPoint]->GetThermalConductivity_ve()  );

    /*--- Calculate the maximum spectral radius from convection ---*/
    Lambda = (fabs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
		if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);
		if (geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Inv(Lambda);

		/*--- Calculate mean viscous quantities ---*/
    Mean_LaminarVisc    = 0.5*(node[iPoint]->GetLaminarViscosity() +
                               node[jPoint]->GetLaminarViscosity()  );
    Mean_ThermalCond    = 0.5*(node[iPoint]->GetThermalConductivity() +
                               node[jPoint]->GetThermalConductivity()  );
    Mean_ThermalCond_ve = 0.5*(node[iPoint]->GetThermalConductivity_ve() +
                               node[jPoint]->GetThermalConductivity_ve()  );
    Mean_Density        = 0.5*(node[iPoint]->GetDensity() +
                               node[jPoint]->GetDensity()  );
    cv = 0.5*(node[iPoint]->GetRhoCv_tr() + node[iPoint]->GetRhoCv_ve() +
              node[jPoint]->GetRhoCv_tr() + node[jPoint]->GetRhoCv_ve()  )/ Mean_Density;

    /*--- Determine the viscous spectral radius and apply it to the control volume ---*/
		Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
		Lambda_2 = (Mean_ThermalCond+Mean_ThermalCond_ve)/cv;
		Lambda   = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
		if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Visc(Lambda);
		if (geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Visc(Lambda);
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
      V_i = node[iPoint]->GetPrimVar();

      /*--- Calculate the required mean values ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Mean_ProjVel      = V_i[VEL_INDEX+iDim]*UnitNormal[iDim];
      Mean_SoundSpeed     = V_i[A_INDEX];
      Mean_Density        = V_i[RHO_INDEX];
      Mean_ThermalCond    = node[iPoint]->GetThermalConductivity();
      Mean_ThermalCond_ve = node[iPoint]->GetThermalConductivity_ve();

      /*--- Calculate the maximum spectral radius from convection ---*/
      Lambda = (fabs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
      if (geometry->node[iPoint]->GetDomain())
        node[iPoint]->AddMax_Lambda_Inv(Lambda);

			/*--- Calculate viscous mean quantities ---*/
      Mean_LaminarVisc    = node[iPoint]->GetLaminarViscosity();
      Mean_ThermalCond    = node[iPoint]->GetThermalConductivity();
      Mean_ThermalCond_ve = node[iPoint]->GetThermalConductivity_ve();
      Mean_Density        = node[iPoint]->GetDensity();
      cv = (node[iPoint]->GetRhoCv_tr() +
            node[iPoint]->GetRhoCv_ve()  ) / Mean_Density;

			Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
			Lambda_2 = (Mean_ThermalCond+Mean_ThermalCond_ve)/cv;
			Lambda   = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;

			if (geometry->node[iPoint]->GetDomain())
        node[iPoint]->AddMax_Lambda_Visc(Lambda);
		}
	}

	/*--- Each element uses their own speed ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		Vol = geometry->node[iPoint]->GetVolume();

    /*--- Calculate local inv. and visc. dTs, take the minimum of the two ---*/
		Local_Delta_Time      = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
		Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
		Local_Delta_Time      = min(Local_Delta_Time, Local_Delta_Time_Visc);
		Global_Delta_Time     = min(Global_Delta_Time, Local_Delta_Time);

    /*--- Store minimum and maximum dt's within the grid for printing ---*/
		Min_Delta_Time        = min(Min_Delta_Time, Local_Delta_Time);
		Max_Delta_Time        = max(Max_Delta_Time, Local_Delta_Time);

    /*--- Set the time step ---*/
		node[iPoint]->SetDelta_Time(Local_Delta_Time);
	}

  /*--- Communicate minimum and maximum time steps ---*/
  #ifndef NO_MPI
  su2double rbuf_time_min, rbuf_time_max, sbuf_time_min, sbuf_time_max;
  sbuf_time_min = Min_Delta_Time;
  sbuf_time_max = Max_Delta_Time;
  #ifdef WINDOWS
  MPI_Reduce(&sbuf_time_min, &rbuf_time_min, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Bcast(&rbuf_time_min, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&sbuf_time_max, &rbuf_time_max, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Bcast(&rbuf_time_max, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  #else
  MPI::COMM_WORLD.Reduce(&sbuf_time_min, &rbuf_time_min, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
  MPI::COMM_WORLD.Bcast(&rbuf_time_min, 1, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Reduce(&sbuf_time_max, &rbuf_time_max, 1, MPI::DOUBLE, MPI::MAX, MASTER_NODE);
  MPI::COMM_WORLD.Bcast(&rbuf_time_max, 1, MPI::DOUBLE, MASTER_NODE);
  MPI::COMM_WORLD.Barrier();
  #endif
  Min_Delta_Time = rbuf_time_min;
  Max_Delta_Time = rbuf_time_max;
  #endif



	/*--- Check if there is any element with only one neighbor...
   a CV that is inside another CV ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		if (geometry->node[iPoint]->GetnPoint() == 1)
			node[iPoint]->SetDelta_Time(Min_Delta_Time);
	}

	/*--- For exact time solution use the minimum delta time of the whole mesh ---*/
	if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
  #ifndef NO_MPI
		su2double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_Time;
  #ifdef WINDOWS
		MPI_Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
		MPI_Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
  #else
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
  #endif
		Global_Delta_Time = rbuf_time;
  #endif
		for(iPoint = 0; iPoint < nPointDomain; iPoint++)
			node[iPoint]->SetDelta_Time(Global_Delta_Time);
	}

	/*--- Recompute the unsteady time step for the dual time stratey
	 if the unsteady CFL is diferent from 0 ---*/
	if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
		Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
  #ifndef NO_MPI
		su2double rbuf_time, sbuf_time;
		sbuf_time = Global_Delta_UnstTimeND;
  #ifdef WINDOWS
		MPI_Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
		MPI_Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
  #else
		MPI::COMM_WORLD.Reduce(&sbuf_time, &rbuf_time, 1, MPI::DOUBLE, MPI::MIN, MASTER_NODE);
		MPI::COMM_WORLD.Bcast(&rbuf_time, 1, MPI::DOUBLE, MASTER_NODE);
		MPI::COMM_WORLD.Barrier();
  #endif
		Global_Delta_UnstTimeND = rbuf_time;
  #endif
		config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
	}

	/*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
	if (dual_time)
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			if (!implicit) {
				Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
				/*--- Check if there is any element with only one neighbor...
				 a CV that is inside another CV ---*/
				if (geometry->node[iPoint]->GetnPoint() == 1) Local_Delta_Time = 0.0;
				node[iPoint]->SetDelta_Time(Local_Delta_Time);
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
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );


  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord() );
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Primitive variables, and gradient ---*/
    numerics->SetConservative   (node[iPoint]->GetSolution(),
                                 node[jPoint]->GetSolution() );
    numerics->SetConsVarGradient(node[iPoint]->GetGradient(),
                                 node[jPoint]->GetGradient() );
    numerics->SetPrimitive      (node[iPoint]->GetPrimVar(),
                                 node[jPoint]->GetPrimVar() );
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                 node[jPoint]->GetGradient_Primitive() );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU  (node[iPoint]->GetdPdU(),   node[jPoint]->GetdPdU());
    numerics->SetdTdU  (node[iPoint]->GetdTdU(),   node[jPoint]->GetdTdU());
    numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node[jPoint]->GetdTvedU());
    numerics->SetEve   (node[iPoint]->GetEve(),    node[jPoint]->GetEve());
    numerics->SetCvve  (node[iPoint]->GetCvve(),   node[jPoint]->GetCvve());

    /*--- Species diffusion coefficients ---*/
    numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(),
                                node[jPoint]->GetDiffusionCoeff() );

    /*--- Laminar viscosity ---*/
    numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(),
                                  node[jPoint]->GetLaminarViscosity() );

    /*--- Thermal conductivity ---*/
    numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(),
                                     node[jPoint]->GetThermalConductivity());

    /*--- Vib-el. thermal conductivity ---*/
    numerics->SetThermalConductivity_ve(node[iPoint]->GetThermalConductivity_ve(),
                                        node[jPoint]->GetThermalConductivity_ve() );

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
  numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  numerics->SetPIndex      ( node[0]->GetPIndex()       );
  numerics->SetTIndex      ( node[0]->GetTIndex()       );
  numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  numerics->SetHIndex      ( node[0]->GetHIndex()       );
  numerics->SetAIndex      ( node[0]->GetAIndex()       );
  numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );

  /*--- loop over interior points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
  //    if (geometry->node[iPoint]->GetSolidBoundary()) {

    /*--- Set conserved & primitive variables  ---*/
    numerics->SetConservative(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
    numerics->SetPrimitive   (node[iPoint]->GetPrimVar(),  node[iPoint]->GetPrimVar() );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU  (node[iPoint]->GetdPdU(),   node[iPoint]->GetdPdU()  );
    numerics->SetdTdU  (node[iPoint]->GetdTdU(),   node[iPoint]->GetdTdU()  );
    numerics->SetdTvedU(node[iPoint]->GetdTvedU(), node[iPoint]->GetdTvedU());
    numerics->SetEve   (node[iPoint]->GetEve(),    node[iPoint]->GetEve()   );
    numerics->SetCvve  (node[iPoint]->GetCvve(),   node[iPoint]->GetCvve()  );

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

    /*--- Compute axisymmetric source terms (if needed) ---*/
    if (config->GetAxisymmetric()) {
      numerics->ComputeAxisymmetric(Residual, Jacobian_i, config);

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
    numerics->ComputeChemistry(Residual, Jacobian_i, config);

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
    numerics->ComputeVibRelaxation(Residual, Jacobian_i, config);


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

  #ifndef NO_MPI
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

void CTNE2NSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {

  unsigned short Boundary, Monitoring, iMarker, iDim, jDim;
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
  su2double dTn, dTven, pnorm, HeatLoad;
  su2double Alpha, Beta, RefLength, RefArea, *Origin;
  su2double factor;

  su2double WallShearStress, WallDistMod, WallDist[3];

  su2double Vel[3], Velocity_Inf[3];

  su2double UnitNormal[3];
  su2double TauElem[3];
  su2double TauTangent[3];
  su2double Tau[3][3];

  /*--- Retrieve index information from CVariable ---*/
  VEL_INDEX = node[0]->GetVelIndex();
  T_INDEX   = node[0]->GetTIndex();
  TVE_INDEX = node[0]->GetTveIndex();

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
    Velocity_Inf[iDim] = node_infty->GetVelocity(iDim);
    RefVel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
	RefDensity  = node_infty->GetDensity();

	factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

	/*-- Initialization --*/
  AllBound_CMx_Visc   = 0.0; AllBound_CMy_Visc   = 0.0; AllBound_CMz_Visc = 0.0;
	AllBound_CFx_Visc   = 0.0; AllBound_CFy_Visc   = 0.0; AllBound_CFz_Visc = 0.0;
	AllBound_CD_Visc = 0.0; AllBound_CL_Visc = 0.0;
	AllBound_HeatFlux_Visc     = 0.0; AllBound_MaxHeatFlux_Visc  = 0.0;
	AllBound_CEff_Visc  = 0.0;

	/*--- Loop over the Navier-Stokes markers ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- Identify boundary information ---*/
    Boundary   = config->GetMarker_All_KindBC(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);

    /*--- Forces initialization at each Marker ---*/
    CD_Visc[iMarker] = 0.0; CL_Visc[iMarker] = 0.0; CEff_Visc[iMarker] = 0.0;
    CMx_Visc[iMarker]   = 0.0; CMy_Visc[iMarker]   = 0.0; CMz_Visc[iMarker]  = 0.0;
    CFx_Visc[iMarker]   = 0.0; CFy_Visc[iMarker]   = 0.0; CFz_Visc[iMarker]  = 0.0;
    Heat_Visc[iMarker]  = 0.0; MaxHeatFlux_Visc[iMarker] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ForceViscous[iDim]  = 0.0;
      MomentViscous[iDim] = 0.0;
    }
    HeatLoad = 0.0;

		if ((Boundary == HEAT_FLUX              ) ||
        (Boundary == HEAT_FLUX_CATALYTIC    ) ||
        (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == ISOTHERMAL_CATALYTIC   ) ||
        (Boundary == ISOTHERMAL_NONCATALYTIC)) {

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
				Grad_PrimVar   = node[iPoint]->GetGradient_Primitive();
        Viscosity      = node[iPoint]->GetLaminarViscosity();
        ThermalCond    = node[iPoint]->GetThermalConductivity();
        ThermalCond_ve = node[iPoint]->GetThermalConductivity_ve();
        Density        = node[iPoint]->GetDensity();

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
          Vel[iDim] = node[iPointNormal]->GetVelocity(iDim);

        for (iDim = 0; iDim < nDim; iDim++)
          WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          WallDistMod += WallDist[iDim]*WallDist[iDim];
        WallDistMod = sqrt(WallDistMod);

				/*--- Compute wall skin friction coefficient, and heat flux on the wall ---*/
				CSkinFriction[iMarker][iVertex] = WallShearStress / (0.5*RefDensity*RefVel2);

				/*--- Compute y+ and non-dimensional velocity ---*/
				FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);

				/*--- Compute heat flux on the wall ---*/
				dTn = 0.0; dTven = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          dTn   += Grad_PrimVar[T_INDEX][iDim]*UnitNormal[iDim];
          dTven += Grad_PrimVar[TVE_INDEX][iDim]*UnitNormal[iDim];
        }

  //        cout << "S: " << node[iPoint]->GetGradient()[0][0] << endl;
  //        cout << "U: " << node[iPoint]->GetSolution(0) << endl;
  //        unsigned short iVar;
  //        for (iVar = 0; iVar < nPrimVar; iVar++) {
  //          for (iDim = 0; iDim < nDim; iDim++)
  //            cout << Grad_PrimVar[iVar][iDim] << "\t";
  //          cout << endl;
  //        }
  //        cin.get();
  //        if (Grad_PrimVar[T_INDEX][0] != 0) {
  //          cout << Grad_PrimVar[T_INDEX][0] << "\t" << Grad_PrimVar[T_INDEX][1] << endl;
  //          cin.get();
  //        }

        HeatFlux[iMarker][iVertex] = ThermalCond*dTn + ThermalCond_ve*dTven;
        Heat_Visc[iMarker] += HeatFlux[iMarker][iVertex]*Area;
        MaxHeatFlux_Visc[iMarker] += pow(HeatFlux[iMarker][iVertex], pnorm)*Area;

				/*--- Compute viscous forces, and moment using the stress tensor ---*/
				if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {

					for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim]*Area*factor;
						ForceViscous[iDim] += Force[iDim];
					}

					if (iDim == 3) {
            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLength;
            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLength;
          }
          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLength;
				}
			}

			/*--- Transform ForceInviscid into CLift and CDrag ---*/
			if  (Monitoring == YES) {

				if (nDim == 2) {
					CD_Visc[iMarker] =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
					CL_Visc[iMarker] = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
					CMz_Visc[iMarker]   = MomentViscous[2];
					CEff_Visc[iMarker]  = CL_Visc[iMarker]/(CD_Visc[iMarker]+EPS);
					CFx_Visc[iMarker]   = ForceViscous[0];
					CFy_Visc[iMarker]   = ForceViscous[1];
					CFz_Visc[iMarker]   = 0.0;
				}

				if (nDim == 3) {
					CD_Visc[iMarker] = ForceViscous[0]*cos(Alpha)*cos(Beta) +
                                ForceViscous[1]*sin(Beta) +
                                ForceViscous[2]*sin(Alpha)*cos(Beta);
					CL_Visc[iMarker] = -ForceViscous[0]*sin(Alpha) +
                                 ForceViscous[2]*cos(Alpha);
          CEff_Visc[iMarker]  = CL_Visc[iMarker]/(CD_Visc[iMarker]+EPS);
					CMx_Visc[iMarker]   = MomentViscous[0];
					CMy_Visc[iMarker]   = MomentViscous[1];
					CMz_Visc[iMarker]   = MomentViscous[2];
					CFx_Visc[iMarker]   = ForceViscous[0];
					CFy_Visc[iMarker]   = ForceViscous[1];
					CFz_Visc[iMarker]   = ForceViscous[2];
				}

				AllBound_CD_Visc       += CD_Visc[iMarker];
				AllBound_CL_Visc       += CL_Visc[iMarker];
        AllBound_CEff_Visc        += CEff_Visc[iMarker];
				AllBound_CMx_Visc         += CMx_Visc[iMarker];
				AllBound_CMy_Visc         += CMy_Visc[iMarker];
				AllBound_CMz_Visc         += CMz_Visc[iMarker];
				AllBound_CFx_Visc         += CFx_Visc[iMarker];
				AllBound_CFy_Visc         += CFy_Visc[iMarker];
				AllBound_CFz_Visc         += CFz_Visc[iMarker];
        AllBound_MaxHeatFlux_Visc += MaxHeatFlux_Visc[iMarker];
        AllBound_HeatFlux_Visc    += Heat_Visc[iMarker];

			}
		}
	}


  #ifndef NO_MPI

  /*--- Add AllBound information using all the nodes ---*/

  su2double MyAllBound_CD_Visc       = AllBound_CD_Visc;
  su2double MyAllBound_CL_Visc       = AllBound_CL_Visc;
  su2double MyAllBound_CEff_Visc        = AllBound_CEff_Visc;
  su2double MyAllBound_CMx_Visc         = AllBound_CMx_Visc;
  su2double MyAllBound_CMy_Visc         = AllBound_CMy_Visc;
  su2double MyAllBound_CMz_Visc         = AllBound_CMz_Visc;
  su2double MyAllBound_CFx_Visc         = AllBound_CFx_Visc;
  su2double MyAllBound_CFy_Visc         = AllBound_CFy_Visc;
  su2double MyAllBound_CFz_Visc         = AllBound_CFz_Visc;
  su2double MyAllBound_HeatFlux_Visc    = AllBound_HeatFlux_Visc;
  su2double MyAllBound_MaxHeatFlux_Visc = AllBound_MaxHeatFlux_Visc;

  AllBound_CD_Visc         = 0.0;
  AllBound_CL_Visc         = 0.0;
  AllBound_CEff_Visc          = 0.0;
  AllBound_CMx_Visc           = 0.0;
  AllBound_CMy_Visc           = 0.0;
  AllBound_CMz_Visc           = 0.0;
  AllBound_CFx_Visc           = 0.0;
  AllBound_CFy_Visc           = 0.0;
  AllBound_CFz_Visc           = 0.0;
  AllBound_HeatFlux_Visc      = 0.0;
  AllBound_MaxHeatFlux_Visc   = 0.0;

  #ifdef WINDOWS
  MPI_Allreduce(&MyAllBound_CD_Visc, &AllBound_CD_Visc, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CL_Visc, &AllBound_CL_Visc, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  MPI_Allreduce(&MyAllBound_CMx_Visc,      &AllBound_CMx_Visc,      1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CMy_Visc,      &AllBound_CMy_Visc,      1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CMz_Visc,      &AllBound_CMz_Visc,      1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFx_Visc,      &AllBound_CFx_Visc,      1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFy_Visc,      &AllBound_CFy_Visc,      1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFz_Visc,      &AllBound_CFz_Visc,      1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_HeatFlux_Visc,     &AllBound_HeatFlux_Visc,     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_MaxHeatFlux_Visc, &AllBound_MaxHeatFlux_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #else
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CD_Visc, &AllBound_CD_Visc, 1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CL_Visc, &AllBound_CL_Visc, 1, MPI::DOUBLE, MPI::SUM);
  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CMx_Visc,         &AllBound_CMx_Visc,         1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CMy_Visc,         &AllBound_CMy_Visc,         1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CMz_Visc,         &AllBound_CMz_Visc,         1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CFx_Visc,         &AllBound_CFx_Visc,         1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CFy_Visc,         &AllBound_CFy_Visc,         1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_CFz_Visc,         &AllBound_CFz_Visc,         1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_HeatFlux_Visc,    &AllBound_HeatFlux_Visc,    1, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MyAllBound_MaxHeatFlux_Visc, &AllBound_MaxHeatFlux_Visc, 1, MPI::DOUBLE, MPI::SUM);
  #endif
  #endif

	Total_CD   += AllBound_CD_Visc;
	Total_CL   += AllBound_CL_Visc;
	Total_CMx     += AllBound_CMx_Visc;
	Total_CMy     += AllBound_CMy_Visc;
	Total_CMz     += AllBound_CMz_Visc;
	Total_CEff     = Total_CL/(Total_CD+EPS);
	Total_CFx     += AllBound_CFx_Visc;
	Total_CFy     += AllBound_CFy_Visc;
	Total_CFz     += AllBound_CFz_Visc;
  Total_Heat     = AllBound_HeatFlux_Visc;
  Total_MaxHeat  = pow(AllBound_MaxHeatFlux_Visc, 1.0/pnorm);
}

void CTNE2NSSolver::BC_Sym_Plane(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker) {

  /*--- Call the Euler wall routine ---*/
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config,
                val_marker);

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
	double Wall_HeatFlux, dTdn, dTvedn, ktr, kve, pcontrol;
	double *Normal, Area;
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
  T_INDEX    = node[0]->GetTIndex();
  TVE_INDEX  = node[0]->GetTveIndex();

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
      GradV  = node[iPoint]->GetGradient_Primitive();
      dTdn   = 0.0;
      dTvedn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
        dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
      }
      ktr = node[iPoint]->GetThermalConductivity();
  //      kve = node[iPoint]->GetThermalConductivity_ve();
  //			Res_Visc[nSpecies+nDim]   += pcontrol*(ktr*dTdn+kve*dTvedn) +
  //                                   Wall_HeatFlux*Area;
  //      Res_Visc[nSpecies+nDim+1] += pcontrol*(kve*dTvedn) +
  //                                   Wall_HeatFlux*Area;
  //
  //			/*--- Apply viscous residual to the linear system ---*/
  //      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Apply the no-slip condition in a strong way ---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
			node[iPoint]->SetVelocity_Old(Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        node[iPoint]->SetVal_ResTruncError_Zero(nSpecies+iDim);
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
  //  RHOS_INDEX = node[0]->GetRhosIndex();
  //  RHO_INDEX  = node[0]->GetRhoIndex();
  //  T_INDEX    = node[0]->GetTIndex();
  //  TVE_INDEX  = node[0]->GetTveIndex();
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
  //		if (geometry->node[iPoint]->GetDomain()) {
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
  //      V = node[iPoint]->GetPrimVar();
  //      GradV  = node[iPoint]->GetGradient_Primitive();
  //
  //      /*--- Rename for convenience ---*/
  //      rho = V[RHO_INDEX];
  //      Ds  = node[iPoint]->GetDiffusionCoeff();
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
  //        eves = node[iPoint]->CalcEve(config, V[TVE_INDEX], iSpecies);
  //        hs   = node[iPoint]->CalcHs(config, V[T_INDEX], eves, iSpecies);
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
	double Wall_HeatFlux, dTdn, dTvedn, ktr, kve, pcontrol;
  su2double rho, Ys, eves, hs;
	double *Normal, Area;
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
  T_INDEX    = node[0]->GetTIndex();
  TVE_INDEX  = node[0]->GetTveIndex();
  RHOS_INDEX = node[0]->GetRhosIndex();
  RHO_INDEX  = node[0]->GetRhoIndex();

  /*--- Allocate arrays ---*/
  dYdn = new su2double[nSpecies];
  GradY = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    GradY[iSpecies] = new su2double[nDim];

  //  /*--- Pass structure of the primitive variable vector to CNumerics ---*/
  //  sour_numerics->SetRhosIndex   ( node[0]->GetRhosIndex()    );
  //  sour_numerics->SetRhoIndex    ( node[0]->GetRhoIndex()     );
  //  sour_numerics->SetPIndex      ( node[0]->GetPIndex()       );
  //  sour_numerics->SetTIndex      ( node[0]->GetTIndex()       );
  //  sour_numerics->SetTveIndex    ( node[0]->GetTveIndex()     );
  //  sour_numerics->SetVelIndex    ( node[0]->GetVelIndex()     );
  //  sour_numerics->SetHIndex      ( node[0]->GetHIndex()       );
  //  sour_numerics->SetAIndex      ( node[0]->GetAIndex()       );
  //  sour_numerics->SetRhoCvtrIndex( node[0]->GetRhoCvtrIndex() );
  //  sour_numerics->SetRhoCvveIndex( node[0]->GetRhoCvveIndex() );

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
			node[iPoint]->SetVelocity_Old(Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        node[iPoint]->SetVal_ResTruncError_Zero(nSpecies+iDim);
      }

      /*--- Get temperature gradient information ---*/
      V = node[iPoint]->GetPrimVar();
      GradV  = node[iPoint]->GetGradient_Primitive();
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
        Ds  = node[iPoint]->GetDiffusionCoeff();

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
          eves = node[iPoint]->CalcEve(config, V[TVE_INDEX], iSpecies);
          hs   = node[iPoint]->CalcHs(config, V[T_INDEX], eves, iSpecies);
          //          Res_Visc[iSpecies] = -rho*Ds[iSpecies]*dYdn[iSpecies] + Ys*SdYdn;
          //          Res_Visc[nSpecies+nDim]   += Res_Visc[iSpecies]*hs;
          //          Res_Visc[nSpecies+nDim+1] += Res_Visc[iSpecies]*eves;
        }
      }

      /*--- Get node thermal conductivity ---*/
      ktr = node[iPoint]->GetThermalConductivity();
      kve = node[iPoint]->GetThermalConductivity_ve();

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
      //      sour_numerics->SetConservative(node[iPoint]->GetSolution(),
      //                                     node[iPoint]->GetSolution() );
      //      sour_numerics->SetPrimitive   (node[iPoint]->GetPrimVar() ,
      //                                     node[iPoint]->GetPrimVar()  );
      //      sour_numerics->SetdPdU        (node[iPoint]->GetdPdU()    ,
      //                                     node[iPoint]->GetdPdU()     );
      //      sour_numerics->SetdTdU        (node[iPoint]->GetdTdU()    ,
      //                                     node[iPoint]->GetdTdU()     );
      //      sour_numerics->SetdTvedU      (node[iPoint]->GetdTvedU()  ,
      //                                     node[iPoint]->GetdTvedU()   );
      //      sour_numerics->SetVolume      (geometry->node[iPoint]->GetVolume());
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

  bool ionization, implicit;
  unsigned short iDim, iVar, jVar;
  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr, kve;
  su2double Ti, Tvei, Tj, Tvej, *dTdU, *dTvedU;
  su2double Twall, dij, theta;
  su2double Area, *Normal, UnitNormal[3];
  su2double *Coord_i, *Coord_j;
  su2double C;

  implicit   = (config->GetKind_TimeIntScheme_TNE2() == EULER_IMPLICIT);
  ionization = config->GetIonization();

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

  RHOS_INDEX    = node[0]->GetRhosIndex();
  T_INDEX       = node[0]->GetTIndex();
  TVE_INDEX     = node[0]->GetTveIndex();
  RHOCVTR_INDEX = node[0]->GetRhoCvtrIndex();
  RHOCVVE_INDEX = node[0]->GetRhoCvveIndex();

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
			node[iPoint]->SetVelocity_Old(Vector);
			for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        node[iPoint]->SetVal_ResTruncError_Zero(nSpecies+iDim);
      }

      /*--- Calculate the gradient of temperature ---*/
      Ti   = node[iPoint]->GetTemperature();
      Tj   = node[jPoint]->GetTemperature();
      Tvei = node[iPoint]->GetTemperature_ve();
      Tvej = node[jPoint]->GetTemperature_ve();

      /*--- Rename variables for convenience ---*/
      ktr     = node[iPoint]->GetThermalConductivity();
      kve     = node[iPoint]->GetThermalConductivity_ve();

      /*--- Apply to the linear system ---*/
      Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                   (ktr*(Twall-Ti) + kve*(Twall-Tvei))*C)*Area/dij;
      Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Twall-Tvei) *C)*Area/dij;
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;

        dTdU   = node[iPoint]->GetdTdU();
        dTvedU = node[iPoint]->GetdTvedU();
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
	double pcontrol;
  su2double rho, *eves, *hs, Ru, *Ms, *xi;
  su2double *dTdU, *dTvedU, *Cvtr, *Cvve;
	double *Normal, Area, dij, UnitNormal[3];
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
  RHOS_INDEX    = node[0]->GetRhosIndex();
  RHO_INDEX     = node[0]->GetRhoIndex();
  T_INDEX       = node[0] ->GetTIndex();

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
  Ru = UNIVERSAL_GAS_CONSTANT;
  Ms = config->GetMolar_Mass();
  xi = config->GetRotationModes();

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
      Vi = node[iPoint]->GetPrimVar();
      Vj = node[jPoint]->GetPrimVar();
      Di = node[iPoint]->GetDiffusionCoeff();
      Dj = node[jPoint]->GetDiffusionCoeff();
      eves = node[iPoint]->GetEve();
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        hs[iSpecies] = node[iPoint]->CalcHs(config, Vi[T_INDEX],
                                            eves[iSpecies], iSpecies);
        Yj[iSpecies] = Vj[RHOS_INDEX+iSpecies]/Vj[RHO_INDEX];
      }
      rho    = Vi[RHO_INDEX];
      dTdU   = node[iPoint]->GetdTdU();
      dTvedU = node[iPoint]->GetdTvedU();


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
        Cvve = node[iPoint]->GetCvve();

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
