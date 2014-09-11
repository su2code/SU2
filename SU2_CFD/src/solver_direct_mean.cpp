/*!
 * \file solution_direct_mean.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.1 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;

  Surface_CLift = NULL; Surface_CDrag = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  /*--- Rotorcraft simulation array initialization ---*/

  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;

  /*--- Supersonic simulation array initialization ---*/

  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;

  /*--- Engine simulation array initialization ---*/

  FanFace_MassFlow = NULL;  FanFace_Pressure = NULL;
  FanFace_Mach = NULL;  FanFace_Area = NULL;
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
  double Density, Velocity2, Pressure, Temperature, dull_val;

  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  double Gas_Constant = config->GetGas_ConstantND();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool adjoint = config->GetAdjoint();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Array initialization ---*/

  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  Surface_CLift = NULL; Surface_CDrag = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;

  ForceInviscid = NULL; MomentInviscid = NULL;
  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;

  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;

  CEquivArea_Inv = NULL;  CNearFieldOF_Inv = NULL;

  FanFace_MassFlow = NULL;  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  FanFace_Pressure = NULL;  FanFace_Mach = NULL;  FanFace_Area = NULL;

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
   Compressible flow, primitive variables nDim+7, (T,vx,vy,vz,P,rho,h,c,lamMu,EddyMu).
   Incompressible flow, primitive variables nDim+5, (P,vx,vy,vz,rho,beta,lamMu,EddyMu).
   FreeSurface Incompressible flow, primitive variables nDim+7, (P,vx,vy,vz,rho,beta,lamMu,EddyMu,LevelSet,Dist).
   ---*/

  nDim = geometry->GetnDim();
  if (incompressible) { nVar = nDim+1; nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  if (freesurface)    { nVar = nDim+2; nPrimVar = nDim+7; nPrimVarGrad = nDim+6; }
  if (compressible)   { nVar = nDim+2;
    nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
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

  Residual      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Residual_i    = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the geometry ---*/

  Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

  /*--- Define some auxiliary vectors related to the primitive solution ---*/

  Primitive   = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the Secondary solution ---*/

  Secondary   = new double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;


  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new double [nPoint];
    jPoint_UndLapl = new double [nPoint];
  }

  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/

  if (roe_turkel) {
    LowMach_Precontioner = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new double[nVar];
  }

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

  } else {
    if (rank == MASTER_NODE) cout << "Explicit scheme. No Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
  }

  /*--- Define some auxiliary vectors for computing flow variable gradients by least squares ---*/

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/

    Smatrix = new double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new double [nDim];

    /*--- c vector := transpose(WA)*(Wb) ---*/

    cvector = new double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      cvector[iVar] = new double [nDim];

  }

  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/

  CharacPrimVar = new double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }

  /*--- Force definition and coefficient arrays for all of the markers ---*/

  CPressure = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Force definition and coefficient arrays for all of the markers ---*/

  CPressureTarget = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressureTarget[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Non-dimensional coefficients ---*/

  ForceInviscid     = new double[nDim];
  MomentInviscid    = new double[3];
  CDrag_Inv         = new double[nMarker];
  CLift_Inv         = new double[nMarker];
  CSideForce_Inv    = new double[nMarker];
  CMx_Inv           = new double[nMarker];
  CMy_Inv           = new double[nMarker];
  CMz_Inv           = new double[nMarker];
  CEff_Inv          = new double[nMarker];
  CFx_Inv           = new double[nMarker];
  CFy_Inv           = new double[nMarker];
  CFz_Inv           = new double[nMarker];
  Surface_CLift_Inv = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv = new double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv   = new double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv   = new double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv   = new double[config->GetnMarker_Monitoring()];
  Surface_CLift     = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag     = new double[config->GetnMarker_Monitoring()];
  Surface_CMx       = new double[config->GetnMarker_Monitoring()];
  Surface_CMy       = new double[config->GetnMarker_Monitoring()];
  Surface_CMz       = new double[config->GetnMarker_Monitoring()];

  /*--- Rotorcraft coefficients ---*/

  CT_Inv           = new double[nMarker];
  CQ_Inv           = new double[nMarker];
  CMerit_Inv       = new double[nMarker];

  /*--- Supersonic coefficients ---*/

  CEquivArea_Inv   = new double[nMarker];
  CNearFieldOF_Inv = new double[nMarker];

  /*--- Engine simulation ---*/

  FanFace_MassFlow  = new double[nMarker];
  Exhaust_MassFlow  = new double[nMarker];
  Exhaust_Area      = new double[nMarker];
  FanFace_Pressure  = new double[nMarker];
  FanFace_Mach      = new double[nMarker];
  FanFace_Area      = new double[nMarker];

  /*--- Init total coefficients ---*/

  Total_CDrag = 0.0;    Total_CLift = 0.0;      Total_CSideForce = 0.0;
  Total_CMx = 0.0;      Total_CMy = 0.0;        Total_CMz = 0.0;
  Total_CEff = 0.0;     Total_CEquivArea = 0.0; Total_CNearFieldOF = 0.0;
  Total_CFx = 0.0;      Total_CFy = 0.0;        Total_CFz = 0.0;
  Total_CT = 0.0;       Total_CQ = 0.0;         Total_CMerit = 0.0;
  Total_MaxHeat = 0.0;  Total_Heat = 0.0;
  Total_CpDiff = 0.0;   Total_HeatFluxDiff = 0.0;

  /*--- Read farfield conditions ---*/

  Density_Inf  = config->GetDensity_FreeStreamND();
  Pressure_Inf = config->GetPressure_FreeStreamND();
  Velocity_Inf = config->GetVelocity_FreeStreamND();
  Energy_Inf   = config->GetEnergy_FreeStreamND();
  Mach_Inf     = config->GetMach();

  /*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    FanFace_MassFlow[iMarker] = 0.0;
    Exhaust_MassFlow[iMarker] = 0.0;
    FanFace_Mach[iMarker] = Mach_Inf;
    FanFace_Pressure[iMarker] = Pressure_Inf;
    FanFace_Area[iMarker] = 0.0;
    Exhaust_Area[iMarker] = 0.0;
  }

  /*--- Initialize the cauchy critera array for fixed CL mode ---*/

  if (config->GetFixed_CL_Mode())
    Cauchy_Serie = new double [config->GetCauchy_Elems()+1];

  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {

    /*--- Restart the solution from the free-stream state ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
  }

  else {

    /*--- Initialize the solution from the restart file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();

    /*--- Modify file name for an unsteady restart ---*/
    if (dual_time) {
      int Unst_RestartIter;
      if (adjoint) {
        Unst_RestartIter = int(config->GetUnst_AdjointIter()) - 1;
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = int(config->GetUnst_RestartIter())-1;
      else
        Unst_RestartIter = int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }

    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
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
      }
      iPoint_Global++;
    }

    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
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
      Pressure    = Gamma_Minus_One*Density*(node[iPoint]->GetSolution(nDim+1)/Density-0.5*Velocity2);
      Temperature = Pressure / ( Gas_Constant * Density);
      if ((Pressure < 0.0) || (Temperature < 0.0)) {
        Solution[0] = Density_Inf;
        for (iDim = 0; iDim < nDim; iDim++)
          Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
        Solution[nDim+1] = Energy_Inf*Density_Inf;
        node[iPoint]->SetSolution(Solution);
        node[iPoint]->SetSolution_Old(Solution);
        counter_local++;
      }
    }
#ifdef HAVE_MPI
    MPI_Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    config->SetNonphysical_Points(counter_global);
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

  /*--- Perform the MPI communication of the solution ---*/
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
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;
  if (Surface_CLift != NULL)    delete [] Surface_CLift;
  if (Surface_CDrag != NULL)    delete [] Surface_CDrag;
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
  if (FanFace_MassFlow != NULL)  delete [] FanFace_MassFlow;
  if (Exhaust_MassFlow != NULL)  delete [] Exhaust_MassFlow;
  if (Exhaust_Area != NULL)      delete [] Exhaust_Area;
  if (FanFace_Pressure != NULL)  delete [] FanFace_Pressure;
  if (FanFace_Mach != NULL)      delete [] FanFace_Mach;
  if (FanFace_Area != NULL)      delete [] FanFace_Area;
  if (iPoint_UndLapl != NULL)       delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL)       delete [] jPoint_UndLapl;
  if (Primitive != NULL)        delete [] Primitive;
  if (Primitive_i != NULL)      delete [] Primitive_i;
  if (Primitive_j != NULL)      delete [] Primitive_j;

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
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  int send_to, receive_from;

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];

      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  int send_to, receive_from;

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
  int send_to, receive_from;

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  int send_to, receive_from;

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new double [nBufferR_Vector];
      Buffer_Send_Lambda = new double[nBufferS_Vector];
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
      MPI_Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      MPI_Sendrecv(Buffer_Send_Neighbor, nBufferS_Vector, MPI_UNSIGNED_SHORT, send_to, 1,
                   Buffer_Receive_Neighbor, nBufferR_Vector, MPI_UNSIGNED_SHORT, receive_from, 1, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  int send_to, receive_from;

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new double [nBufferR_Vector];
      Buffer_Send_Lambda = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  int send_to, receive_from;

  double **Gradient = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new double[nDim];

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  int send_to, receive_from;

  double *Limiter = new double [nVar];

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new double [nBufferR_Vector];
      Buffer_Send_Limit = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  int send_to, receive_from;

  double **Gradient = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Gradient[iVar] = new double[nDim];

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  int send_to, receive_from;

  double *Limiter = new double [nPrimVarGrad];

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nPrimVarGrad;        nBufferR_Vector = nVertexR*nPrimVarGrad;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new double [nBufferR_Vector];
      Buffer_Send_Limit = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Primitive(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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

void CEulerSolver::Set_MPI_Secondary_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  int send_to, receive_from;

  double **Gradient = new double* [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
    Gradient[iVar] = new double[nDim];

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSecondaryVarGrad*nDim;        nBufferR_Vector = nVertexR*nSecondaryVarGrad*nDim;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nSecondaryVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Secondary(iVar, iDim);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
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
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];

        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nSecondaryVarGrad*nVertexR+iVar*nVertexR+iVertex];
          }
        }

        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient_Secondary(iVar, iDim, Gradient[iVar][iDim]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;

    }

  }

  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;

}

void CEulerSolver::Set_MPI_Secondary_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  double *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  int send_to, receive_from;

  double *Limiter = new double [nSecondaryVarGrad];

#ifdef HAVE_MPI
  MPI_Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSecondaryVarGrad;        nBufferR_Vector = nVertexR*nSecondaryVarGrad;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new double [nBufferR_Vector];
      Buffer_Send_Limit = new double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter_Secondary(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          Limiter[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          node[iPoint]->SetLimiter_Secondary(iVar, Limiter[iVar]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;

    }

  }

  delete [] Limiter;

}

void CEulerSolver::SetNondimensionalization(CGeometry *geometry, CConfig *config, unsigned short iMesh) {

  double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0, Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0,
  Velocity_Reynolds = 0.0, Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Tke_FreeStream = 0.0;
  double Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Energy_Ref= 0.0, Froude = 0.0;
  double Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0, Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0, Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;
  unsigned short iDim;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Local variables and memory allocation ---*/

  double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  double Mach             = config->GetMach();
  double Reynolds         = config->GetReynolds();
  bool compressible       = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface        = (config->GetKind_Regime() == FREESURFACE);
  bool unsteady           = (config->GetUnsteady_Simulation() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = config->GetKind_Solver() == RANS;
  bool tkeNeeded          = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));
  bool fs_temperature     = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool ideal_gas       = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );

  if (compressible) {

    /*--- Compute the Free Stream velocity, using the Mach number ---*/

    Pressure_FreeStream = config->GetPressure_FreeStream();
	Density_FreeStream  = config->GetDensity_FreeStream();
	Temperature_FreeStream  = config->GetTemperature_FreeStream();
	switch (config->GetKind_FluidModel()) {

	case STANDARD_AIR:
		if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
		else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(53.3533);

		FluidModel = new CIdealGas(1.4, config->GetGas_Constant());
			if(fs_temperature){
				FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
				Density_FreeStream = FluidModel->GetDensity();
			}else {
				FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
				Temperature_FreeStream = FluidModel->GetTemperature();
			}
			break;

	case IDEAL_GAS:
		FluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
		if(fs_temperature){
			FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
			Density_FreeStream = FluidModel->GetDensity();
		}else {
			FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
			Temperature_FreeStream = FluidModel->GetTemperature();
		}
		break;

	case VW_GAS:
		FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(), config->GetTemperature_Critical());
		if(fs_temperature){
			FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
			Density_FreeStream = FluidModel->GetDensity();
		} else {

			FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
			Temperature_FreeStream = FluidModel->GetTemperature();

		}
		break;

	case PR_GAS:
		FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(), config->GetTemperature_Critical(), config->GetAcentric_Factor());
		if(fs_temperature){
			FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
			Density_FreeStream = FluidModel->GetDensity();
		} else {
			FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
			Temperature_FreeStream = FluidModel->GetTemperature();

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

    ModVel_FreeStream = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
    ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);


    if (viscous) {

    	if (ideal_gas) {

    		/*--- First, check if there is mesh motion. If yes, use the Mach
         number relative to the body to initialize the flow. ---*/

    		if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
    		else Velocity_Reynolds = ModVel_FreeStream;

    		/*--- For viscous flows, pressure will be computed from a density
         that is found from the Reynolds number. The viscosity is computed
         from the dimensional version of Sutherland's law ---*/

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

    		FluidModel->SetViscosityModel(config);

    		Viscosity_FreeStream = FluidModel->GetLaminarViscosity(Temperature_FreeStream, Density_FreeStream);
    		config->SetViscosity_FreeStream(Viscosity_FreeStream);


    		Density_FreeStream   = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds()); config->SetDensity_FreeStream(Density_FreeStream);
    		config->SetDensity_FreeStream(Density_FreeStream);
    		FluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
    		Pressure_FreeStream  = FluidModel->GetPressure();
    		config->SetPressure_FreeStream(Pressure_FreeStream);
    		Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    	} else {
    		FluidModel->SetViscosityModel(config);
    		Viscosity_FreeStream = FluidModel->GetLaminarViscosity(Temperature_FreeStream, Density_FreeStream);
        config->SetViscosity_FreeStream(Viscosity_FreeStream);
    		Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
    	}

      /*--- Turbulence quantities ---*/

      Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
      Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());

    }
    else {

      /*--- For inviscid flow, energy is calculated from the specified
       FreeStream quantities using the proper gas law. ---*/

      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*-- Compute the freestream energy. ---*/

    if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

    /*--- Additional reference values defined by Pref, Tref, Rho_ref. By definition,
     Lref is one because we have converted the grid to meters.---*/

    Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
    Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
    Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
    Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
    Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;                        config->SetForce_Ref(Force_Ref);
    Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
    Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
    Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);         config->SetFroude(Froude);

  }

  else {
    /*--- Incompressible_case */

    /*--- Reference length = 1 (by default)
     Reference density = liquid density or freestream
     Reference viscosity = liquid viscosity or freestream
     Reference velocity = liquid velocity or freestream
     Reference pressure = Reference density * Reference velocity * Reference velocity
     Reynolds number based on the liquid or reference viscosity ---*/

    Pressure_FreeStream = 0.0; config->SetPressure_FreeStream(Pressure_FreeStream);
    Density_FreeStream  = config->GetDensity_FreeStream();
    ModVel_FreeStream = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
    ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

    /*--- Additional reference values defined by Pref, Tref, Rho_ref. By definition,
     Lref is one because we have converted the grid to meters.---*/

    Length_Ref = config->GetLength_Reynolds();                config->SetLength_Ref(Length_Ref);
    Density_Ref = Density_FreeStream;                         config->SetDensity_Ref(Density_Ref);
    Velocity_Ref = ModVel_FreeStream;                         config->SetVelocity_Ref(Velocity_Ref);
    Pressure_Ref = Density_Ref*(Velocity_Ref*Velocity_Ref);   config->SetPressure_Ref(Pressure_Ref);

    if (viscous) {
      Viscosity_FreeStream = config->GetViscosity_FreeStream();
      Reynolds = Density_Ref*Velocity_Ref*Length_Ref / Viscosity_FreeStream;  config->SetReynolds(Reynolds);
      Viscosity_Ref = Viscosity_FreeStream * Reynolds;                        config->SetViscosity_Ref(Viscosity_Ref);
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

  if(viscous){
	  config->SetMu_RefND(config->GetMu_RefND()/Viscosity_Ref);
	  config->SetMu_SND(config->GetMu_SND()/config->GetTemperature_Ref());
	  config->SetMu_Temperature_RefND(config->GetMu_Temperature_RefND()/config->GetTemperature_Ref());
	  FluidModel->SetViscosityModel(config);

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
      cout << "Note: Negative pressure, temperature or density is not allowed!" << endl;
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
      if (viscous) cout << "Reynolds number: " << config->GetReynolds() << ", computed using free-stream values."<<endl;
      cout << "Only dimensional computation, the grid should be dimensional." << endl;
    }

    cout <<"-- Input conditions:"<< endl;

    if (compressible) {
      cout << "Specific gas constant: " << config->GetGas_Constant();
      if (config->GetSystemMeasurements() == SI) cout << " N.m/kg.K." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " lbf.ft/slug.R." << endl;
    }

    if (incompressible || freesurface) {
      cout << "Bulk modulus: " << config->GetBulk_Modulus();
      if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
      cout << "Artificial compressibility factor: " << config->GetArtComp_Factor();
      if (config->GetSystemMeasurements() == SI) cout << " Pa." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " psf." << endl;
    }

    cout << "Free-stream pressure: " << config->GetPressure_FreeStream();
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
      if (turbulent){
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
    }

    if (unsteady) cout << "Reference time: " << config->GetTime_Ref() <<" s." << endl;

    /*--- Print out resulting non-dim values here. ---*/

    cout << "-- Resulting non-dimensional state:" << endl;
    cout << "Mach number (non-dim): " << config->GetMach() << endl;
    if (viscous) {
      cout << "Reynolds number (non-dim): " << config->GetReynolds() <<". Re length: " << config->GetLength_Reynolds();
      if (config->GetSystemMeasurements() == SI) cout << " m." << endl;
      else if (config->GetSystemMeasurements() == US) cout << " in." << endl;
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
      if (turbulent){
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
  double Density, Pressure, yFreeSurface, PressFreeSurface, Froude, yCoord, Velx, Vely, Velz, RhoVelx, RhoVely, RhoVelz, XCoord, YCoord,
  ZCoord, DensityInc, ViscosityInc, Heaviside, LevelSet, lambda, DensityFreeSurface, Area_Children, Area_Parent, LevelSet_Fine, epsilon,
  *Solution_Fine, *Solution, PressRef, yCoordRef;

  unsigned short nDim = geometry[MESH_0]->GetnDim();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool rans = ((config->GetKind_Solver() == RANS) ||
               (config->GetKind_Solver() == ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool aeroelastic = config->GetAeroelastic_Simulation();
  bool gravity     = (config->GetGravityForce() == YES);

  /*--- Set the location and value of the free-surface ---*/

  if (freesurface) {

    for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {

      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {

        /*--- Set initial boundary condition at iter 0 ---*/
        if ((ExtIter == 0) && (!restart)) {

          /*--- Compute the level set value in all the MG levels (basic case, distance to
           the Y/Z plane, and interpolate the solution to the coarse levels ---*/
          if (iMesh == MESH_0) {
            XCoord = geometry[iMesh]->node[iPoint]->GetCoord(0);
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
          DensityFreeSurface = solver_container[iMesh][FLOW_SOL]->GetDensity_Inf();
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

  if (!freesurface && gravity) {

    for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {

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

  /*--- Set subsonic initial condition for Engine intakes ---*/

  if (config->GetEngine_Intake()) {

    /*--- Set initial boundary condition at iteration 0 ---*/
    if ((ExtIter == 0) && (!restart)) {

      double *Coord;
      double Velocity_FreeStream[3] = {0.0, 0.0, 0.0}, Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStream, Density_FreeStream, Pressure_FreeStream, Density_FreeStreamND, Pressure_FreeStreamND, ModVel_FreeStreamND, Energy_FreeStreamND, ModVel_FreeStream, T_ref = 0.0, S = 0.0, Mu_ref = 0.0;

      double Mach = 0.40;
      double Alpha = config->GetAoA()*PI_NUMBER/180.0;
      double Beta  = config->GetAoS()*PI_NUMBER/180.0;
      double Gamma_Minus_One = Gamma - 1.0;
      double Gas_Constant = config->GetGas_ConstantND();
      double Temperature_FreeStream = config->GetTemperature_FreeStream();
      double Mach2Vel_FreeStream = sqrt(Gamma*Gas_Constant*Temperature_FreeStream);

      for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {

        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {

          Velocity_FreeStream[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
          Velocity_FreeStream[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
          Velocity_FreeStream[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;

          ModVel_FreeStream = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            ModVel_FreeStream += Velocity_FreeStream[iDim]*Velocity_FreeStream[iDim];
          ModVel_FreeStream = sqrt(ModVel_FreeStream);

          if (config->GetViscous()) {
            if (config->GetSystemMeasurements() == SI) { T_ref = 273.15; S = 110.4; Mu_ref = 1.716E-5; }
            if (config->GetSystemMeasurements() == US) { T_ref = 518.7; S = 198.72; Mu_ref = 3.62E-7; }
            Viscosity_FreeStream = Mu_ref*(pow(Temperature_FreeStream/T_ref, 1.5) * (T_ref+S)/(Temperature_FreeStream+S));
            Density_FreeStream   = config->GetReynolds()*Viscosity_FreeStream/(ModVel_FreeStream*config->GetLength_Reynolds());
            Pressure_FreeStream  = Density_FreeStream*Gas_Constant*Temperature_FreeStream;
          }
          else {
            Pressure_FreeStream  = config->GetPressure_FreeStreamND();
            Density_FreeStream  = Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);
          }

          Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();
          Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity_FreeStreamND[iDim] = Velocity_FreeStream[iDim]/config->GetVelocity_Ref();

          ModVel_FreeStreamND = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
          ModVel_FreeStreamND = sqrt(ModVel_FreeStreamND);

          Energy_FreeStreamND = Pressure_FreeStreamND/(Density_FreeStreamND*Gamma_Minus_One)+0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

          Coord = geometry[iMesh]->node[iPoint]->GetCoord();

          if (((Coord[0] >= 16.0) && (Coord[0] <= 20.0)) &&
              ((Coord[1] >= 0.0) && (Coord[1] <= 0.7)) &&
              ((Coord[2] >= 2.5) && (Coord[2] <= 4.0))) {

            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(0, Density_FreeStreamND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(iDim+1, Velocity_FreeStreamND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(nVar-1, Energy_FreeStreamND);

            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(0, Density_FreeStreamND);
            for (iDim = 0; iDim < nDim; iDim++)
              solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(iDim+1, Velocity_FreeStreamND[iDim]);
            solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution_Old(nVar-1, Energy_FreeStreamND);
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
    Solution = new double[nVar];
    for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
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
      Solution = new double[nVar_Turb];
      for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
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

  if (dual_time && (ExtIter == 0 || (restart && ExtIter == config->GetUnst_RestartIter()))) {

    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/
    for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
        if (rans) {
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
        }
      }
    }

    if ((restart && ExtIter == config->GetUnst_RestartIter()) &&
        (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {

      /*--- Load an additional restart file for a 2nd-order restart ---*/
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, int(config->GetUnst_RestartIter()-1));

      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, int(config->GetUnst_RestartIter()-1));

      /*--- Push back this new solution to time level N. ---*/
      for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          if (rans) {
            solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }
  }

  if (aeroelastic) {
    /*--- Reset the plunge and pitch value for the new unsteady step. ---*/
    unsigned short iMarker_Monitoring;
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
      config->SetAeroelastic_pitch(iMarker_Monitoring,0.0);
      config->SetAeroelastic_plunge(iMarker_Monitoring,0.0);
    }
  }

}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;
  int rank;

#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
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
  bool compressible     = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible   = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface      = (config->GetKind_Regime() == FREESURFACE);
  bool engine           = ((config->GetnMarker_NacelleInflow() != 0) || (config->GetnMarker_NacelleExhaust() != 0));
  bool fixed_cl         = config->GetFixed_CL_Mode();
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );

  /*--- Compute nacelle inflow and exhaust properties ---*/

  if (engine) { GetNacelle_Properties(geometry, config, iMesh, Output); }

  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/

  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }

  /*--- Compute distance function to zero level set (Set LevelSet and Distance primitive variables)---*/

  if (freesurface) SetFreeSurface_Distance(geometry, config);

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
     FreeSurface Incompressible flow, primitive variables nDim+5, (P,vx,vy,vz,rho,beta,LevelSet,Dist),
     Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/

    if (compressible) {
    	RightSol = node[iPoint]->SetPrimVar_Compressible(FluidModel);
    	node[iPoint]->SetSecondaryVar_Compressible(FluidModel);
    }
    if (incompressible) { RightSol = node[iPoint]->SetPrimVar_Incompressible(Density_Inf, config); }
    if (freesurface) {    RightSol = node[iPoint]->SetPrimVar_FreeSurface(config); }
    if (!RightSol) ErrorCounter++;

    /*--- Initialize the residual vector (except for output of the residuals) ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  /*--- Upwind second order reconstruction ---*/

  if ((second_order && !center) && ((iMesh == MESH_0) || low_fidelity)) {

    /*--- Gradient computation ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS){
    	SetPrimitive_Gradient_GG(geometry, config);
    	if (compressible && !ideal_gas) SetSecondary_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES){
    	SetPrimitive_Gradient_LS(geometry, config);
    	if (compressible && !ideal_gas) SetSecondary_Gradient_LS(geometry, config);
    }


    /*--- Limiter computation ---*/

    if ((limiter) && (iMesh == MESH_0)){
    	SetPrimitive_Limiter(geometry, config);
    	if (compressible && !ideal_gas) SetSecondary_Limiter(geometry, config);
    }

  }

  /*--- Artificial dissipation ---*/

  if (center) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && ((iMesh == MESH_0) || low_fidelity)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Initialize the Jacobian matrices ---*/

  if (implicit) Jacobian.SetValZero();

  /*--- Error message ---*/
#ifdef HAVE_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
  MPI_Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);


}

void CEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) { }

void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {

  double *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Mean_BetaInc2, Lambda, Local_Delta_Time, Mean_DensityInc, Mean_LevelSet,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j, Delta = 0.0, a, b, c, e, f;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  double epsilon = config->GetFreeSurface_Thickness();
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
      double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      double *GridVel_j = geometry->node[jPoint]->GetGridVel();
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
        double *GridVel = geometry->node[iPoint]->GetGridVel();
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
    Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
    Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
    Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
    Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
    node[iPoint]->SetDelta_Time(Local_Delta_Time);
  }

  /*--- Check if there is any element with only one neighbor...
   a CV that is inside another CV ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (geometry->node[iPoint]->GetnPoint() == 1)
      node[iPoint]->SetDelta_Time(Min_Delta_Time);
  }

  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    MPI_Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for(iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }

  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

#ifdef HAVE_MPI
    double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    MPI_Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
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
      Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
      Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
      Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
      Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
    }
  }

}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {

  double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j,
  *V_i, *V_j, *S_i, *S_j, *Limiter_i = NULL, *Limiter_j = NULL, YDistance, GradHidrosPress, sqvel;
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;
  bool neg_density_i = false, neg_density_j = false, neg_pressure_i = false, neg_pressure_j = false;

  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool low_fidelity     = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool second_order     = (((config->GetSpatialOrder_Flow() == SECOND_ORDER) || (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)) && ((iMesh == MESH_0) || low_fidelity));
  bool limiter          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && !low_fidelity);
  bool freesurface      = (config->GetKind_Regime() == FREESURFACE);
  bool compressible     = (config->GetKind_Regime() == COMPRESSIBLE);
  bool grid_movement    = config->GetGrid_Movement();
  bool roe_turkel       = (config->GetKind_Upwind_Flow() == TURKEL);

  for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

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

    /*--- High order reconstruction using MUSCL strategy ---*/

    if (second_order && !freesurface) {

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
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
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

      /*--- Check for non-physical solutions after reconstruction. If found,
       use the cell-average value of the solution. This results in a locally
       first-order approximation, but this is typically only active
       during the start-up of a calculation. ---*/

      if (compressible) {
        neg_pressure_i = (Primitive_i[nDim+1] < 0.0);
        neg_pressure_j = (Primitive_j[nDim+1] < 0.0);
        neg_density_i  = (Primitive_i[nDim+2] < 0.0);
        neg_density_j  = (Primitive_j[nDim+2] < 0.0);
      }

      /*--- If non-physical, use the cell-averaged state. ---*/

      if (neg_density_i || neg_pressure_i) {
        for (iVar = 0; iVar < nVar; iVar++)
          Primitive_i[iVar] = V_i[iVar];
        counter_local++;
      }
      if (neg_density_j || neg_pressure_j) {
        for (iVar = 0; iVar < nVar; iVar++)
          Primitive_j[iVar] = V_j[iVar];
        counter_local++;
      }

      /*--- If compressible, compute 2nd order reconstruction for the secondary variables ---*/

      if (compressible) {

        Gradient_i = node[iPoint]->GetGradient_Secondary();
        Gradient_j = node[jPoint]->GetGradient_Secondary();
        if (limiter) {
          Limiter_i = node[iPoint]->GetLimiter_Secondary();
          Limiter_j = node[jPoint]->GetLimiter_Secondary();
        }

        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
          Project_Grad_i = 0.0; Project_Grad_j = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
          }
          if (limiter) {
            Secondary_i[iVar] = S_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
            Secondary_j[iVar] = S_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
          }
          else {
            Secondary_i[iVar] = S_i[iVar] + Project_Grad_i;
            Secondary_j[iVar] = S_j[iVar] + Project_Grad_j;
          }
        }

      }

      /*--- Set conservative variables with reconstruction ---*/

      numerics->SetPrimitive(Primitive_i, Primitive_j);
      numerics->SetSecondary(Secondary_i, Secondary_j);


    } else {

      /*--- Set conservative variables without reconstruction ---*/

      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);

    }

    /*--- Free surface simulation should include gradient of the hydrostatic pressure ---*/

    if (freesurface) {

      /*--- The zero order reconstruction includes the gradient
       of the hydrostatic pressure constribution ---*/

      YDistance = 0.5*(geometry->node[jPoint]->GetCoord(nDim-1)-geometry->node[iPoint]->GetCoord(nDim-1));
      GradHidrosPress = node[iPoint]->GetDensityInc()/(config->GetFroude()*config->GetFroude());
      Primitive_i[0] = V_i[0] - GradHidrosPress*YDistance;
      GradHidrosPress = node[jPoint]->GetDensityInc()/(config->GetFroude()*config->GetFroude());
      Primitive_j[0] = V_j[0] + GradHidrosPress*YDistance;

      /*--- Copy the rest of primitive variables ---*/

      for (iVar = 1; iVar < nPrimVar; iVar++) {
        Primitive_i[iVar] = V_i[iVar]+EPS;
        Primitive_j[iVar] = V_j[iVar]+EPS;
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

        /*--- Note that the pressure reconstruction always includes the hydrostatic gradient,
         and we should limit only the kinematic contribution ---*/

        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
          Project_Grad_i = 0.0; Project_Grad_j = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
          }
          if (limiter) {
            if (iVar == 0) {
              Primitive_i[iVar] += Limiter_i[iVar]*(V_i[iVar] + Project_Grad_i - Primitive_i[iVar]);
              Primitive_j[iVar] += Limiter_j[iVar]*(V_j[iVar] + Project_Grad_j - Primitive_j[iVar]);
            }
            else {
              Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
              Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
            }
          }
          else {
            Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
            Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
          }

        }

      }

      /*--- Set primitive variables with reconstruction ---*/

      numerics->SetPrimitive(Primitive_i, Primitive_j);

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
#ifdef HAVE_MPI
  MPI_Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
  counter_global = counter_local;
#endif
  if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);

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
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[iPoint]->GetCoord());

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
    double Vol, x_o, x_od, x, z, levelset, DampingFactor;
    double factor = config->GetFreeSurface_Damping_Length();

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

    double Volume, Source;

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
  double *Normal, Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Mean_BetaInc2, Lambda, Mean_DensityInc,
  ProjVel, ProjVel_i, ProjVel_j;
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
      double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      double *GridVel_j = geometry->node[jPoint]->GetGridVel();
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
        double *GridVel = geometry->node[iPoint]->GetGridVel();
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
  double Pressure_i = 0, Pressure_j = 0, *Diff;
  unsigned short iVar;
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);

  Diff = new double[nVar];

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

#ifdef STRUCTURED_GRID

    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);

#else

    bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();

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

#endif

  }

#ifdef STRUCTURED_GRID

  unsigned long Point_Normal = 0, iVertex;
  double Pressure_mirror = 0, *U_mirror;
  unsigned short iMarker;

  U_mirror = new double[nVar];

  /*--- Loop over all boundaries and include an extra contribution
   from a mirror node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

          /*--- Interpolate & compute difference in the conserved variables ---*/
          for (iVar = 0; iVar < nVar; iVar++) {
            U_mirror[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - node[Point_Normal]->GetSolution(iVar);
            Diff[iVar]   = node[iPoint]->GetSolution(iVar) - U_mirror[iVar];
          }

          /*--- Correction for compressible flows ---*/
          if (compressible) {
            Pressure_mirror = 2.0*node[iPoint]->GetPressure() - node[Point_Normal]->GetPressure();
            Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + node[iPoint]->GetPressure()) - (U_mirror[nVar-1] + Pressure_mirror);
          }

          /*--- Subtract contribution at the boundary node only ---*/
          node[iPoint]->SubtractUnd_Lapl(Diff);
        }
      }
    }
  }

  delete [] U_mirror;

#endif

  delete [] Diff;

  /*--- MPI parallelization ---*/
  Set_MPI_Undivided_Laplacian(geometry, config);

}

void CEulerSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {

  unsigned long iEdge, iPoint, jPoint;
  double Pressure_i = 0.0, Pressure_j = 0.0;

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

#ifdef STRUCTURED_GRID

    /*--- Compute numerator and denominator ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i);
      jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j);
    }
    if (geometry->node[jPoint]->GetDomain()) {
      iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j);
      jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j);
    }

#else

    bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();

    /*--- Both points inside the domain, or both on the boundary ---*/
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)){
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    }

    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }

    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }

#endif

  }

#ifdef STRUCTURED_GRID
  unsigned short iMarker;
  unsigned long iVertex, Point_Normal;
  double Press_mirror;

  /*--- Loop over all boundaries and include an extra contribution
   from a mirror node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

          if (compressible) Pressure_i = node[iPoint]->GetPressure();
          if (incompressible || freesurface) Pressure_i = node[iPoint]->GetDensityInc();

          /*--- Interpolate & compute difference in the conserved variables ---*/
          if (compressible) {
            Pressure_i = node[iPoint]->GetPressure();
            Press_mirror = 2.0*Pressure_i - node[Point_Normal]->GetPressure();
          }
          if (incompressible || freesurface) {
            Pressure_i = node[iPoint]->GetDensityInc();
            Press_mirror = 2.0*Pressure_i - node[Point_Normal]->GetDensityInc();
          }

          /*--- Add contribution at the boundary node only ---*/
          iPoint_UndLapl[iPoint] += (Press_mirror - Pressure_i);
          jPoint_UndLapl[iPoint] += (Pressure_i + Press_mirror);

        }
      }
    }
  }

#endif

  /*--- Set pressure switch for each point ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetSensor(fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);

  /*--- MPI parallelization ---*/
  Set_MPI_Dissipation_Switch(geometry, config);

}

void CEulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  double Pressure = 0.0, *Normal = NULL, MomentDist[3], *Coord, Area,
  factor, NFPressOF, RefVel2, RefDensity, RefPressure, Gas_Constant, Mach2Vel, Mach_Motion, UnitNormal[3], Force[3];
  double *Origin = config->GetRefOriginMoment(0);
  string Marker_Tag, Monitoring_Tag;

  bool grid_movement      = config->GetGrid_Movement();
  double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  double RefAreaCoeff     = config->GetRefAreaCoeff();
  double RefLengthMoment  = config->GetRefLengthMoment();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);

  /*--- For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  if (grid_movement) {
    Gas_Constant = config->GetGas_ConstantND();
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  RefDensity  = Density_Inf;
  RefPressure = Pressure_Inf;

  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);

  /*-- Variables initialization ---*/

  Total_CDrag = 0.0;  Total_CLift = 0.0;    Total_CSideForce = 0.0;   Total_CEff = 0.0;
  Total_CMx = 0.0;    Total_CMy = 0.0;      Total_CMz = 0.0;
  Total_CFx = 0.0;    Total_CFy = 0.0;      Total_CFz = 0.0;
  Total_CT = 0.0;     Total_CQ = 0.0;       Total_CMerit = 0.0;
  Total_CNearFieldOF = 0.0;
  Total_Heat = 0.0;  Total_MaxHeat = 0.0;

  AllBound_CDrag_Inv = 0.0;   AllBound_CLift_Inv = 0.0; AllBound_CSideForce_Inv = 0.0;   AllBound_CEff_Inv = 0.0;
  AllBound_CMx_Inv = 0.0;     AllBound_CMy_Inv = 0.0;   AllBound_CMz_Inv = 0.0;
  AllBound_CFx_Inv = 0.0;     AllBound_CFy_Inv = 0.0;   AllBound_CFz_Inv = 0.0;
  AllBound_CT_Inv = 0.0;      AllBound_CQ_Inv = 0.0;    AllBound_CMerit_Inv = 0.0;
  AllBound_CNearFieldOF_Inv = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift_Inv[iMarker_Monitoring] = 0.0;
    Surface_CDrag_Inv[iMarker_Monitoring] = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CLift[iMarker_Monitoring]     = 0.0;
    Surface_CDrag[iMarker_Monitoring]     = 0.0;
    Surface_CMx[iMarker_Monitoring]       = 0.0;
    Surface_CMy[iMarker_Monitoring]       = 0.0;
    Surface_CMz[iMarker_Monitoring]       = 0.0;
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

      CDrag_Inv[iMarker] = 0.0;         CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;  CEff_Inv[iMarker] = 0.0;
      CMx_Inv[iMarker] = 0.0;           CMy_Inv[iMarker] = 0.0;   CMz_Inv[iMarker] = 0.0;
      CFx_Inv[iMarker] = 0.0;           CFy_Inv[iMarker] = 0.0;   CFz_Inv[iMarker] = 0.0;
      CT_Inv[iMarker] = 0.0;            CQ_Inv[iMarker] = 0.0;    CMerit_Inv[iMarker] = 0.0;
      CNearFieldOF_Inv[iMarker] = 0.0;

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
            UnitNormal[iDim] = Normal[iDim]/Area;
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }

          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/

          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf)*Normal[iDim]*factor;
            ForceInviscid[iDim] += Force[iDim];
          }

          /*--- Moment with respect to the reference axis ---*/

          if (iDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLengthMoment;
            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLengthMoment;

        }

      }

      /*--- Project forces and store the non-dimensional coefficients ---*/

      if  (Monitoring == YES) {
        if (nDim == 2) {
          if (Boundary != NEARFIELD_BOUNDARY) {
            CDrag_Inv[iMarker]  =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            CLift_Inv[iMarker]  = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
            CEff_Inv[iMarker]   = CLift_Inv[iMarker] / (CDrag_Inv[iMarker]+config->GetCteViscDrag()+EPS);
            CMz_Inv[iMarker]    = MomentInviscid[2];
            CFx_Inv[iMarker]    = ForceInviscid[0];
            CFy_Inv[iMarker]    = ForceInviscid[1];
            CT_Inv[iMarker]     = -CFx_Inv[iMarker];
            CQ_Inv[iMarker]     = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker] = CT_Inv[iMarker]/CQ_Inv[iMarker];
          }
          else { CNearFieldOF_Inv[iMarker] = NFPressOF; }
        }
        if (nDim == 3) {
          if (Boundary != NEARFIELD_BOUNDARY) {
            CDrag_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
            CLift_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
            CSideForce_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            CEff_Inv[iMarker]       = CLift_Inv[iMarker] / (CDrag_Inv[iMarker]+config->GetCteViscDrag()+EPS);
            CMx_Inv[iMarker]        = MomentInviscid[0];
            CMy_Inv[iMarker]        = MomentInviscid[1];
            CMz_Inv[iMarker]        = MomentInviscid[2];
            CFx_Inv[iMarker]        = ForceInviscid[0];
            CFy_Inv[iMarker]        = ForceInviscid[1];
            CFz_Inv[iMarker]        = ForceInviscid[2];
            CT_Inv[iMarker]         = -CFz_Inv[iMarker];
            CQ_Inv[iMarker]         = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker]     = CT_Inv[iMarker] / CQ_Inv[iMarker];
          }
          else { CNearFieldOF_Inv[iMarker] = NFPressOF; }
        }

        AllBound_CDrag_Inv        += CDrag_Inv[iMarker];
        AllBound_CLift_Inv        += CLift_Inv[iMarker];
        AllBound_CSideForce_Inv   += CSideForce_Inv[iMarker];
        AllBound_CMx_Inv          += CMx_Inv[iMarker];
        AllBound_CMy_Inv          += CMy_Inv[iMarker];
        AllBound_CMz_Inv          += CMz_Inv[iMarker];
        AllBound_CFx_Inv          += CFx_Inv[iMarker];
        AllBound_CFy_Inv          += CFy_Inv[iMarker];
        AllBound_CFz_Inv          += CFz_Inv[iMarker];
        AllBound_CT_Inv           += CT_Inv[iMarker];
        AllBound_CQ_Inv           += CQ_Inv[iMarker];
        AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];

        AllBound_CEff_Inv = AllBound_CLift_Inv / (AllBound_CDrag_Inv + config->GetCteViscDrag() + EPS);
        AllBound_CMerit_Inv = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);

        /*--- Compute the coefficients per surface ---*/

        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CLift_Inv[iMarker_Monitoring] += CLift_Inv[iMarker];
            Surface_CDrag_Inv[iMarker_Monitoring] += CDrag_Inv[iMarker];
            Surface_CMx_Inv[iMarker_Monitoring]   += CMx_Inv[iMarker];
            Surface_CMy_Inv[iMarker_Monitoring]   += CMy_Inv[iMarker];
            Surface_CMz_Inv[iMarker_Monitoring]   += CMz_Inv[iMarker];
          }
        }

      }


    }
  }

#ifdef HAVE_MPI

  /*--- Add AllBound information using all the nodes ---*/

  double MyAllBound_CDrag_Inv        = AllBound_CDrag_Inv;        AllBound_CDrag_Inv = 0.0;
  double MyAllBound_CLift_Inv       = AllBound_CLift_Inv;        AllBound_CLift_Inv = 0.0;
  double MyAllBound_CSideForce_Inv   = AllBound_CSideForce_Inv;   AllBound_CSideForce_Inv = 0.0;
  double MyAllBound_CEff_Inv         = AllBound_CEff_Inv;         AllBound_CEff_Inv = 0.0;
  double MyAllBound_CMx_Inv          = AllBound_CMx_Inv;          AllBound_CMx_Inv = 0.0;
  double MyAllBound_CMy_Inv          = AllBound_CMy_Inv;          AllBound_CMy_Inv = 0.0;
  double MyAllBound_CMz_Inv          = AllBound_CMz_Inv;          AllBound_CMz_Inv = 0.0;
  double MyAllBound_CFx_Inv          = AllBound_CFx_Inv;          AllBound_CFx_Inv = 0.0;
  double MyAllBound_CFy_Inv         = AllBound_CFy_Inv;          AllBound_CFy_Inv = 0.0;
  double MyAllBound_CFz_Inv          = AllBound_CFz_Inv;          AllBound_CFz_Inv = 0.0;
  double MyAllBound_CT_Inv           = AllBound_CT_Inv;           AllBound_CT_Inv = 0.0;
  double MyAllBound_CQ_Inv           = AllBound_CQ_Inv;           AllBound_CQ_Inv = 0.0;
  double MyAllBound_CMerit_Inv       = AllBound_CMerit_Inv;       AllBound_CMerit_Inv = 0.0;
  double MyAllBound_CNearFieldOF_Inv = AllBound_CNearFieldOF_Inv; AllBound_CNearFieldOF_Inv = 0.0;

  MPI_Allreduce(&MyAllBound_CDrag_Inv, &AllBound_CDrag_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CLift_Inv, &AllBound_CLift_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CSideForce_Inv, &AllBound_CSideForce_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Inv = AllBound_CLift_Inv / (AllBound_CDrag_Inv + config->GetCteViscDrag() + EPS);
  MPI_Allreduce(&MyAllBound_CMx_Inv, &AllBound_CMx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CMy_Inv, &AllBound_CMy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CMz_Inv, &AllBound_CMz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFx_Inv, &AllBound_CFx_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFy_Inv, &AllBound_CFy_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFz_Inv, &AllBound_CFz_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CT_Inv, &AllBound_CT_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CQ_Inv, &AllBound_CQ_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Inv = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
  MPI_Allreduce(&MyAllBound_CNearFieldOF_Inv, &AllBound_CNearFieldOF_Inv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /*--- Add the forces on the surfaces using all the nodes ---*/

  double *MySurface_CLift_Inv = NULL;
  double *MySurface_CDrag_Inv = NULL;
  double *MySurface_CMx_Inv   = NULL;
  double *MySurface_CMy_Inv   = NULL;
  double *MySurface_CMz_Inv   = NULL;

  MySurface_CLift_Inv = new double[config->GetnMarker_Monitoring()];
  MySurface_CDrag_Inv = new double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Inv   = new double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Inv   = new double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Inv   = new double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CLift_Inv[iMarker_Monitoring] = Surface_CLift_Inv[iMarker_Monitoring];
    MySurface_CDrag_Inv[iMarker_Monitoring] = Surface_CDrag_Inv[iMarker_Monitoring];
    MySurface_CMx_Inv[iMarker_Monitoring]   = Surface_CMx_Inv[iMarker_Monitoring];
    MySurface_CMy_Inv[iMarker_Monitoring]   = Surface_CMy_Inv[iMarker_Monitoring];
    MySurface_CMz_Inv[iMarker_Monitoring]   = Surface_CMz_Inv[iMarker_Monitoring];
    Surface_CLift_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CDrag_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]     = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]     = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]     = 0.0;
  }

  MPI_Allreduce(MySurface_CLift_Inv, Surface_CLift_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CDrag_Inv, Surface_CDrag_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CMx_Inv, Surface_CMx_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CMy_Inv, Surface_CMy_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CMz_Inv, Surface_CMz_Inv, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete [] MySurface_CLift_Inv;
  delete [] MySurface_CDrag_Inv;
  delete [] MySurface_CMx_Inv;
  delete [] MySurface_CMy_Inv;
  delete [] MySurface_CMz_Inv;

#endif

  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/

  Total_CDrag         = AllBound_CDrag_Inv;
  Total_CLift         = AllBound_CLift_Inv;
  Total_CSideForce    = AllBound_CSideForce_Inv;
  Total_CEff          = Total_CLift / (Total_CDrag + config->GetCteViscDrag() + EPS);
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
    Surface_CLift[iMarker_Monitoring]     = Surface_CLift_Inv[iMarker_Monitoring];
    Surface_CDrag[iMarker_Monitoring]     = Surface_CDrag_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]       = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]       = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]       = Surface_CMz_Inv[iMarker_Monitoring];
  }

}

void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;

  double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
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
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
      }
    }

  }

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);


}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
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
        AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
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
  double Delta, *local_Res_TruncError, Vol;

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

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
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
        node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
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
  double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal;

  /*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
  PrimVar_Vertex = new double [nPrimVarGrad];
  PrimVar_i = new double [nPrimVarGrad];
  PrimVar_j = new double [nPrimVarGrad];

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

void CEulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config) {

  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
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
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0; detR2 = 0.0;

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
  double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Primitive_i, *Primitive_j,
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

  if (config->GetKind_SlopeLimit() == BARTH_JESPERSEN) {

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

  if (config->GetKind_SlopeLimit() == VENKATAKRISHNAN) {

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

void CEulerSolver::SetSecondary_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  double *SecondaryVar_Vertex, *SecondaryVar_i, *SecondaryVar_j, SecondaryVar_Average,
  Partial_Gradient, Partial_Res, *Normal;

  /*--- Gradient Secondary variables compressible (temp, vx, vy, vz, P, rho)
   Gradient Secondary variables incompressible (rho, vx, vy, vz, beta) ---*/
  SecondaryVar_Vertex = new double [nSecondaryVarGrad];
  SecondaryVar_i = new double [nSecondaryVarGrad];
  SecondaryVar_j = new double [nSecondaryVarGrad];

  /*--- Set Gradient_Secondary to zero ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetGradient_SecondaryZero(nSecondaryVarGrad);

  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
      SecondaryVar_i[iVar] = node[iPoint]->GetSecondary(iVar);
      SecondaryVar_j[iVar] = node[jPoint]->GetSecondary(iVar);
    }

    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
      SecondaryVar_Average =  0.5 * ( SecondaryVar_i[iVar] + SecondaryVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = SecondaryVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddGradient_Secondary(iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Secondary(iVar, iDim, Partial_Res);
      }
    }
  }

  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {

        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          SecondaryVar_Vertex[iVar] = node[iPoint]->GetSecondary(iVar);

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = SecondaryVar_Vertex[iVar]*Normal[iDim];
            node[iPoint]->SubtractGradient_Secondary(iVar, iDim, Partial_Res);
          }
      }
    }
  }

  /*--- Update gradient value ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = node[iPoint]->GetGradient_Secondary(iVar,iDim) / (geometry->node[iPoint]->GetVolume());
        node[iPoint]->SetGradient_Secondary(iVar, iDim, Partial_Gradient);
      }
    }
  }

  delete [] SecondaryVar_Vertex;
  delete [] SecondaryVar_i;
  delete [] SecondaryVar_j;

  Set_MPI_Secondary_Gradient(geometry, config);

}

void CEulerSolver::SetSecondary_Gradient_LS(CGeometry *geometry, CConfig *config) {

  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  double *SecondaryVar_i, *SecondaryVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;

  /*--- Loop over points of the grid ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Set the value of the singular ---*/
    singular = false;

    /*--- Get coordinates ---*/

    Coord_i = geometry->node[iPoint]->GetCoord();

    /*--- Get Secondarys from CVariable ---*/

    SecondaryVar_i = node[iPoint]->GetSecondary();

    /*--- Inizialization of variables ---*/

    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        cvector[iVar][iDim] = 0.0;

    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0; detR2 = 0.0;

    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();

      SecondaryVar_j = node[jPoint]->GetSecondary();

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

        for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(SecondaryVar_j[iVar]-SecondaryVar_i[iVar])/weight;

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
    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
        }

        node[iPoint]->SetGradient_Secondary(iVar, iDim, product);
      }
    }

  }

  Set_MPI_Secondary_Gradient(geometry, config);

}

void CEulerSolver::SetSecondary_Limiter(CGeometry *geometry, CConfig *config) {

  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Secondary_i, *Secondary_j,
  dave, LimK, eps2, dm, dp, du, limiter;

  /*--- Initialize solution max and solution min in the entire domain --*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
      node[iPoint]->SetSolution_Max(iVar, -EPS);
      node[iPoint]->SetSolution_Min(iVar, EPS);
    }
  }

  /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Get the conserved variables ---*/
    Secondary_i = node[iPoint]->GetSecondary();
    Secondary_j = node[jPoint]->GetSecondary();

    /*--- Compute the maximum, and minimum values for nodes i & j ---*/
    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
      du = (Secondary_j[iVar] - Secondary_i[iVar]);
      node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
      node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
      node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
      node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
    }
  }

  /*--- Initialize the limiter --*/
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
      node[iPoint]->SetLimiter_Secondary(iVar, 2.0);
    }
  }

  /*--- Venkatakrishnan limiter ---*/

  if (config->GetKind_SlopeLimit() == VENKATAKRISHNAN) {

    /*-- Get limiter parameters from the configuration file ---*/
    dave = config->GetRefElemLength();
    LimK = config->GetLimiterCoeff();
    eps2 = pow((LimK*dave), 3.0);

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient_Secondary();
      Gradient_j = node[jPoint]->GetGradient_Secondary();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();

      for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);

        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if (limiter < node[iPoint]->GetLimiter_Secondary(iVar))
          if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiter_Secondary(iVar, limiter);

        /*-- Repeat for point j on the edge ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);

        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if (limiter < node[jPoint]->GetLimiter_Secondary(iVar))
          if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiter_Secondary(iVar, limiter);
      }
    }
  }

  /*--- Limiter MPI ---*/
  Set_MPI_Secondary_Limiter(geometry, config);

}

void CEulerSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) {
  unsigned short iDim, jDim, iVar, jVar;
  double Beta, local_Mach, Beta2, rho, enthalpy, soundspeed, sq_vel;
  double *U_i = NULL;
  double Beta_min = config->GetminTurkelBeta();
  double Beta_max = config->GetmaxTurkelBeta();


  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();
  Beta 		    = max(Beta_min,min(local_Mach,Beta_max));
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

void CEulerSolver::GetNacelle_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
  unsigned short iDim, iMarker, iVar;
  unsigned long iVertex, iPoint;
  double Pressure, Velocity[3], Velocity2, MassFlow, Density, Energy, Area,
  Mach, SoundSpeed, Flow_Dir[3], alpha;

  unsigned short nMarker_NacelleInflow = config->GetnMarker_NacelleInflow();
  unsigned short nMarker_NacelleExhaust = config->GetnMarker_NacelleExhaust();

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    FanFace_MassFlow[iMarker] = 0.0;
    FanFace_Mach[iMarker] = 0.0;
    FanFace_Pressure[iMarker] = 0.0;
    FanFace_Area[iMarker] = 0.0;

    Exhaust_MassFlow[iMarker] = 0.0;
    Exhaust_Area[iMarker] = 0.0;

    if (config->GetMarker_All_KindBC(iMarker) == NACELLE_INFLOW) {

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

          /*--- Compute the FanFace_MassFlow, FanFace_Pressure, FanFace_Mach, and FanFace_Area ---*/
          FanFace_MassFlow[iMarker] += MassFlow;
          FanFace_Pressure[iMarker] += Pressure*Area;
          FanFace_Mach[iMarker] += Mach*Area;
          FanFace_Area[iMarker] += Area;

        }
      }

    }

    if (config->GetMarker_All_KindBC(iMarker) == NACELLE_EXHAUST) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iPoint]->GetDomain()) {

          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Vector[iDim]*Vector[iDim];
          Area = sqrt (Area);

          MassFlow = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            MassFlow += Vector[iDim]*node[iPoint]->GetSolution(iDim+1);

          /*--- Compute the mass Exhaust_MassFlow ---*/
          Exhaust_MassFlow[iMarker] += MassFlow;
          Exhaust_Area[iMarker] += Area;

        }
      }

    }

  }

  /*--- Copy to the appropriate structure ---*/
  unsigned short iMarker_NacelleInflow, iMarker_NacelleExhaust;

  double *FanFace_MassFlow_Local = new double [nMarker_NacelleInflow];
  double *FanFace_Mach_Local = new double [nMarker_NacelleInflow];
  double *FanFace_Pressure_Local = new double [nMarker_NacelleInflow];
  double *FanFace_Area_Local = new double [nMarker_NacelleInflow];

  double *FanFace_MassFlow_Total = new double [nMarker_NacelleInflow];
  double *FanFace_Mach_Total = new double [nMarker_NacelleInflow];
  double *FanFace_Pressure_Total = new double [nMarker_NacelleInflow];
  double *FanFace_Area_Total = new double [nMarker_NacelleInflow];

  for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
    FanFace_MassFlow_Local[iMarker_NacelleInflow] = 0.0;
    FanFace_Mach_Local[iMarker_NacelleInflow] = 0.0;
    FanFace_Pressure_Local[iMarker_NacelleInflow] = 0.0;
    FanFace_Area_Local[iMarker_NacelleInflow] = 0.0;

    FanFace_MassFlow_Total[iMarker_NacelleInflow] = 0.0;
    FanFace_Mach_Total[iMarker_NacelleInflow] = 0.0;
    FanFace_Pressure_Total[iMarker_NacelleInflow] = 0.0;
    FanFace_Area_Total[iMarker_NacelleInflow] = 0.0;
  }

  double *Exhaust_MassFlow_Local = new double [nMarker_NacelleExhaust];
  double *Exhaust_Area_Local = new double [nMarker_NacelleExhaust];

  double *Exhaust_MassFlow_Total = new double [nMarker_NacelleExhaust];
  double *Exhaust_Area_Total = new double [nMarker_NacelleExhaust];

  for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
    Exhaust_MassFlow_Local[iMarker_NacelleExhaust] = 0.0;
    Exhaust_Area_Local[iMarker_NacelleExhaust] = 0.0;

    Exhaust_MassFlow_Total[iMarker_NacelleExhaust] = 0.0;
    Exhaust_Area_Total[iMarker_NacelleExhaust] = 0.0;
  }

  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == NACELLE_INFLOW) {

      /*--- Loop over all the boundaries with nacelle inflow bc ---*/
      for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {

        /*--- Add the FanFace_MassFlow, FanFace_Mach, FanFace_Pressure and FanFace_Area to the particular boundary ---*/
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_NacelleInflow(iMarker_NacelleInflow)) {
          FanFace_MassFlow_Local[iMarker_NacelleInflow] += FanFace_MassFlow[iMarker];
          FanFace_Mach_Local[iMarker_NacelleInflow] += FanFace_Mach[iMarker];
          FanFace_Pressure_Local[iMarker_NacelleInflow] += FanFace_Pressure[iMarker];
          FanFace_Area_Local[iMarker_NacelleInflow] += FanFace_Area[iMarker];
        }

      }

    }

    if (config->GetMarker_All_KindBC(iMarker) == NACELLE_EXHAUST) {

      /*--- Loop over all the boundaries with nacelle inflow bc ---*/
      for (iMarker_NacelleExhaust= 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {

        /*--- Add the Exhaust_MassFlow, and Exhaust_Area to the particular boundary ---*/
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_NacelleExhaust(iMarker_NacelleExhaust)) {
          Exhaust_MassFlow_Local[iMarker_NacelleExhaust] += Exhaust_MassFlow[iMarker];
          Exhaust_Area_Local[iMarker_NacelleExhaust] += Exhaust_Area[iMarker];
        }

      }

    }

  }

#ifdef HAVE_MPI

  MPI_Allreduce(FanFace_MassFlow_Local, FanFace_MassFlow_Total, nMarker_NacelleInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(FanFace_Mach_Local, FanFace_Mach_Total, nMarker_NacelleInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(FanFace_Pressure_Local, FanFace_Pressure_Total, nMarker_NacelleInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(FanFace_Area_Local, FanFace_Area_Total, nMarker_NacelleInflow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(Exhaust_MassFlow_Local, Exhaust_MassFlow_Total, nMarker_NacelleExhaust, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(Exhaust_Area_Local, Exhaust_Area_Total, nMarker_NacelleExhaust, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

  for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
    FanFace_MassFlow_Total[iMarker_NacelleInflow]   = FanFace_MassFlow_Local[iMarker_NacelleInflow];
    FanFace_Mach_Total[iMarker_NacelleInflow]       = FanFace_Mach_Local[iMarker_NacelleInflow];
    FanFace_Pressure_Total[iMarker_NacelleInflow]   = FanFace_Pressure_Local[iMarker_NacelleInflow];
    FanFace_Area_Total[iMarker_NacelleInflow]       = FanFace_Area_Local[iMarker_NacelleInflow];
  }

  for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
    Exhaust_MassFlow_Total[iMarker_NacelleExhaust]  = Exhaust_MassFlow_Local[iMarker_NacelleExhaust];
    Exhaust_Area_Total[iMarker_NacelleExhaust]      = Exhaust_Area_Local[iMarker_NacelleExhaust];
  }

#endif

  /*--- Compute the value of FanFace_Area_Total, and FanFace_Pressure_Total, and
   set the value in the config structure for future use ---*/
  for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
    if (FanFace_Area_Total[iMarker_NacelleInflow] != 0.0) FanFace_Mach_Total[iMarker_NacelleInflow] /= FanFace_Area_Total[iMarker_NacelleInflow];
    else FanFace_Mach_Total[iMarker_NacelleInflow] = 0.0;
    if (FanFace_Area_Total[iMarker_NacelleInflow] != 0.0) FanFace_Pressure_Total[iMarker_NacelleInflow] /= FanFace_Area_Total[iMarker_NacelleInflow];
    else FanFace_Pressure_Total[iMarker_NacelleInflow] = 0.0;

    if (iMesh == MESH_0) {
      config->SetFanFace_Mach(iMarker_NacelleInflow, FanFace_Mach_Total[iMarker_NacelleInflow]);
      config->SetFanFace_Pressure(iMarker_NacelleInflow, FanFace_Pressure_Total[iMarker_NacelleInflow]);
    }

  }

  bool write_heads = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*20)) == 0));

  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output) {

    cout.precision(4);
    cout.setf(ios::fixed,ios::floatfield);

    cout << endl << "---------------------------- Engine properties --------------------------" << endl;
    for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
      cout << "Nacelle inflow ("<< config->GetMarker_NacelleInflow(iMarker_NacelleInflow)
      << "): MassFlow (kg/s): " << FanFace_MassFlow_Total[iMarker_NacelleInflow] * config->GetDensity_Ref() * config->GetVelocity_Ref()
      << ", Mach: " << FanFace_Mach_Total[iMarker_NacelleInflow]
      << ", Area: " << FanFace_Area_Total[iMarker_NacelleInflow] <<"."<< endl;
    }

    for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
      cout << "Nacelle exhaust ("<< config->GetMarker_NacelleExhaust(iMarker_NacelleExhaust)
      << "): MassFlow (kg/s): " << Exhaust_MassFlow_Total[iMarker_NacelleExhaust] * config->GetDensity_Ref() * config->GetVelocity_Ref()
      << ", Area: " << Exhaust_Area_Total[iMarker_NacelleExhaust] <<"."<< endl;
    }
    cout << "-------------------------------------------------------------------------" << endl;

  }

  /*--- Check the flow orientation in the nacelle inflow ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) == NACELLE_INFLOW) {

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
        if (alpha < 0.0) {

          /*--- Copy the old solution ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->SetSolution(iVar, node[iPoint]->GetSolution_Old(iVar));

        }

      }
    }
  }

  delete [] FanFace_MassFlow_Local;
  delete [] FanFace_Mach_Local;
  delete [] FanFace_Pressure_Local;
  delete [] FanFace_Area_Local;

  delete [] FanFace_MassFlow_Total;
  delete [] FanFace_Mach_Total;
  delete [] FanFace_Pressure_Total;
  delete [] FanFace_Area_Total;

  delete [] Exhaust_MassFlow_Local;
  delete [] Exhaust_Area_Local;

  delete [] Exhaust_MassFlow_Total;
  delete [] Exhaust_Area_Total;

}

void CEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output) {

  unsigned short iDim, iCounter;
  bool Update_AoA = false;
  double Target_CL, AoA_inc, AoA, Eps_Factor = 1e2;
  double DampingFactor = config->GetDamp_Fixed_CL();
  double Beta = config->GetAoS()*PI_NUMBER/180.0;
  double Vel_Infty[3], Vel_Infty_Mag;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Only the fine mesh level should check the convergence criteria ---*/

  if (iMesh == MESH_0) {

    /*--- Initialize the update flag to false in config ---*/

    config->SetUpdate_AoA(false);

    /*--- Initialize the local cauchy criteria on the first iteration ---*/

    if (config->GetExtIter() == 0) {
      Cauchy_Value = 0.0;
      Cauchy_Counter = 0;
      for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
        Cauchy_Serie[iCounter] = 0.0;
      AoA_old = config->GetAoA()*PI_NUMBER/180.0;
    }

    /*--- Check on the level of convergence in the lift coefficient. ---*/

    Old_Func = New_Func;
    New_Func = Total_CLift;
    Cauchy_Func = fabs(New_Func - Old_Func);
    Cauchy_Serie[Cauchy_Counter] = Cauchy_Func;
    Cauchy_Counter++;
    if (Cauchy_Counter == config->GetCauchy_Elems()) Cauchy_Counter = 0;

    Cauchy_Value = 1;
    if (config->GetExtIter() >= config->GetCauchy_Elems()) {
      Cauchy_Value = 0;
      for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
        Cauchy_Value += Cauchy_Serie[iCounter];
    }

    /*--- Check whether we are within two digits of the requested convergence
     epsilon for the cauchy criteria. ---*/

    if (Cauchy_Value >= config->GetCauchy_Eps()*Eps_Factor) Update_AoA = false;
    else Update_AoA = true;

    /*--- Do not apply any convergence criteria if the number
     of iterations is less than a particular value ---*/
    if (config->GetExtIter() < config->GetStartConv_Iter()) {
      Update_AoA = false;

    }
    /*--- Store the update boolean for use on other mesh levels in the MG ---*/

    config->SetUpdate_AoA(Update_AoA);

  } else
    Update_AoA = config->GetUpdate_AoA();

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

  bool write_heads = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*20)) == 0));
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && write_heads && Output) {
    cout.precision(7);
    cout.setf(ios::fixed,ios::floatfield);
    cout << endl << "----------------------------- Fixed CL Mode -----------------------------" << endl;
    cout << " Target CL: " << config->GetTarget_CL();
    cout << ", Current CL: " << Total_CLift;
    cout << ", Current AoA: " << config->GetAoA() << " deg" << endl;
    cout << "-------------------------------------------------------------------------" << endl;
  }

}

void CEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jVar, kVar, jDim;
  unsigned long iPoint, iVertex;
  double Pressure = 0.0, *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3], *NormalArea,
  ProjGridVel = 0.0, turb_ke;

  double Density_b, StaticEnergy_b, Enthalpy_b, *Velocity_b, Kappa_b, Chi_b, ProjVelocity_b, Energy_b, VelMagnitude2_b, Pressure_b;
  double Density_i, *Velocity_i, ProjVelocity_i, Energy_i, VelMagnitude2_i;
  double **Jacobian_b, **DubDu;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  Normal = new double[nDim];
  NormalArea = new double[nDim];
  Velocity_b = new double[nDim];
  Velocity_i = new double[nDim];

  Jacobian_b = new double*[nVar];
  DubDu = new double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_b[iVar] = new double[nVar];
    DubDu[iVar] = new double[nVar];
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
        for(iDim = 0; iDim < nDim; iDim++) {
          Velocity_i[iDim] = node[iPoint]->GetVelocity(iDim);
          ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
          VelMagnitude2_i += Velocity_i[iDim]*Velocity_i[iDim];
        }

        Density_i = node[iPoint]->GetDensity();
        Energy_i = node[iPoint]->GetEnergy();

        /*--- Compute the boundary state b ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          Velocity_b[iDim] = Velocity_i[iDim] - ProjVelocity_i * UnitNormal[iDim];

        if (grid_movement) {
          GridVel = geometry->node[iPoint]->GetGridVel();
          ProjGridVel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel += GridVel[iDim]*UnitNormal[iDim];

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity_b[iDim] += ProjGridVel * UnitNormal[iDim];
        }

        VelMagnitude2_b = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          VelMagnitude2_b += Velocity_b[iDim] * Velocity_b[iDim];

        /*--- Add the kinetic energy correction ---*/
        turb_ke = 0.0;
        if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);

        Density_b = Density_i;
        StaticEnergy_b = Energy_i - 0.5 * VelMagnitude2_i - turb_ke;
        Energy_b = StaticEnergy_b + 0.5 * VelMagnitude2_b;

        FluidModel->SetTDState_rhoe(Density_b, StaticEnergy_b);

        Kappa_b = FluidModel->GetdPde_rho() / Density_b;
        Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;

        Pressure_b = FluidModel->GetPressure();
        Enthalpy_b = Energy_b + Pressure_b/Density_b;

	if (tkeNeeded) {
          Pressure += (2.0/3.0)*node[iPoint]->GetDensity()*turb_ke;
        }

        /*--- Compute the residual ---*/

        numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, NormalArea, Residual);

      }

      /*--- Incompressible solver ---*/

      if (incompressible || freesurface) {

        Pressure = node[iPoint]->GetPressureInc();

        /*--- Add the kinetic energy correction ---*/

        turb_ke = 0.0;
        if (tkeNeeded) {
          turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
          Pressure += (2.0/3.0)*node[iPoint]->GetDensity()*turb_ke;
        }

        /*--- Compute the residual ---*/

        Residual[0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual[iDim+1] = Pressure*UnitNormal[iDim]*Area;

        if (compressible || freesurface) {
          Residual[nVar-1] = 0.0;
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

        if (compressible)  {

          /*--- Compute DubDu ---*/

          for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++)
              DubDu[iVar][jVar]= 0.0;
            DubDu[iVar][iVar]= 1.0;
          }

          for(iDim = 0; iDim<nDim; iDim++)
            for(jDim = 0; jDim<nDim; jDim++)
              DubDu[iDim+1][jDim+1] -= UnitNormal[iDim]*UnitNormal[jDim];

          DubDu[nVar-1][0] += 0.5*ProjVelocity_i*ProjVelocity_i;

          for(iDim = 0; iDim<nDim; iDim++) {
            DubDu[nVar-1][iDim+1] -= ProjVelocity_i*UnitNormal[iDim];
          }

          /*--- Compute flux Jacobian in state b ---*/

          numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, NormalArea, 1, Jacobian_b);

// Check for grid movement, should be already considered since Jacobian b is computed from u_b
//          if (grid_movement) { 
//            Jacobian_b[nVar-1][0] += 0.5*ProjGridVel*ProjGridVel;
//            for (iDim = 0; iDim < nDim; iDim++)
//              Jacobian_b[nVar-1][iDim+1] -= ProjGridVel * UnitNormal[iDim];
//          }


          /*--- Compute numerical flux Jacobian at node i ---*/

          for(iVar = 0; iVar < nVar; iVar++)
            for(jVar = 0; jVar < nVar; jVar++)
              for(kVar = 0; kVar < nVar; kVar++)
                Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];

          /*--- Add the Jacobian to the sparse matrix ---*/

          Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);

        }

        /*--- Incompressible solver ---*/

        if (incompressible || freesurface)  {
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

  double *GridVel;
  double Area, UnitNormal[3];
  double Density, Pressure, Velocity[3], Energy;
  double Density_Bound, Pressure_Bound, Vel_Bound[3];
  double Density_Infty, Pressure_Infty, Vel_Infty[3];
  double SoundSpeed, Entropy, Velocity2, Vn;
  double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound;
  double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  double RiemannPlus, RiemannMinus;
  double *V_infty, *V_domain;

  double Gas_Constant     = config->GetGas_ConstantND();

  bool implicit         = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool grid_movement    = config->GetGrid_Movement();
  bool compressible     = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible   = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface      = (config->GetKind_Regime() == FREESURFACE);
  bool viscous          = config->GetViscous();
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  double *Normal = new double[nDim];

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
        Entropy_Bound    = pow(Density_Bound,Gamma)/Pressure_Bound;

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
        Entropy_Infty    = pow(Density_Infty,Gamma)/Pressure_Infty;

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
  double P_Total, T_Total, P_static, T_static, Rho_static, *Mach, *Flow_Dir, Area, UnitNormal[3];

  double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
  double *Velocity_e, Velocity2_e, VelMag_e, Enthalpy_e, Entropy_e, Energy_e = 0.0, StaticEnthalpy_e, StaticEnergy_e, Density_e = 0.0, Pressure_e;
  double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  double ProjGridVel, ProjVelocity_i, ProjVelocity_b;
  double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_e, *u_i, *u_b;

  double *V_boundary, *V_domain;

  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  double *Normal;

/// Doesn't work for ALE


  Normal = new double[nDim];

  Velocity_i = new double[nDim];
  Velocity_b = new double[nDim];
  Velocity_e = new double[nDim];


  V_boundary = new double[nPrimVar];
  Lambda_i = new double[nVar];
  u_i = new double[nVar];
  u_e = new double[nVar];
  u_b = new double[nVar];
  dw = new double[nVar];

  P_Tensor = new double*[nVar];
  invP_Tensor = new double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new double[nVar];
    invP_Tensor[iVar] = new double[nVar];
  }

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
        V_domain = node[iPoint]->GetPrimitive();

        /* --- Compute the internal state u_i --- */
        Velocity2_i = 0;
        for(iDim=0; iDim < nDim; iDim++)
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

 /*--- Build the external state u_e from boundary data and internal node ---*///
        switch(config->GetKind_Data_Riemann(Marker_Tag))
        {

            case TOTAL_CONDITIONS_PT: case TOTAL_SUPERSONIC_INFLOW:
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

                if(config->GetKind_Data_Riemann(Marker_Tag) == TOTAL_CONDITIONS_PT){
                	Velocity2_e = Velocity2_i;

                	for (iDim = 0; iDim < nDim; iDim++) {
                		Velocity_e[iDim] = sqrt(Velocity2_e)*Flow_Dir[iDim];
                	}
                }else{
                	Velocity2_e = 0.0;
                	for (iDim = 0; iDim < nDim; iDim++) {
						Velocity_e[iDim] = Flow_Dir[iDim]/config->GetVelocity_Ref();
						Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                		}
                }
                StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;

                FluidModel->SetTDState_hs(StaticEnthalpy_e, Entropy_e);

                Density_e = FluidModel->GetDensity();
                StaticEnergy_e = FluidModel->GetStaticEnergy();

                Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;

                if (tkeNeeded) Energy_e += GetTke_Inf();

                Pressure_e = FluidModel->GetPressure();
                Enthalpy_e = Energy_e + Pressure_e/Density_e;

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

                FluidModel->SetTDState_rhoe(Density_e, Energy_e);

                Pressure_e = FluidModel->GetPressure();
                Enthalpy_e = Energy_e + Pressure_e/Density_e;

                break;


            case STATIC_PRESSURE:

                Pressure_e = config->GetRiemann_Var1(Marker_Tag);
                Pressure_e /= config->GetPressure_Ref();

                Density_e = Density_i;

                FluidModel->SetTDState_Prho(Pressure_e, Density_e);

                Velocity2_e = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                  Velocity_e[iDim] = Velocity_i[iDim];
                  Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                }

                Energy_e = FluidModel->GetStaticEnergy() + 0.5*Velocity2_e;
                Enthalpy_e = Energy_e + Pressure_e/Density_e;

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

				Pressure_e = FluidModel->GetPressure();
				Enthalpy_e = Energy_e + Pressure_e/Density_e;

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

				Pressure_e = FluidModel->GetPressure();
				Enthalpy_e = Energy_e + Pressure_e/Density_e;

				break;


            default:
                cout << "Warning! Invalid Riemann input!" << endl;
		exit(1);

                break;

        }


         /*--- Compute P (matrix of right eigenvectors) ---*/
        conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);

        /*--- Compute inverse P (matrix of left eigenvectors)---*/
        conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);

        /*--- Flow eigenvalues ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Lambda_i[iDim] = ProjVelocity_i;
        Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
        Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;

        u_e[0] = Density_e;
        for (iDim = 0; iDim < nDim; iDim++)
          u_e[iDim+1] = Velocity_e[iDim]*Density_e;
        u_e[nVar-1] = Energy_e*Density_e;

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
                if(Lambda_i[jVar] < 0)
                {
                    u_b[iVar] += P_Tensor[iVar][jVar]*dw[jVar];

                }
            }
        }

//       cout << u_e[2]<< " "<< u_i[2] <<" "<< u_b[2] << endl;
//        /*--- Primitive variables, using the derived quantities ---*/
//
//        V_boundary[nDim+2] = u_b[0];
//        double v2 = 0;
//        for (iDim = 0; iDim < nDim; iDim++){
//        	V_boundary[iDim+1] = u_b[iDim]/u_b[0];
//        	v2 += V_boundary[iDim+1]*V_boundary[iDim+1];
//        }
//        double e_b = u_b[nVar-1]/u_b[0] - 0.5*v2;
//        FluidModel->SetTDState_rhoe(u_b[0], e_b);
//        V_boundary[0] = FluidModel->GetTemperature() ;
//        V_boundary[nDim+1] = FluidModel->GetPressure();
//        V_boundary[nDim+3] = u_b[nVar-1]/u_b[0] + V_boundary[nDim+1]/u_b[0];
//        V_boundary[nDim+4] = FluidModel->GetSoundSpeed();
////        cout <<  V_boundary[0] << " "<<  V_boundary[1] << " "<<  V_boundary[2] << " "<<  V_boundary[3] << " "<<  V_boundary[4] << " "<< V_boundary[5] << " "<< Marker_Tag << " "<< endl;
//
//
//        /*--- Set various quantities in the solver class ---*/
//        conv_numerics->SetPrimitive(V_domain, V_boundary);
//
//        if (grid_movement)
//          conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
//
//        /*--- Compute the residual using an upwind scheme ---*/
//        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//        /*--- Update residual value ---*/
//        LinSysRes.AddBlock(iPoint, Residual);
//
//        /*--- Jacobian contribution for implicit integration ---*/
//        if (implicit)
//          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

//        for (iVar = 0; iVar < nVar; iVar++){
//        	cout << u_i[iVar] <<" "<< u_b[iVar]<<" "<< u_e[iVar] << endl;
//        }

        /*--- Compute the thermodynamic state in u_b ---*/
        Density_b = u_b[0];
//        cout << u_b[0] << endl;

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
        Enthalpy_b = Energy_b + Pressure_b/Density_b;

        Kappa_b = FluidModel->GetdPde_rho() / Density_b;
        Chi_b = FluidModel->GetdPdrho_e() - Kappa_b * StaticEnergy_b;

        /*--- Compute the residuals ---*/
        conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);

        if (implicit) {

              Jacobian_b = new double*[nVar];
              DubDu = new double*[nVar];
              for (iVar = 0; iVar < nVar; iVar++)
              {
                  Jacobian_b[iVar] = new double[nVar];
                  DubDu[iVar] = new double[nVar];
              }

               /*--- Initialize DubDu to unit matrix---*/

              for (iVar = 0; iVar < nVar; iVar++)
              {
                  for (jVar = 0; jVar < nVar; jVar++)
                    DubDu[iVar][jVar]= 0;

                  DubDu[iVar][iVar]= 1;
              }

              /*--- Compute DubDu -= RNL---*/
              for(iVar=0; iVar<nVar; iVar++)
              {
                  for(jVar=0; jVar<nVar; jVar++)
                  {
                      for(kVar=0; kVar<nVar; kVar++)
                      {
                          if(Lambda_i[kVar]<0)
                            DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
                      }
                  }
              }


               /*--- Compute flux Jacobian in state b ---*/
//             cout << Enthalpy_b << " " << Chi_b<< " " << Kappa_b <<" "<< endl;
              conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);
/// check ALE




              if (grid_movement)
              {
                Jacobian_b[nVar-1][0] += 0.5*ProjGridVel*ProjGridVel;

                for (iDim = 0; iDim < nDim; iDim++)
                  Jacobian_b[nVar-1][iDim+1] += 0.5 * ProjVelocity_b * UnitNormal[iDim];
              }

               /*--- Compute numerical flux Jacobian at node i ---*/

              for(iVar=0; iVar<nVar; iVar++)
              {
                  for(jVar=0; jVar<nVar; jVar++)
                  {
                      for(kVar=0; kVar<nVar; kVar++)
                      {
                          Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
                      }
                  }
//              cout << Jacobian_i[iVar][0] << " " << Jacobian_i[iVar][1] << " " << Jacobian_i[iVar][2]<< " " << Jacobian_i[iVar][3]<< " " << endl;
              }


//              Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);


              for (iVar = 0; iVar < nVar; iVar++)
              {
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

          V_boundary[nDim+5] = node[iPoint]->GetLaminarViscosity();
          V_boundary[nDim+6] = node[iPoint]->GetEddyViscosity();

        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_boundary);
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
// getchar();
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] Velocity_e;
  delete [] Velocity_b;
  delete [] Velocity_i;

  delete [] V_boundary;
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

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3];
  double *V_inlet, *V_domain;

  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
  double Gas_Constant       = config->GetGas_ConstantND();
  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  double *Normal = new double[nDim];

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

        /*--- The velocity is computed from the infinity values ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim+1] = GetVelocity_Inf(iDim);

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
  double LevelSet, Density_Outlet = 0.0, Pressure, P_Exit, Velocity[3],
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
  Area, UnitNormal[3], Height, yCoordRef, yCoord;
  double *V_outlet, *V_domain;

  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  double Gas_Constant     = config->GetGas_ConstantND();
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement      = config->GetGrid_Movement();
  double FreeSurface_Zero = config->GetFreeSurface_Zero();
  double epsilon          = config->GetFreeSurface_Thickness();
  double RatioDensity     = config->GetRatioDensity();
  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity = (config->GetGravityForce());
  double PressFreeSurface = GetPressure_Inf();
  double Froude           = config->GetFroude();
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  double *Normal = new double[nDim];

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
        Energy     = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
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

          Entropy = Pressure*pow(1.0/Density,Gamma);
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
  double Area, UnitNormal[3];
  double *V_inlet, *V_domain;

  double Density, Pressure, Temperature, Energy, *Velocity, Velocity2;
  double Gas_Constant = config->GetGas_ConstantND();

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool viscous              = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  double *Normal = new double[nDim];

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
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

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

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

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

void CEulerSolver::BC_Nacelle_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double Pressure, P_Fan, Velocity[3], Velocity2, Entropy, Target_FanFace_Mach = 0.0, Density, Energy,
  Riemann, Area, UnitNormal[3], Vn, SoundSpeed, Vn_Exit, P_Fan_inc, P_Fan_old, M_Fan_old;
  double *V_inflow, *V_domain;

  double DampingFactor = config->GetDamp_Nacelle_Inflow();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous              = config->GetViscous();
  double Gas_Constant = config->GetGas_ConstantND();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  double *Normal = new double[nDim];

  /*--- Retrieve the specified target fan face mach in the nacelle. ---*/
  Target_FanFace_Mach = config->GetFanFace_Mach_Target(Marker_Tag);

  /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
  P_Fan_old = config->GetFanFace_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
  M_Fan_old = config->GetFanFace_Mach(Marker_Tag);

  /*--- Compute the Pressure increment ---*/
  P_Fan_inc = ((M_Fan_old/Target_FanFace_Mach) - 1.0) * config->GetPressure_FreeStreamND();

  /*--- Estimate the new fan face pressure ---*/
  P_Fan = (1.0 - DampingFactor)*P_Fan_old + DampingFactor * (P_Fan_old + P_Fan_inc);

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
      Energy     = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Entropy = Pressure*pow(1.0/Density,Gamma);
      Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

      /*--- Compute the new fictious state at the outlet ---*/
      Density    = pow(P_Fan/Entropy,1.0/Gamma);
      Pressure   = P_Fan;
      SoundSpeed = sqrt(Gamma*P_Fan/Density);
      Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
      Velocity2  = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }

      Energy = P_Fan/(Density*Gamma_Minus_One) + 0.5*Velocity2;
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

void CEulerSolver::BC_Nacelle_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double P_Exhaust, T_Exhaust, Velocity[3], Velocity2, H_Exhaust, Temperature, Riemann, Area, UnitNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Exhaust2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
  double *V_exhaust, *V_domain;

  double Gas_Constant = config->GetGas_ConstantND();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous = config->GetViscous();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));

  double *Normal = new double[nDim];

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

      /*--- Retrieve the specified total conditions for this inlet. ---*/
      P_Exhaust  = config->GetNozzle_Ptotal(Marker_Tag);
      T_Exhaust  = config->GetNozzle_Ttotal(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exhaust /= config->GetPressure_Ref();
      T_Exhaust /= config->GetTemperature_Ref();

      /*--- Store primitives and set some variables for clarity. ---*/
      Density = V_domain[nDim+2];
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Energy = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure = V_domain[nDim+1];
      H_Exhaust     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Exhaust;
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
      dd = bb*bb - 4.0*aa*cc;
      dd = sqrt(max(0.0,dd));
      Vel_Mag   = (-bb + dd)/(2.0*aa);
      Vel_Mag   = max(0.0,Vel_Mag);
      Velocity2 = Vel_Mag*Vel_Mag;

      /*--- Compute speed of sound from total speed of sound eqn. ---*/
      SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;

      /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
      Mach2 = Velocity2/SoundSpeed2;
      Mach2 = min(1.0,Mach2);
      Velocity2   = Mach2*SoundSpeed2;
      Vel_Mag     = sqrt(Velocity2);
      SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;

      /*--- Compute new velocity vector at the inlet ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

      /*--- Static temperature from the speed of sound relation ---*/
      Temperature = SoundSpeed2/(Gamma*Gas_Constant);

      /*--- Static pressure using isentropic relation at a point ---*/
      Pressure = P_Exhaust*pow((Temperature/T_Exhaust),Gamma/Gamma_Minus_One);

      /*--- Density at the inlet from the gas law ---*/
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

void CEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {

  /*--- Call the Euler residual --- */

  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);

}

void CEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint, jPoint;
  unsigned short iDim, iVar;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  double *Normal = new double[nDim];
  double *PrimVar_i = new double[nPrimVar];
  double *PrimVar_j = new double[nPrimVar];

#ifndef HAVE_MPI

  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Find the associate pair to the original node ---*/

      jPoint = geometry->vertex[val_marker][iVertex]->GetDonorPoint();

      if (iPoint != jPoint) {

        /*--- Store the solution for both points ---*/

        for (iVar = 0; iVar < nPrimVar; iVar++) {
          PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
          PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
        }

        /*--- Set primitive variables ---*/

        numerics->SetPrimitive(PrimVar_i, PrimVar_j);

        /*--- Set the normal vector ---*/

        geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
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

#else
  int rank, jProcessor;
  MPI_Status status;
  MPI_Request send_req, recv_req;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bool compute;
  double *Buffer_Send_V = new double [nPrimVar];
  double *Buffer_Receive_V = new double [nPrimVar];

  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/

  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Find the associate pair to the original node ---*/

      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];

      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
      else compute = true;

      /*--- We only send the information that belong to other boundary, -1 processor
       means that the boundary condition is not applied ---*/

      if (compute) {

        if (jProcessor != rank) {

          /*--- Copy the primitive variable ---*/

          for (iVar = 0; iVar < nPrimVar; iVar++)
            Buffer_Send_V[iVar] = node[iPoint]->GetPrimitive(iVar);
          MPI_Isend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD,&send_req);
          MPI_Wait(&send_req,&status);
        }

      }

    }
  }

  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Find the associate pair to the original node ---*/

      jPoint = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[0];
      jProcessor = geometry->vertex[val_marker][iVertex]->GetPeriodicPointDomain()[1];

      if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
      else compute = true;

      if (compute) {

        /*--- We only receive the information that belong to other boundary ---*/

        if (jProcessor != rank) {
          MPI_Irecv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &recv_req);
          MPI_Wait(&recv_req,&status);
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

        geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
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

  delete[] Buffer_Send_V;
  delete[] Buffer_Receive_V;

#endif

  /*--- Free locally allocated memory ---*/

  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;

}

void CEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {

  /*--- Call the Interface_Boundary residual --- */

  BC_Interface_Boundary(geometry, solver_container, numerics, config, val_marker);


}

void CEulerSolver::BC_ActDisk_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config) {

  unsigned long iVertex, iPoint, jPoint, Pin = 0, Pout = 0, jProcessor;
  unsigned short iDim, iVar, iMarker;
  int iProcessor;
  double *Coord, radius, R, V_tip, DeltaP_avg, DeltaP_tip;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  unsigned short nMarker_ActDisk_Inlet = config->GetnMarker_ActDisk_Inlet();

  if (nMarker_ActDisk_Inlet != 0) {

    double *Normal = new double[nDim];
    double *PrimVar_out = new double[nPrimVar];
    double *PrimVar_in = new double[nPrimVar];
    double *MeanPrimVar = new double[nPrimVar];
    double *PrimVar_out_ghost = new double[nPrimVar];
    double *PrimVar_in_ghost = new double[nPrimVar];
    double *ActDisk_Jump = new double[nPrimVar];
    double *Buffer_Send_V = new double [nPrimVar];
    double *Buffer_Receive_V = new double [nPrimVar];

#ifndef HAVE_MPI
    iProcessor = MASTER_NODE;
#else
    MPI_Status status;
    MPI_Request send_req, recv_req;
    MPI_Comm_rank(MPI_COMM_WORLD, &iProcessor);
#endif

    /*--- Identify the points that should be sended in a MPI implementation.
     Do the send process, by the moment we are sending each node individually, this must be changed---*/

#ifdef HAVE_MPI

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
          (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {

        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->node[iPoint]->GetDomain()) {

            /*--- Find the associate pair to the original node ---*/

            jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
            jProcessor = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();

            /*--- We only send the information that belong to other boundary, using jPoint as the ID for the message  ---*/

            if (jProcessor != iProcessor) {

              /*--- Copy the primitive variables ---*/

              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Send_V[iVar] = node[iPoint]->GetPrimitive(iVar);

              MPI_Isend(Buffer_Send_V, nPrimVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD,&send_req);
              MPI_Wait(&send_req,&status);

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
        double *origin = config->GetActDisk_Origin(config->GetMarker_All_TagBound(iMarker));     // Center of the rotor
        double R_root = config->GetActDisk_RootRadius(config->GetMarker_All_TagBound(iMarker));
        double R_tip = config->GetActDisk_TipRadius(config->GetMarker_All_TagBound(iMarker));
        double C_T = config->GetActDisk_CT(config->GetMarker_All_TagBound(iMarker));             // Rotor thrust coefficient
        double Omega = config->GetActDisk_Omega(config->GetMarker_All_TagBound(iMarker));        // Revolution per minute

        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->node[iPoint]->GetDomain()) {

            /*--- Find the associate pair to the original node ---*/

            jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
            jProcessor = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();

            /*--- Receive the information, using jPoint as the ID for the message ---*/

            if (jProcessor != iProcessor) {
#ifdef HAVE_MPI
              MPI_Irecv(Buffer_Receive_V, nPrimVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &recv_req);
              MPI_Wait(&recv_req,&status);
#endif
            }
            else {

              /*--- The point is in the same processor... no MPI required ---*/

              for (iVar = 0; iVar < nPrimVar; iVar++)
                Buffer_Receive_V[iVar] = node[jPoint]->GetPrimitive(iVar);

            }

            /*--- Identify the inner and the outer point (based on the normal direction) ---*/

            if (boundary == ACTDISK_INLET) {

              Pin = iPoint; Pout = jPoint;

              for (iVar = 0; iVar < nPrimVar; iVar++) {
                PrimVar_out[iVar] = Buffer_Receive_V[iVar];
                PrimVar_in[iVar] = node[Pin]->GetPrimitive(iVar);
                MeanPrimVar[iVar] = 0.5*(PrimVar_out[iVar] + PrimVar_in[iVar]);
              }

            }
            if (boundary == ACTDISK_OUTLET) {

              Pout = iPoint; Pin = jPoint;

              for (iVar = 0; iVar < nPrimVar; iVar++) {
                PrimVar_out[iVar] = node[Pout]->GetPrimitive(iVar);
                PrimVar_in[iVar] = Buffer_Receive_V[iVar];
                MeanPrimVar[iVar] = 0.5*(PrimVar_out[iVar] + PrimVar_in[iVar]);
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
              radius += (Coord[iDim]-origin[iDim])*(Coord[iDim]-origin[iDim]);
            radius = sqrt(radius);

            /*--- Compute the pressure increment and the jump ---*/

            R = R_tip;
            V_tip = R_tip*Omega*(PI_NUMBER/30.0);
            DeltaP_avg = R*R*V_tip*V_tip*C_T/(R*R-R_root*R_root);
            DeltaP_tip = 3.0*DeltaP_avg*R*(R*R-R_root*R_root)/(2.0*(R*R*R-R_root*R_root*R_root));
            ActDisk_Jump[nDim+1] = -DeltaP_tip*radius/R;

            /*--- Inner point ---*/

            if (iPoint == Pin) {
              for (iVar = 0; iVar < nPrimVar; iVar++)
                PrimVar_in_ghost[iVar] = 2.0*MeanPrimVar[iVar] - PrimVar_in[iVar] - ActDisk_Jump[iVar];
              numerics->SetPrimitive(PrimVar_in, PrimVar_in_ghost);
            }

            /*--- Outer point ---*/

            if (iPoint == Pout) {
              for (iVar = 0; iVar < nPrimVar; iVar++)
                PrimVar_out_ghost[iVar] = 2.0*MeanPrimVar[iVar] - PrimVar_out[iVar] + ActDisk_Jump[iVar];
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

    /*--- Free locally allocated memory ---*/

    delete [] Normal;
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

  double *U_time_nM1, *U_time_n, *U_time_nP1;
  double Volume_nM1, Volume_nP1, TimeStep;
  double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;

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
      for(iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      if ((incompressible || freesurface) && (FlowEq || AdjEq)) Residual[0] = 0.0;
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Compute the GCL component of the source term for node j ---*/

      U_time_n = node[jPoint]->GetSolution_time_n();
      for(iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      if ((incompressible || freesurface) && (FlowEq || AdjEq)) Residual[0] = 0.0;
      LinSysRes.SubtractBlock(jPoint, Residual);

    }

    /*---	Loop over the boundary edges ---*/

    for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

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
        for(iVar = 0; iVar < nVar; iVar++)
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

      for(iVar = 0; iVar < nVar; iVar++) {
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
  //  unsigned short iMarker, iDim;
  //  unsigned long iVertex, iPoint;
  //  double *Coord, VarCoord[3];
  //
  //#ifndef HAVE_MPI
  //  unsigned long iPoint_Donor;
  //  double *CoordDonor, *DisplacementDonor;
  //
  //  for (iMarker = 0; iMarker < flow_config->GetnMarker_All(); iMarker++) {
  //    if (flow_config->GetMarker_All_Moving(iMarker) == YES) {
  //      for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
  //        iPoint = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
  //        iPoint_Donor = flow_geometry[MESH_0]->vertex[iMarker][iVertex]->GetDonorPoint();
  //        Coord = flow_geometry[MESH_0]->node[iPoint]->GetCoord();
  //        CoordDonor = fea_geometry[MESH_0]->node[iPoint_Donor]->GetCoord();
  //        DisplacementDonor = fea_solution[MESH_0][FEA_SOL]->node[iPoint_Donor]->GetSolution();
  //
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          VarCoord[iDim] = (CoordDonor[iDim]+DisplacementDonor[iDim])-Coord[iDim];
  //
  //        flow_geometry[MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
  //      }
  //    }
  //  }
  //  flow_grid_movement->SetVolume_Deformation(flow_geometry[MESH_0], flow_config, true);
  //
  //#else
  //
  //  int rank = MPI::COMM_WORLD.Get_rank(), jProcessor;
  //  double *Buffer_Send_Coord = new double [nDim];
  //  double *Buffer_Receive_Coord = new double [nDim];
  //  unsigned long jPoint;
  //
  //  /*--- Do the send process, by the moment we are sending each
  //   node individually, this must be changed ---*/
  //  for (iMarker = 0; iMarker < fea_config->GetnMarker_All(); iMarker++) {
  //    if (fea_config->GetMarker_All_KindBC(iMarker) == LOAD_BOUNDARY) {
  //      for(iVertex = 0; iVertex < fea_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
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
  //      for(iVertex = 0; iVertex < flow_geometry[MESH_0]->nVertex[iMarker]; iVertex++) {
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
  double Area_Children, Area_Parent, Coord[3], *Solution_Fine, dull_val;
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  string UnstExt, text_line;
  ifstream restart_file;

  string restart_filename = config->GetSolution_FlowFileName();

  /*--- Modify file name for an unsteady restart ---*/
  if (dual_time)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Open the restart file, and throw an error if this fails. ---*/
  restart_file.open(restart_filename.data(), ios::in);
  if (restart_file.fail()) {
    cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
    exit(1);
  }

  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  long *Global2Local = NULL;
  Global2Local = new long[geometry[MESH_0]->GetGlobal_nPointDomain()];
  /*--- First, set all indices to a negative value by default ---*/
  for(iPoint = 0; iPoint < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint++) {
    Global2Local[iPoint] = -1;
  }

  /*--- Now fill array with the transform values only for local points ---*/
  for(iPoint = 0; iPoint < geometry[MESH_0]->GetnPointDomain(); iPoint++) {
    Global2Local[geometry[MESH_0]->node[iPoint]->GetGlobalIndex()] = iPoint;
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
        if (turb_model == SA || turb_model == ML) {
          point_line >> dull_val;
        } else if (turb_model == SST) {
          point_line >> dull_val >> dull_val;
        }

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        double GridVel[3];
        if (nDim == 2) point_line >> GridVel[0] >> GridVel[1];
        if (nDim == 3) point_line >> GridVel[0] >> GridVel[1] >> GridVel[2];
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
  for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
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
    geometry[MESH_0]->SetCG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);

    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/

    for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config,geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config,geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine],config);
    }
  }

}

void CEulerSolver::SetFreeSurface_Distance(CGeometry *geometry, CConfig *config) {
  double *coord = NULL, dist2, *iCoord = NULL, *jCoord = NULL, LevelSet_i, LevelSet_j,
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
  Coord_LevelSet = new double* [nVertex_LevelSet];
  for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
    Coord_LevelSet[iVertex] = new double [nDim];

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
  double *Buffer_Send_Coord = NULL, *Buffer_Receive_Coord = NULL;

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

  MPI_Allreduce(&nLocalVertex_LevelSet, &nGlobalVertex_LevelSet, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex_LevelSet, &MaxLocalVertex_LevelSet, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  nBuffer = MaxLocalVertex_LevelSet*nDim;
  Buffer_Send_Coord = new double [nBuffer];
  Buffer_Receive_Coord = new double [nProcessor*nBuffer];

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

  MPI_Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);

  /*--- Identification of the 0 level set points and coordinates ---*/
  nVertex_LevelSet = 0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
    nVertex_LevelSet += Buffer_Receive_nVertex[iProcessor];

  /*--- Allocate vector of boundary coordinates ---*/
  Coord_LevelSet = new double* [nVertex_LevelSet];
  for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
    Coord_LevelSet[iVertex] = new double [nDim];

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
    double NumberSign = 1.0;
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
      LevelSetDiff_Squared = 0.0; LevelSetDiff = 0.0;

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
      if (config->GetUnsteady_Simulation()){
        if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.dat", int(iExtIter));
        if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.dat", int(iExtIter));
        if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.dat", int(iExtIter));
        if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.dat", int(iExtIter));
        if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.dat", int(iExtIter));
      }
      else {
        sprintf (buffer, ".dat");
      }

      strcat(cstr,buffer);

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

  /*--- Array initialization ---*/
  CDrag_Visc = NULL;
  CLift_Visc = NULL;
  CSideForce_Visc = NULL;
  CEff_Visc = NULL;
  CMx_Visc = NULL;
  CMy_Visc = NULL;
  CMz_Visc = NULL;
  CFx_Visc = NULL;
  CFy_Visc = NULL;
  CFz_Visc = NULL;
  Surface_CLift_Visc = NULL;
  Surface_CDrag_Visc = NULL;
  Surface_CMx_Visc = NULL;
  Surface_CMy_Visc = NULL;
  Surface_CMz_Visc = NULL;
  CMerit_Visc = NULL;
  CT_Visc = NULL;
  CQ_Visc = NULL;
  ForceViscous = NULL;
  MomentViscous = NULL;
  CSkinFriction = NULL;

}

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CEulerSolver() {

  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  double Density, Velocity2, Pressure, Temperature, dull_val;

  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface = (config->GetKind_Regime() == FREESURFACE);
  double Gas_Constant = config->GetGas_ConstantND();
  bool adjoint = config->GetAdjoint();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Array initialization ---*/
  CDrag_Visc = NULL;
  CLift_Visc = NULL;
  CSideForce_Visc = NULL;
  CEff_Visc = NULL;
  CMx_Visc = NULL;
  CMy_Visc = NULL;
  CMz_Visc = NULL;
  CFx_Visc = NULL;
  CFy_Visc = NULL;
  CFz_Visc = NULL;
  Surface_CLift_Visc = NULL;
  Surface_CDrag_Visc = NULL;
  Surface_CMx_Visc = NULL;
  Surface_CMy_Visc = NULL;
  Surface_CMz_Visc = NULL;
  CMerit_Visc = NULL;
  CT_Visc = NULL;
  CQ_Visc = NULL;
  Heat_Visc = NULL;
  MaxHeatFlux_Visc = NULL;
  ForceViscous = NULL;
  MomentViscous = NULL;
  CSkinFriction = NULL;
  Cauchy_Serie = NULL;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure
   Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta,lamMu,EddyMu),
   FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,lamMu,EddyMu, dist),
   Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c,lamMu,EddyMu) ---*/
  nDim = geometry->GetnDim();
  if (incompressible) { nVar = nDim+1; nPrimVar = nDim+5; nPrimVarGrad = nDim+3; }
  if (freesurface)    { nVar = nDim+2; nPrimVar = nDim+7; nPrimVarGrad = nDim+6; }
  if (compressible)   { nVar = nDim+2;
    nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
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
  Residual      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  Point_Max  = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;

  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the geometry ---*/
  Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;

  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  Primitive   = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the Secondary solution ---*/

  Secondary   = new double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar]   = 0.0;
  Secondary_i = new double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_i[iVar] = 0.0;
  Secondary_j = new double[nSecondaryVar]; for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary_j[iVar] = 0.0;

  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new double [nPoint];
    jPoint_UndLapl = new double [nPoint];
  }

  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  if (config->GetKind_Upwind_Flow() == TURKEL) {
    LowMach_Precontioner = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new double[nVar];
  }

  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Jacobians and vector structures for implicit computations ---*/
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    /*--- Point to point Jacobians ---*/
    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }
    /*--- Initialization of the structure of the whole Jacobian ---*/
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }

  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
  }

  /*--- Define some auxiliary vectors for computing flow variable gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new double [nDim];

    /*--- c vector := transpose(WA)*(Wb) ---*/
    cvector = new double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      cvector[iVar] = new double [nDim];
  }

  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  CharacPrimVar = new double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }

  /*--- Inviscid force definition and coefficient in all the markers ---*/
  CPressure = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Inviscid force definition and coefficient in all the markers ---*/
  CPressureTarget = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressureTarget[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Heat tranfer in all the markers ---*/
  HeatFlux = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatFlux[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      HeatFlux[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Heat tranfer in all the markers ---*/
  HeatFluxTarget = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatFluxTarget[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      HeatFluxTarget[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Y plus in all the markers ---*/
  YPlus = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    YPlus[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      YPlus[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Skin friction in all the markers ---*/
  CSkinFriction = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSkinFriction[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CSkinFriction[iMarker][iVertex] = 0.0;
    }
  }

  /*--- Non dimensional coefficients ---*/
  ForceInviscid  = new double[3];
  MomentInviscid = new double[3];
  CDrag_Inv      = new double[nMarker];
  CLift_Inv      = new double[nMarker];
  CSideForce_Inv = new double[nMarker];
  CMx_Inv        = new double[nMarker];
  CMy_Inv        = new double[nMarker];
  CMz_Inv        = new double[nMarker];
  CEff_Inv       = new double[nMarker];
  CFx_Inv        = new double[nMarker];
  CFy_Inv        = new double[nMarker];
  CFz_Inv        = new double[nMarker];

  Surface_CLift_Inv= new double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv= new double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv  = new double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv  = new double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv  = new double[config->GetnMarker_Monitoring()];
  Surface_CLift    = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag    = new double[config->GetnMarker_Monitoring()];
  Surface_CMx      = new double[config->GetnMarker_Monitoring()];
  Surface_CMy      = new double[config->GetnMarker_Monitoring()];
  Surface_CMz      = new double[config->GetnMarker_Monitoring()];

  /*--- Rotational coefficients ---*/
  CMerit_Inv = new double[nMarker];
  CT_Inv     = new double[nMarker];
  CQ_Inv     = new double[nMarker];

  /*--- Supersonic coefficients ---*/
  CEquivArea_Inv   = new double[nMarker];
  CNearFieldOF_Inv = new double[nMarker];

  /*--- Nacelle simulation ---*/
  FanFace_MassFlow  = new double[nMarker];
  Exhaust_MassFlow  = new double[nMarker];
  Exhaust_Area      = new double[nMarker];
  FanFace_Pressure  = new double[nMarker];
  FanFace_Mach      = new double[nMarker];
  FanFace_Area      = new double[nMarker];

  /*--- Init total coefficients ---*/
  Total_CDrag   = 0.0;	Total_CLift       = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy         = 0.0;  Total_CMz          = 0.0;
  Total_CEff    = 0.0;	Total_CEquivArea  = 0.0;  Total_CNearFieldOF = 0.0;
  Total_CFx     = 0.0;	Total_CFy         = 0.0;  Total_CFz          = 0.0;
  Total_CT      = 0.0;	Total_CQ          = 0.0;  Total_CMerit       = 0.0;
  Total_MaxHeat = 0.0;  Total_Heat        = 0.0;
  Total_CpDiff  = 0.0;  Total_HeatFluxDiff    = 0.0;

  ForceViscous    = new double[3];
  MomentViscous   = new double[3];
  CDrag_Visc      = new double[nMarker];
  CLift_Visc      = new double[nMarker];
  CSideForce_Visc = new double[nMarker];
  CMx_Visc        = new double[nMarker];
  CMy_Visc        = new double[nMarker];
  CMz_Visc        = new double[nMarker];
  CEff_Visc       = new double[nMarker];
  CFx_Visc        = new double[nMarker];
  CFy_Visc        = new double[nMarker];
  CFz_Visc        = new double[nMarker];
  CMerit_Visc     = new double[nMarker];
  CT_Visc         = new double[nMarker];
  CQ_Visc         = new double[nMarker];
  Heat_Visc       = new double[nMarker];
  MaxHeatFlux_Visc   = new double[nMarker];

  Surface_CLift_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc = new double[config->GetnMarker_Monitoring()];

  /*--- Read farfield conditions from config ---*/
  Density_Inf   = config->GetDensity_FreeStreamND();
  Pressure_Inf  = config->GetPressure_FreeStreamND();
  Velocity_Inf  = config->GetVelocity_FreeStreamND();
  Energy_Inf    = config->GetEnergy_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Mach_Inf      = config->GetMach();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
  Tke_Inf       = config->GetTke_FreeStreamND();

  /*--- Initializate Fan Face Pressure ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    FanFace_MassFlow[iMarker] = 0.0;
    Exhaust_MassFlow[iMarker] = 0.0;
    Exhaust_Area[iMarker] = 0.0;
    FanFace_Pressure[iMarker] = Pressure_Inf;
    FanFace_Mach[iMarker] = Mach_Inf;
  }

  /*--- Initialize the cauchy critera array for fixed CL mode ---*/

  if (config->GetFixed_CL_Mode())
    Cauchy_Serie = new double [config->GetCauchy_Elems()+1];

  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/

  if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {

    /*--- Restart the solution from the free-stream state ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
  }

  else {

    /*--- Initialize the solution from the restart file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();

    /*--- Modify file name for an unsteady restart ---*/
    if (dual_time) {
      int Unst_RestartIter;
      if (adjoint) {
        Unst_RestartIter = int(config->GetUnst_AdjointIter()) - 1;
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = int(config->GetUnst_RestartIter())-1;
      else
        Unst_RestartIter = int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }

    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
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
      }
      iPoint_Global++;
    }

    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
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
      Pressure    = Gamma_Minus_One*Density*(node[iPoint]->GetSolution(nDim+1)/Density-0.5*Velocity2);
      Temperature = Pressure / ( Gas_Constant * Density);
      if ((Pressure < 0.0) || (Temperature < 0.0)) {
        Solution[0] = Density_Inf;
        for (iDim = 0; iDim < nDim; iDim++)
          Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
        Solution[nDim+1] = Energy_Inf*Density_Inf;
        node[iPoint]->SetSolution(Solution);
        node[iPoint]->SetSolution_Old(Solution);
        counter_local++;
      }
    }
#ifdef HAVE_MPI
    MPI_Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
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
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    space_centered = true;
  else space_centered = false;

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;

  /*--- Perform the MPI communication of the solution ---*/
  Set_MPI_Solution(geometry, config);

}

CNSSolver::~CNSSolver(void) {
  unsigned short iMarker;

  if (CDrag_Visc != NULL)      delete [] CDrag_Visc;
  if (CLift_Visc != NULL)      delete [] CLift_Visc;
  if (CSideForce_Visc != NULL) delete [] CSideForce_Visc;
  if (CMx_Visc != NULL)        delete [] CMx_Visc;
  if (CMy_Visc != NULL)        delete [] CMy_Visc;
  if (CMz_Visc != NULL)        delete [] CMz_Visc;
  if (CFx_Visc != NULL)        delete [] CFx_Visc;
  if (CFy_Visc != NULL)        delete [] CFy_Visc;
  if (CFz_Visc != NULL)        delete [] CFz_Visc;
  if (CEff_Visc != NULL)       delete [] CEff_Visc;
  if (CMerit_Visc != NULL)     delete [] CMerit_Visc;
  if (CT_Visc != NULL)         delete [] CT_Visc;
  if (CQ_Visc != NULL)         delete [] CQ_Visc;
  if (Heat_Visc != NULL)          delete [] Heat_Visc;
  if (MaxHeatFlux_Visc != NULL)       delete [] MaxHeatFlux_Visc;
  if (ForceViscous != NULL)    delete [] ForceViscous;
  if (MomentViscous != NULL)   delete [] MomentViscous;

  if (Surface_CLift_Visc != NULL) delete [] Surface_CLift_Visc;
  if (Surface_CDrag_Visc != NULL) delete [] Surface_CDrag_Visc;
  if (Surface_CMx_Visc != NULL)   delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)   delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)   delete [] Surface_CMz_Visc;

  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }

  if (Cauchy_Serie != NULL)
    delete [] Cauchy_Serie;

}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;
  double eddy_visc = 0.0, turb_ke = 0.0;
  bool RightSol = true;
  int rank;

#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  unsigned long ExtIter     = config->GetExtIter();
  bool adjoint              = config->GetAdjoint();
  unsigned short turb_model = config->GetKind_Turb_Model();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center               = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst           = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow         = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_turb         = ((config->GetSpatialOrder_Turb() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_adjflow      = ((config->GetSpatialOrder_AdjFlow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool compressible         = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible       = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface          = (config->GetKind_Regime() == FREESURFACE);
  bool tkeNeeded            = (turb_model == SST);
  bool fixed_cl             = config->GetFixed_CL_Mode();
  bool engine               = ((config->GetnMarker_NacelleInflow() != 0) || (config->GetnMarker_NacelleExhaust() != 0));
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS );

  /*--- Compute nacelle inflow and exhaust properties ---*/

  if (engine) GetNacelle_Properties(geometry, config, iMesh, Output);

  /*--- Update the angle of attack at the far-field for fixed CL calculations. ---*/

  if (fixed_cl) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }

  /*--- Compute distance function to zero level set ---*/

  if (freesurface) SetFreeSurface_Distance(geometry, config);

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
    }

    /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
     FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
     Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/

    if (compressible) {
    	RightSol = node[iPoint]->SetPrimVar_Compressible(eddy_visc, turb_ke, FluidModel);
    	node[iPoint]->SetSecondaryVar_Compressible(FluidModel);
    }

    if (incompressible) RightSol = node[iPoint]->SetPrimVar_Incompressible(Density_Inf, Viscosity_Inf, eddy_visc, turb_ke, config);
    if (freesurface) RightSol = node[iPoint]->SetPrimVar_FreeSurface(eddy_visc, turb_ke, config);
    if (!RightSol) ErrorCounter++;

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  /*--- Artificial dissipation ---*/

  if (center) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Compute gradient of the primitive variables ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS){
	  SetPrimitive_Gradient_GG(geometry, config);
	  if (compressible && !ideal_gas) SetSecondary_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES){
	  SetPrimitive_Gradient_LS(geometry, config);
	  if (compressible && !ideal_gas) SetSecondary_Gradient_LS(geometry, config);
  }

  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/

  if ((iMesh == MESH_0) && (limiter_flow || limiter_turb || limiter_adjflow)) { SetPrimitive_Limiter(geometry, config); }

  /*--- Initialize the Jacobian matrices ---*/

  if (implicit) Jacobian.SetValZero();

  /*--- Error message ---*/

#ifdef HAVE_MPI
  unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
  MPI_Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);

}

void CNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  double Mean_BetaInc2, *Normal, Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time, Local_Delta_Time_Visc, Mean_DensityInc,
  Global_Delta_Time = 1E6, Mean_LaminarVisc = 0.0, Mean_EddyVisc = 0.0, Mean_Density = 0.0, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  unsigned short iDim, iMarker;
  double ProjVel, ProjVel_i, ProjVel_j;

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
      double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      double *GridVel_j = geometry->node[jPoint]->GetGridVel();
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
        double *GridVel = geometry->node[iPoint]->GetGridVel();
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

  /*--- Each element uses their own speed ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
    Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
    Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
    Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
    Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
    Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
    node[iPoint]->SetDelta_Time(Local_Delta_Time);
  }

  /*--- Check if there is any element with only one neighbor...
   a CV that is inside another CV ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (geometry->node[iPoint]->GetnPoint() == 1)
      node[iPoint]->SetDelta_Time(Min_Delta_Time);
  }

  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    MPI_Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for(iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }

  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

#ifdef HAVE_MPI
    double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    MPI_Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
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

    /*--- Primitive variables, and gradient ---*/

    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());

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
  double Delta, Viscosity = 0.0, **Grad_PrimVar, div_vel, *Normal, MomentDist[3], WallDist[3],
  *Coord, *Coord_Normal, Area, WallShearStress, TauNormal, factor, RefVel2,
  RefDensity, GradTemperature, Density = 0.0, Vel[3], WallDistMod, FrictionVel,
  Mach2Vel, Mach_Motion, *Velocity_Inf, UnitNormal[3], TauElem[3], TauTangent[3], Tau[3][3], Force[3], Cp, thermal_conductivity, MaxNorm = 8.0;
  double *Origin = config->GetRefOriginMoment(0);
  string Marker_Tag, Monitoring_Tag;

  double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  double RefAreaCoeff     = config->GetRefAreaCoeff();
  double RefLengthMoment  = config->GetRefLengthMoment();
  double Gas_Constant     = config->GetGas_ConstantND();
  bool grid_movement      = config->GetGrid_Movement();
  bool compressible       = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface        = (config->GetKind_Regime() == FREESURFACE);
  double Prandtl_Lam      = config->GetPrandtl_Lam();

  /*--- For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/

  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  } else {
    Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  RefDensity  = config->GetDensity_FreeStreamND();
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);

  /*--- Variables initialization ---*/

  AllBound_CDrag_Visc = 0.0;  AllBound_CLift_Visc = 0.0;  AllBound_CSideForce_Visc = 0.0;  AllBound_CEff_Visc = 0.0;
  AllBound_CMx_Visc = 0.0;    AllBound_CMy_Visc = 0.0;    AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc = 0.0;    AllBound_CFy_Visc = 0.0;    AllBound_CFz_Visc = 0.0;
  AllBound_CT_Visc = 0.0;     AllBound_CQ_Visc = 0.0;     AllBound_CMerit_Visc = 0.0;
  AllBound_HeatFlux_Visc = 0.0;      AllBound_MaxHeatFlux_Visc = 0.0;

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift_Visc[iMarker_Monitoring] = 0.0;
    Surface_CDrag_Visc[iMarker_Monitoring] = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]   = 0.0;
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
      CDrag_Visc[iMarker] = 0.0;  CLift_Visc[iMarker] = 0.0; CSideForce_Visc[iMarker] = 0.0;  CEff_Visc[iMarker] = 0.0;
      CMx_Visc[iMarker] = 0.0;    CMy_Visc[iMarker] = 0.0;   CMz_Visc[iMarker] = 0.0;
      CFx_Visc[iMarker] = 0.0;    CFy_Visc[iMarker] = 0.0;   CFz_Visc[iMarker] = 0.0;
      CT_Visc[iMarker] = 0.0;     CQ_Visc[iMarker] = 0.0;    CMerit_Visc[iMarker] = 0.0;
      Heat_Visc[iMarker] = 0.0;      MaxHeatFlux_Visc[iMarker] = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
      MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;

      /*--- Loop over the vertices to compute the forces ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        Coord = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPointNormal]->GetCoord();

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Grad_PrimVar = node[iPoint]->GetGradient_Primitive();
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

        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Delta = 0.0; if (iDim == jDim) Delta = 1.0;
            Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
            TWO3*Viscosity*div_vel*Delta;
          }
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++)
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
        }

        /*--- Compute wall shear stress (using the stress tensor) ---*/

        TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++) TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
        WallShearStress = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallShearStress += TauTangent[iDim]*TauTangent[iDim];
        WallShearStress = sqrt(WallShearStress);

        for (iDim = 0; iDim < nDim; iDim++)
          Vel[iDim] = node[iPointNormal]->GetVelocity(iDim);

        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);

        /*--- Compute wall skin friction coefficient, and heat flux on the wall ---*/

        CSkinFriction[iMarker][iVertex] = WallShearStress / (0.5*RefDensity*RefVel2);

        /*--- Compute y+ and non-dimensional velocity ---*/

        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);

        /*--- Compute total and max heat flux on the wall (compressible solver only) ---*/

        if (compressible) {

          GradTemperature = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GradTemperature += Grad_PrimVar[0][iDim]*(-Normal[iDim]);

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

      /*--- Transform ForceViscous and MomentViscous into non-dimensional coefficient ---*/

      if  (Monitoring == YES) {
        if (nDim == 2) {
          CDrag_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CLift_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker]        = CLift_Visc[iMarker]/(CDrag_Visc[iMarker]+EPS);
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CT_Visc[iMarker]          = -CFx_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker]/CQ_Visc[iMarker];
          MaxHeatFlux_Visc[iMarker] = pow(MaxHeatFlux_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          CDrag_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          CLift_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSideForce_Visc[iMarker]  = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker]        = CLift_Visc[iMarker]/(CDrag_Visc[iMarker]+EPS);
          CMx_Visc[iMarker]         = MomentViscous[0];
          CMy_Visc[iMarker]         = MomentViscous[1];
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CFz_Visc[iMarker]         = ForceViscous[2];
          CT_Visc[iMarker]          = -CFz_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker]/CQ_Visc[iMarker];
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
            Surface_CLift_Visc[iMarker_Monitoring] += CLift_Visc[iMarker];
            Surface_CDrag_Visc[iMarker_Monitoring] += CDrag_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]   += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]   += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]   += CMz_Visc[iMarker];
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

  double MyAllBound_CDrag_Visc        = AllBound_CDrag_Visc;                      AllBound_CDrag_Visc = 0.0;
  double MyAllBound_CLift_Visc        = AllBound_CLift_Visc;                      AllBound_CLift_Visc = 0.0;
  double MyAllBound_CSideForce_Visc   = AllBound_CSideForce_Visc;                 AllBound_CSideForce_Visc = 0.0;
  double MyAllBound_CEff_Visc         = AllBound_CEff_Visc;                       AllBound_CEff_Visc = 0.0;
  double MyAllBound_CMx_Visc          = AllBound_CMx_Visc;                        AllBound_CMx_Visc = 0.0;
  double MyAllBound_CMy_Visc          = AllBound_CMy_Visc;                        AllBound_CMy_Visc = 0.0;
  double MyAllBound_CMz_Visc          = AllBound_CMz_Visc;                        AllBound_CMz_Visc = 0.0;
  double MyAllBound_CFx_Visc          = AllBound_CFx_Visc;                        AllBound_CFx_Visc = 0.0;
  double MyAllBound_CFy_Visc          = AllBound_CFy_Visc;                        AllBound_CFy_Visc = 0.0;
  double MyAllBound_CFz_Visc          = AllBound_CFz_Visc;                        AllBound_CFz_Visc = 0.0;
  double MyAllBound_CT_Visc           = AllBound_CT_Visc;                         AllBound_CT_Visc = 0.0;
  double MyAllBound_CQ_Visc           = AllBound_CQ_Visc;                         AllBound_CQ_Visc = 0.0;
  double MyAllBound_CMerit_Visc       = AllBound_CMerit_Visc;                     AllBound_CMerit_Visc = 0.0;
  double MyAllBound_HeatFlux_Visc     = AllBound_HeatFlux_Visc;                       AllBound_HeatFlux_Visc = 0.0;
  double MyAllBound_MaxHeatFlux_Visc  = pow(AllBound_MaxHeatFlux_Visc, MaxNorm);  AllBound_MaxHeatFlux_Visc = 0.0;

  MPI_Allreduce(&MyAllBound_CDrag_Visc, &AllBound_CDrag_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CLift_Visc, &AllBound_CLift_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CSideForce_Visc, &AllBound_CSideForce_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CEff_Visc = AllBound_CLift_Visc / (AllBound_CDrag_Visc + EPS);
  MPI_Allreduce(&MyAllBound_CMx_Visc, &AllBound_CMx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CMy_Visc, &AllBound_CMy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CMz_Visc, &AllBound_CMz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFx_Visc, &AllBound_CFx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFy_Visc, &AllBound_CFy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CFz_Visc, &AllBound_CFz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CT_Visc, &AllBound_CT_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_CQ_Visc, &AllBound_CQ_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  MPI_Allreduce(&MyAllBound_HeatFlux_Visc, &AllBound_HeatFlux_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyAllBound_MaxHeatFlux_Visc, &AllBound_MaxHeatFlux_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  AllBound_MaxHeatFlux_Visc = pow(AllBound_MaxHeatFlux_Visc, 1.0/MaxNorm);

  /*--- Add the forces on the surfaces using all the nodes ---*/
  double *MySurface_CLift_Visc = NULL;
  double *MySurface_CDrag_Visc = NULL;
  double *MySurface_CMx_Visc   = NULL;
  double *MySurface_CMy_Visc   = NULL;
  double *MySurface_CMz_Visc   = NULL;

  MySurface_CLift_Visc = new double[config->GetnMarker_Monitoring()];
  MySurface_CDrag_Visc = new double[config->GetnMarker_Monitoring()];
  MySurface_CMx_Visc   = new double[config->GetnMarker_Monitoring()];
  MySurface_CMy_Visc   = new double[config->GetnMarker_Monitoring()];
  MySurface_CMz_Visc   = new double[config->GetnMarker_Monitoring()];

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    MySurface_CLift_Visc[iMarker_Monitoring] = Surface_CLift_Visc[iMarker_Monitoring];
    MySurface_CDrag_Visc[iMarker_Monitoring] = Surface_CDrag_Visc[iMarker_Monitoring];
    MySurface_CMx_Visc[iMarker_Monitoring]   = Surface_CMx_Visc[iMarker_Monitoring];
    MySurface_CMy_Visc[iMarker_Monitoring]   = Surface_CMy_Visc[iMarker_Monitoring];
    MySurface_CMz_Visc[iMarker_Monitoring]   = Surface_CMz_Visc[iMarker_Monitoring];
    Surface_CLift_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CDrag_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]     = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]     = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]     = 0.0;
  }

  MPI_Allreduce(MySurface_CLift_Visc, Surface_CLift_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CDrag_Visc, Surface_CDrag_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CMx_Visc, Surface_CMx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CMy_Visc, Surface_CMy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MySurface_CMz_Visc, Surface_CMz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete [] MySurface_CLift_Visc;
  delete [] MySurface_CDrag_Visc;
  delete [] MySurface_CMx_Visc;
  delete [] MySurface_CMy_Visc;
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
  Total_CMerit      += AllBound_CMerit_Visc;
  Total_Heat        = AllBound_HeatFlux_Visc;
  Total_MaxHeat     = AllBound_MaxHeatFlux_Visc;

  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/

  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift[iMarker_Monitoring]     += Surface_CLift_Visc[iMarker_Monitoring];
    Surface_CDrag[iMarker_Monitoring]     += Surface_CDrag_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]       += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]       += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]       += Surface_CMz_Visc[iMarker_Monitoring];
  }

}

void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  /*--- Local variables ---*/
  unsigned short iDim, jDim, iVar, jVar;
  unsigned long iVertex, iPoint, Point_Normal, total_index;

  double Wall_HeatFlux, dist_ij, *Coord_i, *Coord_j, theta2;
  double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  double ProjGridVel, *GridVel, GridVel2, *Normal, Area, Pressure = 0.0;
  double total_viscosity, div_vel, Density, turb_ke, tau_vel[3], UnitNormal[3];
  double laminar_viscosity = 0.0, eddy_viscosity = 0.0, **grad_primvar, tau[3][3];
  double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

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
      if (compressible)   node[iPoint]->SetVelocity_Old(Vector);
      if (incompressible || freesurface) node[iPoint]->SetVelocityInc_Old(Vector);

      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();

      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      Res_Visc[nDim+1] = Wall_HeatFlux * Area;

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
        grad_primvar      = node[iPoint]->GetGradient_Primitive();

        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        else
          turb_ke = 0.0;

        /*--- Divergence of the velocity ---*/
        div_vel = 0.0;
        for (iDim = 0 ; iDim < nDim; iDim++)
          div_vel += grad_primvar[iDim+1][iDim];

        /*--- Compute the viscous stress tensor ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( grad_primvar[jDim+1][iDim]
                                               +grad_primvar[iDim+1][jDim] )
            - TWO3*total_viscosity*div_vel*delta[iDim][jDim]
            - TWO3*Density*turb_ke*delta[iDim][jDim];
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

  double *Normal, *Coord_i, *Coord_j, Area, dist_ij, theta2;
  double Twall, Temperature, dTdn, dTdrho, thermal_conductivity;
  double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  double ProjGridVel, *GridVel, GridVel2, Pressure = 0.0, Density, Vel2, Energy;
  double total_viscosity, div_vel, turb_ke, tau_vel[3], UnitNormal[3];
  double laminar_viscosity, eddy_viscosity, **grad_primvar, tau[3][3];
  double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  double Prandtl_Lam  = config->GetPrandtl_Lam();
  double Prandtl_Turb = config->GetPrandtl_Turb();
  double Gas_Constant = config->GetGas_ConstantND();
  double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool freesurface    = (config->GetKind_Regime() == FREESURFACE);
  bool grid_movement  = config->GetGrid_Movement();

  Point_Normal = 0;

  /*--- Identify the boundary ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature ---*/

  Twall = config->GetIsothermal_Temperature(Marker_Tag);

  /*--- Loop over boundary points ---*/

  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

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
        Energy  = node[iPoint]->GetSolution(nDim+1);
        Temperature = node[iPoint]->GetPrimitive(0);
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
        grad_primvar      = node[iPoint]->GetGradient_Primitive();

        /*--- Turbulent kinetic energy ---*/

        if (config->GetKind_Turb_Model() == SST)
          turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        else
          turb_ke = 0.0;

        /*--- Divergence of the velocity ---*/

        div_vel = 0.0;
        for (iDim = 0 ; iDim < nDim; iDim++)
          div_vel += grad_primvar[iDim+1][iDim];

        /*--- Compute the viscous stress tensor ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( grad_primvar[jDim+1][iDim]
                                               +grad_primvar[iDim+1][jDim] )
            - TWO3*total_viscosity*div_vel*delta[iDim][jDim]
            - TWO3*Density*turb_ke*delta[iDim][jDim];
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
