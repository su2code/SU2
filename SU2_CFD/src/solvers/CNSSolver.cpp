/*!
 * \file CNSSolver.cpp
 * \brief Main subrotuines for solving Finite-Volume Navier-Stokes flow problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CEulerSolver() {

  unsigned long iPoint, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  su2double Density, Velocity2, Pressure, Temperature, StaticEnergy;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart    = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter = 0;
  bool dual_time = (config->GetTime_Marching() == DT_STEPPING_1ST) ||
                   (config->GetTime_Marching() == DT_STEPPING_2ND);
  bool time_stepping = (config->GetTime_Marching() == TIME_STEPPING);

  /*--- A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain ---*/
  dynamic_grid = config->GetDynamic_Grid();

  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  string filename_ = "flow";

  unsigned short direct_diff = config->GetDirectDiff();
  bool rans = (config->GetKind_Turb_Model() != NONE);

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1) &&
      (config->GetFixed_CL_Mode() || config->GetFixed_CM_Mode())) {

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

    Read_SU2_Restart_Metadata(geometry, config, adjoint, filename_);

  }

  /*--- Set the gamma value ---*/

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp).
   ---*/

  nDim = geometry->GetnDim();

  nVar = nDim+2;
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  nSecondaryVar = 8; nSecondaryVarGrad = 2;


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

  /*--- Define some auxiliar vector related with the residual ---*/

  Residual      = new su2double[nVar]();
  Residual_RMS  = new su2double[nVar]();
  Residual_Max  = new su2double[nVar]();
  Residual_i    = new su2double[nVar]();
  Residual_j    = new su2double[nVar]();
  Res_Conv      = new su2double[nVar]();
  Res_Visc      = new su2double[nVar]();
  Res_Sour      = new su2double[nVar]();

  /*--- Define some structures for locating max residuals ---*/

  Point_Max = new unsigned long[nVar];
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim]();
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution = new su2double[nVar]();

  /*--- Define some auxiliary vectors related to the geometry ---*/

  Vector = new su2double[nDim]();

  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
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

  if (config->GetLeastSquaresRequired()) {
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];

    Cvector = new su2double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }

  /*--- Allocates a 2D array with variable "outer" sizes and init to 0. ---*/

  auto Alloc2D = [](unsigned long M, const unsigned long* N, su2double**& X) {
    X = new su2double* [M];
    for(unsigned long i = 0; i < M; ++i)
      X[i] = new su2double [N[i]] ();
  };

  /*--- Allocates a 3D array with variable "middle" sizes and init to 0. ---*/

  auto Alloc3D = [](unsigned long M, const unsigned long* N, unsigned long P, su2double***& X) {
    X = new su2double** [M];
    for(unsigned long i = 0; i < M; ++i) {
      X[i] = new su2double* [N[i]];
      for(unsigned long j = 0; j < N[i]; ++j)
        X[i][j] = new su2double [P] ();
    }
  };

  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/

  Alloc3D(nMarker, nVertex, nPrimVar, CharacPrimVar);

  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/

  Alloc3D(nMarker, nVertex, (rans? nPrimVar+2 : nPrimVar), DonorPrimVar);

  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/

  DonorGlobalIndex = new unsigned long* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    DonorGlobalIndex[iMarker] = new unsigned long [geometry->nVertex[iMarker]]();
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

  /*--- Store the value of the Delta P at the Actuator Disk ---*/

  Alloc2D(nMarker, nVertex, ActDisk_DeltaP);

  /*--- Store the value of the Delta T at the Actuator Disk ---*/

  Alloc2D(nMarker, nVertex, ActDisk_DeltaT);

  /*--- Store the value of the Total Pressure at the inlet BC ---*/

  Alloc2D(nMarker, nVertex, Inlet_Ttotal);

  /*--- Store the value of the Total Temperature at the inlet BC ---*/

  Alloc2D(nMarker, nVertex, Inlet_Ptotal);

  /*--- Store the value of the Flow direction at the inlet BC ---*/

  Alloc3D(nMarker, nVertex, nDim, Inlet_FlowDir);

  /*--- Inviscid force definition and coefficient in all the markers ---*/

  Alloc2D(nMarker, nVertex, CPressure);
  Alloc2D(nMarker, nVertex, CPressureTarget);

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
      CSkinFriction[iMarker][iDim] = new su2double[geometry->nVertex[iMarker]] ();
    }
  }

  /*--- Non dimensional aerodynamic coefficients ---*/

  InvCoeff.allocate(nMarker);
  MntCoeff.allocate(nMarker);
  ViscCoeff.allocate(nMarker);
  SurfaceInvCoeff.allocate(config->GetnMarker_Monitoring());
  SurfaceMntCoeff.allocate(config->GetnMarker_Monitoring());
  SurfaceViscCoeff.allocate(config->GetnMarker_Monitoring());
  SurfaceCoeff.allocate(config->GetnMarker_Monitoring());

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

  /*--- Initialize fan face pressure, fan face mach number, and mass flow rate ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;

    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;

  }
  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++){

    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)]();
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)]();

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++)
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1]();
    }
  }

  /*--- Only initialize when there is a Marker_Fluid_Load
   *--- (this avoids overhead in all other cases while a more permanent structure is being developed) ---*/
  if((config->GetnMarker_Fluid_Load() > 0) && (MGLevel == MESH_0)){

    InitVertexTractionContainer(geometry, config);

    if (config->GetDiscrete_Adjoint())
      InitVertexTractionAdjointContainer(geometry, config);

  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

#ifdef HAVE_OMP
  /*--- On the fine grid get the edge coloring, on coarse grids (which are difficult to color)
   *    setup the reducer strategy, i.e. one loop over edges followed by a point loop to sum
   *    the fluxes for each cell and set the diagonal of the system matrix. ---*/

  const auto& coloring = geometry->GetEdgeColoring();

  if (!coloring.empty()) {
    auto nColor = coloring.getOuterSize();
    EdgeColoring.resize(nColor);

    for(auto iColor = 0ul; iColor < nColor; ++iColor) {
      EdgeColoring[iColor].size = coloring.getNumNonZeros(iColor);
      EdgeColoring[iColor].indices = coloring.innerIdx(iColor);
    }
  }
  ColorGroupSize = geometry->GetEdgeColorGroupSize();

  if (MGLevel != MESH_0) {
    EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
  }

  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), OMP_MAX_SIZE);
#endif

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/

  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    Density = nodes->GetDensity(iPoint);

    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += pow(nodes->GetSolution(iPoint,iDim+1)/Density,2);

    StaticEnergy= nodes->GetEnergy(iPoint) - 0.5*Velocity2;

    GetFluidModel()->SetTDState_rhoe(Density, StaticEnergy);
    Pressure= GetFluidModel()->GetPressure();
    Temperature= GetFluidModel()->GetTemperature();

    /*--- Use the values at the infinity ---*/

    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      nodes->SetSolution(iPoint,Solution);
      nodes->SetSolution_Old(iPoint,Solution);
      counter_local++;
    }

  }

  /*--- Warning message about non-physical points ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  }

  /*--- Initialize the BGS residuals in FSI problems. ---*/
  if (config->GetMultizone_Residual()){
    Residual_BGS      = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]  = 1.0;
    Residual_Max_BGS  = new su2double[nVar];   for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS = new unsigned long[nVar]();
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim]();
    }
  }

  /*--- Define solver parameters needed for execution of destructor ---*/

  space_centered = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  euler_implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  least_squares = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);

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

  /* Store the initial CFL number for all grid points. */

  const su2double CFL = config->GetCFL(MGLevel);
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "C.FLOW";

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

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  /// TODO: Try to start a parallel section here encompassing all "heavy" methods.

  unsigned long ErrorCounter = 0;

  unsigned long InnerIter   = config->GetInnerIter();
  bool cont_adjoint         = config->GetContinuous_Adjoint();
  bool disc_adjoint         = config->GetDiscrete_Adjoint();
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center               = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst           = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow         = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool limiter_turb         = (config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool limiter_adjflow      = (cont_adjoint && (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter()));
  bool fixed_cl             = config->GetFixed_CL_Mode();
  bool engine               = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk        = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool nearfield            = (config->GetnMarker_NearFieldBound() != 0);
  bool van_albada           = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  unsigned short kind_row_dissipation = config->GetKind_RoeLowDiss();
  bool roe_low_dissipation  = (kind_row_dissipation != NO_ROELOWDISS) &&
                              (config->GetKind_Upwind_Flow() == ROE ||
                               config->GetKind_Upwind_Flow() == SLAU ||
                               config->GetKind_Upwind_Flow() == SLAU2);
  bool wall_functions       = config->GetWall_Functions();

  /*--- Update the angle of attack at the far-field for fixed CL calculations (only direct problem). ---*/

  if ((fixed_cl) && (!disc_adjoint) && (!cont_adjoint)) { SetFarfield_AoA(geometry, solver_container, config, iMesh, Output); }

  /*--- Set the primitive variables ---*/

  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Compute the engine properties ---*/

  if (engine) { GetPower_Properties(geometry, config, iMesh, Output); }

  /*--- Compute the actuator disk properties and distortion levels ---*/

  if (actuator_disk) {
    Set_MPI_ActDisk(solver_container, geometry, config);
    SetActDisk_BCThrust(geometry, solver_container, config, iMesh, Output);
  }

  /*--- Compute NearField MPI ---*/

  if (nearfield) { Set_MPI_Nearfield(geometry, config); }

  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Roe Low Dissipation Sensor ---*/

  if (roe_low_dissipation){
    SetRoe_Dissipation(geometry, config);
    if (kind_row_dissipation == FD_DUCROS || kind_row_dissipation == NTS_DUCROS){
      SetUpwind_Ducros_Sensor(geometry, config);
    }
  }

  /*--- Compute gradient for MUSCL reconstruction. ---*/

  if (config->GetReconstructionGradientRequired() && (iMesh == MESH_0)) {
    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetPrimitive_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
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

  /*--- Evaluate the vorticity and strain rate magnitude ---*/

  StrainMag_Max = 0.0; Omega_Max = 0.0;

  SU2_OMP_PARALLEL
  {
    solver_container[FLOW_SOL]->GetNodes()->SetVorticity_StrainMag();

    su2double strainMax = 0.0, omegaMax = 0.0;

    SU2_OMP(for schedule(static,omp_chunk_size) nowait)
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

      su2double StrainMag = solver_container[FLOW_SOL]->GetNodes()->GetStrainMag(iPoint);
      const su2double* Vorticity = solver_container[FLOW_SOL]->GetNodes()->GetVorticity(iPoint);
      su2double Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);

      strainMax = max(strainMax, StrainMag);
      omegaMax = max(omegaMax, Omega);

    }
    SU2_OMP_CRITICAL
    {
      StrainMag_Max = max(StrainMag_Max, strainMax);
      Omega_Max = max(Omega_Max, omegaMax);
    }

  } // end SU2_OMP_PARALLEL

  /*--- Compute the TauWall from the wall functions ---*/

  if (wall_functions)
    SetTauWall_WF(geometry, solver_container, config);

  /*--- Initialize the Jacobian matrices ---*/

  if (implicit && !config->GetDiscrete_Adjoint()) Jacobian.SetValZero();

  /*--- Error message ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    if (iMesh == MESH_0) {

      unsigned long MyErrorCounter = ErrorCounter;
      su2double MyOmega_Max = Omega_Max;
      su2double MyStrainMag_Max = StrainMag_Max;

      SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      config->SetNonphysical_Points(ErrorCounter);
      solver_container[FLOW_SOL]->SetStrainMag_Max(StrainMag_Max);
      solver_container[FLOW_SOL]->SetOmega_Max(Omega_Max);

    }

  }

}

unsigned long CNSSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {

  unsigned long nonPhysicalPoints = 0;

  const unsigned short turb_model = config->GetKind_Turb_Model();
  const bool tkeNeeded = (turb_model == SST) || (turb_model == SST_SUST);

  SU2_OMP_PARALLEL_(for schedule(static,omp_chunk_size) reduction(+:nonPhysicalPoints))
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Retrieve the value of the kinetic energy (if needed). ---*/

    su2double eddy_visc = 0.0, turb_ke = 0.0;

    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);

      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
        su2double DES_LengthScale = solver_container[TURB_SOL]->GetNodes()->GetDES_LengthScale(iPoint);
        nodes->SetDES_LengthScale(iPoint, DES_LengthScale);
      }
    }

    /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/

    bool physical = static_cast<CNSVariable*>(nodes)->SetPrimVar(iPoint, eddy_visc, turb_ke, GetFluidModel());
    nodes->SetSecondaryVar(iPoint, GetFluidModel());

    /*--- Check for non-realizable states for reporting. ---*/

    nonPhysicalPoints += !physical;

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return nonPhysicalPoints;
}

void CNSSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config) {

  const bool implicit  = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  const bool tkeNeeded = (config->GetKind_Turb_Model() == SST) ||
                         (config->GetKind_Turb_Model() == SST_SUST);

  CVariable* turbNodes = nullptr;
  if (tkeNeeded) turbNodes = solver_container[TURB_SOL]->GetNodes();

  /*--- Points, coordinates and normal vector in edge ---*/

  auto iPoint = geometry->edge[iEdge]->GetNode(0);
  auto jPoint = geometry->edge[iEdge]->GetNode(1);

  numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());

  numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

  /*--- Primitive and secondary variables. ---*/

  numerics->SetPrimitive(nodes->GetPrimitive(iPoint),
                         nodes->GetPrimitive(jPoint));

  numerics->SetSecondary(nodes->GetSecondary(iPoint),
                         nodes->GetSecondary(jPoint));

  /*--- Gradients. ---*/

  numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                               nodes->GetGradient_Primitive(jPoint));

  /*--- Turbulent kinetic energy. ---*/

  if (tkeNeeded)
    numerics->SetTurbKineticEnergy(turbNodes->GetSolution(iPoint,0),
                                   turbNodes->GetSolution(jPoint,0));

  /*--- Wall shear stress values (wall functions) ---*/

  numerics->SetTauWall(nodes->GetTauWall(iPoint),
                       nodes->GetTauWall(iPoint));

  /*--- Compute and update residual ---*/

  auto residual = numerics->ComputeResidual(config);

#ifdef HAVE_OMP
  if (MGLevel != MESH_0) {
    EdgeFluxes.SubtractBlock(iEdge, residual);
    if (implicit)
      Jacobian.UpdateBlocksSub(iEdge, residual.jacobian_i, residual.jacobian_j);
  }
  else {
#else
  {
#endif
    LinSysRes.SubtractBlock(iPoint, residual);
    LinSysRes.AddBlock(jPoint, residual);

    if (implicit)
      Jacobian.UpdateBlocksSub(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
  }

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

  su2double Prandtl_Lam = config->GetPrandtl_Lam();
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

        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*delta[iDim][jDim];
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
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);

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


void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                 CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

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

  /*--- Identify the boundary by string name ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config as well as the
        wall function treatment.---*/

  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag)/config->GetHeat_Flux_Ref();

//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != NO_WALL_FUNCTION) {
//    SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
//  }

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- If it is a customizable patch, retrieve the specified wall heat flux. ---*/

      if (config->GetMarker_All_PyCustom(val_marker)) Wall_HeatFlux = geometry->GetCustomBoundaryHeatFlux(val_marker, iVertex);

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

      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/

      Res_Visc[nDim+1] = Wall_HeatFlux * Area;

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

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( Grad_Vel[jDim][iDim]+Grad_Vel[iDim][jDim] ) -
                              TWO3*total_viscosity*div_vel*delta[iDim][jDim];
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

          Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

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

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
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

  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature from config
        as well as the wall function treatment.---*/

  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != NO_WALL_FUNCTION) {
//    SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
//  }

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

  SU2_OMP_PARALLEL_(for schedule(static,omp_chunk_size))
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

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

void CNSSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
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

  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  su2double Temperature_Ref = config->GetTemperature_Ref();

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Identify the boundary ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

//  /*--- Retrieve the specified wall function treatment.---*/
//
//  Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
//  if (Wall_Function != NO_WALL_FUNCTION) {
//      SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
//  }

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

void CNSSolver::SetTauWall_WF(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iDim, jDim, iMarker;
  unsigned long iVertex, iPoint, Point_Normal, counter;

  su2double Area, div_vel, UnitNormal[3], *Normal;
  su2double **grad_primvar, tau[3][3];

  su2double Vel[3] = {0.0, 0.0, 0.0}, VelNormal, VelTang[3], VelTangMod, VelInfMod, WallDist[3], WallDistMod;
  su2double T_Normal, P_Normal;
  su2double Density_Wall, T_Wall, P_Wall, Lam_Visc_Wall, Tau_Wall = 0.0, Tau_Wall_Old = 0.0;
  su2double *Coord, *Coord_Normal;
  su2double diff, Delta;
  su2double U_Tau, U_Plus, Gam, Beta, Phi, Q, Y_Plus_White, Y_Plus;
  su2double TauElem[3], TauNormal, TauTangent[3], WallShearStress;
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;

  unsigned short max_iter = 10;
  su2double tol = 1e-6;

  /*--- Get the freestream velocity magnitude for non-dim. purposes ---*/

  su2double *VelInf = config->GetVelocity_FreeStreamND();
  VelInfMod = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    VelInfMod += VelInf[iDim];
  VelInfMod = sqrt(VelInfMod);

  /*--- Compute the recovery factor ---*/
  // Double-check: laminar or turbulent Pr for this?
  su2double Recovery = pow(config->GetPrandtl_Lam(), (1.0/3.0));

  /*--- Typical constants from boundary layer theory ---*/

  su2double kappa = 0.4;
  su2double B = 5.5;

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

          /*--- Compute mach number ---*/

          // M_Normal = VelTangMod / sqrt(Gamma * Gas_Constant * T_Normal);

          /*--- Compute the wall temperature using the Crocco-Buseman equation ---*/

          //T_Wall = T_Normal * (1.0 + 0.5*Gamma_Minus_One*Recovery*M_Normal*M_Normal);
          T_Wall = T_Normal + Recovery*pow(VelTangMod,2.0)/(2.0*Cp);

          /*--- Extrapolate the pressure from the interior & compute the
           wall density using the equation of state ---*/

          P_Wall = P_Normal;
          Density_Wall = P_Wall/(Gas_Constant*T_Wall);

          /*--- Compute the shear stress at the wall in the regular fashion
           by using the stress tensor on the surface ---*/

          Lam_Visc_Wall = nodes->GetLaminarViscosity(iPoint);
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

          Tau_Wall_Old = WallShearStress;
          counter = 0; diff = 1.0;

          while (diff > tol) {

            /*--- Friction velocity and u+ ---*/

            U_Tau = sqrt(Tau_Wall_Old/Density_Wall);
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

            /*--- Calculate an updated value for the wall shear stress
             using the y+ value, the definition of y+, and the definition of
             the friction velocity. ---*/

            Tau_Wall = (1.0/Density_Wall)*pow(Y_Plus*Lam_Visc_Wall/WallDistMod,2.0);

            /*--- Difference between the old and new Tau. Update old value. ---*/

            diff = fabs(Tau_Wall-Tau_Wall_Old);
            Tau_Wall_Old += 0.25*(Tau_Wall-Tau_Wall_Old);

            counter++;
            if (counter > max_iter) {
              cout << "WARNING: Tau_Wall evaluation has not converged in solver_direct_mean.cpp" << endl;
              cout << Tau_Wall_Old << " " << Tau_Wall << " " << diff << endl;
              break;
            }

          }


          /*--- Store this value for the wall shear stress at the node.  ---*/

          nodes->SetTauWall(iPoint,Tau_Wall);


        }

      }

    }
  }

}
