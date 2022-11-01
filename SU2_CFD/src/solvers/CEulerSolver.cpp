/*!
 * \file CEulerSolver.cpp
 * \brief Main subrotuines for solving Finite-Volume Euler flow problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
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


#include "../../include/solvers/CEulerSolver.hpp"
#include "../../include/variables/CNSVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"

void CEulerSolver::AeroCoeffsArray::allocate(int size) {
    _size = size;
    CD = new su2double[size]; CL = new su2double[size]; CSF = new su2double[size]; CEff = new su2double[size];
    CFx = new su2double[size]; CFy = new su2double[size]; CFz = new su2double[size]; CMx = new su2double[size];
    CMy = new su2double[size]; CMz = new su2double[size]; CoPx = new su2double[size]; CoPy = new su2double[size];
    CoPz = new su2double[size]; CT = new su2double[size]; CQ = new su2double[size]; CMerit = new su2double[size];
    setZero();
}

CEulerSolver::AeroCoeffsArray::~AeroCoeffsArray() {
    delete [] CD; delete [] CL; delete [] CSF; delete [] CEff;
    delete [] CFx; delete [] CFy; delete [] CFz; delete [] CMx;
    delete [] CMy; delete [] CMz; delete [] CoPx; delete [] CoPy;
    delete [] CoPz; delete [] CT; delete [] CQ; delete [] CMerit;
}

void CEulerSolver::AeroCoeffsArray::setZero(int i) {
    CD[i] = CL[i] = CSF[i] = CEff[i] = 0.0;
    CFx[i] = CFy[i] = CFz[i] = CMx[i] = 0.0;
    CMy[i] = CMz[i] = CoPx[i] = CoPy[i] = 0.0;
    CoPz[i] = CT[i] = CQ[i] = CMerit[i] = 0.0;
}

CEulerSolver::CEulerSolver(void) : CSolver() { }

CEulerSolver::CEulerSolver(CGeometry *geometry, CConfig *config,
                           unsigned short iMesh, const bool navier_stokes) : CSolver() {
    
    /*--- Based on the navier_stokes boolean, determine if this constructor is
     *    being called by itself, or by its derived class CNSSolver. ---*/
    string description;
    unsigned short nSecVar;
    if (navier_stokes) {
        description = "Navier-Stokes";
        nSecVar = 8;
    }
    else {
        description = "Euler";
        nSecVar = 2;
    }
    
    unsigned long iPoint, counter_local = 0, counter_global = 0;
    unsigned short iVar, iDim, iMarker, nLineLets;
    su2double StaticEnergy, Density, Velocity2, Pressure, Temperature;
    unsigned short nZone = geometry->GetnZone();
    bool restart = (config->GetRestart() || config->GetRestart_Flow());
    bool rans = (config->GetKind_Turb_Model() != NONE);
    unsigned short direct_diff = config->GetDirectDiff();
    int Unst_RestartIter = 0;
    bool dual_time = (config->GetTime_Marching() == DT_STEPPING_1ST) ||
    (config->GetTime_Marching() == DT_STEPPING_2ND);
    bool time_stepping = (config->GetTime_Marching() == TIME_STEPPING);
    
    /*--- A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain ---*/
    dynamic_grid = config->GetDynamic_Grid();
    
    bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
    string filename_ = "flow";
    
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
    nSecondaryVar = nSecVar; nSecondaryVarGrad = 2;
    
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
    
    Point_Max = new unsigned long[nVar]();
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
    
#ifdef HAVE_OMP
    /*--- Get the edge coloring. If the expected parallel efficiency becomes too low setup the
     *    reducer strategy. Where one loop is performed over edges followed by a point loop to
     *    sum the fluxes for each cell and set the diagonal of the system matrix. ---*/
    
    su2double parallelEff = 1.0;
    const auto& coloring = geometry->GetEdgeColoring(&parallelEff);
    
    /*--- The decision to use the strategy is local to each rank. ---*/
    ReducerStrategy = parallelEff < COLORING_EFF_THRESH;
    
    /*--- When using the reducer force a single color to reduce the color loop overhead. ---*/
    if (ReducerStrategy && (coloring.getOuterSize()>1))
        geometry->SetNaturalEdgeColoring();
    
    if (!coloring.empty()) {
        /*--- If the reducer strategy is used we are not constrained by group
         *    size as we have no other edge loops in the Euler/NS solvers. ---*/
        auto groupSize = ReducerStrategy? 1ul : geometry->GetEdgeColorGroupSize();
        auto nColor = coloring.getOuterSize();
        EdgeColoring.reserve(nColor);
        
        for(auto iColor = 0ul; iColor < nColor; ++iColor)
            EdgeColoring.emplace_back(coloring.innerIdx(iColor), coloring.getNumNonZeros(iColor), groupSize);
    }
    
    /*--- If the reducer strategy is not being forced (by EDGE_COLORING_GROUP_SIZE=0) print some messages. ---*/
    if (config->GetEdgeColoringGroupSize() != 1<<30) {
        
        su2double minEff = 1.0;
        SU2_MPI::Reduce(&parallelEff, &minEff, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
        
        int tmp = ReducerStrategy, numRanksUsingReducer = 0;
        SU2_MPI::Reduce(&tmp, &numRanksUsingReducer, 1, MPI_INT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
        
        if (minEff < COLORING_EFF_THRESH) {
            cout << "WARNING: On " << numRanksUsingReducer << " MPI ranks the coloring efficiency was less than "
            << COLORING_EFF_THRESH << " (min value was " << minEff << ").\n"
            << "         Those ranks will now use a fallback strategy, better performance may be possible\n"
            << "         with a different value of config option EDGE_COLORING_GROUP_SIZE (default 512)." << endl;
        }
    }
    
    if (ReducerStrategy)
        EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);
    
    omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), OMP_MAX_SIZE);
#else
    EdgeColoring[0] = DummyGridColor<>(geometry->GetnEdge());
#endif
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    
    if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
        
        Jacobian_i = new su2double* [nVar];
        Jacobian_j = new su2double* [nVar];
        for (iVar = 0; iVar < nVar; iVar++) {
            Jacobian_i[iVar] = new su2double [nVar];
            Jacobian_j[iVar] = new su2double [nVar];
        }
        
        if (rank == MASTER_NODE)
            cout << "Initialize Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;
        
        Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
        
        if (config->GetKind_Linear_Solver_Prec() == LINELET) {
            nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
            if (rank == MASTER_NODE)
                cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
        }
        
    }
    
    else {
        if (rank == MASTER_NODE)
            cout << "Explicit scheme. No Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;
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
        DonorGlobalIndex[iMarker] = new unsigned long [nVertex[iMarker]]();
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
    
    /*--- Force definition and coefficient arrays for all of the markers ---*/
    
    Alloc2D(nMarker, nVertex, CPressure);
    Alloc2D(nMarker, nVertex, CPressureTarget);
    
    /*--- Non dimensional aerodynamic coefficients ---*/
    
    InvCoeff.allocate(nMarker);
    MntCoeff.allocate(nMarker);
    SurfaceInvCoeff.allocate(config->GetnMarker_Monitoring());
    SurfaceMntCoeff.allocate(config->GetnMarker_Monitoring());
    SurfaceCoeff.allocate(config->GetnMarker_Monitoring());
    
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
    Mach_Inf        = config->GetMach();
    
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
    
    SlidingState       = new su2double*** [nMarker]();
    SlidingStateNodes  = new int*         [nMarker]();
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        
        if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
            
            SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)]();
            SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)]();
            
            for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++)
                SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1]();
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
    
    if (navier_stokes) {
        nodes = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nPoint, nDim, nVar, config);
    } else {
        nodes = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nPoint, nDim, nVar, config);
    }
    SetBaseClassPointerToNodes();
    
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
            cout << "Warning. The original solution contains " << counter_global << " points that are not physical." << endl;
    }
    
    /*--- Initialize the BGS residuals in FSI problems. ---*/
    if (config->GetMultizone_Residual()){
        Residual_BGS     = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar] = 1.0;
        Residual_Max_BGS = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar] = 1.0;
        
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
    
    /*--- Finally, check that the static arrays will be large enough (keep this
     *    check at the bottom to make sure we consider the "final" values). ---*/
    if((nDim > MAXNDIM) || (nPrimVar > MAXNVAR) || (nSecondaryVar > MAXNVAR))
        SU2_MPI::Error("Oops! The CEulerSolver static array sizes are not large enough.",CURRENT_FUNCTION);
}

CEulerSolver::~CEulerSolver(void) {
    
    unsigned short iVar, iMarker, iSpan;
    
    unsigned long iVertex;
    
    /*--- Array deallocation ---*/
    
    delete [] CEquivArea_Inv;
    delete [] CNearFieldOF_Inv;
    
    delete [] Inflow_MassFlow;
    delete [] Exhaust_MassFlow;
    delete [] Exhaust_Area;
    delete [] Inflow_Pressure;
    delete [] Inflow_Mach;
    delete [] Inflow_Area;
    
    delete [] Exhaust_Pressure;
    delete [] Exhaust_Temperature;
    
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
            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
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
    
    if (DonorPrimVar != NULL) {
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                delete [] DonorPrimVar[iMarker][iVertex];
            delete [] DonorPrimVar[iMarker];
        }
        delete [] DonorPrimVar;
    }
    
    if (DonorGlobalIndex != NULL) {
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] DonorGlobalIndex[iMarker];
        delete [] DonorGlobalIndex;
    }
    
    if (ActDisk_DeltaP != NULL) {
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] ActDisk_DeltaP[iMarker];
        delete [] ActDisk_DeltaP;
    }
    
    if (ActDisk_DeltaT != NULL) {
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] ActDisk_DeltaT[iMarker];
        delete [] ActDisk_DeltaT;
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
    
    for(auto& model : FluidModel) delete model;
    
    if(AverageVelocity !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
                delete [] AverageVelocity[iMarker][iSpan];
            delete  [] AverageVelocity[iMarker];
        }
        delete [] AverageVelocity;
    }
    
    if(AverageTurboVelocity !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
                delete [] AverageTurboVelocity[iMarker][iSpan];
            delete  [] AverageTurboVelocity[iMarker];
        }
        delete [] AverageTurboVelocity;
    }
    
    if(OldAverageTurboVelocity != NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
                delete [] OldAverageTurboVelocity[iMarker][iSpan];
            delete  [] OldAverageTurboVelocity[iMarker];
        }
        delete [] OldAverageTurboVelocity;
    }
    
    if(ExtAverageTurboVelocity !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
                delete [] ExtAverageTurboVelocity[iMarker][iSpan];
            delete  [] ExtAverageTurboVelocity[iMarker];
        }
        delete [] ExtAverageTurboVelocity;
    }
    
    if(AverageFlux !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
                delete [] AverageFlux[iMarker][iSpan];
            delete  [] AverageFlux[iMarker];
        }
        delete [] AverageFlux;
    }
    
    if(SpanTotalFlux !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++)
                delete [] SpanTotalFlux[iMarker][iSpan];
            delete  [] SpanTotalFlux[iMarker];
        }
        delete [] SpanTotalFlux;
    }
    
    if(AveragePressure !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] AveragePressure[iMarker];
        delete [] AveragePressure;
    }
    
    if(OldAveragePressure !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] OldAveragePressure[iMarker];
        delete [] OldAveragePressure;
    }
    
    if(RadialEquilibriumPressure !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] RadialEquilibriumPressure[iMarker];
        delete [] RadialEquilibriumPressure;
    }
    
    if(ExtAveragePressure !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] ExtAveragePressure[iMarker];
        delete [] ExtAveragePressure;
    }
    
    if(AverageDensity !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] AverageDensity[iMarker];
        delete [] AverageDensity;
    }
    
    if(OldAverageDensity !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] OldAverageDensity[iMarker];
        delete [] OldAverageDensity;
    }
    
    if(ExtAverageDensity !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] ExtAverageDensity[iMarker];
        delete [] ExtAverageDensity;
    }
    
    if(AverageKine !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] AverageKine[iMarker];
        delete [] AverageKine;
    }
    
    if(AverageOmega !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] AverageOmega[iMarker];
        delete [] AverageOmega;
    }
    
    if(AverageNu !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] AverageNu[iMarker];
        delete [] AverageNu;
    }
    
    if(ExtAverageKine !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] ExtAverageKine[iMarker];
        delete [] ExtAverageKine;
    }
    
    if(ExtAverageOmega !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] ExtAverageOmega[iMarker];
        delete [] ExtAverageOmega;
    }
    
    if(ExtAverageNu !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            delete [] ExtAverageNu[iMarker];
        delete [] ExtAverageNu;
    }
    
    if(TurboVelocityIn !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
            for (iSpan = 0; iSpan < nSpanMax + 1; iSpan++){
                delete [] TurboVelocityIn[iMarker][iSpan];
            }
            delete [] TurboVelocityIn[iMarker];
        }
        delete [] TurboVelocityIn;
    }
    
    if(TurboVelocityOut !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
            for (iSpan = 0; iSpan < nSpanMax + 1; iSpan++){
                delete [] TurboVelocityOut[iMarker][iSpan];
            }
            delete [] TurboVelocityOut[iMarker];
        }
        delete [] TurboVelocityOut;
    }
    
    if(DensityIn !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] DensityIn[iMarker];
        delete [] DensityIn;
    }
    
    if(PressureIn !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] PressureIn[iMarker];
        delete [] PressureIn;
    }
    
    if(DensityOut !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] DensityOut[iMarker];
        delete [] DensityOut;
    }
    
    if(PressureOut !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] PressureOut[iMarker];
        delete [] PressureOut;
    }
    
    if(KineIn !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] KineIn[iMarker];
        delete [] KineIn;
    }
    
    if(OmegaIn !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] OmegaIn[iMarker];
        delete [] OmegaIn;
    }
    
    if(NuIn !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] NuIn[iMarker];
        delete [] NuIn;
    }
    
    if(KineOut !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] KineOut[iMarker];
        delete [] KineOut;
    }
    
    if(OmegaOut !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] OmegaOut[iMarker];
        delete [] OmegaOut;
    }
    
    if(NuOut !=NULL){
        for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++)
            delete  [] NuOut[iMarker];
        delete [] NuOut;
    }
    
    if(CkInflow !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
                delete [] CkInflow[iMarker][iSpan];
            delete  [] CkInflow[iMarker];
        }
        delete [] CkInflow;
    }
    
    if(CkOutflow1 !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
                delete [] CkOutflow1[iMarker][iSpan];
            delete  [] CkOutflow1[iMarker];
        }
        delete [] CkOutflow1;
    }
    
    if(CkOutflow2 !=NULL){
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for(iSpan = 0; iSpan <nSpanWiseSections ; iSpan++)
                delete [] CkOutflow2[iMarker][iSpan];
            delete  [] CkOutflow2[iMarker];
        }
        delete [] CkOutflow2;
    }
    
    delete nodes;
}

void CEulerSolver::InitTurboContainers(CGeometry *geometry, CConfig *config){
    unsigned short iMarker, iSpan, iVar;
    nSpanMax    = config->GetnSpanMaxAllZones();
    
    
    /*--- Initialize quantities for the average process for internal flow ---*/
    
    nSpanWiseSections = config->GetnSpanWiseSections();
    
    AverageVelocity                       = new su2double** [nMarker];
    AverageTurboVelocity                  = new su2double** [nMarker];
    OldAverageTurboVelocity               = new su2double** [nMarker];
    ExtAverageTurboVelocity               = new su2double** [nMarker];
    AverageFlux                           = new su2double** [nMarker];
    SpanTotalFlux                         = new su2double** [nMarker];
    AveragePressure                       = new su2double* [nMarker];
    OldAveragePressure                    = new su2double* [nMarker];
    RadialEquilibriumPressure             = new su2double* [nMarker];
    ExtAveragePressure                    = new su2double* [nMarker];
    AverageDensity                        = new su2double* [nMarker];
    OldAverageDensity                     = new su2double* [nMarker];
    ExtAverageDensity                     = new su2double* [nMarker];
    AverageNu                             = new su2double* [nMarker];
    AverageKine                           = new su2double* [nMarker];
    AverageOmega                          = new su2double* [nMarker];
    ExtAverageNu                          = new su2double* [nMarker];
    ExtAverageKine                        = new su2double* [nMarker];
    ExtAverageOmega                       = new su2double* [nMarker];
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        AverageVelocity[iMarker]                = new su2double* [nSpanWiseSections + 1];
        AverageTurboVelocity[iMarker]           = new su2double* [nSpanWiseSections + 1];
        OldAverageTurboVelocity[iMarker]        = new su2double* [nSpanWiseSections + 1];
        ExtAverageTurboVelocity[iMarker]        = new su2double* [nSpanWiseSections + 1];
        AverageFlux[iMarker]                    = new su2double* [nSpanWiseSections + 1];
        SpanTotalFlux[iMarker]                  = new su2double* [nSpanWiseSections + 1];
        AveragePressure[iMarker]                = new su2double [nSpanWiseSections + 1]();
        OldAveragePressure[iMarker]             = new su2double [nSpanWiseSections + 1]();
        RadialEquilibriumPressure[iMarker]      = new su2double [nSpanWiseSections + 1]();
        ExtAveragePressure[iMarker]             = new su2double [nSpanWiseSections + 1]();
        AverageDensity[iMarker]                 = new su2double [nSpanWiseSections + 1]();
        OldAverageDensity[iMarker]              = new su2double [nSpanWiseSections + 1]();
        ExtAverageDensity[iMarker]              = new su2double [nSpanWiseSections + 1]();
        AverageNu[iMarker]                      = new su2double [nSpanWiseSections + 1]();
        AverageKine[iMarker]                    = new su2double [nSpanWiseSections + 1]();
        AverageOmega[iMarker]                   = new su2double [nSpanWiseSections + 1]();
        ExtAverageNu[iMarker]                   = new su2double [nSpanWiseSections + 1]();
        ExtAverageKine[iMarker]                 = new su2double [nSpanWiseSections + 1]();
        ExtAverageOmega[iMarker]                = new su2double [nSpanWiseSections + 1]();
        
        for(iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++){
            AverageVelocity[iMarker][iSpan]           = new su2double [nDim]();
            AverageTurboVelocity[iMarker][iSpan]      = new su2double [nDim]();
            OldAverageTurboVelocity[iMarker][iSpan]   = new su2double [nDim]();
            ExtAverageTurboVelocity[iMarker][iSpan]   = new su2double [nDim]();
            AverageFlux[iMarker][iSpan]               = new su2double [nVar]();
            SpanTotalFlux[iMarker][iSpan]             = new su2double [nVar]();
        }
    }
    
    /*--- Initialize primitive quantities for turboperformace ---*/
    
    nMarkerTurboPerf = config->GetnMarker_TurboPerformance();
    
    DensityIn                     = new su2double*[nMarkerTurboPerf];
    PressureIn                    = new su2double*[nMarkerTurboPerf];
    TurboVelocityIn               = new su2double**[nMarkerTurboPerf];
    DensityOut                    = new su2double*[nMarkerTurboPerf];
    PressureOut                   = new su2double*[nMarkerTurboPerf];
    TurboVelocityOut              = new su2double**[nMarkerTurboPerf];
    KineIn                        = new su2double*[nMarkerTurboPerf];
    OmegaIn                       = new su2double*[nMarkerTurboPerf];
    NuIn                          = new su2double*[nMarkerTurboPerf];
    KineOut                       = new su2double*[nMarkerTurboPerf];
    OmegaOut                      = new su2double*[nMarkerTurboPerf];
    NuOut                         = new su2double*[nMarkerTurboPerf];
    
    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++) {
        DensityIn[iMarker]          = new su2double [nSpanMax + 1]();
        PressureIn[iMarker]         = new su2double [nSpanMax + 1]();
        TurboVelocityIn[iMarker]    = new su2double*[nSpanMax + 1]();
        DensityOut[iMarker]         = new su2double [nSpanMax + 1]();
        PressureOut[iMarker]        = new su2double [nSpanMax + 1]();
        TurboVelocityOut[iMarker]   = new su2double*[nSpanMax + 1]();
        KineIn[iMarker]             = new su2double [nSpanMax + 1]();
        OmegaIn[iMarker]            = new su2double [nSpanMax + 1]();
        NuIn[iMarker]               = new su2double [nSpanMax + 1]();
        KineOut[iMarker]            = new su2double [nSpanMax + 1]();
        OmegaOut[iMarker]           = new su2double [nSpanMax + 1]();
        NuOut[iMarker]              = new su2double [nSpanMax + 1]();
        
        for (iSpan = 0; iSpan < nSpanMax + 1; iSpan++) {
            TurboVelocityIn[iMarker][iSpan]    = new su2double[nDim]();
            TurboVelocityOut [iMarker][iSpan]  = new su2double[nDim]();
        }
    }
    
    
    /*--- Initialize quantities for NR BC ---*/
    
    if(config->GetBoolGiles()){
        
        CkInflow= new complex<su2double>** [nMarker];
        CkOutflow1= new complex<su2double>** [nMarker];
        CkOutflow2= new complex<su2double>** [nMarker];
        
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            CkInflow[iMarker]= new complex<su2double>*[nSpanWiseSections];
            CkOutflow1[iMarker]= new complex<su2double>*[nSpanWiseSections];
            CkOutflow2[iMarker]= new complex<su2double>*[nSpanWiseSections];
            
            for(iSpan = 0; iSpan < nSpanWiseSections; iSpan++) {
                CkInflow[iMarker][iSpan]= new complex<su2double>[2*geometry->GetnFreqSpanMax(INFLOW)+1];
                CkOutflow1[iMarker][iSpan]= new complex<su2double>[2*geometry->GetnFreqSpanMax(OUTFLOW)+1];
                CkOutflow2[iMarker][iSpan]= new complex<su2double>[2*geometry->GetnFreqSpanMax(OUTFLOW)+1];
                
                for (iVar = 0; iVar <2*geometry->GetnFreqSpanMax(INFLOW)+1; iVar++) {
                    CkInflow[iMarker][iSpan][iVar]= complex<su2double>(0.0,0.0);
                }
                
                for (iVar = 0; iVar <2*geometry->GetnFreqSpanMax(OUTFLOW)+1; iVar++) {
                    CkOutflow1[iMarker][iSpan][iVar]= complex<su2double>(0.0,0.0);
                    CkOutflow2[iMarker][iSpan][iVar]= complex<su2double>(0.0,0.0);
                }
            }
        }
    }
    
    
}

void CEulerSolver::Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config) {
    
    unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
    Buffer_Send_nPointTotal = 0;
    long iGlobalIndex, iGlobal;
    unsigned short iVar, iMarker, jMarker;
    long nDomain = 0, iDomain, jDomain;
    //bool ActDisk_Perimeter;
    bool rans = (config->GetKind_Turb_Model() != NONE) && (solver_container[TURB_SOL] != nullptr);
    
    unsigned short nPrimVar_ = nPrimVar;
    if (rans) nPrimVar_ += 2; // Add two extra variables for the turbulence.
    
#ifdef HAVE_MPI
    
    /*--- MPI status and request arrays for non-blocking communications ---*/
    
    SU2_MPI::Status status;
    
#endif
    
    /*--- Define buffer vector interior domain ---*/
    
    su2double *Buffer_Send_PrimVar = NULL;
    long      *Buffer_Send_Data    = NULL;
    
    unsigned long *nPointTotal_s = new unsigned long[size];
    unsigned long *nPointTotal_r = new unsigned long[size];
    su2double *iPrimVar = new su2double [nPrimVar_];
    
    unsigned long Buffer_Size_PrimVar = 0;
    unsigned long Buffer_Size_Data    = 0;
    
    unsigned long PointTotal_Counter = 0;
    
    /*--- Allocate the memory that we only need if we have MPI support ---*/
    
    su2double *Buffer_Receive_PrimVar = NULL;
    long      *Buffer_Receive_Data    = NULL;
    
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
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    //ActDisk_Perimeter = geometry->vertex[iMarker][iVertex]->GetActDisk_Perimeter();
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
                    //          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain()) && (!ActDisk_Perimeter)) {
                    if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
                        Buffer_Send_nPointTotal++;
                    }
                }
            }
        }
        
        /*--- Store the counts on a partition by partition basis. ---*/
        
        nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
        
        /*--- Total counts for allocating send buffers below ---*/
        
        Buffer_Size_PrimVar += nPointTotal_s[iDomain]*(nPrimVar_);
        Buffer_Size_Data += nPointTotal_s[iDomain]*(3);
        
    }
    
    /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
    
    Buffer_Send_PrimVar = new su2double[Buffer_Size_PrimVar];
    Buffer_Send_Data    = new long[Buffer_Size_Data];
    
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
            
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    
                    jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
                    //ActDisk_Perimeter = geometry->vertex[iMarker][iVertex]->GetActDisk_Perimeter();
                    
                    //          if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain()) && (!ActDisk_Perimeter)) {
                    if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
                        
                        for (iVar = 0; iVar < nPrimVar; iVar++) {
                            Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+iVar] = nodes->GetPrimitive(iPoint,iVar);
                        }
                        if (rans) {
                            Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+nPrimVar] = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
                            Buffer_Send_PrimVar[(nPrimVar_)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = 0.0;
                        }
                        
                        iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
                        jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
                        jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
                        
                        Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(0)]  = iGlobalIndex;
                        Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(1)] = jVertex;
                        Buffer_Send_Data[(3)*(PointTotal_Counter+iPointTotal)+(2)]  = jMarker;
                        
                        iPointTotal++;
                        
                    }
                    
                }
                
            }
            
        }
        
        /*--- Send the buffers with the geometrical information ---*/
        
        if (iDomain != rank) {
            
#ifdef HAVE_MPI
            
            /*--- Communicate the coordinates, global index, colors, and element
             date to iDomain with non-blocking sends. ---*/
            
            SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar_)],
                           nPointTotal_s[iDomain]*(nPrimVar_), MPI_DOUBLE, iDomain,
                           iDomain,  MPI_COMM_WORLD);
            
            SU2_MPI::Bsend(&Buffer_Send_Data[PointTotal_Counter*(3)],
                           nPointTotal_s[iDomain]*(3), MPI_LONG, iDomain,
                           iDomain+nDomain,  MPI_COMM_WORLD);
            
#endif
            
        }
        
        else {
            
            /*--- Allocate local memory for the local recv of the elements ---*/
            
            Buffer_Receive_PrimVar            = new su2double[nPointTotal_s[iDomain]*(nPrimVar_)];
            Buffer_Receive_Data               = new long[nPointTotal_s[iDomain]*(3)];
            
            for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar_); iter++)
                Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar_)+iter];
            
            for (iter = 0; iter < nPointTotal_s[iDomain]*(3); iter++)
                Buffer_Receive_Data[iter] = Buffer_Send_Data[PointTotal_Counter*(3)+iter];
            
            
            /*--- Recv the point data from ourselves (same procedure as above) ---*/
            
            for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
                
                for (iVar = 0; iVar < nPrimVar_; iVar++)
                    iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar_)+iVar];
                
                iGlobal       =  Buffer_Receive_Data[iPoint*(3)+(0)];
                iVertex      =  Buffer_Receive_Data[iPoint*(3)+(1)];
                iMarker      = Buffer_Receive_Data[iPoint*(3)+(2)];
                
                for (iVar = 0; iVar < nPrimVar_; iVar++)
                    SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
                
                SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
                
            }
            
            /*--- Delete memory for recv the point stuff ---*/
            
            delete [] Buffer_Receive_PrimVar;
            delete [] Buffer_Receive_Data;
            
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
            
            Buffer_Receive_PrimVar            = new su2double [nPointTotal_r[iDomain]*(nPrimVar_)];
            Buffer_Receive_Data               = new long [nPointTotal_r[iDomain]*(3)];
            
            /*--- Receive the buffers with the coords, global index, and colors ---*/
            
            SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar_) , MPI_DOUBLE,
                          iDomain, rank, MPI_COMM_WORLD, &status);
            
            SU2_MPI::Recv(Buffer_Receive_Data, nPointTotal_r[iDomain]*(3) , MPI_LONG,
                          iDomain, rank+nDomain, MPI_COMM_WORLD, &status);
            
            /*--- Loop over all of the points that we have recv'd and store the
             coords, global index vertex and markers ---*/
            
            for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
                
                iGlobal      = Buffer_Receive_Data[iPoint*(3)+(0)];
                iVertex      = Buffer_Receive_Data[iPoint*(3)+(1)];
                iMarker      = Buffer_Receive_Data[iPoint*(3)+(2)];
                
                for (iVar = 0; iVar < nPrimVar_; iVar++)
                    iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar_)+iVar];
                
                for (iVar = 0; iVar < nPrimVar_; iVar++) {
                    SetDonorPrimVar(iMarker, iVertex, iVar, iPrimVar[iVar]);
                }
                
                SetDonorGlobalIndex(iMarker, iVertex, iGlobal);
                
            }
            
            /*--- Delete memory for recv the point stuff ---*/
            
            delete [] Buffer_Receive_PrimVar;
            delete [] Buffer_Receive_Data;
            
#endif
            
        }
        
    }
    
    /*--- Wait for the non-blocking sends to complete. ---*/
    
#ifdef HAVE_MPI
    
    SU2_MPI::Barrier(MPI_COMM_WORLD);
    
#endif
    
    /*--- Free all of the memory used for communicating points and elements ---*/
    
    delete[] Buffer_Send_PrimVar;
    delete[] Buffer_Send_Data;
    
    /*--- Release all of the temporary memory ---*/
    
    delete [] nPointTotal_s;
    delete [] nPointTotal_r;
    delete [] iPrimVar;
    
}

void CEulerSolver::Set_MPI_Nearfield(CGeometry *geometry, CConfig *config) {
    
    unsigned long iter,  iPoint, iVertex, jVertex, iPointTotal,
    Buffer_Send_nPointTotal = 0;
    long iGlobalIndex, iGlobal;
    unsigned short iVar, iMarker, jMarker;
    long nDomain = 0, iDomain, jDomain;
    
#ifdef HAVE_MPI
    
    /*--- MPI status and request arrays for non-blocking communications ---*/
    
    SU2_MPI::Status status, status_;
    
    
#endif
    
    /*--- Define buffer vector interior domain ---*/
    
    su2double        *Buffer_Send_PrimVar          = NULL;
    
    unsigned long *nPointTotal_s = new unsigned long[size];
    unsigned long *nPointTotal_r = new unsigned long[size];
    su2double     *iPrimVar      = new su2double [nPrimVar];
    
    unsigned long Buffer_Size_PrimVar = 0;
    unsigned long PointTotal_Counter = 0;
    
    /*--- Allocate the memory that we only need if we have MPI support ---*/
    
    su2double* Buffer_Receive_PrimVar = NULL;
    
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
            if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
                    if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
                        Buffer_Send_nPointTotal++;
                    }
                }
            }
        }
        
        /*--- Store the counts on a partition by partition basis. ---*/
        
        nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
        
        /*--- Total counts for allocating send buffers below ---*/
        
        Buffer_Size_PrimVar += nPointTotal_s[iDomain]*(nPrimVar+3);
        
    }
    
    /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
    
    Buffer_Send_PrimVar = new su2double[Buffer_Size_PrimVar];
    
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
            
            if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
                
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    jDomain = geometry->vertex[iMarker][iVertex]->GetDonorProcessor();
                    if ((iDomain == jDomain) && (geometry->node[iPoint]->GetDomain())) {
                        iGlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
                        jVertex = geometry->vertex[iMarker][iVertex]->GetDonorVertex();
                        jMarker = geometry->vertex[iMarker][iVertex]->GetDonorMarker();
                        for (iVar = 0; iVar < nPrimVar; iVar++) {
                            Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+iVar] = nodes->GetPrimitive(iPoint,iVar);
                        }
                        Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+0)] = su2double(iGlobalIndex);
                        Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+1)] = su2double(jVertex);
                        Buffer_Send_PrimVar[(nPrimVar+3)*(PointTotal_Counter+iPointTotal)+(nPrimVar+2)] = su2double(jMarker);
                        
                        iPointTotal++;
                        
                    }
                    
                }
                
            }
            
        }
        
        /*--- Send the buffers with the geometrical information ---*/
        
        if (iDomain != rank) {
            
#ifdef HAVE_MPI
            
            /*--- Communicate the coordinates, global index, colors, and element
             date to iDomain with non-blocking sends. ---*/
            
            SU2_MPI::Bsend(&Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)],
                           nPointTotal_s[iDomain]*(nPrimVar+3), MPI_DOUBLE, iDomain,
                           iDomain,  MPI_COMM_WORLD);
            
#endif
            
        }
        
        else {
            
            /*--- Allocate local memory for the local recv of the elements ---*/
            
            Buffer_Receive_PrimVar = new su2double[nPointTotal_s[iDomain]*(nPrimVar+3)];
            
            for (iter = 0; iter < nPointTotal_s[iDomain]*(nPrimVar+3); iter++)
                Buffer_Receive_PrimVar[iter] = Buffer_Send_PrimVar[PointTotal_Counter*(nPrimVar+3)+iter];
            
            /*--- Recv the point data from ourselves (same procedure as above) ---*/
            
            for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
                
                iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
                iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
                iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
                for (iVar = 0; iVar < nPrimVar; iVar++)
                    iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
                
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
            
            Buffer_Receive_PrimVar = new su2double [nPointTotal_r[iDomain]*(nPrimVar+3)];
            
            /*--- Receive the buffers with the coords, global index, and colors ---*/
            
            SU2_MPI::Recv(Buffer_Receive_PrimVar, nPointTotal_r[iDomain]*(nPrimVar+3) , MPI_DOUBLE,
                          iDomain, rank, MPI_COMM_WORLD, &status_);
            
            /*--- Loop over all of the points that we have recv'd and store the
             coords, global index vertex and markers ---*/
            
            for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {
                
                iGlobal      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+0)]);
                iVertex      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+1)]);
                iMarker      = SU2_TYPE::Int(Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+(nPrimVar+2)]);
                for (iVar = 0; iVar < nPrimVar; iVar++)
                    iPrimVar[iVar] = Buffer_Receive_PrimVar[iPoint*(nPrimVar+3)+iVar];
                
                if (iVertex < 0.0) cout <<" Negative iVertex (receive)" << endl;
                if (iMarker < 0.0) cout <<" Negative iMarker (receive)" << endl;
                
                if (iMarker > nMarker) cout << "ERROR" <<  endl;
                if (iVertex > geometry->nVertex[iMarker]) cout << "ERROR" <<  endl;
                
                for (iVar = 0; iVar < nPrimVar; iVar++)
                    SetDonorPrimVar(iMarker, iVertex, iVar,  iPrimVar[iVar]);
                
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
    
}

void CEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {
    
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
    Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0, TgammaR = 0.0, Heat_Flux_Ref = 0.0;
    
    unsigned short iDim;
    
    /*--- Local variables ---*/
    
    su2double Alpha         = config->GetAoA()*PI_NUMBER/180.0;
    su2double Beta          = config->GetAoS()*PI_NUMBER/180.0;
    su2double Mach          = config->GetMach();
    su2double Reynolds      = config->GetReynolds();
    bool unsteady           = (config->GetTime_Marching() != NO);
    bool viscous            = config->GetViscous();
    bool gravity            = config->GetGravityForce();
    bool turbulent          = (config->GetKind_Turb_Model() != NONE);
    bool tkeNeeded          = (turbulent && (config->GetKind_Turb_Model() == SST || config->GetKind_Turb_Model() == SST_SUST));
    bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
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
    
    /*--- Compute the Free Stream velocity, using the Mach number ---*/
    
    Pressure_FreeStream = config->GetPressure_FreeStream();
    Density_FreeStream  = config->GetDensity_FreeStream();
    Temperature_FreeStream = config->GetTemperature_FreeStream();
    
    CFluidModel* auxFluidModel = nullptr;
    
    switch (config->GetKind_FluidModel()) {
            
        case STANDARD_AIR:
            
            switch (config->GetSystemMeasurements()) {
                case SI: config->SetGas_Constant(287.058); break;
                case US: config->SetGas_Constant(1716.49); break;
            }
            
            auxFluidModel = new CIdealGas(1.4, config->GetGas_Constant());
            
            if (free_stream_temp && aeroelastic) {
                Temperature_FreeStream = TgammaR / (config->GetGas_Constant()*1.4);
                config->SetTemperature_FreeStream(Temperature_FreeStream);
            }
            break;
            
        case IDEAL_GAS:
            
            auxFluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
            break;
            
        case VW_GAS:
            
            auxFluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                                config->GetPressure_Critical(), config->GetTemperature_Critical());
            break;
            
        case PR_GAS:
            
            auxFluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                              config->GetTemperature_Critical(), config->GetAcentric_Factor());
            break;
            
        default:
            SU2_MPI::Error("Unknown fluid model.", CURRENT_FUNCTION);
            break;
    }
    
    if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
    }
    else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
    }
    
    Mach2Vel_FreeStream = auxFluidModel->GetSoundSpeed();
    
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
            
            if (dynamic_grid) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
            else Velocity_Reynolds = ModVel_FreeStream;
            
            /*--- For viscous flows, pressure will be computed from a density
             that is found from the Reynolds number. The viscosity is computed
             from the dimensional version of Sutherland's law or the constant
             viscosity, depending on the input option.---*/
            
            auxFluidModel->SetLaminarViscosityModel(config);
            
            Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
            config->SetViscosity_FreeStream(Viscosity_FreeStream);
            
            Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
            config->SetDensity_FreeStream(Density_FreeStream);
            auxFluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
            Pressure_FreeStream = auxFluidModel->GetPressure();
            config->SetPressure_FreeStream(Pressure_FreeStream);
            Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
            
        }
        
        /*--- Thermodynamics quantities based initialization ---*/
        
        else {
            
            auxFluidModel->SetLaminarViscosityModel(config);
            Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
            config->SetViscosity_FreeStream(Viscosity_FreeStream);
            Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
            
        }
        
        /*--- Turbulence kinetic energy ---*/
        
        Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
        
    }
    else {
        
        /*--- For inviscid flow, energy is calculated from the specified
         FreeStream quantities using the proper gas law. ---*/
        
        Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
        
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
    
    /*--- Auxilary (dimensional) FluidModel no longer needed. ---*/
    delete auxFluidModel;
    
    /*--- Set viscosity ND constants before defining the visc. model of the fluid models. ---*/
    
    if (viscous) {
        /*--- Constant viscosity model. ---*/
        config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);
        
        /*--- Sutherland's model. ---*/
        config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
        config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
        config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());
        
        /*--- Constant thermal conductivity model. ---*/
        config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);
    }
    
    /*--- Create one final fluid model object per OpenMP thread to be able to use them in parallel.
     *    GetFluidModel() should be used to automatically access the "right" object of each thread. ---*/
    
    assert(FluidModel.empty() && "Potential memory leak!");
    FluidModel.resize(omp_get_max_threads());
    
    SU2_OMP_PARALLEL
    {
        const int thread = omp_get_thread_num();
        
        switch (config->GetKind_FluidModel()) {
                
            case STANDARD_AIR:
                FluidModel[thread] = new CIdealGas(1.4, Gas_ConstantND);
                break;
                
            case IDEAL_GAS:
                FluidModel[thread] = new CIdealGas(Gamma, Gas_ConstantND);
                break;
                
            case VW_GAS:
                FluidModel[thread] = new CVanDerWaalsGas(Gamma, Gas_ConstantND,
                                                         config->GetPressure_Critical() / config->GetPressure_Ref(),
                                                         config->GetTemperature_Critical() / config->GetTemperature_Ref());
                break;
                
            case PR_GAS:
                FluidModel[thread] = new CPengRobinson(Gamma, Gas_ConstantND,
                                                       config->GetPressure_Critical() / config->GetPressure_Ref(),
                                                       config->GetTemperature_Critical() / config->GetTemperature_Ref(),
                                                       config->GetAcentric_Factor());
                break;
        }
        
        GetFluidModel()->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
        if (viscous) {
            GetFluidModel()->SetLaminarViscosityModel(config);
            GetFluidModel()->SetThermalConductivityModel(config);
        }
        
    } // end SU2_OMP_PARALLEL
    
    Energy_FreeStreamND = GetFluidModel()->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
    
    if (tkeNeeded) Energy_FreeStreamND += Tke_FreeStreamND;
    
    config->SetEnergy_FreeStreamND(Energy_FreeStreamND);
    
    Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);
    
    Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
    Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);
    
    /*--- Write output to the console if this is the master node and first domain ---*/
    
    if ((rank == MASTER_NODE) && (MGLevel == MESH_0)) {
        
        cout.precision(6);
        
        if (viscous) {
            cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
            cout << "based on the free-stream temperature and a density computed" << endl;
            cout << "from the Reynolds number." << endl;
        } else {
            cout << "Inviscid flow: Computing density based on free-stream" << endl;
            cout << "temperature and pressure using the ideal gas law." << endl;
        }
        
        if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
        else cout << "Force coefficients computed using free-stream values." << endl;
        
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
        
        if (viscous) {
            
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
                    
            }
        } else {
            ModelTable << "-" << "-";
        }
        
        if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
        NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
        NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << Gamma;
        Unit.str("");
        
        switch(config->GetKind_FluidModel()){
            case STANDARD_AIR:
                ModelTable << "STANDARD_AIR";
                break;
            case IDEAL_GAS:
                ModelTable << "IDEAL_GAS";
                break;
            case VW_GAS:
                ModelTable << "VW_GAS";
                break;
            case PR_GAS:
                ModelTable << "PR_GAS";
                break;
        }
        
        if (config->GetKind_FluidModel() == VW_GAS || config->GetKind_FluidModel() == PR_GAS){
            NonDimTable << "Critical Pressure" << config->GetPressure_Critical() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_Critical() /config->GetPressure_Ref();
            Unit.str("");
            Unit << "K";
            NonDimTable << "Critical Temperature" << config->GetTemperature_Critical() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_Critical() /config->GetTemperature_Ref();
            Unit.str("");
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
        NonDimTable << "Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
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

void CEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {
    
    const bool restart = (config->GetRestart() || config->GetRestart_Flow());
    const bool rans = (config->GetKind_Turb_Model() != NONE);
    const bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                            (config->GetTime_Marching() == DT_STEPPING_2ND));
    const bool SubsonicEngine = config->GetSubsonicEngine();
    
    /*--- Start OpenMP parallel region. ---*/
    
    SU2_OMP_PARALLEL
    {
        
        unsigned long iPoint;
        unsigned short iMesh, iDim;
        su2double X0[3] = {0.0,0.0,0.0}, X1[3] = {0.0,0.0,0.0}, X2[3] = {0.0,0.0,0.0},
        X1_X0[3] = {0.0,0.0,0.0}, X2_X0[3] = {0.0,0.0,0.0}, X2_X1[3] = {0.0,0.0,0.0},
        CP[3] = {0.0,0.0,0.0}, Distance, DotCheck, Radius;
        
        /*--- Check if a verification solution is to be computed. ---*/
        if ((VerificationSolution) && (TimeIter == 0) && !restart) {
            
            /*--- Loop over the multigrid levels. ---*/
            for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
                
                /*--- Loop over all grid points. ---*/
                SU2_OMP_FOR_STAT(omp_chunk_size)
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
        
        /*--- Set subsonic initial condition for engine intakes ---*/
        
        if (SubsonicEngine) {
            
            /*--- Set initial boundary condition at iteration 0 ---*/
            
            if ((TimeIter == 0) && (!restart)) {
                
                su2double Velocity_Cyl[3] = {0.0, 0.0, 0.0}, Velocity_CylND[3] = {0.0, 0.0, 0.0}, Viscosity_Cyl,
                Density_Cyl, Density_CylND, Pressure_CylND, ModVel_Cyl, ModVel_CylND, Energy_CylND,
                T_ref = 0.0, S = 0.0, Mu_ref = 0.0;
                const su2double *Coord, *SubsonicEngine_Cyl, *SubsonicEngine_Values;
                
                SubsonicEngine_Values = config->GetSubsonicEngine_Values();
                su2double Mach_Cyl        = SubsonicEngine_Values[0];
                su2double Alpha_Cyl       = SubsonicEngine_Values[1];
                su2double Beta_Cyl        = SubsonicEngine_Values[2];
                su2double Pressure_Cyl    = SubsonicEngine_Values[3];
                su2double Temperature_Cyl = SubsonicEngine_Values[4];
                
                su2double Alpha = Alpha_Cyl*PI_NUMBER/180.0;
                su2double Beta  = Beta_Cyl*PI_NUMBER/180.0;
                
                su2double Gamma_Minus_One = Gamma - 1.0;
                su2double Gas_Constant = config->GetGas_Constant();
                
                su2double Mach2Vel_Cyl = sqrt(Gamma*Gas_Constant*Temperature_Cyl);
                
                for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
                    
                    SU2_OMP_FOR_STAT(omp_chunk_size)
                    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
                        
                        Velocity_Cyl[0] = cos(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;
                        Velocity_Cyl[1] = sin(Beta)*Mach_Cyl*Mach2Vel_Cyl;
                        Velocity_Cyl[2] = sin(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;
                        
                        ModVel_Cyl = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            ModVel_Cyl += Velocity_Cyl[iDim]*Velocity_Cyl[iDim];
                        }
                        ModVel_Cyl = sqrt(ModVel_Cyl);
                        
                        if (config->GetViscous()) {
                            if (config->GetSystemMeasurements() == SI) { T_ref = 273.15; S = 110.4; Mu_ref = 1.716E-5; }
                            if (config->GetSystemMeasurements() == US) {
                                T_ref = (273.15 - 273.15) * 1.8 + 491.67;
                                S = (110.4 - 273.15) * 1.8 + 491.67;
                                Mu_ref = 1.716E-5/47.88025898;
                            }
                            Viscosity_Cyl = Mu_ref*(pow(Temperature_Cyl/T_ref, 1.5) * (T_ref+S)/(Temperature_Cyl+S));
                            Density_Cyl   = config->GetReynolds()*Viscosity_Cyl/(ModVel_Cyl*config->GetLength_Reynolds());
                            Pressure_Cyl  = Density_Cyl*Gas_Constant*Temperature_Cyl;
                        }
                        else {
                            Density_Cyl = Pressure_Cyl/(Gas_Constant*Temperature_Cyl);
                        }
                        
                        Density_CylND  = Density_Cyl/config->GetDensity_Ref();
                        Pressure_CylND = Pressure_Cyl/config->GetPressure_Ref();
                        
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Velocity_CylND[iDim] = Velocity_Cyl[iDim]/config->GetVelocity_Ref();
                        }
                        
                        ModVel_CylND = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            ModVel_CylND += Velocity_CylND[iDim]*Velocity_CylND[iDim];
                        }
                        ModVel_CylND = sqrt(ModVel_CylND);
                        
                        Energy_CylND = Pressure_CylND/(Density_CylND*Gamma_Minus_One)+0.5*ModVel_CylND*ModVel_CylND;
                        
                        Coord = geometry[iMesh]->node[iPoint]->GetCoord();
                        
                        SubsonicEngine_Cyl = config->GetSubsonicEngine_Cyl();
                        
                        X0[0] = Coord[0];               X0[1] = Coord[1];               X0[2] = Coord[2];
                        X1[0] = SubsonicEngine_Cyl[0];  X1[1] = SubsonicEngine_Cyl[1];  X1[2] = SubsonicEngine_Cyl[2];
                        X2[0] = SubsonicEngine_Cyl[3];  X2[1] = SubsonicEngine_Cyl[4];  X2[2] = SubsonicEngine_Cyl[5];
                        Radius = SubsonicEngine_Cyl[6];
                        
                        for (iDim = 0; iDim < nDim; iDim++) {
                            X2_X1[iDim]= X1[iDim] - X2[iDim];
                            X1_X0[iDim]= X0[iDim] - X1[iDim];
                            X2_X0[iDim]= X0[iDim] - X2[iDim];
                        }
                        
                        CP[0] = (X2_X1[1]*X1_X0[2] - X2_X1[2]*X1_X0[1]);
                        CP[1] = (X2_X1[2]*X1_X0[0] - X2_X1[0]*X1_X0[2]);
                        CP[2] = (X2_X1[0]*X1_X0[1] - X2_X1[1]*X1_X0[0]);
                        
                        Distance = sqrt((CP[0]*CP[0]+CP[1]*CP[1]+CP[2]*CP[2])/(X2_X1[0]*X2_X1[0]+X2_X1[1]*X2_X1[1]+X2_X1[2]*X2_X1[2]));
                        
                        DotCheck = -(X1_X0[0]*X2_X1[0]+X1_X0[1]*X2_X1[1]+X1_X0[2]*X2_X1[2]);
                        if (DotCheck < 0.0) Distance = sqrt(X1_X0[0]*X1_X0[0]+X1_X0[1]*X1_X0[1]+X1_X0[2]*X1_X0[2]);
                        
                        DotCheck = (X2_X0[0]*X2_X1[0]+X2_X0[1]*X2_X1[1]+X2_X0[2]*X2_X1[2]);
                        if (DotCheck < 0.0) Distance = sqrt(X2_X0[0]*X2_X0[0]+X2_X0[1]*X2_X0[1]+X2_X0[2]*X2_X0[2]);
                        
                        if (Distance < Radius) {
                            
                            solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint, 0, Density_CylND);
                            for (iDim = 0; iDim < nDim; iDim++)
                                solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint, iDim+1, Density_CylND*Velocity_CylND[iDim]);
                            solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint, nVar-1, Density_CylND*Energy_CylND);
                            
                            solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution_Old(iPoint, 0, Density_CylND);
                            for (iDim = 0; iDim < nDim; iDim++)
                                solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution_Old(iPoint, iDim+1, Density_CylND*Velocity_CylND[iDim]);
                            solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution_Old(iPoint, nVar-1, Density_CylND*Energy_CylND);
                            
                        }
                        
                    }
                    
                    /*--- Set the MPI communication ---*/
                    
                    SU2_OMP_MASTER
                    {
                        solver_container[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
                        solver_container[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
                        
                        solver_container[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION_OLD);
                        solver_container[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION_OLD);
                    }
                    SU2_OMP_BARRIER
                    
                }
                
            }
            
        }
        
        /*--- Make sure that the solution is well initialized for unsteady
         calculations with dual time-stepping (load additional restarts for 2nd-order). ---*/
        
        if (dual_time && ((TimeIter == 0) || (restart && (TimeIter == config->GetRestart_Iter()))) ) {
            
            /*--- Push back the initial condition to previous solution containers
             for a 1st-order restart or when simply intitializing to freestream. ---*/
            
            for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
                solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
                solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n1();
                if (rans) {
                    solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
                    solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1();
                }
            }
            
            if (restart && (TimeIter == config->GetRestart_Iter()) && (config->GetTime_Marching() == DT_STEPPING_2ND)) {
                
                /*--- Load an additional restart file for a 2nd-order restart ---*/
                
                solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, config->GetRestart_Iter()-1, true);
                
                /*--- Load an additional restart file for the turbulence model ---*/
                if (rans)
                    solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, config->GetRestart_Iter()-1, false);
                
                /*--- Push back this new solution to time level N. ---*/
                
                for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
                    solver_container[iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
                    if (rans) {
                        solver_container[iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
                    }
                }
            }
        }
        
    } // end SU2_OMP_PARALLEL
    
}

void CEulerSolver::CommonPreprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                       unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
    
    bool cont_adjoint     = config->GetContinuous_Adjoint();
    bool disc_adjoint     = config->GetDiscrete_Adjoint();
    bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
    bool center_jst       = (config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0);
    bool engine           = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
    bool actuator_disk    = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
    bool nearfield        = (config->GetnMarker_NearFieldBound() != 0);
    bool fixed_cl         = config->GetFixed_CL_Mode();
    unsigned short kind_row_dissipation = config->GetKind_RoeLowDiss();
    bool roe_low_dissipation  = (kind_row_dissipation != NO_ROELOWDISS) &&
    (config->GetKind_Upwind_Flow() == ROE ||
     config->GetKind_Upwind_Flow() == SLAU ||
     config->GetKind_Upwind_Flow() == SLAU2);
    bool outlet           = ((config->GetnMarker_Outlet() != 0));

    /*--- Update the angle of attack at the far-field for fixed CL calculations (only direct problem). ---*/
    
    if (fixed_cl && !disc_adjoint && !cont_adjoint) {
        SU2_OMP_MASTER
        SetFarfield_AoA(geometry, solver_container, config, iMesh, Output);
        SU2_OMP_BARRIER
    }
    
    /*--- Set the primitive variables ---*/
    
    SU2_OMP_MASTER
    ErrorCounter = 0;
    SU2_OMP_BARRIER
    
    SU2_OMP_ATOMIC
    ErrorCounter += SetPrimitive_Variables(solver_container, config, Output);
    
    if ((iMesh == MESH_0) && (config->GetComm_Level() == COMM_FULL)) {
        SU2_OMP_BARRIER
        SU2_OMP_MASTER
        {
            unsigned long tmp = ErrorCounter;
            SU2_MPI::Allreduce(&tmp, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            config->SetNonphysical_Points(ErrorCounter);
        }
        SU2_OMP_BARRIER
    }
    
    /*--- Compute the engine properties ---*/
    
    if (engine) {
        SU2_OMP_MASTER
        GetPower_Properties(geometry, config, iMesh, Output);
        SU2_OMP_BARRIER
    }
    
    /*--- Compute properties needed for mass flow BCs. ---*/

    if (outlet) GetOutlet_Properties(geometry, config, iMesh, Output);

    /*--- Compute the actuator disk properties and distortion levels ---*/
    
    if (actuator_disk) {
        SU2_OMP_MASTER
        {
            Set_MPI_ActDisk(solver_container, geometry, config);
            GetPower_Properties(geometry, config, iMesh, Output);
            SetActDisk_BCThrust(geometry, solver_container, config, iMesh, Output);
        }
        SU2_OMP_BARRIER
    }
    
    /*--- Compute NearField MPI ---*/
    
    if (nearfield) {
        SU2_OMP_MASTER
        Set_MPI_Nearfield(geometry, config);
        SU2_OMP_BARRIER
    }
    
    /*--- Artificial dissipation ---*/
    
    if (center && !Output) {
        SetMax_Eigenvalue(geometry, config);
        if (center_jst)
            SetUndivided_Laplacian_And_Centered_Dissipation_Sensor(geometry, config);
    }
    
    /*--- Roe Low Dissipation Sensor ---*/
    
    if (roe_low_dissipation) {
        SetRoe_Dissipation(geometry, config);
        if (kind_row_dissipation == FD_DUCROS || kind_row_dissipation == NTS_DUCROS){
            SetUpwind_Ducros_Sensor(geometry, config);
        }
    }
    
    /*--- Initialize the Jacobian matrix and residual, not needed for the reducer strategy
     *    as we set blocks (including diagonal ones) and completely overwrite. ---*/
    
    if(!ReducerStrategy && !Output) {
        LinSysRes.SetValZero();
        if (implicit && !disc_adjoint) Jacobian.SetValZero();
        else {SU2_OMP_BARRIER} // because of "nowait" in LinSysRes
    }
    
}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                 unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
    
    unsigned long InnerIter = config->GetInnerIter();
    bool cont_adjoint     = config->GetContinuous_Adjoint();
    bool muscl            = (config->GetMUSCL_Flow() || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
    bool limiter          = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
    bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
    bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
    
    /*--- Common preprocessing steps. ---*/
    
    CommonPreprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);
    
    /*--- Upwind second order reconstruction ---*/
    
    if ((muscl && !center) && (iMesh == MESH_0) && !Output) {
        
        /*--- Gradient computation for MUSCL reconstruction. ---*/
        
        switch (config->GetKind_Gradient_Method_Recon()) {
            case GREEN_GAUSS:
                SetPrimitive_Gradient_GG(geometry, config, true); break;
            case LEAST_SQUARES:
            case WEIGHTED_LEAST_SQUARES:
                SetPrimitive_Gradient_LS(geometry, config, true); break;
            default: break;
        }
        
        /*--- Limiter computation ---*/
        
        if (limiter && (iMesh == MESH_0) && !Output && !van_albada)
            SetPrimitive_Limiter(geometry, config);
    }
    
}

unsigned long CEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
    
    /*--- Number of non-physical points, local to the thread, needs
     *    further reduction if function is called in parallel ---*/
    unsigned long nonPhysicalPoints = 0;
    
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {
        
        /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/
        
        bool physical = nodes->SetPrimVar(iPoint, GetFluidModel());
        nodes->SetSecondaryVar(iPoint, GetFluidModel());
        
        /* Check for non-realizable states for reporting. */
        
        if (!physical) nonPhysicalPoints++;
    }
    
    return nonPhysicalPoints;
}

void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {
    
    const bool viscous       = config->GetViscous();
    const bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    const bool time_stepping = (config->GetTime_Marching() == TIME_STEPPING);
    const bool dual_time     = (config->GetTime_Marching() == DT_STEPPING_1ST) ||
    (config->GetTime_Marching() == DT_STEPPING_2ND);
    const su2double K_v = 0.25;
    
    /*--- Init thread-shared variables to compute min/max values.
     *    Critical sections are used for this instead of reduction
     *    clauses for compatibility with OpenMP 2.0 (Windows...). ---*/
    
    SU2_OMP_MASTER
    {
        Min_Delta_Time = 1e30;
        Max_Delta_Time = 0.0;
        Global_Delta_UnstTimeND = 1e30;
    }
    SU2_OMP_BARRIER
    
    const su2double *Normal = nullptr;
    su2double Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Local_Delta_Time_Visc;
    su2double Mean_LaminarVisc, Mean_EddyVisc, Mean_Density, Lambda_1, Lambda_2;
    unsigned long iEdge, iVertex, iPoint, jPoint;
    unsigned short iDim, iMarker;
    
    /*--- Loop domain points. ---*/
    
    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
        
        auto node_i = geometry->node[iPoint];
        
        /*--- Set maximum eigenvalues to zero. ---*/
        
        nodes->SetMax_Lambda_Inv(iPoint,0.0);
        
        if (viscous)
            nodes->SetMax_Lambda_Visc(iPoint,0.0);
        
        /*--- Loop over the neighbors of point i. ---*/
        
        for (unsigned short iNeigh = 0; iNeigh < node_i->GetnPoint(); ++iNeigh)
        {
            jPoint = node_i->GetPoint(iNeigh);
            auto node_j = geometry->node[jPoint];
            
            iEdge = node_i->GetEdge(iNeigh);
            Normal = geometry->edges->GetNormal(iEdge);
            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += pow(Normal[iDim],2); Area = sqrt(Area);
            
            /*--- Mean Values ---*/
            
            Mean_ProjVel = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
            Mean_SoundSpeed = 0.5 * (nodes->GetSoundSpeed(iPoint) + nodes->GetSoundSpeed(jPoint)) * Area;
            
            /*--- Adjustment for grid movement ---*/
            
            if (dynamic_grid) {
                const su2double *GridVel_i = node_i->GetGridVel();
                const su2double *GridVel_j = node_j->GetGridVel();
                
                for (iDim = 0; iDim < nDim; iDim++)
                    Mean_ProjVel -= 0.5 * (GridVel_i[iDim] + GridVel_j[iDim]) * Normal[iDim];
            }
            
            /*--- Inviscid contribution ---*/
            
            Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
            nodes->AddMax_Lambda_Inv(iPoint,Lambda);
            
            /*--- Viscous contribution ---*/
            
            if (!viscous) continue;
            
            Mean_LaminarVisc = 0.5*(nodes->GetLaminarViscosity(iPoint) + nodes->GetLaminarViscosity(jPoint));
            Mean_EddyVisc    = 0.5*(nodes->GetEddyViscosity(iPoint) + nodes->GetEddyViscosity(jPoint));
            Mean_Density     = 0.5*(nodes->GetDensity(iPoint) + nodes->GetDensity(jPoint));
            
            Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
            //TODO (REAL_GAS) removing Gamma it cannot work with FLUIDPROP
            Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
            
            Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
            nodes->AddMax_Lambda_Visc(iPoint, Lambda);
        }
        
    }
    
    /*--- Loop boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
            (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
            
            SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                
                /*--- Point identification, Normal vector and area ---*/
                
                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                
                if (!geometry->node[iPoint]->GetDomain()) continue;
                
                Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
                
                /*--- Mean Values ---*/
                
                Mean_ProjVel = nodes->GetProjVel(iPoint,Normal);
                Mean_SoundSpeed = nodes->GetSoundSpeed(iPoint) * Area;
                
                /*--- Adjustment for grid movement ---*/
                
                if (dynamic_grid) {
                    const su2double *GridVel = geometry->node[iPoint]->GetGridVel();
                    
                    for (iDim = 0; iDim < nDim; iDim++)
                        Mean_ProjVel -= GridVel[iDim]*Normal[iDim];
                }
                
                /*--- Inviscid contribution ---*/
                
                Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
                nodes->AddMax_Lambda_Inv(iPoint,Lambda);
                
                /*--- Viscous contribution ---*/
                
                if (!viscous) continue;
                
                Mean_LaminarVisc = nodes->GetLaminarViscosity(iPoint);
                Mean_EddyVisc    = nodes->GetEddyViscosity(iPoint);
                Mean_Density     = nodes->GetDensity(iPoint);
                
                Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
                Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
                Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
                
                nodes->AddMax_Lambda_Visc(iPoint, Lambda);
                
            }
        }
    }
    
    /*--- Each element uses their own speed, steady state simulation. ---*/
    {
        /*--- Thread-local variables for min/max reduction. ---*/
        su2double minDt = 1e30, maxDt = 0.0;
        
        SU2_OMP(for schedule(static,omp_chunk_size) nowait)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            Vol = geometry->node[iPoint]->GetVolume();
            
            if (Vol != 0.0) {
                Local_Delta_Time = nodes->GetLocalCFL(iPoint)*Vol / nodes->GetMax_Lambda_Inv(iPoint);
                
                if(viscous) {
                    Local_Delta_Time_Visc = nodes->GetLocalCFL(iPoint)*K_v*Vol*Vol/ nodes->GetMax_Lambda_Visc(iPoint);
                    Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
                }
                
                minDt = min(minDt, Local_Delta_Time);
                maxDt = max(maxDt, Local_Delta_Time);
                
                nodes->SetDelta_Time(iPoint, min(Local_Delta_Time, config->GetMax_DeltaTime()));
            }
            else {
                nodes->SetDelta_Time(iPoint,0.0);
            }
        }
        /*--- Min/max over threads. ---*/
        SU2_OMP_CRITICAL
        {
            Min_Delta_Time = min(Min_Delta_Time, minDt);
            Max_Delta_Time = max(Max_Delta_Time, maxDt);
            Global_Delta_Time = Min_Delta_Time;
        }
        SU2_OMP_BARRIER
    }
    
    /*--- Compute the min/max dt (in parallel, now over mpi ranks). ---*/
    
    SU2_OMP_MASTER
    if (config->GetComm_Level() == COMM_FULL) {
        su2double rbuf_time;
        SU2_MPI::Allreduce(&Min_Delta_Time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        Min_Delta_Time = rbuf_time;
        
        SU2_MPI::Allreduce(&Max_Delta_Time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        Max_Delta_Time = rbuf_time;
    }
    SU2_OMP_BARRIER
    
    /*--- For exact time solution use the minimum delta time of the whole mesh. ---*/
    if (time_stepping) {
        
        /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step,
         *    otherwise it computes the time step based on the unsteady CFL. ---*/
        
        SU2_OMP_MASTER
        {
            if (config->GetUnst_CFL() == 0.0) {
                Global_Delta_Time = config->GetDelta_UnstTime();
            }
            else {
                Global_Delta_Time = Min_Delta_Time;
            }
            Max_Delta_Time = Global_Delta_Time;
            
            config->SetDelta_UnstTimeND(Global_Delta_Time);
        }
        SU2_OMP_BARRIER
        
        /*--- Sets the regular CFL equal to the unsteady CFL. ---*/
        
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            nodes->SetLocalCFL(iPoint, config->GetUnst_CFL());
            nodes->SetDelta_Time(iPoint, Global_Delta_Time);
        }
        
    }
    
    /*--- Recompute the unsteady time step for the dual time strategy if the unsteady CFL is diferent from 0. ---*/
    
    if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
        
        /*--- Thread-local variable for reduction. ---*/
        su2double glbDtND = 1e30;
        
        SU2_OMP(for schedule(static,omp_chunk_size) nowait)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            glbDtND = min(glbDtND, config->GetUnst_CFL()*Global_Delta_Time / nodes->GetLocalCFL(iPoint));
        }
        SU2_OMP_CRITICAL
        Global_Delta_UnstTimeND = min(Global_Delta_UnstTimeND, glbDtND);
        SU2_OMP_BARRIER
        
        SU2_OMP_MASTER
        {
            SU2_MPI::Allreduce(&Global_Delta_UnstTimeND, &glbDtND, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            Global_Delta_UnstTimeND = glbDtND;
            
            config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
        }
        SU2_OMP_BARRIER
    }
    
    /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
    
    if (dual_time && !implicit) {
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), nodes->GetDelta_Time(iPoint));
            nodes->SetDelta_Time(iPoint, Local_Delta_Time);
        }
    }
    
}

void CEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
    
    const bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    const bool jst_scheme = (config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0);
    
    /*--- Pick one numerics object per thread. ---*/
    CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];
    
    /*--- Loop over edge colors. ---*/
    for (auto color : EdgeColoring)
    {
        /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
        SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
        for(auto k = 0ul; k < color.size; ++k) {
            
            auto iEdge = color.indices[k];
            
            /*--- Points in edge, set normal vectors, and number of neighbors ---*/
            
            auto iPoint = geometry->edges->GetNode(iEdge,0);
            auto jPoint = geometry->edges->GetNode(iEdge,1);
            
            numerics->SetNormal(geometry->edges->GetNormal(iEdge));
            numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
            
            /*--- Set primitive variables w/o reconstruction ---*/
            
            numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(jPoint));
            
            /*--- Set the largest convective eigenvalue ---*/
            
            numerics->SetLambda(nodes->GetLambda(iPoint), nodes->GetLambda(jPoint));
            
            /*--- Set undivided laplacian an pressure based sensor ---*/
            
            if (jst_scheme) {
                numerics->SetUndivided_Laplacian(nodes->GetUndivided_Laplacian(iPoint),
                                                 nodes->GetUndivided_Laplacian(jPoint));
                numerics->SetSensor(nodes->GetSensor(iPoint),
                                    nodes->GetSensor(jPoint));
            }
            
            /*--- Grid movement ---*/
            
            if (dynamic_grid) {
                numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
            }
            
            /*--- Compute residuals, and Jacobians ---*/
            
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Update convective and artificial dissipation residuals. ---*/
            
            if (ReducerStrategy) {
                EdgeFluxes.SetBlock(iEdge, residual);
                if (implicit)
                    Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
            }
            else {
                LinSysRes.AddBlock(iPoint, residual);
                LinSysRes.SubtractBlock(jPoint, residual);
                if (implicit)
                    Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
            }
            
            /*--- Viscous contribution. ---*/
            
            Viscous_Residual(iEdge, geometry, solver_container,
                             numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);
        }
    } // end color loop
    
    if (ReducerStrategy) {
        SumEdgeFluxes(geometry);
        if (implicit)
            Jacobian.SetDiagonalAsColumnSum();
    }
    
}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                   CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
    
    const auto InnerIter        = config->GetInnerIter();
    const bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    const bool ideal_gas        = (config->GetKind_FluidModel() == STANDARD_AIR) ||
    (config->GetKind_FluidModel() == IDEAL_GAS);
    
    const bool roe_turkel       = (config->GetKind_Upwind_Flow() == TURKEL);
    const bool low_mach_corr    = config->Low_Mach_Correction();
    const auto kind_dissipation = config->GetKind_RoeLowDiss();
    
    const bool muscl            = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
    const bool limiter          = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) &&
    (InnerIter <= config->GetLimiterIter());
    const bool van_albada       = (config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE);
    
    /*--- Non-physical counter. ---*/
    unsigned long counter_local = 0;
    SU2_OMP_MASTER
    ErrorCounter = 0;
    
    /*--- Pick one numerics object per thread. ---*/
    CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];
    
    /*--- Static arrays of MUSCL-reconstructed primitives and secondaries (thread safety). ---*/
    su2double Primitive_i[MAXNVAR] = {0.0}, Primitive_j[MAXNVAR] = {0.0};
    su2double Secondary_i[MAXNVAR] = {0.0}, Secondary_j[MAXNVAR] = {0.0};
    
    /*--- Loop over edge colors. ---*/
    for (auto color : EdgeColoring)
    {
        /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
        SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
        for(auto k = 0ul; k < color.size; ++k) {
            
            auto iEdge = color.indices[k];
            
            unsigned short iDim, iVar;
            
            /*--- Points in edge and normal vectors ---*/
            
            auto iPoint = geometry->edges->GetNode(iEdge,0);
            auto jPoint = geometry->edges->GetNode(iEdge,1);
            
            numerics->SetNormal(geometry->edges->GetNormal(iEdge));
            
            auto Coord_i = geometry->node[iPoint]->GetCoord();
            auto Coord_j = geometry->node[jPoint]->GetCoord();
            
            /*--- Roe Turkel preconditioning ---*/
            
            if (roe_turkel) {
                su2double sqvel = 0.0;
                for (iDim = 0; iDim < nDim; iDim ++)
                    sqvel += pow(config->GetVelocity_FreeStream()[iDim], 2);
                numerics->SetVelocity2_Inf(sqvel);
            }
            
            /*--- Grid movement ---*/
            
            if (dynamic_grid) {
                numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                     geometry->node[jPoint]->GetGridVel());
            }
            
            /*--- Get primitive and secondary variables ---*/
            
            auto V_i = nodes->GetPrimitive(iPoint); auto V_j = nodes->GetPrimitive(jPoint);
            auto S_i = nodes->GetSecondary(iPoint); auto S_j = nodes->GetSecondary(jPoint);
            
            /*--- Set them with or without high order reconstruction using MUSCL strategy. ---*/
            
            if (!muscl) {
                
                numerics->SetPrimitive(V_i, V_j);
                numerics->SetSecondary(S_i, S_j);
                
            }
            else {
                /*--- Reconstruction ---*/
                
                su2double Vector_ij[MAXNDIM] = {0.0};
                for (iDim = 0; iDim < nDim; iDim++) {
                    Vector_ij[iDim] = 0.5*(Coord_j[iDim] - Coord_i[iDim]);
                }
                
                auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
                auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);
                
                su2double *Limiter_i = nullptr, *Limiter_j = nullptr;
                
                if (limiter) {
                    Limiter_i = nodes->GetLimiter_Primitive(iPoint);
                    Limiter_j = nodes->GetLimiter_Primitive(jPoint);
                }
                
                for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                    
                    su2double Project_Grad_i = 0.0;
                    su2double Project_Grad_j = 0.0;
                    
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Project_Grad_i += Vector_ij[iDim]*Gradient_i[iVar][iDim];
                        Project_Grad_j -= Vector_ij[iDim]*Gradient_j[iVar][iDim];
                    }
                    
                    if (limiter) {
                        if (van_albada) {
                            su2double V_ij = V_j[iVar] - V_i[iVar];
                            Limiter_i[iVar] = V_ij*( 2.0*Project_Grad_i + V_ij) / (4*pow(Project_Grad_i, 2) + pow(V_ij, 2) + EPS);
                            Limiter_j[iVar] = V_ij*(-2.0*Project_Grad_j + V_ij) / (4*pow(Project_Grad_j, 2) + pow(V_ij, 2) + EPS);
                        }
                        Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
                        Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
                    }
                    else {
                        Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
                        Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
                    }
                    
                }
                
                /*--- Recompute the reconstructed quantities in a thermodynamically consistent way. ---*/
                
                if (!ideal_gas || low_mach_corr) {
                    ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_i, Secondary_i);
                    ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_j, Secondary_j);
                }
                
                /*--- Low-Mach number correction. ---*/
                
                if (low_mach_corr) {
                    LowMachPrimitiveCorrection(GetFluidModel(), nDim, Primitive_i, Primitive_j);
                }
                
                /*--- Check for non-physical solutions after reconstruction. If found, use the
                 cell-average value of the solution. This is a locally 1st order approximation,
                 which is typically only active during the start-up of a calculation. ---*/
                
                bool neg_pres_or_rho_i = (Primitive_i[nDim+1] < 0.0) || (Primitive_i[nDim+2] < 0.0);
                bool neg_pres_or_rho_j = (Primitive_j[nDim+1] < 0.0) || (Primitive_j[nDim+2] < 0.0);
                
                su2double R = sqrt(fabs(Primitive_j[nDim+2]/Primitive_i[nDim+2]));
                su2double sq_vel = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    su2double RoeVelocity = (R*Primitive_j[iDim+1]+Primitive_i[iDim+1])/(R+1);
                    sq_vel += pow(RoeVelocity, 2);
                }
                su2double RoeEnthalpy = (R*Primitive_j[nDim+3]+Primitive_i[nDim+3])/(R+1);
                
                bool neg_sound_speed = ((Gamma-1)*(RoeEnthalpy-0.5*sq_vel) < 0.0);
                
                bool bad_i = neg_sound_speed || neg_pres_or_rho_i;
                bool bad_j = neg_sound_speed || neg_pres_or_rho_j;
                
                nodes->SetNon_Physical(iPoint, bad_i);
                nodes->SetNon_Physical(jPoint, bad_j);
                
                /*--- Get updated state, in case the point recovered after the set. ---*/
                bad_i = nodes->GetNon_Physical(iPoint);
                bad_j = nodes->GetNon_Physical(jPoint);
                
                counter_local += bad_i+bad_j;
                
                numerics->SetPrimitive(bad_i? V_i : Primitive_i,  bad_j? V_j : Primitive_j);
                numerics->SetSecondary(bad_i? S_i : Secondary_i,  bad_j? S_j : Secondary_j);
                
            }
            
            /*--- Roe Low Dissipation Scheme ---*/
            
            if (kind_dissipation != NO_ROELOWDISS) {
                
                numerics->SetDissipation(nodes->GetRoe_Dissipation(iPoint),
                                         nodes->GetRoe_Dissipation(jPoint));
                
                if (kind_dissipation == FD_DUCROS || kind_dissipation == NTS_DUCROS){
                    numerics->SetSensor(nodes->GetSensor(iPoint),
                                        nodes->GetSensor(jPoint));
                }
                if (kind_dissipation == NTS || kind_dissipation == NTS_DUCROS){
                    numerics->SetCoord(Coord_i, Coord_j);
                }
            }
            
            /*--- Compute the residual ---*/
            
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Set the final value of the Roe dissipation coefficient ---*/
            
            if ((kind_dissipation != NO_ROELOWDISS) && (MGLevel != MESH_0)) {
                nodes->SetRoe_Dissipation(iPoint,numerics->GetDissipation());
                nodes->SetRoe_Dissipation(jPoint,numerics->GetDissipation());
            }
            
            /*--- Update residual value ---*/
            
            if (ReducerStrategy) {
                EdgeFluxes.SetBlock(iEdge, residual);
                if (implicit)
                    Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
            }
            else {
                LinSysRes.AddBlock(iPoint, residual);
                LinSysRes.SubtractBlock(jPoint, residual);
                
                /*--- Set implicit computation ---*/
                if (implicit)
                    Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
            }
            
            /*--- Viscous contribution. ---*/
            
            Viscous_Residual(iEdge, geometry, solver_container,
                             numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);
        }
    } // end color loop
    
    if (ReducerStrategy) {
        SumEdgeFluxes(geometry);
        if (implicit)
            Jacobian.SetDiagonalAsColumnSum();
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
            SU2_MPI::Reduce(&counter_local, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
            config->SetNonphysical_Reconstr(ErrorCounter);
        }
        SU2_OMP_BARRIER
    }
    
}

void CEulerSolver::SumEdgeFluxes(CGeometry* geometry) {
    
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
        
        LinSysRes.SetBlock_Zero(iPoint);
        
        for (unsigned short iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
            
            auto iEdge = geometry->node[iPoint]->GetEdge(iNeigh);
            
            if (iPoint == geometry->edges->GetNode(iEdge,0))
                LinSysRes.AddBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
            else
                LinSysRes.SubtractBlock(iPoint, EdgeFluxes.GetBlock(iEdge));
        }
    }
    
}

void CEulerSolver::ComputeConsistentExtrapolation(CFluidModel *fluidModel, unsigned short nDim,
                                                  su2double *primitive, su2double *secondary) {
    
    su2double density = primitive[nDim+2];
    su2double pressure = primitive[nDim+1];
    
    su2double velocity2 = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
        velocity2 += pow(primitive[iDim+1], 2);
    
    fluidModel->SetTDState_Prho(pressure, density);
    
    primitive[0] = fluidModel->GetTemperature();
    primitive[nDim+3] = fluidModel->GetStaticEnergy() + primitive[nDim+1]/primitive[nDim+2] + 0.5*velocity2;
    primitive[nDim+4] = fluidModel->GetSoundSpeed();
    secondary[0] = fluidModel->GetdPdrho_e();
    secondary[1] = fluidModel->GetdPde_rho();
    
}

void CEulerSolver::LowMachPrimitiveCorrection(CFluidModel *fluidModel, unsigned short nDim,
                                              su2double *primitive_i, su2double *primitive_j) {
    unsigned short iDim;
    
    su2double velocity2_i = 0.0;
    su2double velocity2_j = 0.0;
    
    for (iDim = 0; iDim < nDim; iDim++) {
        velocity2_i += pow(primitive_i[iDim+1], 2);
        velocity2_j += pow(primitive_j[iDim+1], 2);
    }
    su2double mach_i = sqrt(velocity2_i)/primitive_i[nDim+4];
    su2double mach_j = sqrt(velocity2_j)/primitive_j[nDim+4];
    
    su2double z = min(max(mach_i,mach_j),1.0);
    velocity2_i = 0.0;
    velocity2_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
        su2double vel_i_corr = ( primitive_i[iDim+1] + primitive_j[iDim+1] )/2.0
        + z * ( primitive_i[iDim+1] - primitive_j[iDim+1] )/2.0;
        su2double vel_j_corr = ( primitive_i[iDim+1] + primitive_j[iDim+1] )/2.0
        + z * ( primitive_j[iDim+1] - primitive_i[iDim+1] )/2.0;
        
        velocity2_i += pow(vel_i_corr, 2);
        velocity2_j += pow(vel_j_corr, 2);
        
        primitive_i[iDim+1] = vel_i_corr;
        primitive_j[iDim+1] = vel_j_corr;
    }
    
    fluidModel->SetEnergy_Prho(primitive_i[nDim+1], primitive_i[nDim+2]);
    primitive_i[nDim+3]= fluidModel->GetStaticEnergy() + primitive_i[nDim+1]/primitive_i[nDim+2] + 0.5*velocity2_i;
    
    fluidModel->SetEnergy_Prho(primitive_j[nDim+1], primitive_j[nDim+2]);
    primitive_j[nDim+3]= fluidModel->GetStaticEnergy() + primitive_j[nDim+1]/primitive_j[nDim+2] + 0.5*velocity2_j;
    
}

void CEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                   CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
    
    const bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    const bool rotating_frame   = config->GetRotating_Frame();
    const bool axisymmetric     = config->GetAxisymmetric();
    const bool gravity          = (config->GetGravityForce() == YES);
    const bool harmonic_balance = (config->GetTime_Marching() == HARMONIC_BALANCE);
    const bool windgust         = config->GetWind_Gust();
    const bool body_force       = config->GetBody_Force();
    
    /*--- Pick one numerics object per thread. ---*/
    CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];
    
    unsigned short iVar;
    unsigned long iPoint;
    
    if (body_force) {
        
        /*--- Loop over all points ---*/
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Load the conservative variables ---*/
            numerics->SetConservative(nodes->GetSolution(iPoint),
                                      nodes->GetSolution(iPoint));
            
            /*--- Load the volume of the dual mesh cell ---*/
            numerics->SetVolume(geometry->node[iPoint]->GetVolume());
            
            /*--- Compute the rotating frame source residual ---*/
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Add the source residual to the total ---*/
            LinSysRes.AddBlock(iPoint, residual);
            
        }
    }
    
    if (rotating_frame) {
        
        /*--- Include the residual contribution from GCL due to the static
         mesh movement that is set for rotating frame. ---*/
        
        SU2_OMP_MASTER
        SetRotatingFrame_GCL(geometry, config);
        SU2_OMP_BARRIER
        
        /*--- Loop over all points ---*/
        SU2_OMP_FOR_DYN(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Load the conservative variables ---*/
            numerics->SetConservative(nodes->GetSolution(iPoint),
                                      nodes->GetSolution(iPoint));
            
            /*--- Load the volume of the dual mesh cell ---*/
            numerics->SetVolume(geometry->node[iPoint]->GetVolume());
            
            /*--- Compute the rotating frame source residual ---*/
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Add the source residual to the total ---*/
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Add the implicit Jacobian contribution ---*/
            if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
        }
    }
    
    if (axisymmetric) {
        
        /*--- loop over points ---*/
        SU2_OMP_FOR_DYN(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Set solution  ---*/
            numerics->SetConservative(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));
            
            /*--- Set control volume ---*/
            numerics->SetVolume(geometry->node[iPoint]->GetVolume());
            
            /*--- Set y coordinate ---*/
            numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
            
            /*--- Compute Source term Residual ---*/
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Add Residual ---*/
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Implicit part ---*/
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
        }
    }
    
    if (gravity) {
        
        /*--- loop over points ---*/
        SU2_OMP_FOR_DYN(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Set solution  ---*/
            numerics->SetConservative(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));
            
            /*--- Set control volume ---*/
            numerics->SetVolume(geometry->node[iPoint]->GetVolume());
            
            /*--- Compute Source term Residual ---*/
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Add Residual ---*/
            LinSysRes.AddBlock(iPoint, residual);
            
        }
        
    }
    
    if (harmonic_balance) {
        
        /*--- loop over points ---*/
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Get control volume ---*/
            su2double Volume = geometry->node[iPoint]->GetVolume();
            
            /*--- Get stored time spectral source term and add to residual ---*/
            for (iVar = 0; iVar < nVar; iVar++) {
                LinSysRes(iPoint,iVar) += Volume * nodes->GetHarmonicBalance_Source(iPoint,iVar);
            }
        }
    }
    
    if (windgust) {
        
        /*--- Loop over all points ---*/
        SU2_OMP_FOR_DYN(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Load the wind gust ---*/
            numerics->SetWindGust(nodes->GetWindGust(iPoint), nodes->GetWindGust(iPoint));
            
            /*--- Load the wind gust derivatives ---*/
            numerics->SetWindGustDer(nodes->GetWindGustDer(iPoint), nodes->GetWindGustDer(iPoint));
            
            /*--- Load the primitive variables ---*/
            numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(iPoint));
            
            /*--- Load the volume of the dual mesh cell ---*/
            numerics->SetVolume(geometry->node[iPoint]->GetVolume());
            
            /*--- Compute the rotating frame source residual ---*/
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Add the source residual to the total ---*/
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Add the implicit Jacobian contribution ---*/
            if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
        }
    }
    
    /*--- Check if a verification solution is to be computed. ---*/
    
    if ( VerificationSolution ) {
        if ( VerificationSolution->IsManufacturedSolution() ) {
            
            /*--- Get the physical time. ---*/
            su2double time = 0.0;
            if (config->GetTime_Marching()) time = config->GetPhysicalTime();
            
            /*--- Loop over points ---*/
            SU2_OMP_FOR_DYN(omp_chunk_size)
            for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
                
                /*--- Get control volume size. ---*/
                su2double Volume = geometry->node[iPoint]->GetVolume();
                
                /*--- Get the current point coordinates. ---*/
                const su2double *coor = geometry->node[iPoint]->GetCoord();
                
                /*--- Get the MMS source term. ---*/
                vector<su2double> sourceMan(nVar,0.0);
                VerificationSolution->GetMMSSourceTerm(coor, time, sourceMan.data());
                
                /*--- Compute the residual for this control volume and subtract. ---*/
                for (iVar = 0; iVar < nVar; iVar++) {
                    LinSysRes(iPoint,iVar) -= sourceMan[iVar]*Volume;
                }
            }
        }
    }
    
}

void CEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
    
    /* This method should be used to call any new source terms for a particular problem*/
    /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */
    
    /* Next we describe how to get access to some important quanties for this method */
    /* Access to all points in the current geometric mesh by saying: nPointDomain */
    /* Get the vector of conservative variables at some point iPoint = nodes->GetSolution(iPoint) */
    /* Get the volume (or area in 2D) associated with iPoint = nodes->GetVolume(iPoint) */
    /* Get the vector of geometric coordinates of point iPoint = nodes->GetCoord(iPoint) */
    
}

void CEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
    
    /*--- Loop domain points. ---*/
    
    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
        
        /*--- Set eigenvalues to zero. ---*/
        nodes->SetLambda(iPoint,0.0);
        
        /*--- Loop over the neighbors of point i. ---*/
        for (unsigned short iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh)
        {
            auto jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
            
            auto iEdge = geometry->node[iPoint]->GetEdge(iNeigh);
            auto Normal = geometry->edges->GetNormal(iEdge);
            su2double Area = 0.0;
            for (unsigned short iDim = 0; iDim < nDim; iDim++) Area += pow(Normal[iDim],2);
            Area = sqrt(Area);
            
            /*--- Mean Values ---*/
            
            su2double Mean_ProjVel = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
            su2double Mean_SoundSpeed = 0.5 * (nodes->GetSoundSpeed(iPoint) + nodes->GetSoundSpeed(jPoint)) * Area;
            
            /*--- Adjustment for grid movement ---*/
            
            if (dynamic_grid) {
                const su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
                const su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
                
                for (unsigned short iDim = 0; iDim < nDim; iDim++)
                    Mean_ProjVel -= 0.5 * (GridVel_i[iDim] + GridVel_j[iDim]) * Normal[iDim];
            }
            
            /*--- Inviscid contribution ---*/
            
            su2double Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
            nodes->AddLambda(iPoint, Lambda);
        }
        
    }
    
    /*--- Loop boundary edges ---*/
    
    for (unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
            (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
            
            SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
            for (unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                
                /*--- Point identification, Normal vector and area ---*/
                
                auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                su2double Area = 0.0;
                for (unsigned short iDim = 0; iDim < nDim; iDim++)
                    Area += pow(Normal[iDim],2);
                Area = sqrt(Area);
                
                /*--- Mean Values ---*/
                
                su2double Mean_ProjVel = nodes->GetProjVel(iPoint,Normal);
                su2double Mean_SoundSpeed = nodes->GetSoundSpeed(iPoint) * Area;
                
                /*--- Adjustment for grid movement ---*/
                
                if (dynamic_grid) {
                    auto GridVel = geometry->node[iPoint]->GetGridVel();
                    for (unsigned short iDim = 0; iDim < nDim; iDim++)
                        Mean_ProjVel -= GridVel[iDim]*Normal[iDim];
                }
                
                /*--- Inviscid contribution ---*/
                
                su2double Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
                if (geometry->node[iPoint]->GetDomain()) {
                    nodes->AddLambda(iPoint,Lambda);
                }
            }
        }
    }
    
    SU2_OMP_MASTER
    {
        /*--- Correct the eigenvalue values across any periodic boundaries. ---*/
        
        for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
            InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
            CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
        }
        
        /*--- MPI parallelization ---*/
        
        InitiateComms(geometry, config, MAX_EIGENVALUE);
        CompleteComms(geometry, config, MAX_EIGENVALUE);
    }
    SU2_OMP_BARRIER
    
}

void CEulerSolver::SetUndivided_Laplacian_And_Centered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) {
    
    /*--- We can access memory more efficiently if there are no periodic boundaries. ---*/
    
    const bool isPeriodic = (config->GetnMarker_Periodic() > 0);
    
    /*--- Loop domain points. ---*/
    
    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
        
        const bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
        const su2double Pressure_i = nodes->GetPressure(iPoint);
        
        /*--- Initialize. ---*/
        for (unsigned short iVar = 0; iVar < nVar; iVar++)
            nodes->SetUnd_Lapl(iPoint, iVar, 0.0);
        
        iPoint_UndLapl[iPoint] = 0.0;
        jPoint_UndLapl[iPoint] = 0.0;
        
        /*--- Loop over the neighbors of point i. ---*/
        for (unsigned short iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh)
        {
            auto jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
            bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
            
            /*--- If iPoint is boundary it only takes contributions from other boundary points. ---*/
            if (boundary_i && !boundary_j) continue;
            
            /*--- Add solution differences, with correction for compressible flows which use the enthalpy. ---*/
            
            for (unsigned short iVar = 0; iVar < nVar; iVar++)
                nodes->AddUnd_Lapl(iPoint, iVar, nodes->GetSolution(jPoint,iVar)-nodes->GetSolution(iPoint,iVar));
            
            su2double Pressure_j = nodes->GetPressure(jPoint);
            nodes->AddUnd_Lapl(iPoint, nVar-1, Pressure_j-Pressure_i);
            
            /*--- Dissipation sensor, add pressure difference and pressure sum. ---*/
            iPoint_UndLapl[iPoint] += Pressure_j - Pressure_i;
            jPoint_UndLapl[iPoint] += Pressure_j + Pressure_i;
        }
        
        if (!isPeriodic)
            nodes->SetSensor(iPoint, fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
    }
    
    if (isPeriodic) {
        /*--- Correct the Laplacian and sensor values across any periodic boundaries. ---*/
        
        SU2_OMP_MASTER
        {
            for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
                InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
                CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
                
                InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
                CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
            }
        }
        SU2_OMP_BARRIER
        
        /*--- Set final pressure switch for each point ---*/
        
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
            nodes->SetSensor(iPoint, fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
    }
    
    SU2_OMP_MASTER
    {
        /*--- MPI parallelization ---*/
        
        InitiateComms(geometry, config, UNDIVIDED_LAPLACIAN);
        CompleteComms(geometry, config, UNDIVIDED_LAPLACIAN);
        
        InitiateComms(geometry, config, SENSOR);
        CompleteComms(geometry, config, SENSOR);
    }
    SU2_OMP_BARRIER
    
}

void CEulerSolver::SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config){
    
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*---- Ducros sensor for iPoint and its neighbor points to avoid lower dissipation near shocks. ---*/
        
        su2double Ducros_i = 0.0;
        auto nNeigh = geometry->node[iPoint]->GetnNeighbor();
        
        for (unsigned short iNeigh = 0; iNeigh <= nNeigh; iNeigh++) {
            
            auto jPoint = iPoint;
            if (iNeigh < nNeigh) jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
            
            /*---- Dilatation for jPoint ---*/
            
            su2double uixi=0.0;
            for(unsigned short iDim = 0; iDim < nDim; iDim++){
                uixi += nodes->GetGradient_Primitive(jPoint,iDim+1, iDim);
            }
            
            /*--- Compute norm of vorticity ---*/
            
            const su2double* Vorticity = nodes->GetVorticity(jPoint);
            su2double Omega = 0.0;
            for (unsigned short iDim = 0; iDim < nDim; iDim++) {
                Omega += pow(Vorticity[iDim], 2);
            }
            Omega = sqrt(Omega);
            
            su2double Ducros_j = 0.0;
            
            if (config->GetKind_RoeLowDiss() == FD_DUCROS) {
                Ducros_j = -uixi / (fabs(uixi) + Omega + 1e-20);
            }
            else if (config->GetKind_RoeLowDiss() == NTS_DUCROS) {
                Ducros_j = pow(uixi,2.0) /(pow(uixi,2.0)+ pow(Omega,2.0) + 1e-20);
            }
            Ducros_i = max(Ducros_i, Ducros_j);
        }
        
        nodes->SetSensor(iPoint, Ducros_i);
    }
    
    SU2_OMP_MASTER
    {
        InitiateComms(geometry, config, SENSOR);
        CompleteComms(geometry, config, SENSOR);
    }
    SU2_OMP_BARRIER
    
}

void CEulerSolver::Pressure_Forces(CGeometry *geometry, CConfig *config) {
    
    unsigned long iVertex, iPoint;
    unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
    su2double Pressure = 0.0, factor, NFPressOF, RefVel2,
    RefTemp, RefDensity, RefPressure, Mach2Vel, Mach_Motion;
    const su2double *Normal = nullptr, *Coord = nullptr;
    string Marker_Tag, Monitoring_Tag;
    su2double AxiFactor;
    
    su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
    su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
    su2double RefArea = config->GetRefArea();
    su2double RefLength = config->GetRefLength();
    su2double Gas_Constant = config->GetGas_ConstantND();
    const su2double *Origin = nullptr;
    if (config->GetnMarker_Monitoring() != 0){
        Origin = config->GetRefOriginMoment(0);
    }
    bool axisymmetric = config->GetAxisymmetric();
    
    /*--- Evaluate reference values for non-dimensionalization.
     For dynamic meshes, use the motion Mach number as a reference value
     for computing the force coefficients. Otherwise, use the freestream
     values, which is the standard convention. ---*/
    
    RefTemp = Temperature_Inf;
    RefDensity = Density_Inf;
    RefPressure = Pressure_Inf;
    if (dynamic_grid) {
        Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
        Mach_Motion = config->GetMach_Motion();
        RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    }
    else {
        RefVel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
            RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
    
    factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
    
    /*-- Variables initialization ---*/
    
    TotalCoeff.setZero();
    
    Total_CNearFieldOF = 0.0; Total_Heat = 0.0;  Total_MaxHeat = 0.0;
    
    AllBoundInvCoeff.setZero();
    
    AllBound_CNearFieldOF_Inv = 0.0;
    
    SurfaceInvCoeff.setZero();
    SurfaceCoeff.setZero();
    
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
            
            InvCoeff.setZero(iMarker);
            
            CNearFieldOF_Inv[iMarker] = 0.0;
            
            su2double ForceInviscid[MAXNDIM] = {0.0}, MomentInviscid[MAXNDIM] = {0.0};
            su2double MomentX_Force[MAXNDIM] = {0.0}, MomentY_Force[MAXNDIM] = {0.0}, MomentZ_Force[MAXNDIM] = {0.0};
            
            NFPressOF = 0.0;
            
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
                    
                    /*--- Quadratic objective function for the near-field.
                     This uses the infinity pressure regardless of Mach number. ---*/
                    
                    NFPressOF += 0.5*(Pressure - Pressure_Inf)*(Pressure - Pressure_Inf)*Normal[nDim-1];
                    
                    su2double MomentDist[MAXNDIM] = {0.0};
                    for (iDim = 0; iDim < nDim; iDim++) {
                        MomentDist[iDim] = Coord[iDim] - Origin[iDim];
                    }
                    
                    /*--- Axisymmetric simulations ---*/
                    
                    if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
                    else AxiFactor = 1.0;
                    
                    /*--- Force computation, note the minus sign due to the
                     orientation of the normal (outward) ---*/
                    
                    su2double Force[MAXNDIM] = {0.0};
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
                        InvCoeff.CD[iMarker]     =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
                        InvCoeff.CL[iMarker]     = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
                        InvCoeff.CEff[iMarker]   = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker]+EPS);
                        InvCoeff.CMz[iMarker]    = MomentInviscid[2];
                        InvCoeff.CoPx[iMarker]   = MomentZ_Force[1];
                        InvCoeff.CoPy[iMarker]   = -MomentZ_Force[0];
                        InvCoeff.CFx[iMarker]    = ForceInviscid[0];
                        InvCoeff.CFy[iMarker]    = ForceInviscid[1];
                        InvCoeff.CT[iMarker]     = -InvCoeff.CFx[iMarker];
                        InvCoeff.CQ[iMarker]     = -InvCoeff.CMz[iMarker];
                        InvCoeff.CMerit[iMarker] = InvCoeff.CT[iMarker] / (InvCoeff.CQ[iMarker] + EPS);
                    }
                    if (nDim == 3) {
                        InvCoeff.CD[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
                        InvCoeff.CL[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
                        InvCoeff.CSF[iMarker]     = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
                        InvCoeff.CEff[iMarker]    = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker] + EPS);
                        InvCoeff.CMx[iMarker]     = MomentInviscid[0];
                        InvCoeff.CMy[iMarker]     = MomentInviscid[1];
                        InvCoeff.CMz[iMarker]     = MomentInviscid[2];
                        InvCoeff.CoPx[iMarker]    = -MomentY_Force[0];
                        InvCoeff.CoPz[iMarker]    = MomentY_Force[2];
                        InvCoeff.CFx[iMarker]     = ForceInviscid[0];
                        InvCoeff.CFy[iMarker]     = ForceInviscid[1];
                        InvCoeff.CFz[iMarker]     = ForceInviscid[2];
                        InvCoeff.CT[iMarker]      = -InvCoeff.CFz[iMarker];
                        InvCoeff.CQ[iMarker]      = -InvCoeff.CMz[iMarker];
                        InvCoeff.CMerit[iMarker]  = InvCoeff.CT[iMarker] / (InvCoeff.CQ[iMarker] + EPS);
                    }
                    
                    AllBoundInvCoeff.CD           += InvCoeff.CD[iMarker];
                    AllBoundInvCoeff.CL           += InvCoeff.CL[iMarker];
                    AllBoundInvCoeff.CSF          += InvCoeff.CSF[iMarker];
                    AllBoundInvCoeff.CEff          = AllBoundInvCoeff.CL / (AllBoundInvCoeff.CD + EPS);
                    AllBoundInvCoeff.CMx          += InvCoeff.CMx[iMarker];
                    AllBoundInvCoeff.CMy          += InvCoeff.CMy[iMarker];
                    AllBoundInvCoeff.CMz          += InvCoeff.CMz[iMarker];
                    AllBoundInvCoeff.CoPx         += InvCoeff.CoPx[iMarker];
                    AllBoundInvCoeff.CoPy         += InvCoeff.CoPy[iMarker];
                    AllBoundInvCoeff.CoPz         += InvCoeff.CoPz[iMarker];
                    AllBoundInvCoeff.CFx          += InvCoeff.CFx[iMarker];
                    AllBoundInvCoeff.CFy          += InvCoeff.CFy[iMarker];
                    AllBoundInvCoeff.CFz          += InvCoeff.CFz[iMarker];
                    AllBoundInvCoeff.CT           += InvCoeff.CT[iMarker];
                    AllBoundInvCoeff.CQ           += InvCoeff.CQ[iMarker];
                    AllBoundInvCoeff.CMerit        = AllBoundInvCoeff.CT / (AllBoundInvCoeff.CQ + EPS);
                    
                    /*--- Compute the coefficients per surface ---*/
                    
                    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
                        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                        if (Marker_Tag == Monitoring_Tag) {
                            SurfaceInvCoeff.CL[iMarker_Monitoring]      += InvCoeff.CL[iMarker];
                            SurfaceInvCoeff.CD[iMarker_Monitoring]      += InvCoeff.CD[iMarker];
                            SurfaceInvCoeff.CSF[iMarker_Monitoring]     += InvCoeff.CSF[iMarker];
                            SurfaceInvCoeff.CEff[iMarker_Monitoring]     = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker] + EPS);
                            SurfaceInvCoeff.CFx[iMarker_Monitoring]     += InvCoeff.CFx[iMarker];
                            SurfaceInvCoeff.CFy[iMarker_Monitoring]     += InvCoeff.CFy[iMarker];
                            SurfaceInvCoeff.CFz[iMarker_Monitoring]     += InvCoeff.CFz[iMarker];
                            SurfaceInvCoeff.CMx[iMarker_Monitoring]     += InvCoeff.CMx[iMarker];
                            SurfaceInvCoeff.CMy[iMarker_Monitoring]     += InvCoeff.CMy[iMarker];
                            SurfaceInvCoeff.CMz[iMarker_Monitoring]     += InvCoeff.CMz[iMarker];
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
    
    if (config->GetComm_Level() == COMM_FULL) {
        
        auto Allreduce = [](su2double x) {
            su2double tmp = x; x = 0.0;
            SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return x;
        };
        AllBoundInvCoeff.CD = Allreduce(AllBoundInvCoeff.CD);
        AllBoundInvCoeff.CL = Allreduce(AllBoundInvCoeff.CL);
        AllBoundInvCoeff.CSF = Allreduce(AllBoundInvCoeff.CSF);
        AllBoundInvCoeff.CEff = AllBoundInvCoeff.CL / (AllBoundInvCoeff.CD + EPS);
        
        AllBoundInvCoeff.CMx = Allreduce(AllBoundInvCoeff.CMx);
        AllBoundInvCoeff.CMy = Allreduce(AllBoundInvCoeff.CMy);
        AllBoundInvCoeff.CMz = Allreduce(AllBoundInvCoeff.CMz);
        
        AllBoundInvCoeff.CoPx = Allreduce(AllBoundInvCoeff.CoPx);
        AllBoundInvCoeff.CoPy = Allreduce(AllBoundInvCoeff.CoPy);
        AllBoundInvCoeff.CoPz = Allreduce(AllBoundInvCoeff.CoPz);
        
        AllBoundInvCoeff.CFx = Allreduce(AllBoundInvCoeff.CFx);
        AllBoundInvCoeff.CFy = Allreduce(AllBoundInvCoeff.CFy);
        AllBoundInvCoeff.CFz = Allreduce(AllBoundInvCoeff.CFz);
        
        AllBoundInvCoeff.CT = Allreduce(AllBoundInvCoeff.CT);
        AllBoundInvCoeff.CQ = Allreduce(AllBoundInvCoeff.CQ);
        AllBoundInvCoeff.CMerit = AllBoundInvCoeff.CT / (AllBoundInvCoeff.CQ + EPS);
        AllBound_CNearFieldOF_Inv = Allreduce(AllBound_CNearFieldOF_Inv);
        
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
        
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CL);
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CD);
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CSF);
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; iMarker_Monitoring++)
            SurfaceInvCoeff.CEff[iMarker_Monitoring] = SurfaceInvCoeff.CL[iMarker_Monitoring] / (SurfaceInvCoeff.CD[iMarker_Monitoring] + EPS);
        
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFx);
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFy);
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFz);
        
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMx);
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMy);
        Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMz);
        
        delete [] buffer;
        
    }
    
#endif
    
    /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
    
    TotalCoeff.CD            = AllBoundInvCoeff.CD;
    TotalCoeff.CL            = AllBoundInvCoeff.CL;
    TotalCoeff.CSF           = AllBoundInvCoeff.CSF;
    TotalCoeff.CEff          = TotalCoeff.CL / (TotalCoeff.CD + EPS);
    TotalCoeff.CFx           = AllBoundInvCoeff.CFx;
    TotalCoeff.CFy           = AllBoundInvCoeff.CFy;
    TotalCoeff.CFz           = AllBoundInvCoeff.CFz;
    TotalCoeff.CMx           = AllBoundInvCoeff.CMx;
    TotalCoeff.CMy           = AllBoundInvCoeff.CMy;
    TotalCoeff.CMz           = AllBoundInvCoeff.CMz;
    TotalCoeff.CoPx          = AllBoundInvCoeff.CoPx;
    TotalCoeff.CoPy          = AllBoundInvCoeff.CoPy;
    TotalCoeff.CoPz          = AllBoundInvCoeff.CoPz;
    TotalCoeff.CT            = AllBoundInvCoeff.CT;
    TotalCoeff.CQ            = AllBoundInvCoeff.CQ;
    TotalCoeff.CMerit        = TotalCoeff.CT / (TotalCoeff.CQ + EPS);
    Total_CNearFieldOF       = AllBound_CNearFieldOF_Inv;
    
    /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
    
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        SurfaceCoeff.CL[iMarker_Monitoring]      = SurfaceInvCoeff.CL[iMarker_Monitoring];
        SurfaceCoeff.CD[iMarker_Monitoring]      = SurfaceInvCoeff.CD[iMarker_Monitoring];
        SurfaceCoeff.CSF[iMarker_Monitoring]     = SurfaceInvCoeff.CSF[iMarker_Monitoring];
        SurfaceCoeff.CEff[iMarker_Monitoring]    = SurfaceInvCoeff.CL[iMarker_Monitoring] / (SurfaceInvCoeff.CD[iMarker_Monitoring] + EPS);
        SurfaceCoeff.CFx[iMarker_Monitoring]     = SurfaceInvCoeff.CFx[iMarker_Monitoring];
        SurfaceCoeff.CFy[iMarker_Monitoring]     = SurfaceInvCoeff.CFy[iMarker_Monitoring];
        SurfaceCoeff.CFz[iMarker_Monitoring]     = SurfaceInvCoeff.CFz[iMarker_Monitoring];
        SurfaceCoeff.CMx[iMarker_Monitoring]     = SurfaceInvCoeff.CMx[iMarker_Monitoring];
        SurfaceCoeff.CMy[iMarker_Monitoring]     = SurfaceInvCoeff.CMy[iMarker_Monitoring];
        SurfaceCoeff.CMz[iMarker_Monitoring]     = SurfaceInvCoeff.CMz[iMarker_Monitoring];
    }
    
}

void CEulerSolver::Momentum_Forces(CGeometry *geometry, CConfig *config) {
    
    unsigned long iVertex, iPoint;
    unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
    su2double Area, factor, RefVel2, RefTemp, RefDensity,  Mach2Vel, Mach_Motion, MassFlow, Density;
    const su2double *Normal = nullptr, *Coord = nullptr;
    string Marker_Tag, Monitoring_Tag;
    su2double AxiFactor;
    
    su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
    su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
    su2double RefArea = config->GetRefArea();
    su2double RefLength = config->GetRefLength();
    su2double Gas_Constant = config->GetGas_ConstantND();
    const su2double *Origin = nullptr;
    if (config->GetnMarker_Monitoring() != 0){
        Origin = config->GetRefOriginMoment(0);
    }
    bool axisymmetric          = config->GetAxisymmetric();
    
    /*--- Evaluate reference values for non-dimensionalization.
     For dynamic meshes, use the motion Mach number as a reference value
     for computing the force coefficients. Otherwise, use the freestream values,
     which is the standard convention. ---*/
    
    RefTemp     = Temperature_Inf;
    RefDensity  = Density_Inf;
    if (dynamic_grid) {
        Mach2Vel = sqrt(Gamma*Gas_Constant*RefTemp);
        Mach_Motion = config->GetMach_Motion();
        RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    }
    else {
        RefVel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
            RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    }
    
    factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
    
    /*-- Variables initialization ---*/
    
    AllBoundMntCoeff.setZero();
    SurfaceMntCoeff.setZero();
    
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
            
            MntCoeff.setZero(iMarker);
            
            su2double ForceMomentum[MAXNDIM] = {0.0}, MomentMomentum[MAXNDIM] = {0.0};
            su2double MomentX_Force[3] = {0.0}, MomentY_Force[3] = {0.0}, MomentZ_Force[3] = {0.0};
            
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
                    su2double Velocity[MAXNDIM] = {0.0}, MomentDist[MAXNDIM] = {0.0};
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
                    
                    su2double Force[MAXNDIM] = {0.0};
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
                    MntCoeff.CD[iMarker]     =  ForceMomentum[0]*cos(Alpha) + ForceMomentum[1]*sin(Alpha);
                    MntCoeff.CL[iMarker]     = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[1]*cos(Alpha);
                    MntCoeff.CEff[iMarker]   = MntCoeff.CL[iMarker] / (MntCoeff.CD[iMarker]+EPS);
                    MntCoeff.CFx[iMarker]    = ForceMomentum[0];
                    MntCoeff.CFy[iMarker]    = ForceMomentum[1];
                    MntCoeff.CMz[iMarker]    = MomentMomentum[2];
                    MntCoeff.CoPx[iMarker]   = MomentZ_Force[1];
                    MntCoeff.CoPy[iMarker]   = -MomentZ_Force[0];
                    MntCoeff.CT[iMarker]     = -MntCoeff.CFx[iMarker];
                    MntCoeff.CQ[iMarker]     = -MntCoeff.CMz[iMarker];
                    MntCoeff.CMerit[iMarker] = MntCoeff.CT[iMarker] / (MntCoeff.CQ[iMarker] + EPS);
                }
                if (nDim == 3) {
                    MntCoeff.CD[iMarker]         =  ForceMomentum[0]*cos(Alpha)*cos(Beta) + ForceMomentum[1]*sin(Beta) + ForceMomentum[2]*sin(Alpha)*cos(Beta);
                    MntCoeff.CL[iMarker]         = -ForceMomentum[0]*sin(Alpha) + ForceMomentum[2]*cos(Alpha);
                    MntCoeff.CSF[iMarker]        = -ForceMomentum[0]*sin(Beta)*cos(Alpha) + ForceMomentum[1]*cos(Beta) - ForceMomentum[2]*sin(Beta)*sin(Alpha);
                    MntCoeff.CEff[iMarker]       = MntCoeff.CL[iMarker] / (MntCoeff.CD[iMarker] + EPS);
                    MntCoeff.CFx[iMarker]        = ForceMomentum[0];
                    MntCoeff.CFy[iMarker]        = ForceMomentum[1];
                    MntCoeff.CFz[iMarker]        = ForceMomentum[2];
                    MntCoeff.CMx[iMarker]        = MomentMomentum[0];
                    MntCoeff.CMy[iMarker]        = MomentMomentum[1];
                    MntCoeff.CMz[iMarker]        = MomentMomentum[2];
                    MntCoeff.CoPx[iMarker]       = -MomentY_Force[0];
                    MntCoeff.CoPz[iMarker]       =  MomentY_Force[2];
                    MntCoeff.CT[iMarker]         = -MntCoeff.CFz[iMarker];
                    MntCoeff.CQ[iMarker]         = -MntCoeff.CMz[iMarker];
                    MntCoeff.CMerit[iMarker]     = MntCoeff.CT[iMarker] / (MntCoeff.CQ[iMarker] + EPS);
                }
                
                AllBoundMntCoeff.CD           += MntCoeff.CD[iMarker];
                AllBoundMntCoeff.CL           += MntCoeff.CL[iMarker];
                AllBoundMntCoeff.CSF          += MntCoeff.CSF[iMarker];
                AllBoundMntCoeff.CEff          = AllBoundMntCoeff.CL / (AllBoundMntCoeff.CD + EPS);
                AllBoundMntCoeff.CFx          += MntCoeff.CFx[iMarker];
                AllBoundMntCoeff.CFy          += MntCoeff.CFy[iMarker];
                AllBoundMntCoeff.CFz          += MntCoeff.CFz[iMarker];
                AllBoundMntCoeff.CMx          += MntCoeff.CMx[iMarker];
                AllBoundMntCoeff.CMy          += MntCoeff.CMy[iMarker];
                AllBoundMntCoeff.CMx          += MntCoeff.CMz[iMarker];
                AllBoundMntCoeff.CoPx         += MntCoeff.CoPx[iMarker];
                AllBoundMntCoeff.CoPy         += MntCoeff.CoPy[iMarker];
                AllBoundMntCoeff.CoPz         += MntCoeff.CoPz[iMarker];
                AllBoundMntCoeff.CT           += MntCoeff.CT[iMarker];
                AllBoundMntCoeff.CQ           += MntCoeff.CQ[iMarker];
                AllBoundMntCoeff.CMerit       += AllBoundMntCoeff.CT / (AllBoundMntCoeff.CQ + EPS);
                
                /*--- Compute the coefficients per surface ---*/
                
                for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                    Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
                    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                    if (Marker_Tag == Monitoring_Tag) {
                        SurfaceMntCoeff.CL[iMarker_Monitoring]      += MntCoeff.CL[iMarker];
                        SurfaceMntCoeff.CD[iMarker_Monitoring]      += MntCoeff.CD[iMarker];
                        SurfaceMntCoeff.CSF[iMarker_Monitoring]     += MntCoeff.CSF[iMarker];
                        SurfaceMntCoeff.CEff[iMarker_Monitoring]     = MntCoeff.CL[iMarker] / (MntCoeff.CD[iMarker] + EPS);
                        SurfaceMntCoeff.CFx[iMarker_Monitoring]     += MntCoeff.CFx[iMarker];
                        SurfaceMntCoeff.CFy[iMarker_Monitoring]     += MntCoeff.CFy[iMarker];
                        SurfaceMntCoeff.CFz[iMarker_Monitoring]     += MntCoeff.CFz[iMarker];
                        SurfaceMntCoeff.CMx[iMarker_Monitoring]     += MntCoeff.CMx[iMarker];
                        SurfaceMntCoeff.CMy[iMarker_Monitoring]     += MntCoeff.CMy[iMarker];
                        SurfaceMntCoeff.CMz[iMarker_Monitoring]     += MntCoeff.CMz[iMarker];
                    }
                }
                
            }
            
        }
    }
    
#ifdef HAVE_MPI
    
    /*--- Add AllBound information using all the nodes ---*/
    
    if (config->GetComm_Level() == COMM_FULL) {
        
        auto Allreduce = [](su2double x) {
            su2double tmp = x; x = 0.0;
            SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return x;
        };
        
        AllBoundMntCoeff.CD = Allreduce(AllBoundMntCoeff.CD);
        AllBoundMntCoeff.CL = Allreduce(AllBoundMntCoeff.CL);
        AllBoundMntCoeff.CSF = Allreduce(AllBoundMntCoeff.CSF);
        AllBoundMntCoeff.CEff = AllBoundMntCoeff.CL / (AllBoundMntCoeff.CD + EPS);
        
        AllBoundMntCoeff.CFx = Allreduce(AllBoundMntCoeff.CFx);
        AllBoundMntCoeff.CFy = Allreduce(AllBoundMntCoeff.CFy);
        AllBoundMntCoeff.CFz = Allreduce(AllBoundMntCoeff.CFz);
        
        AllBoundMntCoeff.CMx = Allreduce(AllBoundMntCoeff.CMx);
        AllBoundMntCoeff.CMy = Allreduce(AllBoundMntCoeff.CMy);
        AllBoundMntCoeff.CMz = Allreduce(AllBoundMntCoeff.CMz);
        
        AllBoundMntCoeff.CoPx = Allreduce(AllBoundMntCoeff.CoPx);
        AllBoundMntCoeff.CoPy = Allreduce(AllBoundMntCoeff.CoPy);
        AllBoundMntCoeff.CoPz = Allreduce(AllBoundMntCoeff.CoPz);
        
        AllBoundMntCoeff.CT = Allreduce(AllBoundMntCoeff.CT);
        AllBoundMntCoeff.CQ = Allreduce(AllBoundMntCoeff.CQ);
        AllBoundMntCoeff.CMerit = AllBoundMntCoeff.CT / (AllBoundMntCoeff.CQ + EPS);
        
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
        
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CL);
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CD);
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CSF);
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; iMarker_Monitoring++)
            SurfaceMntCoeff.CEff[iMarker_Monitoring] = SurfaceMntCoeff.CL[iMarker_Monitoring] / (SurfaceMntCoeff.CD[iMarker_Monitoring] + EPS);
        
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CFx);
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CFy);
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CFz);
        
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CMx);
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CMy);
        Allreduce_inplace(nMarkerMon, SurfaceMntCoeff.CMz);
        
        delete [] buffer;
        
    }
    
#endif
    
    /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
    
    TotalCoeff.CD            += AllBoundMntCoeff.CD;
    TotalCoeff.CL            += AllBoundMntCoeff.CL;
    TotalCoeff.CSF           += AllBoundMntCoeff.CSF;
    TotalCoeff.CEff          = TotalCoeff.CL / (TotalCoeff.CD + EPS);
    TotalCoeff.CFx           += AllBoundMntCoeff.CFx;
    TotalCoeff.CFy           += AllBoundMntCoeff.CFy;
    TotalCoeff.CFz           += AllBoundMntCoeff.CFz;
    TotalCoeff.CMx           += AllBoundMntCoeff.CMx;
    TotalCoeff.CMy           += AllBoundMntCoeff.CMy;
    TotalCoeff.CMz           += AllBoundMntCoeff.CMz;
    TotalCoeff.CoPx          += AllBoundMntCoeff.CoPx;
    TotalCoeff.CoPy          += AllBoundMntCoeff.CoPy;
    TotalCoeff.CoPz          += AllBoundMntCoeff.CoPz;
    TotalCoeff.CT            += AllBoundMntCoeff.CT;
    TotalCoeff.CQ            += AllBoundMntCoeff.CQ;
    TotalCoeff.CMerit        = TotalCoeff.CT / (TotalCoeff.CQ + EPS);
    
    /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
    
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        SurfaceCoeff.CL[iMarker_Monitoring]         += SurfaceMntCoeff.CL[iMarker_Monitoring];
        SurfaceCoeff.CD[iMarker_Monitoring]         += SurfaceMntCoeff.CD[iMarker_Monitoring];
        SurfaceCoeff.CSF[iMarker_Monitoring]        += SurfaceMntCoeff.CSF[iMarker_Monitoring];
        SurfaceCoeff.CEff[iMarker_Monitoring]       += SurfaceMntCoeff.CL[iMarker_Monitoring] / (SurfaceMntCoeff.CD[iMarker_Monitoring] + EPS);
        SurfaceCoeff.CFx[iMarker_Monitoring]        += SurfaceMntCoeff.CFx[iMarker_Monitoring];
        SurfaceCoeff.CFy[iMarker_Monitoring]        += SurfaceMntCoeff.CFy[iMarker_Monitoring];
        SurfaceCoeff.CFz[iMarker_Monitoring]        += SurfaceMntCoeff.CFz[iMarker_Monitoring];
        SurfaceCoeff.CMx[iMarker_Monitoring]        += SurfaceMntCoeff.CMx[iMarker_Monitoring];
        SurfaceCoeff.CMy[iMarker_Monitoring]        += SurfaceMntCoeff.CMy[iMarker_Monitoring];
        SurfaceCoeff.CMz[iMarker_Monitoring]        += SurfaceMntCoeff.CMz[iMarker_Monitoring];
    }
    
}

template<ENUM_TIME_INT IntegrationType>
void CEulerSolver::Explicit_Iteration(CGeometry *geometry, CSolver **solver_container,
                                      CConfig *config, unsigned short iRKStep) {
    
    static_assert(IntegrationType == CLASSICAL_RK4_EXPLICIT ||
                  IntegrationType == RUNGE_KUTTA_EXPLICIT ||
                  IntegrationType == EULER_EXPLICIT, "");
    
    const bool adjoint = config->GetContinuous_Adjoint();
    
    const su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
    
    /*--- Hard-coded classical RK4 coefficients. Will be added to config. ---*/
    const su2double RK_FuncCoeff[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    const su2double RK_TimeCoeff[] = {0.5, 0.5, 1.0, 1.0};
    
    /*--- Set shared residual variables to 0 and declare
     *    local ones for current thread to work on. ---*/
    
    SU2_OMP_MASTER
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
    }
    SU2_OMP_BARRIER
    
    su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
    const su2double* coordMax[MAXNVAR] = {nullptr};
    unsigned long idxMax[MAXNVAR] = {0};
    
    /*--- Update the solution and residuals ---*/
    
    SU2_OMP(for schedule(static,omp_chunk_size) nowait)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
        su2double Vol = geometry->node[iPoint]->GetVolume() + geometry->node[iPoint]->GetPeriodicVolume();
        su2double Delta = nodes->GetDelta_Time(iPoint) / Vol;
        
        const su2double* Res_TruncError = nodes->GetResTruncError(iPoint);
        const su2double* Residual = LinSysRes.GetBlock(iPoint);
        
        if (!adjoint) {
            for (unsigned short iVar = 0; iVar < nVar; iVar++) {
                
                su2double Res = Residual[iVar] + Res_TruncError[iVar];
                
                /*--- "Static" switch which should be optimized at compile time. ---*/
                switch(IntegrationType) {
                        
                    case EULER_EXPLICIT:
                        nodes->AddSolution(iPoint,iVar, -Res*Delta);
                        break;
                        
                    case RUNGE_KUTTA_EXPLICIT:
                        nodes->AddSolution(iPoint, iVar, -Res*Delta*RK_AlphaCoeff);
                        break;
                        
                    case CLASSICAL_RK4_EXPLICIT:
                    {
                        su2double tmp_time = -1.0*RK_TimeCoeff[iRKStep]*Delta;
                        su2double tmp_func = -1.0*RK_FuncCoeff[iRKStep]*Delta;
                        
                        if (iRKStep < 3) {
                            /* Base Solution Update */
                            nodes->AddSolution(iPoint,iVar, tmp_time*Res);
                            
                            /* New Solution Update */
                            nodes->AddSolution_New(iPoint,iVar, tmp_func*Res);
                        } else {
                            nodes->SetSolution(iPoint, iVar, nodes->GetSolution_New(iPoint, iVar) + tmp_func*Res);
                        }
                    }
                        break;
                }
                
                /*--- Update residual information for current thread. ---*/
                resRMS[iVar] += Res*Res;
                if (fabs(Res) > resMax[iVar]) {
                    resMax[iVar] = fabs(Res);
                    idxMax[iVar] = iPoint;
                    coordMax[iVar] = geometry->node[iPoint]->GetCoord();
                }
            }
        }
    }
    if (!adjoint) {
        /*--- Reduce residual information over all threads in this rank. ---*/
        SU2_OMP_CRITICAL
        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
            AddRes_RMS(iVar, resRMS[iVar]);
            AddRes_Max(iVar, resMax[iVar], geometry->node[idxMax[iVar]]->GetGlobalIndex(), coordMax[iVar]);
        }
    }
    SU2_OMP_BARRIER
    
    SU2_OMP_MASTER
    {
        /*--- MPI solution ---*/
        
        InitiateComms(geometry, config, SOLUTION);
        CompleteComms(geometry, config, SOLUTION);
        
        /*--- Compute the root mean square residual ---*/
        
        SetResidual_RMS(geometry, config);
        
        /*--- For verification cases, compute the global error metrics. ---*/
        
        ComputeVerificationError(geometry, config);
    }
    SU2_OMP_BARRIER
    
}

void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
    
    Explicit_Iteration<RUNGE_KUTTA_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CEulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                          CConfig *config, unsigned short iRKStep) {
    
    Explicit_Iteration<CLASSICAL_RK4_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
    
    Explicit_Iteration<EULER_EXPLICIT>(geometry, solver_container, config, 0);
}

void CEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
    
    const bool adjoint = config->GetContinuous_Adjoint();
    const bool roe_turkel = config->GetKind_Upwind_Flow() == TURKEL;
    const bool low_mach_prec = config->Low_Mach_Preconditioning();
    
    /*--- Local matrix for preconditioning. ---*/
    su2double** LowMachPrec = nullptr;
    if (roe_turkel || low_mach_prec) {
        LowMachPrec = new su2double* [nVar];
        for(unsigned short iVar = 0; iVar < nVar; ++iVar)
            LowMachPrec[iVar] = new su2double [nVar];
    }
    
    /*--- Set shared residual variables to 0 and declare
     *    local ones for current thread to work on. ---*/
    
    SU2_OMP_MASTER
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
    }
    SU2_OMP_BARRIER
    
    su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
    const su2double* coordMax[MAXNVAR] = {nullptr};
    unsigned long idxMax[MAXNVAR] = {0};
    
    /*--- Build implicit system ---*/
    
    SU2_OMP(for schedule(static,omp_chunk_size) nowait)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
        /*--- Read the residual ---*/
        
        su2double* local_Res_TruncError = nodes->GetResTruncError(iPoint);
        
        /*--- Read the volume ---*/
        
        su2double Vol = geometry->node[iPoint]->GetVolume() + geometry->node[iPoint]->GetPeriodicVolume();
        
        /*--- Modify matrix diagonal to assure diagonal dominance ---*/
        
        if (nodes->GetDelta_Time(iPoint) != 0.0) {
            
            su2double Delta = Vol / nodes->GetDelta_Time(iPoint);
            
            if (roe_turkel || low_mach_prec) {
                SetPreconditioner(config, iPoint, Delta, LowMachPrec);
                Jacobian.AddBlock2Diag(iPoint, LowMachPrec);
            }
            else {
                Jacobian.AddVal2Diag(iPoint, Delta);
            }
        }
        else {
            Jacobian.SetVal2Diag(iPoint, 1.0);
            for (unsigned short iVar = 0; iVar < nVar; iVar++) {
                LinSysRes(iPoint,iVar) = 0.0;
                local_Res_TruncError[iVar] = 0.0;
            }
        }
        
        /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
        
        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
            unsigned long total_index = iPoint*nVar + iVar;
            LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
            LinSysSol[total_index] = 0.0;
            
            su2double Res = fabs(LinSysRes[total_index]);
            resRMS[iVar] += Res*Res;
            if (Res > resMax[iVar]) {
                resMax[iVar] = Res;
                idxMax[iVar] = iPoint;
                coordMax[iVar] = geometry->node[iPoint]->GetCoord();
            }
        }
    }
    SU2_OMP_CRITICAL
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        AddRes_RMS(iVar, resRMS[iVar]);
        AddRes_Max(iVar, resMax[iVar], geometry->node[idxMax[iVar]]->GetGlobalIndex(), coordMax[iVar]);
    }
    
    /*--- Initialize residual and solution at the ghost points ---*/
    
    SU2_OMP(sections nowait)
    {
        SU2_OMP(section)
        for (unsigned long iPoint = nPointDomain; iPoint < nPoint; iPoint++)
            LinSysRes.SetBlock_Zero(iPoint);
        
        SU2_OMP(section)
        for (unsigned long iPoint = nPointDomain; iPoint < nPoint; iPoint++)
            LinSysSol.SetBlock_Zero(iPoint);
    }
    
    /*--- Free local preconditioner. ---*/
    if (LowMachPrec) {
        for(unsigned short iVar = 0; iVar < nVar; ++iVar)
            delete [] LowMachPrec[iVar];
        delete [] LowMachPrec;
    }
    
    /*--- Solve or smooth the linear system. ---*/
    
    auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
    SU2_OMP_MASTER
    {
        SetIterLinSolver(iter);
        SetResLinSolver(System.GetResidual());
    }
    SU2_OMP_BARRIER
    
    
    ComputeUnderRelaxationFactor(solver_container, config);
    
    /*--- Update solution (system written in terms of increments) ---*/
    
    if (!adjoint) {
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
            for (unsigned short iVar = 0; iVar < nVar; iVar++) {
                nodes->AddSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint)*LinSysSol[iPoint*nVar+iVar]);
            }
        }
    }
    
    SU2_OMP_MASTER
    {
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
    SU2_OMP_BARRIER
    
}

void CEulerSolver::ComputeUnderRelaxationFactor(CSolver **solver_container, CConfig *config) {
    
    /* Loop over the solution update given by relaxing the linear
     system for this nonlinear iteration. */
    
    const su2double allowableRatio = 0.2;
    
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
        su2double localUnderRelaxation = 1.0;
        
        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
            
            /* We impose a limit on the maximum percentage that the
             density and energy can change over a nonlinear iteration. */
            
            if ((iVar == 0) || (iVar == nVar-1)) {
                const unsigned long index = iPoint*nVar + iVar;
                su2double ratio = fabs(LinSysSol[index])/(nodes->GetSolution(iPoint, iVar)+EPS);
                if (ratio > allowableRatio) {
                    localUnderRelaxation = min(allowableRatio/ratio, localUnderRelaxation);
                }
            }
        }
        
        /* In case of turbulence, take the min of the under-relaxation factor
         between the mean flow and the turb model. */
        
        if (config->GetKind_Turb_Model() != NONE)
            localUnderRelaxation = min(localUnderRelaxation, solver_container[TURB_SOL]->GetNodes()->GetUnderRelaxation(iPoint));
        
        /* Threshold the relaxation factor in the event that there is
         a very small value. This helps avoid catastrophic crashes due
         to non-realizable states by canceling the update. */
        
        if (localUnderRelaxation < 1e-10) localUnderRelaxation = 0.0;
        
        /* Store the under-relaxation factor for this point. */
        
        nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);
        
    }
    
}

void CEulerSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction) {
    
    const auto& primitives = nodes->GetPrimitive();
    auto& gradient = reconstruction? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
    
    computeGradientsGreenGauss(this, PRIMITIVE_GRADIENT, PERIODIC_PRIM_GG, *geometry,
                               *config, primitives, 0, nPrimVarGrad, gradient);
}

void CEulerSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction) {
    
    /*--- Set a flag for unweighted or weighted least-squares. ---*/
    bool weighted;
    
    if (reconstruction)
        weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
    else
        weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);
    
    const auto& primitives = nodes->GetPrimitive();
    auto& rmatrix = nodes->GetRmatrix();
    auto& gradient = reconstruction? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
    PERIODIC_QUANTITIES kindPeriodicComm = weighted? PERIODIC_PRIM_LS : PERIODIC_PRIM_ULS;
    
    computeGradientsLeastSquares(this, PRIMITIVE_GRADIENT, kindPeriodicComm, *geometry, *config,
                                 weighted, primitives, 0, nPrimVarGrad, gradient, rmatrix);
}

void CEulerSolver::SetPrimitive_Limiter(CGeometry *geometry, CConfig *config) {
    
    auto kindLimiter = static_cast<ENUM_LIMITER>(config->GetKind_SlopeLimit_Flow());
    const auto& primitives = nodes->GetPrimitive();
    const auto& gradient = nodes->GetGradient_Reconstruction();
    auto& primMin = nodes->GetSolution_Min();
    auto& primMax = nodes->GetSolution_Max();
    auto& limiter = nodes->GetLimiter_Primitive();
    
    computeLimiters(kindLimiter, this, PRIMITIVE_LIMITER, PERIODIC_LIM_PRIM_1, PERIODIC_LIM_PRIM_2,
                    *geometry, *config, 0, nPrimVarGrad, primitives, gradient, primMin, primMax, limiter);
}

void CEulerSolver::SetPreconditioner(const CConfig *config, unsigned long iPoint,
                                     su2double delta, su2double** preconditioner) const {
    
    unsigned short iDim, jDim, iVar, jVar;
    su2double local_Mach, rho, enthalpy, soundspeed, sq_vel;
    su2double *U_i = NULL;
    su2double Beta_max = config->GetmaxTurkelBeta();
    su2double Mach_infty2, Mach_lim2, aux, parameter;
    
    /*--- Variables to calculate the preconditioner parameter Beta ---*/
    local_Mach = sqrt(nodes->GetVelocity2(iPoint))/nodes->GetSoundSpeed(iPoint);
    
    /*--- Weiss and Smith Preconditioning---*/
    Mach_infty2 = pow(config->GetMach(),2.0);
    Mach_lim2 = pow(0.00001,2.0);
    aux = max(pow(local_Mach,2.0),Mach_lim2);
    parameter = min(1.0, max(aux,Beta_max*Mach_infty2));
    
    U_i = nodes->GetSolution(iPoint);
    
    rho = U_i[0];
    enthalpy = nodes->GetEnthalpy(iPoint);
    soundspeed = nodes->GetSoundSpeed(iPoint);
    sq_vel = nodes->GetVelocity2(iPoint);
    
    /*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
    preconditioner[0][0] = 0.5*sq_vel;
    preconditioner[0][nVar-1] = 1.0;
    for (iDim = 0; iDim < nDim; iDim ++)
        preconditioner[0][1+iDim] = -1.0*U_i[iDim+1]/rho;
    
    for (iDim = 0; iDim < nDim; iDim ++) {
        preconditioner[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
        preconditioner[iDim+1][nVar-1] = U_i[iDim+1]/rho;
        for (jDim = 0; jDim < nDim; jDim ++) {
            preconditioner[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
        }
    }
    
    preconditioner[nVar-1][0] = 0.5*sq_vel*enthalpy;
    preconditioner[nVar-1][nVar-1] = enthalpy;
    for (iDim = 0; iDim < nDim; iDim ++)
        preconditioner[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;
    
    
    for (iVar = 0; iVar < nVar; iVar ++ ) {
        for (jVar = 0; jVar < nVar; jVar ++ ) {
            preconditioner[iVar][jVar] = (parameter - 1.0) * ((Gamma-1.0)/(soundspeed*soundspeed))*preconditioner[iVar][jVar];
            if (iVar == jVar)
                preconditioner[iVar][iVar] += 1.0;
            
            preconditioner[iVar][jVar] *= delta;
        }
    }
    
}

void CEulerSolver::GetPower_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {
    
    unsigned short iDim, iMarker, jMarker;
    unsigned long iVertex, iPoint;
    su2double  *V_inlet = NULL, *V_outlet = NULL, Pressure, Temperature, Velocity[3], Vn,
    Velocity2, Density, Area, SoundSpeed, TotalPressure, Vel_Infty2, RamDrag,
    TotalTemperature, VelocityJet,
    Vel_Infty, MaxPressure, MinPressure, MFR, InfVel2;
    unsigned short iMarker_Inlet, iMarker_Outlet, nMarker_Inlet, nMarker_Outlet;
    string Inlet_TagBound, Outlet_TagBound;
    su2double DeltaPress = 0.0, DeltaTemp = 0.0, TotalPressRatio = 0.0, TotalTempRatio = 0.0, StaticPressRatio = 0.0, StaticTempRatio = 0.0,
    NetThrust = 0.0, GrossThrust = 0.0, Power = 0.0, MassFlow = 0.0, Mach = 0.0, Force = 0.0;
    bool ReverseFlow, Engine = false, Pair = true;
    
    su2double Gas_Constant = config->GetGas_ConstantND();
    su2double Cp = Gas_Constant*Gamma / (Gamma-1.0);
    su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
    su2double Beta = config->GetAoS()*PI_NUMBER/180.0;
    bool write_heads = ((((config->GetInnerIter() % (config->GetWrt_Con_Freq()*40)) == 0) && (config->GetInnerIter()!= 0)) || (config->GetInnerIter() == 1));
    bool Evaluate_BC = ((((config->GetInnerIter() % (config->GetWrt_Con_Freq()*40)) == 0)) || (config->GetInnerIter() == 1) || (config->GetDiscrete_Adjoint()));
    
    if ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0)) Engine = true;
    if ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0)) Engine = false;
    if ((config->GetnMarker_EngineInflow()) != (config->GetnMarker_EngineExhaust())) Pair = false;
    
    if (Engine) { nMarker_Inlet  = config->GetnMarker_EngineInflow(); nMarker_Outlet = config->GetnMarker_EngineExhaust(); }
    else  { nMarker_Inlet = config->GetnMarker_ActDiskInlet(); nMarker_Outlet  = config->GetnMarker_ActDiskOutlet(); }
    
    /*--- Evaluate the MPI for the actuator disk IO ---*/
    
    if (Evaluate_BC) {
        
        su2double *Inlet_MassFlow         = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_ReverseMassFlow  = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_Pressure         = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_Mach             = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_MaxPressure      = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_MinPressure      = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_TotalPressure    = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_Temperature      = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_TotalTemperature = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_Area             = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_RamDrag          = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_Force            = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_Power            = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_XCG              = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_YCG              = new su2double [config->GetnMarker_All()]();
        su2double *Inlet_ZCG              = new su2double [config->GetnMarker_All()]();
        
        su2double *Outlet_MassFlow         = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_Pressure         = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_TotalPressure    = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_Temperature      = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_TotalTemperature = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_Area             = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_GrossThrust      = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_Force            = new su2double [config->GetnMarker_All()]();
        su2double *Outlet_Power            = new su2double [config->GetnMarker_All()]();
        
        /*--- Comute MassFlow, average temp, press, etc. ---*/
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
            
            MinPressure = +1E10; MaxPressure = -1E10;
            
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) {
                
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    
                    if (geometry->node[iPoint]->GetDomain()) {
                        
                        V_inlet = nodes->GetPrimitive(iPoint);
                        
                        geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
                        
                        Temperature = V_inlet[0];
                        Pressure = V_inlet[nDim+1];
                        
                        Density = V_inlet[nDim+2];
                        SoundSpeed = sqrt(Gamma*Pressure/Density);
                        
                        Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vel_Infty2 =0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Area += Vector[iDim]*Vector[iDim];
                            Velocity[iDim] = V_inlet[iDim+1];
                            Velocity2 += Velocity[iDim]*Velocity[iDim];
                            Vel_Infty2 += GetVelocity_Inf(iDim)*GetVelocity_Inf(iDim);
                            MassFlow -= Vector[iDim]*Velocity[iDim]*Density;
                        }
                        
                        Vn = 0.0; ReverseFlow = false;
                        for (iDim = 0; iDim < nDim; iDim++) {  Vn -= Velocity[iDim]*Vector[iDim]/Area; }
                        if (Vn < 0.0) { ReverseFlow = true; }
                        
                        Vel_Infty = sqrt (Vel_Infty2);
                        Area = sqrt (Area);
                        Mach = sqrt(Velocity2)/SoundSpeed;
                        TotalPressure = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));
                        TotalTemperature = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
                        MinPressure = min(MinPressure, TotalPressure);
                        MaxPressure = max(MaxPressure, TotalPressure);
                        
                        RamDrag = MassFlow * Vel_Infty;
                        
                        Inlet_MassFlow[iMarker]         += MassFlow;
                        Inlet_Pressure[iMarker]         += Pressure*MassFlow;
                        Inlet_Mach[iMarker]             += Mach*MassFlow;
                        Inlet_MinPressure[iMarker]       = min (MinPressure, Inlet_MinPressure[iMarker]);
                        Inlet_MaxPressure[iMarker]       = max(MaxPressure, Inlet_MaxPressure[iMarker]);
                        Inlet_TotalPressure[iMarker]    += TotalPressure*MassFlow;
                        Inlet_Temperature[iMarker]      += Temperature*MassFlow;
                        Inlet_TotalTemperature[iMarker] += TotalTemperature*MassFlow;
                        Inlet_Area[iMarker]             += Area;
                        Inlet_RamDrag[iMarker]          += RamDrag;
                        Inlet_Power[iMarker] += MassFlow*Cp*TotalTemperature;
                        if (ReverseFlow) Inlet_ReverseMassFlow[iMarker]  += MassFlow;
                        
                        su2double Inlet_ForceX = -(Pressure - Pressure_Inf)*Vector[0] + MassFlow*Velocity[0];
                        su2double Inlet_ForceY = -(Pressure - Pressure_Inf)*Vector[1] + MassFlow*Velocity[1];
                        su2double Inlet_ForceZ = 0.0;
                        if (nDim == 3) Inlet_ForceZ = -(Pressure - Pressure_Inf)*Vector[2] + MassFlow*Velocity[2];
                        Inlet_Force[iMarker] +=  Inlet_ForceX*cos(Alpha)*cos(Beta) + Inlet_ForceY*sin(Beta) +Inlet_ForceZ*sin(Alpha)*cos(Beta);
                        
                        Inlet_XCG[iMarker] += geometry->node[iPoint]->GetCoord(0)*Area;
                        Inlet_YCG[iMarker] += geometry->node[iPoint]->GetCoord(1)*Area;
                        if (nDim == 3) Inlet_ZCG[iMarker] += geometry->node[iPoint]->GetCoord(2)*Area;
                        
                    }
                }
                
            }
            
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) ||
                (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST)) {
                
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    
                    if (geometry->node[iPoint]->GetDomain()) {
                        
                        V_outlet = nodes->GetPrimitive(iPoint);
                        
                        geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
                        
                        Temperature = V_outlet[0];
                        Pressure = V_outlet[nDim+1];
                        
                        Density = V_outlet[nDim+2];
                        SoundSpeed  = sqrt(Gamma*Pressure/Density);
                        
                        Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vel_Infty2 = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Area += Vector[iDim]*Vector[iDim];
                            Velocity[iDim] = V_outlet[iDim+1];
                            Velocity2 += Velocity[iDim]*Velocity[iDim];
                            Vel_Infty2 += GetVelocity_Inf(iDim)*GetVelocity_Inf(iDim);
                            MassFlow += Vector[iDim]*Velocity[iDim]*Density;
                        }
                        
                        Vel_Infty = sqrt (Vel_Infty2);
                        Area = sqrt (Area);
                        Mach = sqrt(Velocity2)/SoundSpeed;
                        TotalPressure = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
                        TotalTemperature = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
                        VelocityJet = sqrt(Velocity2) ;
                        
                        GrossThrust = MassFlow * VelocityJet;
                        
                        Outlet_MassFlow[iMarker] += MassFlow;
                        Outlet_Pressure[iMarker] += Pressure*MassFlow;
                        Outlet_TotalPressure[iMarker] += TotalPressure*MassFlow;
                        Outlet_Temperature[iMarker] += Temperature*MassFlow;
                        Outlet_TotalTemperature[iMarker] += TotalTemperature*MassFlow;
                        Outlet_Area[iMarker] += Area;
                        Outlet_GrossThrust[iMarker] += GrossThrust;
                        Outlet_Power[iMarker] += MassFlow*Cp*TotalTemperature;
                        
                        su2double Outlet_ForceX = -(Pressure - Pressure_Inf)*Vector[0] -MassFlow*Velocity[0];
                        su2double Outlet_ForceY = -(Pressure - Pressure_Inf)*Vector[1] -MassFlow*Velocity[1];
                        su2double Outlet_ForceZ = 0.0;
                        if (nDim == 3) Outlet_ForceZ = -(Pressure - Pressure_Inf)*Vector[2] -MassFlow*Velocity[2];
                        
                        if (nDim == 2) Outlet_Force[iMarker] +=  Outlet_ForceX*cos(Alpha) + Outlet_ForceY*sin(Alpha);
                        if (nDim == 3) Outlet_Force[iMarker] +=  Outlet_ForceX*cos(Alpha)*cos(Beta) + Outlet_ForceY*sin(Beta) + Outlet_ForceZ*sin(Alpha)*cos(Beta);
                        
                    }
                }
                
            }
            
        }
        
        /*--- Copy to the appropriate structure ---*/
        
        su2double *Inlet_MassFlow_Local             = new su2double [nMarker_Inlet]();
        su2double *Inlet_ReverseMassFlow_Local      = new su2double [nMarker_Inlet]();
        su2double *Inlet_Temperature_Local          = new su2double [nMarker_Inlet]();
        su2double *Inlet_TotalTemperature_Local     = new su2double [nMarker_Inlet]();
        su2double *Inlet_Pressure_Local             = new su2double [nMarker_Inlet]();
        su2double *Inlet_Mach_Local                 = new su2double [nMarker_Inlet]();
        su2double *Inlet_MinPressure_Local          = new su2double [nMarker_Inlet]();
        su2double *Inlet_MaxPressure_Local          = new su2double [nMarker_Inlet]();
        su2double *Inlet_Power_Local                = new su2double [nMarker_Inlet]();
        su2double *Inlet_TotalPressure_Local        = new su2double [nMarker_Inlet]();
        su2double *Inlet_RamDrag_Local              = new su2double [nMarker_Inlet]();
        su2double *Inlet_Force_Local                = new su2double [nMarker_Inlet]();
        su2double *Inlet_Area_Local                 = new su2double [nMarker_Inlet]();
        su2double *Inlet_XCG_Local                  = new su2double [nMarker_Inlet]();
        su2double *Inlet_YCG_Local                  = new su2double [nMarker_Inlet]();
        su2double *Inlet_ZCG_Local                  = new su2double [nMarker_Inlet]();
        
        su2double *Inlet_MassFlow_Total             = new su2double [nMarker_Inlet]();
        su2double *Inlet_ReverseMassFlow_Total      = new su2double [nMarker_Inlet]();
        su2double *Inlet_Pressure_Total             = new su2double [nMarker_Inlet]();
        su2double *Inlet_Mach_Total                 = new su2double [nMarker_Inlet]();
        su2double *Inlet_MinPressure_Total          = new su2double [nMarker_Inlet]();
        su2double *Inlet_MaxPressure_Total          = new su2double [nMarker_Inlet]();
        su2double *Inlet_Power_Total                = new su2double [nMarker_Inlet]();
        su2double *Inlet_TotalPressure_Total        = new su2double [nMarker_Inlet]();
        su2double *Inlet_Temperature_Total          = new su2double [nMarker_Inlet]();
        su2double *Inlet_TotalTemperature_Total     = new su2double [nMarker_Inlet]();
        su2double *Inlet_RamDrag_Total              = new su2double [nMarker_Inlet]();
        su2double *Inlet_Force_Total                = new su2double [nMarker_Inlet]();
        su2double *Inlet_Area_Total                 = new su2double [nMarker_Inlet]();
        su2double *Inlet_XCG_Total                  = new su2double [nMarker_Inlet]();
        su2double *Inlet_YCG_Total                  = new su2double [nMarker_Inlet]();
        su2double *Inlet_ZCG_Total                  = new su2double [nMarker_Inlet]();
        
        su2double *Outlet_MassFlow_Local            = new su2double [nMarker_Outlet]();
        su2double *Outlet_Pressure_Local            = new su2double [nMarker_Outlet]();
        su2double *Outlet_TotalPressure_Local       = new su2double [nMarker_Outlet]();
        su2double *Outlet_Temperature_Local         = new su2double [nMarker_Outlet]();
        su2double *Outlet_TotalTemperature_Local    = new su2double [nMarker_Outlet]();
        su2double *Outlet_GrossThrust_Local         = new su2double [nMarker_Outlet]();
        su2double *Outlet_Force_Local               = new su2double [nMarker_Outlet]();
        su2double *Outlet_Power_Local               = new su2double [nMarker_Outlet]();
        su2double *Outlet_Area_Local                = new su2double [nMarker_Outlet]();
        
        su2double *Outlet_MassFlow_Total            = new su2double [nMarker_Outlet]();
        su2double *Outlet_Pressure_Total            = new su2double [nMarker_Outlet]();
        su2double *Outlet_TotalPressure_Total       = new su2double [nMarker_Outlet]();
        su2double *Outlet_Temperature_Total         = new su2double [nMarker_Outlet]();
        su2double *Outlet_TotalTemperature_Total    = new su2double [nMarker_Outlet]();
        su2double *Outlet_GrossThrust_Total         = new su2double [nMarker_Outlet]();
        su2double *Outlet_Force_Total               = new su2double [nMarker_Outlet]();
        su2double *Outlet_Power_Total               = new su2double [nMarker_Outlet]();
        su2double *Outlet_Area_Total                = new su2double [nMarker_Outlet]();
        
        /*--- Copy the values to the local array for MPI ---*/
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
            
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||  (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) {
                for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
                    
                    if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
                    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW) Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
                    
                    if (config->GetMarker_All_TagBound(iMarker) == Inlet_TagBound) {
                        Inlet_MassFlow_Local[iMarker_Inlet]             += Inlet_MassFlow[iMarker];
                        Inlet_ReverseMassFlow_Local[iMarker_Inlet]      += Inlet_ReverseMassFlow[iMarker];
                        Inlet_Pressure_Local[iMarker_Inlet]             += Inlet_Pressure[iMarker];
                        Inlet_Mach_Local[iMarker_Inlet]                 += Inlet_Mach[iMarker];
                        Inlet_MinPressure_Local[iMarker_Inlet]          += Inlet_MinPressure[iMarker];
                        Inlet_MaxPressure_Local[iMarker_Inlet]          += Inlet_MaxPressure[iMarker];
                        Inlet_TotalPressure_Local[iMarker_Inlet]        += Inlet_TotalPressure[iMarker];
                        Inlet_Temperature_Local[iMarker_Inlet]          += Inlet_Temperature[iMarker];
                        Inlet_TotalTemperature_Local[iMarker_Inlet]     += Inlet_TotalTemperature[iMarker];
                        Inlet_RamDrag_Local[iMarker_Inlet]              += Inlet_RamDrag[iMarker];
                        Inlet_Force_Local[iMarker_Inlet]                += Inlet_Force[iMarker];
                        Inlet_Power_Local[iMarker_Inlet]                += Inlet_Power[iMarker];
                        Inlet_Area_Local[iMarker_Inlet]                 += Inlet_Area[iMarker];
                        Inlet_XCG_Local[iMarker_Inlet]                  += Inlet_XCG[iMarker];
                        Inlet_YCG_Local[iMarker_Inlet]                  += Inlet_YCG[iMarker];
                        if (nDim == 3) Inlet_ZCG_Local[iMarker_Inlet]   += Inlet_ZCG[iMarker];
                    }
                    
                }
            }
            
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) || (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST)) {
                for (iMarker_Outlet= 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
                    
                    if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET) Outlet_TagBound = config->GetMarker_ActDiskOutlet_TagBound(iMarker_Outlet);
                    if (config->GetMarker_All_KindBC(iMarker) == ENGINE_EXHAUST) Outlet_TagBound = config->GetMarker_EngineExhaust_TagBound(iMarker_Outlet);
                    
                    if (config->GetMarker_All_TagBound(iMarker) == Outlet_TagBound) {
                        Outlet_MassFlow_Local[iMarker_Outlet]               += Outlet_MassFlow[iMarker];
                        Outlet_Pressure_Local[iMarker_Outlet]               += Outlet_Pressure[iMarker];
                        Outlet_TotalPressure_Local[iMarker_Outlet]          += Outlet_TotalPressure[iMarker];
                        Outlet_Temperature_Local[iMarker_Outlet]            += Outlet_Temperature[iMarker];
                        Outlet_TotalTemperature_Local[iMarker_Outlet]       += Outlet_TotalTemperature[iMarker];
                        Outlet_GrossThrust_Local[iMarker_Outlet]            += Outlet_GrossThrust[iMarker];
                        Outlet_Force_Local[iMarker_Outlet]                  += Outlet_Force[iMarker];
                        Outlet_Power_Local[iMarker_Outlet]                  += Outlet_Power[iMarker];
                        Outlet_Area_Local[iMarker_Outlet]                   += Outlet_Area[iMarker];
                    }
                    
                }
            }
            
        }
        
        /*--- Correct the min max values for the MPI ---*/
        
        bool ActDisk  = false;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                (config->GetMarker_All_KindBC(iMarker) == ENGINE_INFLOW)) { ActDisk  = true; break; }
        }
        
        if (!ActDisk) {
            for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
                Inlet_MinPressure_Local[iMarker_Inlet]      =  1E10;
                Inlet_MaxPressure_Local[iMarker_Inlet]      = -1E10;
            }
        }
        
        /*--- All the ranks to compute the total value ---*/
        
        SU2_MPI::Allreduce(Inlet_MassFlow_Local, Inlet_MassFlow_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_ReverseMassFlow_Local, Inlet_ReverseMassFlow_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_Pressure_Local, Inlet_Pressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_Mach_Local, Inlet_Mach_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_MinPressure_Local, Inlet_MinPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_MaxPressure_Local, Inlet_MaxPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_TotalPressure_Local, Inlet_TotalPressure_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_Temperature_Local, Inlet_Temperature_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_TotalTemperature_Local, Inlet_TotalTemperature_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_RamDrag_Local, Inlet_RamDrag_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_Force_Local, Inlet_Force_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_Power_Local, Inlet_Power_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_Area_Local, Inlet_Area_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_XCG_Local, Inlet_XCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Inlet_YCG_Local, Inlet_YCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (nDim == 3) SU2_MPI::Allreduce(Inlet_ZCG_Local, Inlet_ZCG_Total, nMarker_Inlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        SU2_MPI::Allreduce(Outlet_MassFlow_Local, Outlet_MassFlow_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_Pressure_Local, Outlet_Pressure_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_TotalPressure_Local, Outlet_TotalPressure_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_Temperature_Local, Outlet_Temperature_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_TotalTemperature_Local, Outlet_TotalTemperature_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_GrossThrust_Local, Outlet_GrossThrust_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_Force_Local, Outlet_Force_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_Power_Local, Outlet_Power_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(Outlet_Area_Local, Outlet_Area_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        /*--- Compute the value of the average surface temperature and pressure and
         set the value in the config structure for future use ---*/
        
        for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
            if (Inlet_Area_Total[iMarker_Inlet] != 0.0) {
                Inlet_Pressure_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
                Inlet_Mach_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
                Inlet_TotalPressure_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
                Inlet_Temperature_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
                Inlet_TotalTemperature_Total[iMarker_Inlet] /= Inlet_MassFlow_Total[iMarker_Inlet];
                Inlet_XCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
                Inlet_YCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
                if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet] /= Inlet_Area_Total[iMarker_Inlet];
            }
            else {
                Inlet_Pressure_Total[iMarker_Inlet] = 0.0;
                Inlet_Mach_Total[iMarker_Inlet] = 0.0;
                Inlet_TotalPressure_Total[iMarker_Inlet] = 0.0;
                Inlet_Temperature_Total[iMarker_Inlet] = 0.0;
                Inlet_TotalTemperature_Total[iMarker_Inlet] = 0.0;
                Inlet_XCG_Total[iMarker_Inlet] = 0.0;
                Inlet_YCG_Total[iMarker_Inlet] = 0.0;
                if (nDim == 3) Inlet_ZCG_Total[iMarker_Inlet] = 0.0;
            }
            
            if (iMesh == MESH_0) {
                
                if (Engine) {
                    config->SetInflow_MassFlow(iMarker_Inlet, Inlet_MassFlow_Total[iMarker_Inlet]);
                    config->SetInflow_ReverseMassFlow(iMarker_Inlet, Inlet_ReverseMassFlow_Total[iMarker_Inlet]);
                    config->SetInflow_Pressure(iMarker_Inlet, Inlet_Pressure_Total[iMarker_Inlet]);
                    config->SetInflow_TotalPressure(iMarker_Inlet, Inlet_TotalPressure_Total[iMarker_Inlet]);
                    config->SetInflow_Temperature(iMarker_Inlet, Inlet_Temperature_Total[iMarker_Inlet]);
                    config->SetInflow_TotalTemperature(iMarker_Inlet, Inlet_TotalTemperature_Total[iMarker_Inlet]);
                    config->SetInflow_RamDrag(iMarker_Inlet, Inlet_RamDrag_Total[iMarker_Inlet]);
                    config->SetInflow_Force(iMarker_Inlet, Inlet_Force_Total[iMarker_Inlet]);
                    config->SetInflow_Power(iMarker_Inlet, Inlet_Power_Total[iMarker_Inlet]);
                }
                else {
                    config->SetActDiskInlet_MassFlow(iMarker_Inlet, Inlet_MassFlow_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_ReverseMassFlow(iMarker_Inlet, Inlet_ReverseMassFlow_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_Pressure(iMarker_Inlet, Inlet_Pressure_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_TotalPressure(iMarker_Inlet, Inlet_TotalPressure_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_Temperature(iMarker_Inlet, Inlet_Temperature_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_TotalTemperature(iMarker_Inlet, Inlet_TotalTemperature_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_RamDrag(iMarker_Inlet, Inlet_RamDrag_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_Force(iMarker_Inlet, Inlet_Force_Total[iMarker_Inlet]);
                    config->SetActDiskInlet_Power(iMarker_Inlet, Inlet_Power_Total[iMarker_Inlet]);
                }
                
            }
            
        }
        
        for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
            if (Outlet_Area_Total[iMarker_Outlet] != 0.0) {
                Outlet_Pressure_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
                Outlet_TotalPressure_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
                Outlet_Temperature_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
                Outlet_TotalTemperature_Total[iMarker_Outlet] /= Outlet_MassFlow_Total[iMarker_Outlet];
            }
            else {
                Outlet_Pressure_Total[iMarker_Outlet] = 0.0;
                Outlet_TotalPressure_Total[iMarker_Outlet] = 0.0;
                Outlet_Temperature_Total[iMarker_Outlet] = 0.0;
                Outlet_TotalTemperature_Total[iMarker_Outlet] = 0.0;
            }
            
            if (iMesh == MESH_0) {
                
                if (Engine) {
                    config->SetExhaust_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
                    config->SetExhaust_Pressure(iMarker_Outlet, Outlet_Pressure_Total[iMarker_Outlet]);
                    config->SetExhaust_TotalPressure(iMarker_Outlet, Outlet_TotalPressure_Total[iMarker_Outlet]);
                    config->SetExhaust_Temperature(iMarker_Outlet, Outlet_Temperature_Total[iMarker_Outlet]);
                    config->SetExhaust_TotalTemperature(iMarker_Outlet, Outlet_TotalTemperature_Total[iMarker_Outlet]);
                    config->SetExhaust_GrossThrust(iMarker_Outlet, Outlet_GrossThrust_Total[iMarker_Outlet]);
                    config->SetExhaust_Force(iMarker_Outlet, Outlet_Force_Total[iMarker_Outlet]);
                    config->SetExhaust_Power(iMarker_Outlet, Outlet_Power_Total[iMarker_Outlet]);
                }
                else {
                    config->SetActDiskOutlet_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
                    config->SetActDiskOutlet_Pressure(iMarker_Outlet, Outlet_Pressure_Total[iMarker_Outlet]);
                    config->SetActDiskOutlet_TotalPressure(iMarker_Outlet, Outlet_TotalPressure_Total[iMarker_Outlet]);
                    config->SetActDiskOutlet_Temperature(iMarker_Outlet, Outlet_Temperature_Total[iMarker_Outlet]);
                    config->SetActDiskOutlet_TotalTemperature(iMarker_Outlet, Outlet_TotalTemperature_Total[iMarker_Outlet]);
                    config->SetActDiskOutlet_GrossThrust(iMarker_Outlet, Outlet_GrossThrust_Total[iMarker_Outlet]);
                    config->SetActDiskOutlet_Force(iMarker_Outlet, Outlet_Force_Total[iMarker_Outlet]);
                    config->SetActDiskOutlet_Power(iMarker_Outlet, Outlet_Power_Total[iMarker_Outlet]);
                }
                
            }
            
        }
        
        
        if (Pair) {
            
            /*--- Store delta pressure, temperature, thrust, and area ---*/
            
            for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
                
                if (Engine) {
                    Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
                    jMarker = config->GetMarker_CfgFile_EngineExhaust(Inlet_TagBound);
                    Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
                }
                else {
                    Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
                    jMarker = config->GetMarker_CfgFile_ActDiskOutlet(Inlet_TagBound);
                    Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
                }
                
                
                su2double DeltaPress = 0.0, DeltaTemp = 0.0, NetThrust = 0.0, GrossThrust = 0.0, TotalPressRatio = 0.0, TotalTempRatio = 0.0, StaticPressRatio = 0.0, StaticTempRatio = 0.0;
                
                if (Engine) {
                    DeltaPress   = config->GetExhaust_Pressure(Outlet_TagBound) - config->GetInflow_Pressure(Inlet_TagBound);
                    DeltaTemp         = config->GetExhaust_Temperature(Outlet_TagBound) - config->GetInflow_Temperature(Inlet_TagBound);
                    NetThrust    = config->GetExhaust_GrossThrust(Outlet_TagBound) - config->GetInflow_RamDrag(Inlet_TagBound);
                    GrossThrust   = config->GetExhaust_GrossThrust(Outlet_TagBound);
                    TotalPressRatio   = config->GetExhaust_TotalPressure(Outlet_TagBound)/config->GetInflow_TotalPressure(Inlet_TagBound);
                    TotalTempRatio    = config->GetExhaust_TotalTemperature(Outlet_TagBound)/config->GetInflow_TotalTemperature(Inlet_TagBound);
                    StaticPressRatio   = config->GetExhaust_Pressure(Outlet_TagBound)/config->GetInflow_Pressure(Inlet_TagBound);
                    StaticTempRatio    = config->GetExhaust_Temperature(Outlet_TagBound)/config->GetInflow_Temperature(Inlet_TagBound);
                    Force = config->GetInflow_Force(Inlet_TagBound) + config->GetExhaust_Force(Outlet_TagBound);
                    Power = config->GetExhaust_Power(Outlet_TagBound) - config->GetInflow_Power(Inlet_TagBound);
                }
                else {
                    DeltaPress   = config->GetActDiskOutlet_Pressure(Outlet_TagBound) - config->GetActDiskInlet_Pressure(Inlet_TagBound);
                    DeltaTemp         = config->GetActDiskOutlet_Temperature(Outlet_TagBound) - config->GetActDiskInlet_Temperature(Inlet_TagBound);
                    NetThrust    = config->GetActDiskOutlet_GrossThrust(Outlet_TagBound) - config->GetActDiskInlet_RamDrag(Inlet_TagBound);
                    GrossThrust   = config->GetActDiskOutlet_GrossThrust(Outlet_TagBound);
                    TotalPressRatio   = config->GetActDiskOutlet_TotalPressure(Outlet_TagBound)/config->GetActDiskInlet_TotalPressure(Inlet_TagBound);
                    TotalTempRatio    = config->GetActDiskOutlet_TotalTemperature(Outlet_TagBound)/config->GetActDiskInlet_TotalTemperature(Inlet_TagBound);
                    StaticPressRatio   = config->GetActDiskOutlet_Pressure(Outlet_TagBound)/config->GetActDiskInlet_Pressure(Inlet_TagBound);
                    StaticTempRatio    = config->GetActDiskOutlet_Temperature(Outlet_TagBound)/config->GetActDiskInlet_Temperature(Inlet_TagBound);
                    Force = config->GetActDiskInlet_Force(Inlet_TagBound) + config->GetActDiskOutlet_Force(Outlet_TagBound);
                    Power =  config->GetActDiskOutlet_Power(Outlet_TagBound) - config->GetActDiskInlet_Power(Inlet_TagBound);
                    MassFlow =  config->GetActDiskInlet_MassFlow(Inlet_TagBound);
                }
                
                Mach        = Inlet_Mach_Total[iMarker_Inlet];
                Area              = Inlet_Area_Total[iMarker_Inlet];
                
                if (Engine) {
                    config->SetEngine_Mach(iMarker_Inlet, Mach);
                    config->SetEngine_Force(iMarker_Inlet, Force);
                    config->SetEngine_Power(iMarker_Inlet, Power);
                    config->SetEngine_NetThrust(iMarker_Inlet, NetThrust);
                    config->SetEngine_GrossThrust(iMarker_Inlet, GrossThrust);
                    config->SetEngine_Area(iMarker_Inlet, Area);
                }
                else {
                    config->SetActDisk_DeltaPress(iMarker_Inlet, DeltaPress);
                    config->SetActDisk_DeltaTemp(iMarker_Inlet, DeltaTemp);
                    config->SetActDisk_Mach(iMarker_Inlet, Mach);
                    config->SetActDisk_Force(iMarker_Inlet, Force);
                    config->SetActDisk_Power(iMarker_Inlet, Power);
                    config->SetActDisk_MassFlow(iMarker_Inlet, MassFlow);
                    config->SetActDisk_TotalPressRatio(iMarker_Inlet, TotalPressRatio);
                    config->SetActDisk_TotalTempRatio(iMarker_Inlet, TotalTempRatio);
                    config->SetActDisk_StaticPressRatio(iMarker_Inlet, StaticPressRatio);
                    config->SetActDisk_StaticTempRatio(iMarker_Inlet, StaticTempRatio);
                    config->SetActDisk_NetThrust(iMarker_Inlet, NetThrust);
                    config->SetActDisk_GrossThrust(iMarker_Inlet, GrossThrust);
                    config->SetActDisk_Area(iMarker_Inlet, Area);
                }
                
            }
            
            /*--- Screen output using the values already stored in the config container ---*/
            
            if ((rank == MASTER_NODE) && (iMesh == MESH_0) ) {
                
                cout.precision(5);
                cout.setf(ios::fixed, ios::floatfield);
                
                if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
                    if (Engine) cout << endl   << "---------------------------- Engine properties --------------------------" << endl;
                    else cout << endl   << "------------------------ Actuator Disk properties -----------------------" << endl;
                }
                
                for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
                    
                    if (Engine) {
                        Inlet_TagBound = config->GetMarker_EngineInflow_TagBound(iMarker_Inlet);
                        jMarker = config->GetMarker_CfgFile_EngineExhaust(Inlet_TagBound);
                        Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
                    }
                    else {
                        Inlet_TagBound = config->GetMarker_ActDiskInlet_TagBound(iMarker_Inlet);
                        jMarker = config->GetMarker_CfgFile_ActDiskOutlet(Inlet_TagBound);
                        Outlet_TagBound = config->GetMarker_CfgFile_TagBound(jMarker);
                    }
                    
                    
                    if (Engine) {
                        NetThrust             =  config->GetEngine_NetThrust(iMarker_Inlet);
                        GrossThrust   = config->GetEngine_GrossThrust(iMarker_Inlet);
                        Power               = config->GetEngine_Power(iMarker_Inlet);
                        Mach                  = config->GetEngine_Mach(iMarker_Inlet);
                        Force               = config->GetEngine_Force(iMarker_Inlet);
                    }
                    else {
                        DeltaPress      = config->GetActDisk_DeltaPress(iMarker_Inlet);
                        DeltaTemp         = config->GetActDisk_DeltaTemp(iMarker_Inlet);
                        TotalPressRatio       = config->GetActDisk_TotalPressRatio(iMarker_Inlet);
                        TotalTempRatio        = config->GetActDisk_TotalTempRatio(iMarker_Inlet);
                        StaticPressRatio      = config->GetActDisk_StaticPressRatio(iMarker_Inlet);
                        StaticTempRatio           = config->GetActDisk_StaticTempRatio(iMarker_Inlet);
                        NetThrust             =  config->GetActDisk_NetThrust(iMarker_Inlet);
                        GrossThrust   = config->GetActDisk_GrossThrust(iMarker_Inlet);
                        Power               = config->GetActDisk_Power(iMarker_Inlet);
                        Mach                  = config->GetActDisk_Mach(iMarker_Inlet);
                        Force               = config->GetActDisk_Force(iMarker_Inlet);
                    }
                    
                    su2double Mach_Inf                  = config->GetMach();
                    su2double Pressure_Inf  = config->GetPressure_FreeStreamND();
                    
                    su2double TotalPressure_Inf  = Pressure_Inf * pow( 1.0 + Mach_Inf * Mach_Inf * 0.5 * (Gamma - 1.0), Gamma     / (Gamma - 1.0));
                    
                    su2double MinPressure = Inlet_MinPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
                    su2double MaxPressure = Inlet_MaxPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
                    su2double AvePressure = Inlet_TotalPressure_Total[iMarker_Inlet]/TotalPressure_Inf;
                    
                    su2double RefDensity  = Density_Inf;
                    su2double RefArea     = config->GetRefArea();
                    su2double RefVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
                    
                    su2double Factor = (0.5*RefDensity*RefArea*RefVel2);
                    su2double Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
                    su2double DmT = GetTotal_CD() * Factor;
                    
                    //          su2double ModDmT = 0.0;
                    //          if (nDim == 2) ModDmT = sqrt(GetTotal_CFx()*GetTotal_CFx() +
                    //                                       GetTotal_CFy()*GetTotal_CFy());
                    //
                    //          if (nDim == 3) ModDmT = sqrt(GetTotal_CFx()*GetTotal_CFx() +
                    //                                       GetTotal_CFy()*GetTotal_CFy() +
                    //                                       GetTotal_CFz()*GetTotal_CFz());
                    //
                    //          DmTVector[0] = GetTotal_CFx()/ModDmT;
                    //          DmTVector[1] = GetTotal_CFy()/ModDmT;
                    //          if (nDim == 3)  DmTVector[2] = GetTotal_CFz()/ModDmT;
                    
                    /*--- Set the aero drag ---*/
                    
                    su2double Aero_Drag = DmT - Force;
                    su2double Aero_CD = Aero_Drag / Factor;
                    
                    SetTotal_AeroCD(Aero_CD);
                    
                    /*--- Set the solid surface drag ---*/
                    
                    su2double Solid_Drag = DmT - Force;
                    su2double Solid_CD = Solid_Drag / Factor;
                    
                    SetTotal_SolidCD(Solid_CD);
                    
                    /*--- Set the net thrust value---*/
                    
                    su2double CT = NetThrust / Factor;
                    
                    SetTotal_NetThrust(CT);
                    
                    /*--- Set the total power ---*/
                    
                    su2double PowerHP = Power * Ref *  config->GetVelocity_Ref() / 550.0;
                    
                    SetTotal_Power(PowerHP);
                    
                    /*--- Set the total ReverseFlow ---*/
                    
                    su2double ReverseFlow;
                    if (Engine) ReverseFlow = fabs(config->GetInflow_ReverseMassFlow(iMarker_Inlet)  / config->GetInflow_MassFlow(Inlet_TagBound));
                    else ReverseFlow = fabs(config->GetActDisk_ReverseMassFlow(iMarker_Inlet)  / config->GetActDiskInlet_MassFlow(Inlet_TagBound));
                    
                    SetTotal_ReverseFlow(ReverseFlow);
                    
                    /*--- Set the total mass flow ratio ---*/
                    
                    InfVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) InfVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
                    if (Engine) MFR =fabs(config->GetInflow_MassFlow(Inlet_TagBound)) / (Density_Inf * sqrt(InfVel2) * config->GetHighlite_Area());
                    else MFR = fabs(config->GetActDiskInlet_MassFlow(Inlet_TagBound)) / (Density_Inf * sqrt(InfVel2) * config->GetHighlite_Area());
                    SetTotal_MFR(MFR);
                    
                    /*--- Evaluate shaft power and adiabatic efficiency (average) ---*/
                    
                    su2double Pstatic1, P1, P2, T1, T2;
                    if (Engine) {
                        Pstatic1 = config->GetInflow_Pressure(Inlet_TagBound);
                        P1 = config->GetInflow_TotalPressure(Inlet_TagBound);
                        P2 = config->GetExhaust_TotalPressure(Outlet_TagBound);
                        T1 = config->GetInflow_TotalTemperature(Inlet_TagBound);
                        T2 = config->GetExhaust_TotalTemperature(Outlet_TagBound);
                    }
                    else {
                        Pstatic1 = config->GetActDiskInlet_Pressure(Inlet_TagBound);
                        P1 = config->GetActDiskInlet_TotalPressure(Inlet_TagBound);
                        P2 = config->GetActDiskOutlet_TotalPressure(Outlet_TagBound);
                        T1 = config->GetActDiskInlet_TotalTemperature(Inlet_TagBound);
                        T2 = config->GetActDiskOutlet_TotalTemperature(Outlet_TagBound);
                    }
                    
                    /*-- Set the propulsive efficiency ---*/
                    
                    su2double mu_prop = fabs(DmT)*sqrt(RefVel2)/Power;
                    SetTotal_Prop_Eff(mu_prop);
                    
                    /*-- Set the bypass propulsive efficiency ---*/
                    
                    su2double mu_bypass_prop = NetThrust*sqrt(RefVel2)/Power;
                    SetTotal_ByPassProp_Eff(mu_bypass_prop);
                    
                    /*-- Set the fan adiabatic efficiency ---*/
                    
                    su2double mu_isentropic = 0.0;
                    if ((P2/P1) > 0.0) mu_isentropic =    (T1/(T2-T1))*(pow((P2/P1),(Gamma-1.0)/Gamma)-1.0);
                    SetTotal_Adiab_Eff(mu_isentropic);
                    
                    /*-- Set the polytropic efficiency ---*/
                    
                    su2double poly_coeff = 1.0/(1.0-log(T2/T1)/log(P2/P1));
                    su2double mu_polytropic = ((Gamma-1.0)/Gamma)/((poly_coeff-1.0)/poly_coeff);
                    SetTotal_Poly_Eff(mu_polytropic);
                    
                    if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
                        
                        if (iMarker_Inlet > 0) cout << endl;
                        
                        /*--- Geometry defintion ---*/
                        
                        if (Engine) cout <<"Engine surfaces: " << Inlet_TagBound << ", " << Outlet_TagBound << "." << endl;
                        else cout <<"Actuator disk surfaces: " << Inlet_TagBound << ", " << Outlet_TagBound << "." << endl;
                        
                        if (nDim == 2) {
                            if (config->GetSystemMeasurements() == SI)
                                cout <<"CG (m): (" << Inlet_XCG_Total[iMarker_Inlet] <<", " << Inlet_YCG_Total[iMarker_Inlet] << "). Length (m): " << Inlet_Area_Total[iMarker_Inlet] << "." << endl;
                            else if (config->GetSystemMeasurements() == US)
                                cout <<"CG (in): (" << Inlet_XCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_YCG_Total[iMarker_Inlet]*12.0 << "). Length (in): " << Inlet_Area_Total[iMarker_Inlet]*12.0 << "." << endl;
                            cout << endl;
                        }
                        
                        if (nDim ==3) {
                            if (config->GetSystemMeasurements() == SI)
                                cout <<"CG (m): (" << Inlet_XCG_Total[iMarker_Inlet] <<", " << Inlet_YCG_Total[iMarker_Inlet] <<", " << Inlet_ZCG_Total[iMarker_Inlet] << "). Area (m^2): " << Inlet_Area_Total[iMarker_Inlet] << ". Radius (m): " << sqrt(Inlet_Area_Total[iMarker_Inlet]/PI_NUMBER) << "." << endl;
                            else if (config->GetSystemMeasurements() == US)
                                cout <<"CG (in): (" << Inlet_XCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_YCG_Total[iMarker_Inlet]*12.0 <<", " << Inlet_ZCG_Total[iMarker_Inlet]*12.0 << "). Area (in^2): " << Inlet_Area_Total[iMarker_Inlet]*12.0*12.0 << "." << endl;
                            cout << endl;
                        }
                        
                        
                        /*--- Flow field descritption  ---*/
                        
                        if (config->GetSystemMeasurements() == SI) {
                            cout << setprecision(2) << "Inlet Ave. P (Pa): " << Pstatic1*config->GetPressure_Ref() << setprecision(3) <<  ". Inlet Ave. Mach: " << Mach << "." << endl;
                            cout << setprecision(2) << "Outlet Ave. PT (Pa): " << P2*config->GetPressure_Ref() << ". Outlet Ave. TT (K): " << T2*config->GetTemperature_Ref() << "." << endl;
                        }
                        else if (config->GetSystemMeasurements() == US) {
                            cout << setprecision(2) << "Inlet Ave. P (psf): " << Pstatic1*config->GetPressure_Ref() << setprecision(3) <<  ". Inlet Ave. Mach: " << Mach << "." << endl;
                            cout << setprecision(2) << "Outlet Ave. PT (psf): " << P2*config->GetPressure_Ref() << ". Outlet Ave. TT (R): " << T2*config->GetTemperature_Ref() << "." << endl;
                        }
                        
                        cout << "Inlet min. PT/PTinf: " << MinPressure << ". Inlet max. PT/PTinf: " << MaxPressure << ". Inlet Ave. PT/PTinf: " << AvePressure << endl;
                        
                        su2double InfVel2, Inlet_MassFlow, Outlet_MassFlow;
                        
                        if (Engine) Inlet_MassFlow = fabs(config->GetInflow_MassFlow(Inlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
                        else Inlet_MassFlow = fabs(config->GetActDiskInlet_MassFlow(Inlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
                        
                        if (config->GetSystemMeasurements() == SI) { cout << "Inlet mass flow (kg/s): "; cout << setprecision(2) << Inlet_MassFlow; }
                        else if (config->GetSystemMeasurements() == US) { cout << "Inlet mass flow (lbs/s): "; cout << setprecision(2) << Inlet_MassFlow * 32.174; }
                        
                        if (Engine) Outlet_MassFlow = fabs(config->GetExhaust_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
                        else Outlet_MassFlow = fabs(config->GetActDiskOutlet_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
                        
                        //          if (config->GetSystemMeasurements() == SI) { cout << ". Outlet mass flow (kg/s): "; cout << setprecision(2) << Outlet_MassFlow; }
                        //          else if (config->GetSystemMeasurements() == US) { cout << ". Outlet mass flow (lbs/s): "; cout << setprecision(2) << Outlet_MassFlow * 32.174; }
                        
                        if (Inlet_MassFlow > Outlet_MassFlow) cout << ". I/O diff.: " << setprecision(2) << 100.0*fabs(1.0-(Outlet_MassFlow/Inlet_MassFlow)) << "%";
                        else cout << ". I/O diff.: " << setprecision(2) << -100.0*fabs(1.0-(Inlet_MassFlow/Outlet_MassFlow)) << "%";
                        
                        InfVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) InfVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
                        cout << setprecision(2) << ". MFR: " << MFR << "." << endl;
                        
                        if (!Engine) {
                            
                            cout << setprecision(3) << "PT in/out ratio: " << TotalPressRatio << ". TT in/out ratio: " << TotalTempRatio  <<  "." << endl;
                            
                            if (config->GetActDisk_Jump() == VARIABLES_JUMP) {
                                if (config->GetSystemMeasurements() == SI) cout << setprecision(3) << "P in/out jump (Pa): ";
                                else if (config->GetSystemMeasurements() == US) cout << setprecision(3) << "P in/out jump (psf): ";
                                cout << setprecision(3) << DeltaPress * config->GetPressure_Ref();
                                if (config->GetSystemMeasurements() == SI) cout << setprecision(3) << ". T in/out jump (K): ";
                                else if (config->GetSystemMeasurements() == US) cout << setprecision(3) << ". T in/out jump (R): ";
                                cout << setprecision(3) << DeltaTemp * config->GetTemperature_Ref() <<"."<< endl;
                            }
                            else  if (config->GetActDisk_Jump() == RATIO) {
                                cout << setprecision(3) << "P in/out ratio: ";
                                cout << setprecision(3) << StaticPressRatio;
                                cout  << setprecision(3) <<". T in/out ratio: ";
                                cout << setprecision(3) << StaticTempRatio <<"."<< endl;
                            }
                        }
                        
                        cout << setprecision(1) << "\nProp. eff. (D-T.V/Shaft P): " << 100*mu_prop << "%. By-pass prop. eff. (NetT.V/Shaft P): " << 100*mu_bypass_prop <<   "%." << endl;
                        cout << setprecision(1) << "Fan adiabatic eff.: " << 100*mu_isentropic  << "%. Fan poly. eff.: " << 100*mu_polytropic << "%. Poly coeff. (n): " << setprecision(4) << poly_coeff << "." << endl;
                        
                        cout << endl;
                        
                        
                        /*--- Forces descritption  ---*/
                        
                        if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << "Ram Drag (N): ";
                        else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << "Ram Drag (lbf): ";
                        cout << (GrossThrust-NetThrust) * Ref;
                        
                        if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << ". Gross Thrust (N): ";
                        else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << ". Gross Thrust (lbf): ";
                        cout << -GrossThrust * Ref  << "." << endl;
                        
                        if (config->GetSystemMeasurements() == SI) cout << setprecision(1) << "Open surfaces Thurst (N): ";
                        else if (config->GetSystemMeasurements() == US) cout  << setprecision(1) << "Open surfaces Thrust (lbf): ";
                        cout<< setprecision(1) << Force * Ref << ". Open surfaces CT: " << setprecision(5) << -Force / Factor << "." << endl;
                        
                        if (config->GetSystemMeasurements() == SI) cout << "Solid surfaces Drag (N): ";
                        else if (config->GetSystemMeasurements() == US) cout << "Solid surfaces Drag (lbf): ";
                        cout << setprecision(1) << Solid_Drag * Ref << ". Solid surfaces CD: " << setprecision(5) << Solid_CD << "." << endl;
                        
                        if (config->GetSystemMeasurements() == SI) cout << setprecision(1) <<"Net Thrust (N): ";
                        else if (config->GetSystemMeasurements() == US) cout << setprecision(1) << "Net Thrust (lbf): ";
                        cout << setprecision(5) << -NetThrust * Ref  << ". Net CT: " << CT;
                        
                        if (config->GetSystemMeasurements() == SI) {
                            cout << ". Power (W): ";
                            cout << setprecision(1) << Power * Ref *  config->GetVelocity_Ref()  << "." << endl;
                        }
                        else if (config->GetSystemMeasurements() == US) {
                            cout << ". Power (HP): ";
                            cout << setprecision(1) << Power * Ref *  config->GetVelocity_Ref() / 550.0 << "." << endl;
                        }
                        
                    }
                    
                }
                
                if (write_heads && Output && !config->GetDiscrete_Adjoint()) cout << "-------------------------------------------------------------------------" << endl << endl;
                
            }
            
        }
        
        delete [] Outlet_MassFlow_Local;
        delete [] Outlet_Temperature_Local;
        delete [] Outlet_TotalTemperature_Local;
        delete [] Outlet_Pressure_Local;
        delete [] Outlet_TotalPressure_Local;
        delete [] Outlet_Area_Local;
        delete [] Outlet_GrossThrust_Local;
        delete [] Outlet_Force_Local;
        delete [] Outlet_Power_Local;
        
        delete [] Outlet_MassFlow_Total;
        delete [] Outlet_Temperature_Total;
        delete [] Outlet_TotalTemperature_Total;
        delete [] Outlet_Pressure_Total;
        delete [] Outlet_TotalPressure_Total;
        delete [] Outlet_Area_Total;
        delete [] Outlet_GrossThrust_Total;
        delete [] Outlet_Force_Total;
        delete [] Outlet_Power_Total;
        
        delete [] Inlet_MassFlow_Local;
        delete [] Inlet_ReverseMassFlow_Local;
        delete [] Inlet_Temperature_Local;
        delete [] Inlet_TotalTemperature_Local;
        delete [] Inlet_Pressure_Local;
        delete [] Inlet_Mach_Local;
        delete [] Inlet_MinPressure_Local;
        delete [] Inlet_MaxPressure_Local;
        delete [] Inlet_TotalPressure_Local;
        delete [] Inlet_Area_Local;
        delete [] Inlet_RamDrag_Local;
        delete [] Inlet_Force_Local;
        delete [] Inlet_Power_Local;
        delete [] Inlet_XCG_Local;
        delete [] Inlet_YCG_Local;
        delete [] Inlet_ZCG_Local;
        
        delete [] Inlet_MassFlow_Total;
        delete [] Inlet_ReverseMassFlow_Total;
        delete [] Inlet_Temperature_Total;
        delete [] Inlet_TotalTemperature_Total;
        delete [] Inlet_Pressure_Total;
        delete [] Inlet_Mach_Total;
        delete [] Inlet_MinPressure_Total;
        delete [] Inlet_MaxPressure_Total;
        delete [] Inlet_TotalPressure_Total;
        delete [] Inlet_Area_Total;
        delete [] Inlet_RamDrag_Total;
        delete [] Inlet_Force_Total;
        delete [] Inlet_Power_Total;
        delete [] Inlet_XCG_Total;
        delete [] Inlet_YCG_Total;
        delete [] Inlet_ZCG_Total;
        
        delete [] Inlet_MassFlow;
        delete [] Inlet_Mach;
        delete [] Inlet_MinPressure;
        delete [] Inlet_MaxPressure;
        delete [] Inlet_ReverseMassFlow;
        delete [] Inlet_Pressure;
        delete [] Inlet_TotalPressure;
        delete [] Inlet_Temperature;
        delete [] Inlet_TotalTemperature;
        delete [] Inlet_Area;
        delete [] Inlet_RamDrag;
        delete [] Inlet_Force;
        delete [] Inlet_Power;
        delete [] Inlet_XCG;
        delete [] Inlet_YCG;
        delete [] Inlet_ZCG;
        
        delete [] Outlet_MassFlow;
        delete [] Outlet_Pressure;
        delete [] Outlet_TotalPressure;
        delete [] Outlet_Temperature;
        delete [] Outlet_TotalTemperature;
        delete [] Outlet_Area;
        delete [] Outlet_GrossThrust;
        delete [] Outlet_Force;
        delete [] Outlet_Power;
        
    }
    
}

void CEulerSolver::SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                                       CConfig *config, unsigned short iMesh, bool Output) {
    
    su2double Massflow = 0.0 , Target_Massflow = 0.0, DragMinusThrust = 0.0 ,
    Target_DragMinusThrust = 0.0, Target_NetThrust = 0.0, BCThrust = 0.0, BCThrust_inc = 0.0;
    unsigned short iDim, iMarker;
    unsigned long iVertex, iPoint;
    su2double  *V_inlet = NULL, Pressure,
    Density, T0_Ti, ATerm, BTerm, LHS, RHS, RHS_PDelta, RHS_MDelta, F, DF_DLa, CTerm_, DTerm_,
    ETerm, La, La_old, TotalArea, To_Ti, DeltaT, Po_Pi, DeltaP, Area, Velocity_Normal,
    SoundSpeed2, Force_Normal,
    RefDensity, RefArea, RefVel2, Factor, Ref;
    unsigned short iter;
    string Marker_Tag;
    su2double Target_Force, Force, Target_Power, Power, NetThrust, BCThrust_old, Initial_BCThrust;
    bool ActDisk_Info;
    su2double MyBCThrust, BCThrust_Init;
    
    su2double dNetThrust_dBCThrust        = config->GetdNetThrust_dBCThrust();
    unsigned short Kind_ActDisk           = config->GetKind_ActDisk();
    bool ratio                            = (config->GetActDisk_Jump() == RATIO);
    unsigned long Update_BCThrust         = config->GetUpdate_BCThrust();
    unsigned long Iter_Fixed_NetThrust    = config->GetIter_Fixed_NetThrust();
    unsigned long InnerIter                 = config->GetInnerIter();
    bool Update_BCThrust_Bool             = false;
    bool restart                          = (config->GetRestart() || config->GetRestart_Flow());
    su2double Fan_Poly_Eff                = config->GetFan_Poly_Eff();
    su2double PolyCoeff                   = 1.0/(1.0-((Gamma-1.0)/Gamma)/Fan_Poly_Eff);
    
    RefDensity   = Density_Inf;
    RefArea = config->GetRefArea();
    RefVel2 = 0.0;  for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
    
    Factor = (0.5*RefDensity*RefArea*RefVel2);
    Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
    
    /*--- Delta P and delta T are inputs ---*/
    
    if (Kind_ActDisk == VARIABLES_JUMP) {
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
            
            if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
                
                Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                
                if (ratio) {
                    if (config->GetMach()  < 0.5) {
                        DeltaP       = config->GetActDisk_PressJump(Marker_Tag, 0);
                        DeltaT       = config->GetActDisk_TempJump(Marker_Tag, 0);
                    }
                    else {
                        DeltaP       = config->GetActDisk_PressJump(Marker_Tag, 1);
                        DeltaT       = config->GetActDisk_TempJump(Marker_Tag, 1);
                    }
                }
                else {
                    if (config->GetMach()  < 0.5) {
                        DeltaP       = max(0.0, config->GetActDisk_PressJump(Marker_Tag, 0) / config->GetPressure_Ref());
                        DeltaT       = max(0.0, config->GetActDisk_TempJump(Marker_Tag, 0) / config->GetTemperature_Ref());
                    }
                    else {
                        DeltaP       = max(0.0, config->GetActDisk_PressJump(Marker_Tag, 1) / config->GetPressure_Ref());
                        DeltaT       = max(0.0, config->GetActDisk_TempJump(Marker_Tag, 1) / config->GetTemperature_Ref());
                    }
                }
                
                /*--- Set the Delta P, Delta T values at each discrete point (uniform distribution)  ---*/
                
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    SetActDisk_DeltaP(iMarker, iVertex, DeltaP);
                    SetActDisk_DeltaT(iMarker, iVertex, DeltaT);
                }
                
            }
        }
    }
    
    /*--- Iteration using BCThrust ---*/
    
    else {
        
        if (InnerIter == 0) BCThrust_Counter = 0;
        
        /*--- Only the fine mesh level should check the convergence criteria ---*/
        
        if ((iMesh == MESH_0) && Output) {
            
            /*--- Initialize the update flag to false ---*/
            
            Update_BCThrust_Bool = false;
            
            /*--- Reevaluate BCThrust at a fix number of iterations ---*/
            
            if ((InnerIter % Iter_Fixed_NetThrust == 0) && (InnerIter != 0)) {
                BCThrust_Counter++;
                if ((BCThrust_Counter != 0) &&
                    (BCThrust_Counter != 1) &&
                    (BCThrust_Counter != Update_BCThrust) &&
                    (BCThrust_Counter != Update_BCThrust + 2) &&
                    (BCThrust_Counter != Update_BCThrust + 4) ) Update_BCThrust_Bool = true;
                else Update_BCThrust_Bool = false;
            }
            
            /*--- Store the update boolean for use on other mesh levels in the MG ---*/
            
            config->SetUpdate_BCThrust_Bool(Update_BCThrust_Bool);
            
        }
        
        else {
            Update_BCThrust_Bool = config->GetUpdate_BCThrust_Bool();
        }
        
        
        /*--- If it is the first iteration, set the BCThrust to a meaning full target value,
         * this can be done at an initialization level, for the time being it is OK here ---*/
        
        if (InnerIter == 0) {
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                    (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
                    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                    
                    if (Kind_ActDisk == NET_THRUST) {
                        if (restart)
                            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
                        else {
                            if (config->GetMach() < 0.5) Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
                            else Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
                        }
                        config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
                        config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
                    }
                    
                    if (Kind_ActDisk == BC_THRUST) {
                        if (restart)
                            Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
                        else {
                            if (config->GetMach() < 0.5) Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
                            else Initial_BCThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
                        }
                        config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
                        config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
                    }
                    
                    if (Kind_ActDisk == POWER) {
                        Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
                        config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
                        config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
                    }
                    
                    if (Kind_ActDisk == DRAG_MINUS_THRUST) {
                        Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
                        config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
                        config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
                    }
                    
                    if (Kind_ActDisk == MASSFLOW) {
                        Initial_BCThrust = config->GetInitial_BCThrust() / Ref;
                        config->SetActDisk_BCThrust(Marker_Tag, Initial_BCThrust);
                        config->SetActDisk_BCThrust_Old(Marker_Tag, Initial_BCThrust);
                    }
                    
                }
            }
        }
        
        /*--- Typical iteration to set the value of BC Thrust at each actuator disk ---*/
        
        if (Update_BCThrust_Bool && Output) {
            
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                
                if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                    (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
                    
                    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                    
                    if (Kind_ActDisk == NET_THRUST) {
                        
                        if (config->GetMach() < 0.5) Target_NetThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
                        else Target_NetThrust = fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
                        NetThrust    = config->GetActDisk_NetThrust(Marker_Tag);
                        BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
                        BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_NetThrust - NetThrust);
                        
                        if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
                        else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
                        
                        if (iMesh == MESH_0) {
                            config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
                            BCThrust_Init = BCThrust*Ref;
                            config->SetInitial_BCThrust(BCThrust_Init);
                        }
                        
                    }
                    
                    if (Kind_ActDisk == BC_THRUST) {
                        
                        if (config->GetMach() < 0.5) Target_Force =  fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / Ref);
                        else Target_Force =  fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / Ref);
                        Force        = -config->GetActDisk_Force(Marker_Tag);
                        BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
                        BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_Force - Force);
                        
                        if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
                        else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
                        
                        if (iMesh == MESH_0) {
                            config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
                            BCThrust_Init = BCThrust*Ref;
                            config->SetInitial_BCThrust(BCThrust_Init);
                        }
                        
                    }
                    
                    if (Kind_ActDisk == POWER) {
                        
                        if (config->GetMach() < 0.5) Target_Power =  fabs( config->GetActDisk_PressJump(Marker_Tag, 0) / (Ref * config->GetVelocity_Ref() /  550.0));
                        else Target_Power =  fabs( config->GetActDisk_PressJump(Marker_Tag, 1) / (Ref * config->GetVelocity_Ref() /  550.0));
                        Power        = config->GetActDisk_Power(Marker_Tag);
                        BCThrust_old = config->GetActDisk_BCThrust_Old(Marker_Tag);
                        BCThrust_inc = (1.0/dNetThrust_dBCThrust)*(Target_Power - Power);
                        
                        if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
                        else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
                        
                        if (iMesh == MESH_0) {
                            config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
                            BCThrust_Init = BCThrust*Ref;
                            config->SetInitial_BCThrust(BCThrust_Init);
                        }
                        
                    }
                    
                    if (Kind_ActDisk == DRAG_MINUS_THRUST) {
                        
                        if (config->GetMach() < 0.5) Target_DragMinusThrust  =  -fabs(config->GetActDisk_PressJump(Marker_Tag, 0)) * Factor;
                        else Target_DragMinusThrust  =  -fabs(config->GetActDisk_PressJump(Marker_Tag, 1)) * Factor;
                        DragMinusThrust = GetTotal_CD() * Factor;
                        BCThrust_old    = config->GetActDisk_BCThrust_Old(Marker_Tag);
                        BCThrust_inc    = -(1.0/dNetThrust_dBCThrust)*(Target_DragMinusThrust - DragMinusThrust);
                        
                        if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
                        else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
                        
                        if (iMesh == MESH_0) {
                            config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
                            BCThrust_Init = BCThrust*Ref;
                            config->SetInitial_BCThrust(BCThrust_Init);
                        }
                        
                    }
                    
                    if (Kind_ActDisk == MASSFLOW) {
                        
                        if (config->GetMach() < 0.5) {
                            Target_Massflow  =  fabs(config->GetActDisk_PressJump(Marker_Tag, 0) / (config->GetDensity_Ref() * config->GetVelocity_Ref()));
                            if (config->GetSystemMeasurements() == US) Target_Massflow /= 32.174;
                        }
                        else {
                            Target_Massflow  =  fabs(config->GetActDisk_PressJump(Marker_Tag, 1) / (config->GetDensity_Ref() * config->GetVelocity_Ref()));
                            if (config->GetSystemMeasurements() == US) Target_Massflow /= 32.174;
                        }
                        
                        Massflow = config->GetActDisk_MassFlow(Marker_Tag);
                        BCThrust_old    = config->GetActDisk_BCThrust_Old(Marker_Tag);
                        BCThrust_inc    = (1.0/dNetThrust_dBCThrust)*(Target_Massflow - Massflow);
                        if (iMesh == MESH_0) BCThrust = max(0.0,(BCThrust_old + BCThrust_inc));
                        else BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
                        if (iMesh == MESH_0) {
                            config->SetActDisk_BCThrust(Marker_Tag, BCThrust);
                            BCThrust_Init = BCThrust*Ref;
                            config->SetInitial_BCThrust(BCThrust_Init);
                        }
                        
                    }
                    
                }
                
            }
            
            /*--- After a complete update of BC_Thrust
             update the value of BC Thrust (old) for future iterations ---*/
            
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                    (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
                    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                    if ((Kind_ActDisk == NET_THRUST) || (Kind_ActDisk == BC_THRUST) ||
                        (Kind_ActDisk == POWER) || (Kind_ActDisk == DRAG_MINUS_THRUST) ||
                        (Kind_ActDisk == MASSFLOW)) {
                        BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
                        config->SetActDisk_BCThrust_Old(Marker_Tag, BCThrust);
                    }
                }
            }
            
        }
        
        /*--- Evaluate the pressure jump at each node using the total thrust ---*/
        
        if ((Update_BCThrust_Bool && Output) || (InnerIter == 0)) {
            
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                
                if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
                    (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
                    
                    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                    RefDensity  = Density_Inf;
                    RefArea = config->GetRefArea();
                    RefVel2 = 0.0; for (iDim = 0; iDim < nDim; iDim++) RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
                    
                    Factor = (0.5*RefDensity*RefArea*RefVel2);
                    Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;
                    BCThrust = config->GetActDisk_BCThrust(Marker_Tag);
                    
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        
                        if (geometry->node[iPoint]->GetDomain()) {
                            
                            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
                            
                            if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) {
                                for (iDim = 0; iDim < nDim; iDim++) { Vector[iDim] = -Vector[iDim]; }
                            }
                            
                            Area = 0.0;
                            for (iDim = 0; iDim < nDim; iDim++) { Area += Vector[iDim]*Vector[iDim]; }
                            Area = sqrt (Area);
                            
                            /*--- Use the inlet state to compute the Pressure and Temperature jumps ---*/
                            
                            if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET)
                                V_inlet = nodes->GetPrimitive(iPoint);
                            if (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)
                                V_inlet = GetDonorPrimVar(iMarker, iVertex);
                            
                            Density       = V_inlet[nDim+2];
                            Pressure      = V_inlet[nDim+1];
                            SoundSpeed2   = Pressure*Gamma/Density;
                            TotalArea     = config->GetActDisk_Area(Marker_Tag);
                            Force_Normal  = Area*(BCThrust/TotalArea);
                            
                            
                            Velocity_Normal = 0.0;
                            for (iDim = 0; iDim < nDim; iDim++) {
                                Velocity_Normal += V_inlet[iDim+1]*Vector[iDim]/Area;
                            }
                            
                            
                            if (Velocity_Normal > EPS) {
                                
                                /*--- Ratio of the total temperature to the temperature at the inflow ---*/
                                
                                T0_Ti = 1.0 + ((Gamma-1.0)/SoundSpeed2)*(0.5*Velocity_Normal*Velocity_Normal + Force_Normal/(Density*Area));
                                
                                
                                ATerm = 2.0*T0_Ti/(Gamma+1.0);
                                BTerm = 0.5*(Gamma+1.0)/(Gamma-1.0);
                                LHS = fabs(Velocity_Normal)/(sqrt(SoundSpeed2)*pow(ATerm,BTerm));
                                
                                CTerm_ = (PolyCoeff-1.0)/(PolyCoeff+1.0);
                                DTerm_ = 1.0/(PolyCoeff-1.0);
                                
                                La = EPS; La_old = EPS;
                                
                                for (iter = 0; iter < 100; iter++) {
                                    
                                    ETerm = ((1.0-CTerm_*La*La)/(1.0-CTerm_+EPS));
                                    
                                    RHS = La*pow(ETerm, DTerm_);
                                    
                                    ETerm = ((1.0-CTerm_*(La+1E-6)*(La+1E-6))/(1.0-CTerm_+EPS));
                                    RHS_PDelta = (La+1E-6)*pow(ETerm, DTerm_);
                                    
                                    ETerm = ((1.0-CTerm_*(La-1E-6)*(La-1E-6))/(1.0-CTerm_+EPS));
                                    RHS_MDelta = (La-1E-6)*pow(ETerm, DTerm_);
                                    
                                    /*--- Objective function and finitte differences derivative ---*/
                                    
                                    F = RHS - LHS;
                                    DF_DLa = (RHS_PDelta - RHS_MDelta)/2E-6;
                                    
                                    /*--- Newton's step ---*/
                                    
                                    La_old = La;
                                    La = La_old - 0.75*(F/DF_DLa);
                                    
                                    if (fabs(F) < 1E-10) break;
                                    
                                }
                                
                                if (iter == 99) cout << "The laval number evaluation is not converging." << endl;
                                
                                /*--- Laval is bounded ---*/
                                
                                La = min(La, sqrt(6.0));  La = max(La, 0.0);
                                
                                To_Ti = max(1.0, T0_Ti*(1.0-CTerm_*La*La));
                                SetActDisk_DeltaT(iMarker, iVertex, To_Ti);
                                
                                Po_Pi = max(1.0, pow(To_Ti, PolyCoeff*DTerm_));
                                SetActDisk_DeltaP(iMarker, iVertex, Po_Pi);
                                
                            }
                            else {
                                SetActDisk_DeltaT(iMarker, iVertex, 1.0);
                                SetActDisk_DeltaP(iMarker, iVertex, 1.0);
                            }
                            
                        }
                        
                    }
                }
                
            }
        }
    }
    
    /*--- Broadcast some information to the master node ---*/
    
    ActDisk_Info = false;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
            (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
            ActDisk_Info = true;
        }
    }
    if (!ActDisk_Info) config->SetInitial_BCThrust(0.0);
    
    MyBCThrust = config->GetInitial_BCThrust();
    SU2_MPI::Allreduce(&MyBCThrust, &BCThrust, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    config->SetInitial_BCThrust(BCThrust);
    
}

void CEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output) {
    
    su2double AoA = 0.0, Vel_Infty[3], Vel_Infty_Mag;
    unsigned short iDim;
    unsigned long InnerIter = config->GetInnerIter();
    su2double Beta = config->GetAoS();
    
    /* --- Initialize values at first iteration --- */
    
    if (InnerIter == 0) {
        Total_CD_Prev = 0.0;
        Total_CL_Prev = 0.0;
        Total_CMx_Prev = 0.0;
        Total_CMy_Prev = 0.0;
        Total_CMz_Prev = 0.0;
        AoA_Prev = config->GetAoA();
        dCL_dAlpha = config->GetdCL_dAlpha();
        AoA_inc = 0.0;
    }
    
    /*--- Retrieve the AoA (degrees) ---*/
    
    AoA = config->GetAoA();
    
    /* --- Set new AoA if needed --- */
    
    if (fabs(AoA_inc) > 0.0 && Output) {
        
        /* --- Update *_Prev values with current coefficients --- */
        
        SetCoefficient_Gradients(config);
        
        Total_CD_Prev = TotalCoeff.CD;
        Total_CL_Prev = TotalCoeff.CL;
        Total_CMx_Prev = TotalCoeff.CMx;
        Total_CMy_Prev = TotalCoeff.CMy;
        Total_CMz_Prev = TotalCoeff.CMz;
        AoA_Prev = AoA;
        
        /*--- Compute a new value for AoA on the fine mesh only (degrees)---*/
        
        if (iMesh == MESH_0) AoA = AoA + AoA_inc;
        else { AoA = config->GetAoA(); }
        
        /*--- Only the fine mesh stores the updated values for AoA in config ---*/
        
        if (iMesh == MESH_0) {
            config->SetAoA(AoA);
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
            Vel_Infty[0] = cos(AoA*PI_NUMBER/180.0)*Vel_Infty_Mag;
            Vel_Infty[1] = sin(AoA*PI_NUMBER/180.0)*Vel_Infty_Mag;
        }
        if (nDim == 3) {
            Vel_Infty[0] = cos(AoA*PI_NUMBER/180.0)*cos(Beta*PI_NUMBER/180.0)*Vel_Infty_Mag;
            Vel_Infty[1] = sin(Beta)*Vel_Infty_Mag;
            Vel_Infty[2] = sin(AoA*PI_NUMBER/180.0)*cos(Beta*PI_NUMBER/180.0)*Vel_Infty_Mag;
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
        
    }
}

bool CEulerSolver::FixedCL_Convergence(CConfig* config, bool convergence) {
    su2double Target_CL = config->GetTarget_CL();
    unsigned long curr_iter = config->GetInnerIter();
    unsigned long Iter_dCL_dAlpha = config->GetIter_dCL_dAlpha();
    bool fixed_cl_conv = false;
    AoA_inc = 0.0;
    
    
    /*--- if in Fixed CL mode, before finite differencing --- */
    
    if (!Start_AoA_FD){
        if (convergence){
            
            /* --- C_L and solution are converged, start finite differencing --- */
            
            if (fabs(TotalCoeff.CL-Target_CL) < (config->GetCauchy_Eps()/2)) {
                
                /* --- If no finite differencing required --- */
                
                if (Iter_dCL_dAlpha == 0){
                    fixed_cl_conv = true;
                    return fixed_cl_conv;
                }
                
                /* --- Else, set up finite differencing routine ---*/
                
                Iter_Update_AoA = curr_iter;
                Start_AoA_FD = true;
                fixed_cl_conv = false;
                AoA_inc = 0.001;
            }
            
            /* --- C_L is not converged to target value and some iterations
             have passed since last update, so update AoA --- */
            
            else if ((curr_iter - Iter_Update_AoA) > config->GetStartConv_Iter()){
                Iter_Update_AoA = curr_iter;
                fixed_cl_conv = false;
                if (fabs(TotalCoeff.CL-Target_CL) > (config->GetCauchy_Eps()/2)) {
                    AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - TotalCoeff.CL);
                }
            }
        }
        
        /* --- If the iteration limit between AoA updates is met, so update AoA --- */
        
        else if ((curr_iter - Iter_Update_AoA) == config->GetUpdate_AoA_Iter_Limit()) {
            Iter_Update_AoA = curr_iter;
            fixed_cl_conv = false;
            if (fabs(TotalCoeff.CL-Target_CL) > (config->GetCauchy_Eps()/2)) {
                AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - TotalCoeff.CL);
            }
        }
        
        /* --- If the total iteration limit is reached, start finite differencing --- */
        
        if (curr_iter == config->GetnInner_Iter() - Iter_dCL_dAlpha){
            if (Iter_dCL_dAlpha == 0){
                End_AoA_FD = true;
            }
            Iter_Update_AoA = curr_iter;
            Start_AoA_FD = true;
            fixed_cl_conv = false;
            AoA_inc = 0.001;
        }
    }
    
    /* --- If Finite Difference Mode has ended, end simulation --- */
    
    if (End_AoA_FD){
        //fixed_cl_conv = true;
        return true;
    }
    
    /* --- If starting Finite Difference Mode --- */
    
    if (Start_AoA_FD){
        
        /* --- Disable history writing --- */
        
        config->SetHistory_Wrt_Freq(2, 0);
        
        /* --- End Finite Difference Mode if iteration limit is reached, so simualtion is converged --- */
        
        End_AoA_FD = ((curr_iter - Iter_Update_AoA - 2) == Iter_dCL_dAlpha ||
                      curr_iter == config->GetnInner_Iter()- 2 );
        
        if (convergence && (curr_iter - Iter_Update_AoA) > config->GetStartConv_Iter())
            End_AoA_FD = true;
        
        
        /* --- If Finite Difference mode is ending, reset AoA and calculate Coefficient Gradients --- */
        
        if (End_AoA_FD){
            SetCoefficient_Gradients(config);
            config->SetAoA(AoA_Prev);
        }
    }
    
    return fixed_cl_conv;
    
}

void CEulerSolver::SetCoefficient_Gradients(CConfig *config){
    su2double dCL_dAlpha_, dCD_dCL_, dCMx_dCL_, dCMy_dCL_, dCMz_dCL_;
    su2double AoA = config->GetAoA();
    
    if (AoA != AoA_Prev) {
        /* --- Calculate gradients of coefficients w.r.t. CL --- */
        
        dCL_dAlpha_ = (TotalCoeff.CL-Total_CL_Prev)/(AoA - AoA_Prev);
        dCD_dCL_    = (TotalCoeff.CD-Total_CD_Prev)/(TotalCoeff.CL-Total_CL_Prev);
        dCMx_dCL_   = (TotalCoeff.CMx-Total_CMx_Prev)/(TotalCoeff.CL-Total_CL_Prev);
        dCMy_dCL_   = (TotalCoeff.CMy-Total_CMy_Prev)/(TotalCoeff.CL-Total_CL_Prev);
        dCMz_dCL_   = (TotalCoeff.CMz-Total_CMz_Prev)/(TotalCoeff.CL-Total_CL_Prev);
        
        /*--- Set the value of the  dOF/dCL in the config file ---*/
        
        config->SetdCD_dCL(dCD_dCL_);
        config->SetdCMx_dCL(dCMx_dCL_);
        config->SetdCMy_dCL(dCMy_dCL_);
        config->SetdCMz_dCL(dCMz_dCL_);
        config->SetdCL_dAlpha(dCL_dAlpha_);
    }
}

void CEulerSolver::SetInletAtVertex(su2double *val_inlet,
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

su2double CEulerSolver::GetInletAtVertex(su2double *val_inlet,
                                         unsigned long val_inlet_point,
                                         unsigned short val_kind_marker,
                                         string val_marker,
                                         CGeometry *geometry,
                                         CConfig *config) const {
    
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

void CEulerSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {
    
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
        
        string Marker_Tag   = config->GetMarker_All_TagBound(iMarker);
        su2double p_total   = config->GetInlet_Ptotal(Marker_Tag);
        su2double t_total   = config->GetInlet_Ttotal(Marker_Tag);
        su2double* flow_dir = config->GetInlet_FlowDir(Marker_Tag);
        
        for (unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
            Inlet_Ttotal[iMarker][iVertex] = t_total;
            Inlet_Ptotal[iMarker][iVertex] = p_total;
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
                Inlet_FlowDir[iMarker][iVertex][iDim] = flow_dir[iDim];
        }
        
    } else {
        
        /*--- For now, non-inlets just get set to zero. In the future, we
         can do more customization for other boundary types here. ---*/
        
        for (unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
            Inlet_Ttotal[iMarker][iVertex] = 0.0;
            Inlet_Ptotal[iMarker][iVertex] = 0.0;
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
                Inlet_FlowDir[iMarker][iVertex][iDim] = 0.0;
        }
    }
    
}

void CEulerSolver::UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config){
    
    unsigned short nMGlevel, iMarker;
    
    // TODO: Update the fluid boundary conditions for MG
    nMGlevel = config->GetnMGLevels();
    if (nMGlevel > 1) {
        for (iMarker=0; iMarker < nMarker; iMarker++) {
            bool isCustomizable = config->GetMarker_All_PyCustom(iMarker);
            bool isInlet = (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
            if (isCustomizable && isInlet)
                SU2_MPI::Error("Custom inlet BCs are not currently compatible with multigrid.", CURRENT_FUNCTION);
        }
    }
}

void CEulerSolver::Evaluate_ObjFunc(CConfig *config) {
    
    unsigned short iMarker_Monitoring, Kind_ObjFunc;
    su2double Weight_ObjFunc;
    
    Total_ComboObj = 0.0;
    
    /*--- Loop over all monitored markers, add to the 'combo' objective ---*/
    
    for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        
        Weight_ObjFunc = config->GetWeight_ObjFunc(iMarker_Monitoring);
        Kind_ObjFunc = config->GetKind_ObjFunc(iMarker_Monitoring);
        
        switch(Kind_ObjFunc) {
            case DRAG_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*(SurfaceCoeff.CD[iMarker_Monitoring]);
                if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCD_dCL()*(SurfaceCoeff.CL[iMarker_Monitoring]);
                if (config->GetFixed_CM_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCD_dCMy()*(SurfaceCoeff.CMy[iMarker_Monitoring]);
                break;
            case LIFT_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*(SurfaceCoeff.CL[iMarker_Monitoring]);
                break;
            case SIDEFORCE_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*(SurfaceCoeff.CSF[iMarker_Monitoring]);
                break;
            case EFFICIENCY:
                Total_ComboObj+=Weight_ObjFunc*(SurfaceCoeff.CEff[iMarker_Monitoring]);
                break;
            case MOMENT_X_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*(SurfaceCoeff.CMx[iMarker_Monitoring]);
                if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMx_dCL()*(SurfaceCoeff.CL[iMarker_Monitoring]);
                break;
            case MOMENT_Y_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*(SurfaceCoeff.CMy[iMarker_Monitoring]);
                if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMy_dCL()*(SurfaceCoeff.CL[iMarker_Monitoring]);
                break;
            case MOMENT_Z_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*(SurfaceCoeff.CMz[iMarker_Monitoring]);
                if (config->GetFixed_CL_Mode()) Total_ComboObj -= Weight_ObjFunc*config->GetdCMz_dCL()*(SurfaceCoeff.CL[iMarker_Monitoring]);
                break;
            case FORCE_X_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*SurfaceCoeff.CFx[iMarker_Monitoring];
                break;
            case FORCE_Y_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*SurfaceCoeff.CFy[iMarker_Monitoring];
                break;
            case FORCE_Z_COEFFICIENT:
                Total_ComboObj+=Weight_ObjFunc*SurfaceCoeff.CFz[iMarker_Monitoring];
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
            Total_ComboObj+=Weight_ObjFunc*TotalCoeff.CT;
            break;
        case TORQUE_COEFFICIENT:
            Total_ComboObj+=Weight_ObjFunc*TotalCoeff.CQ;
            break;
        case FIGURE_OF_MERIT:
            Total_ComboObj+=Weight_ObjFunc*TotalCoeff.CMerit;
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
        case SURFACE_DC60_DISTORTION:
            Total_ComboObj+=Weight_ObjFunc*config->GetSurface_DC60Distortion(0);
            break;
        case AVG_NORMAL_VEL:
            Total_ComboObj+=Weight_ObjFunc*config->GetSurface_NormalVelocity(0);
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


void CEulerSolver::BC_Euler_Wall(CGeometry      *geometry,
                                 CSolver        **solver_container,
                                 CNumerics      *conv_numerics,
                                 CNumerics      *visc_numerics,
                                 CConfig        *config,
                                 unsigned short val_marker) {
    
    /*--- Call the equivalent symmetry plane boundary condition. ---*/
    BC_Sym_Plane(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    
}


void CEulerSolver::BC_Sym_Plane(CGeometry      *geometry,
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
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                          geometry->node[iPoint]->GetGridVel());
            
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
            conv_numerics->SetSecondary(nodes->GetSecondary(iPoint),
                                        nodes->GetSecondary(iPoint));
            
            /*--- Compute the residual using an upwind scheme. ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            /*--- Update residual value ---*/
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration. ---*/
            if (implicit) {
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
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
                visc_numerics->SetSecondary(nodes->GetSecondary(iPoint),
                                            nodes->GetSecondary(iPoint));
                
                /*--- For viscous Fluxes also the gradients of the primitives need to be determined.
                 1. The gradients of scalars are mirrored along the sym plane just as velocity for the primitives
                 2. The gradients of the velocity components need more attention, i.e. the gradient of the
                 normal velocity in tangential direction is mirrored and the gradient of the tangential velocity in
                 normal direction is mirrored. ---*/
                
                /*--- Get gradients of primitives of boundary cell ---*/
                for (iVar = 0; iVar < nPrimVarGrad; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        Grad_Reflected[iVar][iDim] = nodes->GetGradient_Primitive(iPoint, iVar, iDim);
                
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
                if (config->GetKind_Turb_Model() == SST)
                    visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                                        solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
                
                /*--- Compute and update residual. Note that the viscous shear stress tensor is computed in the
                 following routine based upon the velocity-component gradients. ---*/
                auto residual = visc_numerics->ComputeResidual(config);
                
                LinSysRes.SubtractBlock(iPoint, residual);
                
                /*--- Jacobian contribution for implicit integration. ---*/
                if (implicit)
                    Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
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
    bool viscous        = config->GetViscous();
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    
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
            Pressure_Bound   = nodes->GetPressure(iPoint);
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
            
            
            
            /*--- Set various quantities in the numerics class ---*/
            
            conv_numerics->SetPrimitive(V_domain, V_infty);
            
            if (dynamic_grid) {
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                          geometry->node[iPoint]->GetGridVel());
            }
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            /*--- Update residual value ---*/
            
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Convective Jacobian contribution for implicit integration ---*/
            
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            /*--- Viscous residual contribution ---*/
            
            if (viscous) {
                
                /*--- Set laminar and eddy viscosity at the infinity ---*/
                
                V_infty[nDim+5] = nodes->GetLaminarViscosity(iPoint);
                V_infty[nDim+6] = nodes->GetEddyViscosity(iPoint);
                
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
                
                /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
                
                visc_numerics->SetTauWall(-1.0, -1.0);
                
                /*--- Compute and update viscous residual ---*/
                
                auto residual = visc_numerics->ComputeResidual(config);
                LinSysRes.SubtractBlock(iPoint, residual);
                
                /*--- Viscous Jacobian contribution for implicit integration ---*/
                
                if (implicit)
                    Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
                
            }
            
        }
    }
    
    /*--- Free locally allocated memory ---*/
    delete [] Normal;
    
}

void CEulerSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container,
                              CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) {
    
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
    string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
    bool viscous              = config->GetViscous();
    bool gravity = (config->GetGravityForce());
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    
    su2double *Normal, *FlowDirMix, TangVelocity, NormalVelocity;
    Normal = new su2double[nDim];
    
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
            V_domain = nodes->GetPrimitive(iPoint);
            
            /*--- Compute the internal state u_i ---*/
            Velocity2_i = 0;
            for (iDim=0; iDim < nDim; iDim++)
            {
                Velocity_i[iDim] = nodes->GetVelocity(iPoint,iDim);
                Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
            }
            
            
            Density_i = nodes->GetDensity(iPoint);
            
            Energy_i = nodes->GetEnergy(iPoint);
            StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
            
            GetFluidModel()->SetTDState_rhoe(Density_i, StaticEnergy_i);
            
            Pressure_i = GetFluidModel()->GetPressure();
            Enthalpy_i = Energy_i + Pressure_i/Density_i;
            
            SoundSpeed_i = GetFluidModel()->GetSoundSpeed();
            
            Kappa_i = GetFluidModel()->GetdPde_rho() / Density_i;
            Chi_i = GetFluidModel()->GetdPdrho_e() - Kappa_i * StaticEnergy_i;
            
            ProjVelocity_i = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
                ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
            
            /*--- Build the external state u_e from boundary data and internal node ---*/
            
            switch(config->GetKind_Data_Riemann(Marker_Tag))
            {
                case TOTAL_CONDITIONS_PT:
                    
                    /*--- Retrieve the specified total conditions for this boundary. ---*/
                    if (gravity) P_Total = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;/// check in which case is true (only freesurface?)
                    else P_Total  = config->GetRiemann_Var1(Marker_Tag);
                    T_Total  = config->GetRiemann_Var2(Marker_Tag);
                    Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);
                    
                    /*--- Non-dim. the inputs if necessary. ---*/
                    P_Total /= config->GetPressure_Ref();
                    T_Total /= config->GetTemperature_Ref();
                    
                    /* --- Computes the total state --- */
                    GetFluidModel()->SetTDState_PT(P_Total, T_Total);
                    Enthalpy_e = GetFluidModel()->GetStaticEnergy()+ GetFluidModel()->GetPressure()/GetFluidModel()->GetDensity();
                    Entropy_e = GetFluidModel()->GetEntropy();
                    
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
                    }
                    StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
                    GetFluidModel()->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
                    Density_e = GetFluidModel()->GetDensity();
                    StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
                    Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
                    if (tkeNeeded) Energy_e += GetTke_Inf();
                    break;
                    
                case STATIC_SUPERSONIC_INFLOW_PT:
                    
                    /*--- Retrieve the specified total conditions for this boundary. ---*/
                    if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;/// check in which case is true (only freesurface?)
                    else P_static  = config->GetRiemann_Var1(Marker_Tag);
                    T_static  = config->GetRiemann_Var2(Marker_Tag);
                    Mach = config->GetRiemann_FlowDir(Marker_Tag);
                    
                    /*--- Non-dim. the inputs if necessary. ---*/
                    P_static /= config->GetPressure_Ref();
                    T_static /= config->GetTemperature_Ref();
                    
                    /* --- Computes the total state --- */
                    GetFluidModel()->SetTDState_PT(P_static, T_static);
                    
                    /* --- Compute the boundary state u_e --- */
                    Velocity2_e = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Velocity_e[iDim] = Mach[iDim]*GetFluidModel()->GetSoundSpeed();
                        Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                    }
                    Density_e = GetFluidModel()->GetDensity();
                    StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
                    Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
                    if (tkeNeeded) Energy_e += GetTke_Inf();
                    break;
                    
                case STATIC_SUPERSONIC_INFLOW_PD:
                    
                    /*--- Retrieve the specified total conditions for this boundary. ---*/
                    
                    if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;/// check in which case is true (only freesurface?)
                    else P_static  = config->GetRiemann_Var1(Marker_Tag);
                    Rho_static  = config->GetRiemann_Var2(Marker_Tag);
                    Mach = config->GetRiemann_FlowDir(Marker_Tag);
                    
                    /*--- Non-dim. the inputs if necessary. ---*/
                    P_static /= config->GetPressure_Ref();
                    Rho_static /= config->GetDensity_Ref();
                    
                    /* --- Computes the total state --- */
                    GetFluidModel()->SetTDState_Prho(P_static, Rho_static);
                    
                    /* --- Compute the boundary state u_e --- */
                    Velocity2_e = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Velocity_e[iDim] = Mach[iDim]*GetFluidModel()->GetSoundSpeed();
                        Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                    }
                    Density_e = GetFluidModel()->GetDensity();
                    StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
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
                    
                case STATIC_PRESSURE:
                    
                    /*--- Retrieve the staic pressure for this boundary. ---*/
                    Pressure_e = config->GetRiemann_Var1(Marker_Tag);
                    Pressure_e /= config->GetPressure_Ref();
                    Density_e = Density_i;
                    
                    /* --- Compute the boundary state u_e --- */
                    GetFluidModel()->SetTDState_Prho(Pressure_e, Density_e);
                    Velocity2_e = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Velocity_e[iDim] = Velocity_i[iDim];
                        Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                    }
                    Energy_e = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_e;
                    break;
                    
                default:
                    SU2_MPI::Error("Invalid Riemann input!", CURRENT_FUNCTION);
                    break;
            }
            
            /*--- Compute P (matrix of right eigenvectors) ---*/
            conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);
            
            /*--- Compute inverse P (matrix of left eigenvectors)---*/
            conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);
            
            /*--- eigenvalues contribution due to grid motion ---*/
            if (dynamic_grid) {
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
            
            /*--- Compute the boundary state u_e ---*/
            u_e[0] = Density_e;
            for (iDim = 0; iDim < nDim; iDim++)
                u_e[iDim+1] = Velocity_e[iDim]*Density_e;
            u_e[nVar-1] = Energy_e*Density_e;
            
            /*--- Compute the boundary state u_i ---*/
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
            GetFluidModel()->SetTDState_rhoe(Density_b, StaticEnergy_b);
            Pressure_b = GetFluidModel()->GetPressure();
            Temperature_b = GetFluidModel()->GetTemperature();
            Enthalpy_b = Energy_b + Pressure_b/Density_b;
            Kappa_b = GetFluidModel()->GetdPde_rho() / Density_b;
            Chi_b = GetFluidModel()->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
            
            /*--- Compute the residuals ---*/
            conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);
            
            /*--- Residual contribution due to grid motion ---*/
            if (dynamic_grid) {
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
                if (dynamic_grid)
                {
                    gridVel = geometry->node[iPoint]->GetGridVel();
                    su2double projVelocity = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++)
                        projVelocity +=  gridVel[iDim]*Normal[iDim];
                    for (iVar = 0; iVar < nVar; iVar++) {
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
                Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
            
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
                
                V_boundary[nDim+5] = GetFluidModel()->GetLaminarViscosity();
                V_boundary[nDim+6] = nodes->GetEddyViscosity(iPoint);
                V_boundary[nDim+7] = GetFluidModel()->GetThermalConductivity();
                V_boundary[nDim+8] = GetFluidModel()->GetCp();
                
                /*--- Set the normal vector and the coordinates ---*/
                
                visc_numerics->SetNormal(Normal);
                visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
                
                /*--- Primitive variables, and gradient ---*/
                
                visc_numerics->SetPrimitive(V_domain, V_boundary);
                visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
                
                /*--- Secondary variables ---*/
                
                S_domain = nodes->GetSecondary(iPoint);
                
                /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
                
                S_boundary[0]= GetFluidModel()->GetdPdrho_e();
                S_boundary[1]= GetFluidModel()->GetdPde_rho();
                
                S_boundary[2]= GetFluidModel()->GetdTdrho_e();
                S_boundary[3]= GetFluidModel()->GetdTde_rho();
                
                /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
                
                S_boundary[4]= GetFluidModel()->Getdmudrho_T();
                S_boundary[5]= GetFluidModel()->GetdmudT_rho();
                
                S_boundary[6]= GetFluidModel()->Getdktdrho_T();
                S_boundary[7]= GetFluidModel()->GetdktdT_rho();
                
                visc_numerics->SetSecondary(S_domain, S_boundary);
                
                /*--- Turbulent kinetic energy ---*/
                
                if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
                    visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                                        solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
                
                /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
                
                visc_numerics->SetTauWall(-1.0, -1.0);
                
                /*--- Compute and update residual ---*/
                
                auto residual = visc_numerics->ComputeResidual(config);
                LinSysRes.SubtractBlock(iPoint, residual);
                
                /*--- Jacobian contribution for implicit integration ---*/
                
                if (implicit)
                    Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
                
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


void CEulerSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container,
                                   CNumerics *conv_numerics, CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker) {
    
    unsigned short iDim, iVar, jVar, kVar, iSpan;
    unsigned long iPoint, Point_Normal, oldVertex, iVertex;
    su2double P_Total, T_Total, *Flow_Dir;
    su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b, Temperature_b;
    su2double *Velocity_e, Velocity2_e, Enthalpy_e, Entropy_e, Energy_e = 0.0, StaticEnthalpy_e, StaticEnergy_e, Density_e = 0.0, Pressure_e;
    su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
    su2double ProjVelocity_i;
    su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **DubDu, *dw, *u_e, *u_i, *u_b;
    su2double *gridVel;
    su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
    su2double AverageEnthalpy, AverageEntropy;
    unsigned short  iZone  = config->GetiZone();
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
    unsigned short nSpanWiseSections = geometry->GetnSpanWiseSections(config->GetMarker_All_TurbomachineryFlag(val_marker));
    bool viscous = config->GetViscous();
    bool gravity = (config->GetGravityForce());
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    
    su2double *Normal, *turboNormal, *UnitNormal, *FlowDirMix, FlowDirMixMag, *turboVelocity;
    Normal = new su2double[nDim];
    turboNormal = new su2double[nDim];
    UnitNormal = new su2double[nDim];
    
    Velocity_i = new su2double[nDim];
    Velocity_b = new su2double[nDim];
    Velocity_e = new su2double[nDim];
    turboVelocity = new su2double[nDim];
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
    for (iSpan= 0; iSpan < nSpanWiseSections; iSpan++){
        for (iVertex = 0; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {
            
            /*--- using the other vertex information for retrieving some information ---*/
            oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();
            V_boundary= GetCharacPrimVar(val_marker, oldVertex);
            
            /*--- Index of the closest interior node ---*/
            Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();
            
            /*--- Normal vector for this vertex (negate for outward convention),
             *    this normal is scaled with the area of the face of the element  ---*/
            geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            conv_numerics->SetNormal(Normal);
            
            /*--- find the node related to the vertex ---*/
            iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();
            if (geometry->node[iPoint]->GetDomain()){
                
                /*--- Normalize Normal vector for this vertex (already for outward convention) ---*/
                geometry->turbovertex[val_marker][iSpan][iVertex]->GetNormal(UnitNormal);
                geometry->turbovertex[val_marker][iSpan][iVertex]->GetTurboNormal(turboNormal);
                
                /*--- Retrieve solution at this boundary node ---*/
                V_domain = nodes->GetPrimitive(iPoint);
                
                /* --- Compute the internal state u_i --- */
                Velocity2_i = 0;
                for (iDim=0; iDim < nDim; iDim++)
                {
                    Velocity_i[iDim] = nodes->GetVelocity(iPoint,iDim);
                    Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
                }
                
                Density_i = nodes->GetDensity(iPoint);
                
                Energy_i = nodes->GetEnergy(iPoint);
                StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
                
                GetFluidModel()->SetTDState_rhoe(Density_i, StaticEnergy_i);
                
                Pressure_i = GetFluidModel()->GetPressure();
                Enthalpy_i = Energy_i + Pressure_i/Density_i;
                
                SoundSpeed_i = GetFluidModel()->GetSoundSpeed();
                
                Kappa_i = GetFluidModel()->GetdPde_rho() / Density_i;
                Chi_i = GetFluidModel()->GetdPdrho_e() - Kappa_i * StaticEnergy_i;
                
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
                        if (gravity) P_Total = config->GetRiemann_Var1(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;/// check in which case is true (only freesurface?)
                        else P_Total  = config->GetRiemann_Var1(Marker_Tag);
                        T_Total  = config->GetRiemann_Var2(Marker_Tag);
                        Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);
                        
                        /*--- Non-dim. the inputs if necessary. ---*/
                        P_Total /= config->GetPressure_Ref();
                        T_Total /= config->GetTemperature_Ref();
                        
                        /* --- Computes the total state --- */
                        GetFluidModel()->SetTDState_PT(P_Total, T_Total);
                        Enthalpy_e = GetFluidModel()->GetStaticEnergy()+ GetFluidModel()->GetPressure()/GetFluidModel()->GetDensity();
                        Entropy_e = GetFluidModel()->GetEntropy();
                        
                        /* --- Compute the boundary state u_e --- */
                        Velocity2_e = Velocity2_i;
                        for (iDim = 0; iDim < nDim; iDim++)
                            turboVelocity[iDim] = sqrt(Velocity2_e)*Flow_Dir[iDim];
                        ComputeBackVelocity(turboVelocity,turboNormal, Velocity_e, config->GetMarker_All_TurbomachineryFlag(val_marker),config->GetKind_TurboMachinery(iZone));
                        StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
                        GetFluidModel()->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
                        Density_e = GetFluidModel()->GetDensity();
                        StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
                        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
                        if (tkeNeeded) Energy_e += GetTke_Inf();
                        break;
                        
                    case MIXING_IN:
                        
                        /* --- compute total averaged quantities ---*/
                        GetFluidModel()->SetTDState_Prho(ExtAveragePressure[val_marker][iSpan], ExtAverageDensity[val_marker][iSpan]);
                        AverageEnthalpy = GetFluidModel()->GetStaticEnergy() + ExtAveragePressure[val_marker][iSpan]/ExtAverageDensity[val_marker][iSpan];
                        AverageEntropy  = GetFluidModel()->GetEntropy();
                        
                        FlowDirMixMag = 0;
                        for (iDim = 0; iDim < nDim; iDim++)
                            FlowDirMixMag += ExtAverageTurboVelocity[val_marker][iSpan][iDim]*ExtAverageTurboVelocity[val_marker][iSpan][iDim];
                        for (iDim = 0; iDim < nDim; iDim++){
                            FlowDirMix[iDim] = ExtAverageTurboVelocity[val_marker][iSpan][iDim]/sqrt(FlowDirMixMag);
                        }
                        
                        
                        /* --- Computes the total state --- */
                        Enthalpy_e = AverageEnthalpy;
                        Entropy_e = AverageEntropy;
                        
                        /* --- Compute the boundary state u_e --- */
                        Velocity2_e = Velocity2_i;
                        for (iDim = 0; iDim < nDim; iDim++){
                            turboVelocity[iDim] = sqrt(Velocity2_e)*FlowDirMix[iDim];
                            
                        }
                        ComputeBackVelocity(turboVelocity,turboNormal, Velocity_e, config->GetMarker_All_TurbomachineryFlag(val_marker),config->GetKind_TurboMachinery(iZone));
                        
                        StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
                        GetFluidModel()->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
                        Density_e = GetFluidModel()->GetDensity();
                        StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
                        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
                        // if (tkeNeeded) Energy_e += GetTke_Inf();
                        break;
                        
                        
                    case MIXING_OUT:
                        
                        /*--- Retrieve the staic pressure for this boundary. ---*/
                        Pressure_e = ExtAveragePressure[val_marker][iSpan];
                        Density_e = Density_i;
                        
                        /* --- Compute the boundary state u_e --- */
                        GetFluidModel()->SetTDState_Prho(Pressure_e, Density_e);
                        Velocity2_e = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Velocity_e[iDim] = Velocity_i[iDim];
                            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                        }
                        Energy_e = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_e;
                        break;
                        
                    case STATIC_PRESSURE:
                        
                        /*--- Retrieve the staic pressure for this boundary. ---*/
                        Pressure_e = config->GetRiemann_Var1(Marker_Tag);
                        Pressure_e /= config->GetPressure_Ref();
                        Density_e = Density_i;
                        
                        /* --- Compute the boundary state u_e --- */
                        GetFluidModel()->SetTDState_Prho(Pressure_e, Density_e);
                        Velocity2_e = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Velocity_e[iDim] = Velocity_i[iDim];
                            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                        }
                        Energy_e = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_e;
                        break;
                        
                        
                    case RADIAL_EQUILIBRIUM:
                        
                        /*--- Retrieve the staic pressure for this boundary. ---*/
                        Pressure_e = RadialEquilibriumPressure[val_marker][iSpan];
                        Density_e = Density_i;
                        
                        /* --- Compute the boundary state u_e --- */
                        GetFluidModel()->SetTDState_Prho(Pressure_e, Density_e);
                        Velocity2_e = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Velocity_e[iDim] = Velocity_i[iDim];
                            Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
                        }
                        Energy_e = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_e;
                        break;
                        
                    default:
                        SU2_MPI::Error("Invalid Riemann input!", CURRENT_FUNCTION);
                        break;
                }
                
                /*--- Compute P (matrix of right eigenvectors) ---*/
                conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);
                
                /*--- Compute inverse P (matrix of left eigenvectors)---*/
                conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);
                
                /*--- eigenvalues contribution due to grid motion ---*/
                if (dynamic_grid){
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
                GetFluidModel()->SetTDState_rhoe(Density_b, StaticEnergy_b);
                Pressure_b = GetFluidModel()->GetPressure();
                Temperature_b = GetFluidModel()->GetTemperature();
                Enthalpy_b = Energy_b + Pressure_b/Density_b;
                Kappa_b = GetFluidModel()->GetdPde_rho() / Density_b;
                Chi_b = GetFluidModel()->GetdPdrho_e() - Kappa_b * StaticEnergy_b;
                
                /*--- Compute the residuals ---*/
                conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);
                
                /*--- Residual contribution due to grid motion ---*/
                if (dynamic_grid) {
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
                    if (dynamic_grid)
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
                    Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
                
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
                    
                    V_boundary[nDim+5] = GetFluidModel()->GetLaminarViscosity();
                    V_boundary[nDim+6] = nodes->GetEddyViscosity(iPoint);
                    V_boundary[nDim+7] = GetFluidModel()->GetThermalConductivity();
                    V_boundary[nDim+8] = GetFluidModel()->GetCp();
                    
                    /*--- Set the normal vector and the coordinates ---*/
                    
                    visc_numerics->SetNormal(Normal);
                    visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
                    
                    /*--- Primitive variables, and gradient ---*/
                    
                    visc_numerics->SetPrimitive(V_domain, V_boundary);
                    visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
                    
                    /*--- Secondary variables ---*/
                    
                    S_domain = nodes->GetSecondary(iPoint);
                    
                    /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
                    
                    S_boundary[0]= GetFluidModel()->GetdPdrho_e();
                    S_boundary[1]= GetFluidModel()->GetdPde_rho();
                    
                    S_boundary[2]= GetFluidModel()->GetdTdrho_e();
                    S_boundary[3]= GetFluidModel()->GetdTde_rho();
                    
                    /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
                    
                    S_boundary[4]= GetFluidModel()->Getdmudrho_T();
                    S_boundary[5]= GetFluidModel()->GetdmudT_rho();
                    
                    S_boundary[6]= GetFluidModel()->Getdktdrho_T();
                    S_boundary[7]= GetFluidModel()->GetdktdT_rho();
                    
                    visc_numerics->SetSecondary(S_domain, S_boundary);
                    
                    /*--- Turbulent kinetic energy ---*/
                    
                    if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
                        visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                                            solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
                    
                    /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
                    
                    visc_numerics->SetTauWall(-1.0, -1.0);
                    
                    /*--- Compute and update residual ---*/
                    
                    auto residual = visc_numerics->ComputeResidual(config);
                    LinSysRes.SubtractBlock(iPoint, residual);
                    
                    /*--- Jacobian contribution for implicit integration ---*/
                    
                    if (implicit)
                        Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
                    
                }
            }
        }
    }
    
    /*--- Free locally allocated memory ---*/
    delete [] Normal;
    delete [] UnitNormal;
    delete [] turboNormal;
    delete [] turboVelocity;
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

void CEulerSolver::PreprocessBC_Giles(CGeometry *geometry, CConfig *config, CNumerics *conv_numerics, unsigned short marker_flag) {
    /* Implementation of Fuorier Transformations for non-regfelcting BC will come soon */
    su2double cj_inf,cj_out1, cj_out2, Density_i, Pressure_i, *turboNormal, *turboVelocity, *Velocity_i, AverageSoundSpeed;
    su2double *deltaprim, *cj, TwoPiThetaFreq_Pitch, pitch, theta, deltaTheta;
    unsigned short iMarker, iSpan, iMarkerTP, iDim;
    unsigned long  iPoint, kend_max, k, iVertex;
    long freq;
    unsigned short  iZone     = config->GetiZone();
    unsigned short nSpanWiseSections = geometry->GetnSpanWiseSections(marker_flag);
    turboNormal   = new su2double[nDim];
    turboVelocity = new su2double[nDim];
    Velocity_i    = new su2double[nDim];
    deltaprim     = new su2double[nVar];
    cj            = new su2double[nVar];
    complex<su2double> I, cktemp_inf,cktemp_out1, cktemp_out2, expArg;
    I = complex<su2double>(0.0,1.0);
    
#ifdef HAVE_MPI
    su2double MyIm_inf, MyRe_inf, Im_inf, Re_inf, MyIm_out1, MyRe_out1, Im_out1, Re_out1, MyIm_out2, MyRe_out2, Im_out2, Re_out2;
#endif
    
    kend_max = geometry->GetnFreqSpanMax(marker_flag);
    for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
        for(k=0; k < 2*kend_max+1; k++){
            freq = k - kend_max;
            cktemp_inf = complex<su2double>(0.0,0.0);
            cktemp_out1 = complex<su2double>(0.0,0.0);
            cktemp_out2 = complex<su2double>(0.0,0.0);
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
                for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
                    if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                            for (iVertex = 0; iVertex < geometry->GetnVertexSpan(iMarker,iSpan); iVertex++) {
                                
                                /*--- find the node related to the vertex ---*/
                                iPoint = geometry->turbovertex[iMarker][iSpan][iVertex]->GetNode();
                                
                                geometry->turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(turboNormal);
                                /*--- Compute the internal state _i ---*/
                                
                                Pressure_i = nodes->GetPressure(iPoint);
                                Density_i = nodes->GetDensity(iPoint);
                                for (iDim = 0; iDim < nDim; iDim++)
                                {
                                    Velocity_i[iDim] = nodes->GetVelocity(iPoint,iDim);
                                }
                                ComputeTurboVelocity(Velocity_i, turboNormal, turboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));
                                
                                if(nDim ==2){
                                    deltaprim[0] = Density_i - AverageDensity[iMarker][iSpan];
                                    deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[iMarker][iSpan][0];
                                    deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[iMarker][iSpan][1];
                                    deltaprim[3] = Pressure_i - AveragePressure[iMarker][iSpan];
                                }
                                else{
                                    //Here 3d
                                    deltaprim[0] = Density_i - AverageDensity[iMarker][iSpan];
                                    deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[iMarker][iSpan][0];
                                    deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[iMarker][iSpan][1];
                                    deltaprim[3] = turboVelocity[2] - AverageTurboVelocity[iMarker][iSpan][2]; //New char
                                    deltaprim[4] = Pressure_i - AveragePressure[iMarker][iSpan];
                                }
                                
                                GetFluidModel()->SetTDState_Prho(AveragePressure[iMarker][iSpan], AverageDensity[iMarker][iSpan]);
                                AverageSoundSpeed = GetFluidModel()->GetSoundSpeed();
                                conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[iMarker][iSpan], deltaprim, cj);
                                
                                /*-----this is only valid 2D ----*/
                                if(nDim ==2){
                                    cj_out1 = cj[1];
                                    cj_out2 = cj[2];
                                    cj_inf  = cj[3];
                                }
                                else{
                                    //Here 3D
                                    cj_out1 = cj[1];
                                    cj_out2 = cj[3];
                                    cj_inf  = cj[4];
                                }
                                pitch      = geometry->GetMaxAngularCoord(iMarker, iSpan) - geometry->GetMinAngularCoord(iMarker,iSpan);
                                theta      = geometry->turbovertex[iMarker][iSpan][iVertex]->GetRelAngularCoord();
                                deltaTheta = geometry->turbovertex[iMarker][iSpan][iVertex]->GetDeltaAngularCoord();
                                TwoPiThetaFreq_Pitch = 2*PI_NUMBER*freq*theta/pitch;
                                
                                expArg = complex<su2double>(cos(TwoPiThetaFreq_Pitch)) - I*complex<su2double>(sin(TwoPiThetaFreq_Pitch));
                                if (freq != 0){
                                    cktemp_out1 +=  cj_out1*expArg*deltaTheta/pitch;
                                    cktemp_out2 +=  cj_out2*expArg*deltaTheta/pitch;
                                    cktemp_inf  +=  cj_inf*expArg*deltaTheta/pitch;
                                }
                                else{
                                    cktemp_inf += complex<su2double>(0.0,0.0);
                                    cktemp_out1 += complex<su2double>(0.0,0.0);
                                    cktemp_out2 += complex<su2double>(0.0,0.0);
                                }
                            }
                            
                        }
                    }
                }
            }
            
#ifdef HAVE_MPI
            MyRe_inf = cktemp_inf.real(); Re_inf = 0.0;
            MyIm_inf = cktemp_inf.imag(); Im_inf = 0.0;
            cktemp_inf = complex<su2double>(0.0,0.0);
            
            MyRe_out1 = cktemp_out1.real(); Re_out1 = 0.0;
            MyIm_out1 = cktemp_out1.imag(); Im_out1 = 0.0;
            cktemp_out1 = complex<su2double>(0.0,0.0);
            
            MyRe_out2 = cktemp_out2.real(); Re_out2 = 0.0;
            MyIm_out2 = cktemp_out2.imag(); Im_out2 = 0.0;
            cktemp_out2 = complex<su2double>(0.0,0.0);
            
            
            SU2_MPI::Allreduce(&MyRe_inf, &Re_inf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            SU2_MPI::Allreduce(&MyIm_inf, &Im_inf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            SU2_MPI::Allreduce(&MyRe_out1, &Re_out1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            SU2_MPI::Allreduce(&MyIm_out1, &Im_out1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            SU2_MPI::Allreduce(&MyRe_out2, &Re_out2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            SU2_MPI::Allreduce(&MyIm_out2, &Im_out2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            
            cktemp_inf = complex<su2double>(Re_inf,Im_inf);
            cktemp_out1 = complex<su2double>(Re_out1,Im_out1);
            cktemp_out2 = complex<su2double>(Re_out2,Im_out2);
            
#endif
            
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
                for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
                    if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                            /*-----this is only valid 2D ----*/
                            if (marker_flag == INFLOW){
                                CkInflow[iMarker][iSpan][k]= cktemp_inf;
                            }else{
                                CkOutflow1[iMarker][iSpan][k]=cktemp_out1;
                                CkOutflow2[iMarker][iSpan][k]=cktemp_out2;
                            }
                        }
                    }
                }
            }
        }
    }
    
    delete [] turboVelocity;
    delete [] turboNormal;
    delete [] Velocity_i;
    delete [] deltaprim;
    delete [] cj;
    
}

void CEulerSolver::BC_Giles(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                            CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
    
    unsigned short iDim, iVar, jVar, iSpan;
    unsigned long  iPoint, Point_Normal, oldVertex, k, kend, kend_max, iVertex;
    su2double  *UnitNormal, *turboVelocity, *turboNormal;
    
    su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, Density_b, Pressure_b, Temperature_b;
    su2double *Velocity_i, Velocity2_i, Energy_i, StaticEnergy_i, Density_i, Pressure_i;
    su2double Pressure_e;
    su2double *V_boundary, *V_domain, *S_boundary, *S_domain;
    unsigned short  iZone     = config->GetiZone();
    bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
    bool viscous              = config->GetViscous();
    unsigned short nSpanWiseSections = geometry->GetnSpanWiseSections(config->GetMarker_All_TurbomachineryFlag(val_marker));
    su2double relfacAvgCfg       = config->GetGiles_RelaxFactorAverage(Marker_Tag);
    su2double relfacFouCfg       = config->GetGiles_RelaxFactorFourier(Marker_Tag);
    su2double *Normal;
    su2double TwoPiThetaFreq_Pitch, pitch,theta;
    const su2double *SpanWiseValues = NULL;
    su2double spanPercent, extrarelfacAvg = 0.0, deltaSpan = 0.0, relfacAvg, relfacFou, coeffrelfacAvg = 0.0;
    unsigned short Turbo_Flag;
    
    Normal                = new su2double[nDim];
    turboNormal           = new su2double[nDim];
    UnitNormal            = new su2double[nDim];
    turboVelocity         = new su2double[nDim];
    Velocity_i            = new su2double[nDim];
    Velocity_b            = new su2double[nDim];
    
    
    su2double AverageSoundSpeed, *AverageTurboMach, AverageEntropy, AverageEnthalpy;
    AverageTurboMach = new su2double[nDim];
    S_boundary       = new su2double[8];
    
    su2double  AvgMach , *cj, GilesBeta, *delta_c, **R_Matrix, *deltaprim, **R_c_inv,**R_c, alphaIn_BC, gammaIn_BC = 0,
    P_Total, T_Total, *FlowDir, Enthalpy_BC, Entropy_BC, *R, *c_avg,*dcjs, Beta_inf2, c2js_Re, c3js_Re, cOutjs_Re, avgVel2 =0.0;
    
    long freq;
    
    delta_c       = new su2double[nVar];
    deltaprim     = new su2double[nVar];
    cj            = new su2double[nVar];
    R_Matrix      = new su2double*[nVar];
    R_c           = new su2double*[nVar-1];
    R_c_inv       = new su2double*[nVar-1];
    R             = new su2double[nVar-1];
    c_avg         = new su2double[nVar];
    dcjs          = new su2double[nVar];
    
    for (iVar = 0; iVar < nVar; iVar++)
    {
        R_Matrix[iVar] = new su2double[nVar];
        c_avg[iVar]    =  0.0;
        dcjs[iVar]     =  0.0;
    }
    for (iVar = 0; iVar < nVar-1; iVar++)
    {
        R_c[iVar] = new su2double[nVar-1];
        R_c_inv[iVar] = new su2double[nVar-1];
    }
    
    
    complex<su2double> I, c2ks, c2js, c3ks, c3js, c4ks, c4js, cOutks, cOutjs, Beta_inf;
    I = complex<su2double>(0.0,1.0);
    
    /*--- Compute coeff for under relaxation of Avg and Fourier Coefficient for hub and shroud---*/
    if (nDim == 3){
        extrarelfacAvg  = config->GetExtraRelFacGiles(0);
        spanPercent     = config->GetExtraRelFacGiles(1);
        Turbo_Flag      = config->GetMarker_All_TurbomachineryFlag(val_marker);
        SpanWiseValues  = geometry->GetSpanWiseValue(Turbo_Flag);
        deltaSpan       = SpanWiseValues[nSpanWiseSections-1]*spanPercent;
        coeffrelfacAvg  = (relfacAvgCfg - extrarelfacAvg)/deltaSpan;
    }
    
    for (iSpan= 0; iSpan < nSpanWiseSections ; iSpan++){
        /*--- Compute under relaxation for the Hub and Shroud Avg and Fourier Coefficient---*/
        if(nDim == 3){
            if(SpanWiseValues[iSpan] <= SpanWiseValues[0] + deltaSpan){
                relfacAvg = extrarelfacAvg + coeffrelfacAvg*(SpanWiseValues[iSpan] - SpanWiseValues[0]);
                relfacFou = 0.0;
            }
            else if(SpanWiseValues[iSpan] >= SpanWiseValues[nSpanWiseSections -1] - deltaSpan){
                relfacAvg = extrarelfacAvg - coeffrelfacAvg*(SpanWiseValues[iSpan] - SpanWiseValues[nSpanWiseSections -1]);
                relfacFou = 0.0;
            }
            else{
                relfacAvg = relfacAvgCfg;
                relfacFou = relfacFouCfg;
            }
        }
        else{
            {
                relfacAvg = relfacAvgCfg;
                relfacFou = relfacFouCfg;
            }
        }
        
        GetFluidModel()->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
        AverageSoundSpeed = GetFluidModel()->GetSoundSpeed();
        AverageTurboMach[0] = AverageTurboVelocity[val_marker][iSpan][0]/AverageSoundSpeed;
        AverageTurboMach[1] = AverageTurboVelocity[val_marker][iSpan][1]/AverageSoundSpeed;
        
        if(dynamic_grid){
            AverageTurboMach[1] -= geometry->GetAverageTangGridVel(val_marker,iSpan)/AverageSoundSpeed;
        }
        
        AvgMach = AverageTurboMach[0]*AverageTurboMach[0] + AverageTurboMach[1]*AverageTurboMach[1];
        
        kend     = geometry->GetnFreqSpan(val_marker, iSpan);
        kend_max = geometry->GetnFreqSpanMax(config->GetMarker_All_TurbomachineryFlag(val_marker));
        conv_numerics->GetRMatrix(AverageSoundSpeed, AverageDensity[val_marker][iSpan], R_Matrix);
        
        switch(config->GetKind_Data_Giles(Marker_Tag)){
                
            case TOTAL_CONDITIONS_PT:
                
                /*--- Retrieve the specified total conditions for this inlet. ---*/
                P_Total  = config->GetGiles_Var1(Marker_Tag);
                T_Total  = config->GetGiles_Var2(Marker_Tag);
                FlowDir = config->GetGiles_FlowDir(Marker_Tag);
                alphaIn_BC = atan(FlowDir[1]/FlowDir[0]);
                
                gammaIn_BC = 0;
                if (nDim == 3){
                    gammaIn_BC = FlowDir[2]; //atan(FlowDir[2]/FlowDir[0]);
                }
                
                /*--- Non-dim. the inputs---*/
                P_Total /= config->GetPressure_Ref();
                T_Total /= config->GetTemperature_Ref();
                
                /* --- Computes the total state --- */
                GetFluidModel()->SetTDState_PT(P_Total, T_Total);
                Enthalpy_BC = GetFluidModel()->GetStaticEnergy()+ GetFluidModel()->GetPressure()/GetFluidModel()->GetDensity();
                Entropy_BC = GetFluidModel()->GetEntropy();
                
                
                /* --- Computes the inverse matrix R_c --- */
                conv_numerics->ComputeResJacobianGiles(GetFluidModel(), AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan],
                                                       AverageTurboVelocity[val_marker][iSpan], alphaIn_BC, gammaIn_BC, R_c, R_c_inv);
                
                GetFluidModel()->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
                AverageEnthalpy = GetFluidModel()->GetStaticEnergy() + AveragePressure[val_marker][iSpan]/AverageDensity[val_marker][iSpan];
                AverageEntropy  = GetFluidModel()->GetEntropy();
                
                avgVel2 = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) avgVel2 += AverageVelocity[val_marker][iSpan][iDim]*AverageVelocity[val_marker][iSpan][iDim];
                if (nDim == 2){
                    R[0] = -(AverageEntropy - Entropy_BC);
                    R[1] = -(AverageTurboVelocity[val_marker][iSpan][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][iSpan][0]);
                    R[2] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);
                }
                
                else{
                    R[0] = -(AverageEntropy - Entropy_BC);
                    R[1] = -(AverageTurboVelocity[val_marker][iSpan][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][iSpan][0]);
                    R[2] = -(AverageTurboVelocity[val_marker][iSpan][2] - tan(gammaIn_BC)*AverageTurboVelocity[val_marker][iSpan][0]);
                    R[3] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);
                    
                }
                /* --- Compute the avg component  c_avg = R_c^-1 * R --- */
                for (iVar = 0; iVar < nVar-1; iVar++){
                    c_avg[iVar] = 0.0;
                    for (jVar = 0; jVar < nVar-1; jVar++){
                        c_avg[iVar] += R_c_inv[iVar][jVar]*R[jVar];
                    }
                }
                break;
                
            case TOTAL_CONDITIONS_PT_1D:
                
                /*--- Retrieve the specified total conditions for this inlet. ---*/
                P_Total  = config->GetGiles_Var1(Marker_Tag);
                T_Total  = config->GetGiles_Var2(Marker_Tag);
                FlowDir = config->GetGiles_FlowDir(Marker_Tag);
                alphaIn_BC = atan(FlowDir[1]/FlowDir[0]);
                
                gammaIn_BC = 0;
                if (nDim == 3){
                    // Review definition of angle
                    gammaIn_BC = FlowDir[2]; //atan(FlowDir[2]/FlowDir[0]);
                }
                
                /*--- Non-dim. the inputs---*/
                P_Total /= config->GetPressure_Ref();
                T_Total /= config->GetTemperature_Ref();
                
                /* --- Computes the total state --- */
                GetFluidModel()->SetTDState_PT(P_Total, T_Total);
                Enthalpy_BC = GetFluidModel()->GetStaticEnergy()+ GetFluidModel()->GetPressure()/GetFluidModel()->GetDensity();
                Entropy_BC = GetFluidModel()->GetEntropy();
                
                
                /* --- Computes the inverse matrix R_c --- */
                conv_numerics->ComputeResJacobianGiles(GetFluidModel(), AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan],
                                                       AverageTurboVelocity[val_marker][iSpan], alphaIn_BC, gammaIn_BC, R_c, R_c_inv);
                
                GetFluidModel()->SetTDState_Prho(AveragePressure[val_marker][nSpanWiseSections], AverageDensity[val_marker][nSpanWiseSections]);
                AverageEnthalpy = GetFluidModel()->GetStaticEnergy() + AveragePressure[val_marker][nSpanWiseSections]/AverageDensity[val_marker][nSpanWiseSections];
                AverageEntropy  = GetFluidModel()->GetEntropy();
                
                
                avgVel2 = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) avgVel2 += AverageVelocity[val_marker][iSpan][iDim]*AverageVelocity[val_marker][iSpan][iDim];
                if (nDim == 2){
                    R[0] = -(AverageEntropy - Entropy_BC);
                    R[1] = -(AverageTurboVelocity[val_marker][nSpanWiseSections][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][nSpanWiseSections][0]);
                    R[2] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);
                }
                
                else{
                    R[0] = -(AverageEntropy - Entropy_BC);
                    R[1] = -(AverageTurboVelocity[val_marker][nSpanWiseSections][1] - tan(alphaIn_BC)*AverageTurboVelocity[val_marker][nSpanWiseSections][0]);
                    R[2] = -(AverageTurboVelocity[val_marker][nSpanWiseSections][2] - tan(gammaIn_BC)*AverageTurboVelocity[val_marker][nSpanWiseSections][0]);
                    R[3] = -(AverageEnthalpy + 0.5*avgVel2 - Enthalpy_BC);
                    
                }
                /* --- Compute the avg component  c_avg = R_c^-1 * R --- */
                for (iVar = 0; iVar < nVar-1; iVar++){
                    c_avg[iVar] = 0.0;
                    for (jVar = 0; jVar < nVar-1; jVar++){
                        c_avg[iVar] += R_c_inv[iVar][jVar]*R[jVar];
                    }
                }
                break;
                
            case MIXING_IN: case MIXING_OUT:
                
                /* --- Compute average jump of primitive at the mixing-plane interface--- */
                deltaprim[0] = ExtAverageDensity[val_marker][iSpan] - AverageDensity[val_marker][iSpan];
                deltaprim[1] = ExtAverageTurboVelocity[val_marker][iSpan][0] - AverageTurboVelocity[val_marker][iSpan][0];
                deltaprim[2] = ExtAverageTurboVelocity[val_marker][iSpan][1] - AverageTurboVelocity[val_marker][iSpan][1];
                if (nDim == 2){
                    deltaprim[3] = ExtAveragePressure[val_marker][iSpan] - AveragePressure[val_marker][iSpan];
                }
                else
                {
                    deltaprim[3] = ExtAverageTurboVelocity[val_marker][iSpan][2] - AverageTurboVelocity[val_marker][iSpan][2];
                    deltaprim[4] = ExtAveragePressure[val_marker][iSpan] - AveragePressure[val_marker][iSpan];
                }
                
                
                /* --- Compute average jump of charachteristic variable at the mixing-plane interface--- */
                GetFluidModel()->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
                AverageSoundSpeed = GetFluidModel()->GetSoundSpeed();
                conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[val_marker][iSpan], deltaprim, c_avg);
                break;
                
            case MIXING_IN_1D: case MIXING_OUT_1D:
                
                /* --- Compute average jump of primitive at the mixing-plane interface--- */
                deltaprim[0] = ExtAverageDensity[val_marker][nSpanWiseSections] - AverageDensity[val_marker][nSpanWiseSections];
                deltaprim[1] = ExtAverageTurboVelocity[val_marker][nSpanWiseSections][0] - AverageTurboVelocity[val_marker][nSpanWiseSections][0];
                deltaprim[2] = ExtAverageTurboVelocity[val_marker][nSpanWiseSections][1] - AverageTurboVelocity[val_marker][nSpanWiseSections][1];
                if (nDim == 2){
                    deltaprim[3] = ExtAveragePressure[val_marker][nSpanWiseSections] - AveragePressure[val_marker][nSpanWiseSections];
                }
                else
                {
                    deltaprim[3] = ExtAverageTurboVelocity[val_marker][nSpanWiseSections][2] - AverageTurboVelocity[val_marker][nSpanWiseSections][2];
                    deltaprim[4] = ExtAveragePressure[val_marker][nSpanWiseSections] - AveragePressure[val_marker][nSpanWiseSections];
                }
                
                /* --- Compute average jump of charachteristic variable at the mixing-plane interface--- */
                GetFluidModel()->SetTDState_Prho(AveragePressure[val_marker][nSpanWiseSections], AverageDensity[val_marker][nSpanWiseSections]);
                AverageSoundSpeed = GetFluidModel()->GetSoundSpeed();
                conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[val_marker][nSpanWiseSections], deltaprim, c_avg);
                break;
                
                
            case STATIC_PRESSURE:
                Pressure_e = config->GetGiles_Var1(Marker_Tag);
                Pressure_e /= config->GetPressure_Ref();
                
                /* --- Compute avg characteristic jump  --- */
                if (nDim == 2){
                    c_avg[3] = -2.0*(AveragePressure[val_marker][iSpan]-Pressure_e);
                }
                else
                {
                    c_avg[4] = -2.0*(AveragePressure[val_marker][iSpan]-Pressure_e);
                }
                break;
                
            case STATIC_PRESSURE_1D:
                Pressure_e = config->GetGiles_Var1(Marker_Tag);
                Pressure_e /= config->GetPressure_Ref();
                
                /* --- Compute avg characteristic jump  --- */
                if (nDim == 2){
                    c_avg[3] = -2.0*(AveragePressure[val_marker][nSpanWiseSections]-Pressure_e);
                }
                else
                {
                    c_avg[4] = -2.0*(AveragePressure[val_marker][nSpanWiseSections]-Pressure_e);
                }
                break;
                
            case RADIAL_EQUILIBRIUM:
                Pressure_e = RadialEquilibriumPressure[val_marker][iSpan];
                
                /* --- Compute avg characteristic jump  --- */
                c_avg[4] = -2.0*(AveragePressure[val_marker][iSpan]-Pressure_e);
                
                break;
                
        }
        
        /*--- Loop over all the vertices on this boundary marker ---*/
        
        for (iVertex = 0; iVertex < geometry->GetnVertexSpan(val_marker,iSpan); iVertex++) {
            
            /*--- using the other vertex information for retrieving some information ---*/
            oldVertex = geometry->turbovertex[val_marker][iSpan][iVertex]->GetOldVertex();
            V_boundary= GetCharacPrimVar(val_marker, oldVertex);
            
            /*--- Index of the closest interior node ---*/
            Point_Normal = geometry->vertex[val_marker][oldVertex]->GetNormal_Neighbor();
            
            /*--- Normal vector for this vertex (negate for outward convention),
             *    this normal is scaled with the area of the face of the element  ---*/
            geometry->vertex[val_marker][oldVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            conv_numerics->SetNormal(Normal);
            
            /*--- find the node related to the vertex ---*/
            iPoint = geometry->turbovertex[val_marker][iSpan][iVertex]->GetNode();
            
            /*--- Normalize Normal vector for this vertex (already for outward convention) ---*/
            geometry->turbovertex[val_marker][iSpan][iVertex]->GetNormal(UnitNormal);
            geometry->turbovertex[val_marker][iSpan][iVertex]->GetTurboNormal(turboNormal);
            
            /*--- Retrieve solution at this boundary node ---*/
            V_domain = nodes->GetPrimitive(iPoint);
            
            /*--- Retrieve domain Secondary variables ---*/
            S_domain = nodes->GetSecondary(iPoint);
            
            
            /*--- Compute the internal state u_i ---*/
            Velocity2_i = 0;
            for (iDim = 0; iDim < nDim; iDim++)
            {
                Velocity_i[iDim] = nodes->GetVelocity(iPoint,iDim);
                Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
            }
            
            
            Density_i = nodes->GetDensity(iPoint);
            
            Energy_i = nodes->GetEnergy(iPoint);
            StaticEnergy_i = Energy_i - 0.5*Velocity2_i;
            
            GetFluidModel()->SetTDState_rhoe(Density_i, StaticEnergy_i);
            
            Pressure_i = GetFluidModel()->GetPressure();
            
            ComputeTurboVelocity(Velocity_i, turboNormal, turboVelocity, config->GetMarker_All_TurbomachineryFlag(val_marker),config->GetKind_TurboMachinery(iZone));
            if (nDim == 2){
                deltaprim[0] = Density_i - AverageDensity[val_marker][iSpan];
                deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[val_marker][iSpan][0];
                deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[val_marker][iSpan][1];
                deltaprim[3] = Pressure_i - AveragePressure[val_marker][iSpan];
            }
            else{
                deltaprim[0] = Density_i - AverageDensity[val_marker][iSpan];
                deltaprim[1] = turboVelocity[0] - AverageTurboVelocity[val_marker][iSpan][0];
                deltaprim[2] = turboVelocity[1] - AverageTurboVelocity[val_marker][iSpan][1];
                deltaprim[3] = turboVelocity[2] - AverageTurboVelocity[val_marker][iSpan][2];
                deltaprim[4] = Pressure_i - AveragePressure[val_marker][iSpan];
            }
            
            GetFluidModel()->SetTDState_Prho(AveragePressure[val_marker][iSpan], AverageDensity[val_marker][iSpan]);
            AverageSoundSpeed = GetFluidModel()->GetSoundSpeed();
            conv_numerics->GetCharJump(AverageSoundSpeed, AverageDensity[val_marker][iSpan], deltaprim, cj);
            
            
            pitch      = geometry->GetMaxAngularCoord(val_marker, iSpan) - geometry->GetMinAngularCoord(val_marker,iSpan);
            theta      = geometry->turbovertex[val_marker][iSpan][iVertex]->GetRelAngularCoord();
            
            switch(config->GetKind_Data_Giles(Marker_Tag))
            {
                    
                    //Done, generilize for 3D case
                    //TODO(turbo), generilize for Inlet and Outlet in for backflow treatment
                    
                case TOTAL_CONDITIONS_PT: case MIXING_IN:case TOTAL_CONDITIONS_PT_1D: case MIXING_IN_1D:
                    if(config->GetSpatialFourier()){
                        if (AvgMach <= 1.0){
                            Beta_inf= I*complex<su2double>(sqrt(1.0 - AvgMach));
                            c2js = complex<su2double>(0.0,0.0);
                            c3js = complex<su2double>(0.0,0.0);
                            for(k=0; k < 2*kend_max+1; k++){
                                freq = k - kend_max;
                                if(freq >= (long)(-kend) && freq <= (long)(kend) && AverageTurboMach[0] > config->GetAverageMachLimit()){
                                    TwoPiThetaFreq_Pitch = 2*PI_NUMBER*freq*theta/pitch;
                                    
                                    c2ks = -CkInflow[val_marker][iSpan][k]*complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>( 1.0 + AverageTurboMach[0]);
                                    c3ks =  CkInflow[val_marker][iSpan][k]*complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>( 1.0 + AverageTurboMach[0]);
                                    c3ks *= complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>( 1.0 + AverageTurboMach[0]);
                                    c2js += c2ks*(complex<su2double>(cos(TwoPiThetaFreq_Pitch))+I*complex<su2double>(sin(TwoPiThetaFreq_Pitch)));
                                    c3js += c3ks*(complex<su2double>(cos(TwoPiThetaFreq_Pitch))+I*complex<su2double>(sin(TwoPiThetaFreq_Pitch)));
                                }
                                else{
                                    c2js += complex<su2double>(0.0,0.0);
                                    c3js += complex<su2double>(0.0,0.0);
                                }
                            }
                            c2js_Re = c2js.real();
                            c3js_Re = c3js.real();
                            
                            if (nDim == 2){
                                dcjs[0] = 0.0     - cj[0];
                                dcjs[1] = c2js_Re - cj[1];
                                dcjs[2] = c3js_Re - cj[2];
                            }else{
                                dcjs[0] = 0.0     - cj[0];
                                dcjs[1] = c2js_Re - cj[1];
                                dcjs[2] = 0.0     - cj[2];
                                dcjs[3] = c3js_Re - cj[3];
                            }
                            
                        }else{
                            if (AverageTurboVelocity[val_marker][iSpan][1] >= 0.0){
                                Beta_inf2= -sqrt(AvgMach - 1.0);
                            }else{
                                Beta_inf2= sqrt(AvgMach-1.0);
                            }
                            if (nDim == 2){
                                c2js_Re = -cj[3]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
                                c3js_Re = cj[3]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
                                c3js_Re *= (Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
                            }else{
                                c2js_Re = -cj[4]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
                                c3js_Re = cj[4]*(Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
                                c3js_Re *= (Beta_inf2 + AverageTurboMach[1])/( 1.0 + AverageTurboMach[0]);
                            }
                            
                            
                            if (nDim == 2){
                                dcjs[0] = 0.0     - cj[0];
                                dcjs[1] = c2js_Re - cj[1];
                                dcjs[2] = c3js_Re - cj[2];
                            }else{
                                dcjs[0] = 0.0     - cj[0];
                                dcjs[1] = c2js_Re - cj[1];
                                dcjs[2] = 0.0     - cj[2];
                                dcjs[3] = c3js_Re - cj[3];
                            }
                        }
                    }
                    else{
                        if (nDim == 2){
                            dcjs[0] = 0.0;
                            dcjs[1] = 0.0;
                            dcjs[2] = 0.0;
                        }else{
                            dcjs[0] = 0.0;
                            dcjs[1] = 0.0;
                            dcjs[2] = 0.0;
                            dcjs[3] = 0.0;
                        }
                    }
                    /* --- Impose Inlet BC Reflecting--- */
                    delta_c[0] = relfacAvg*c_avg[0] + relfacFou*dcjs[0];
                    delta_c[1] = relfacAvg*c_avg[1] + relfacFou*dcjs[1];
                    delta_c[2] = relfacAvg*c_avg[2] + relfacFou*dcjs[2];
                    if (nDim == 2){
                        delta_c[3] = cj[3];
                    }else{
                        delta_c[3] = relfacAvg*c_avg[3] + relfacFou*dcjs[3];
                        delta_c[4] = cj[4];
                    }
                    break;
                    
                    
                case STATIC_PRESSURE:case STATIC_PRESSURE_1D:case MIXING_OUT:case RADIAL_EQUILIBRIUM:case MIXING_OUT_1D:
                    
                    /* --- implementation of Giles BC---*/
                    if(config->GetSpatialFourier()){
                        if (AvgMach > 1.0){
                            /* --- supersonic Giles implementation ---*/
                            if (AverageTurboVelocity[val_marker][iSpan][1] >= 0.0){
                                GilesBeta= -sqrt(AvgMach - 1.0);
                                
                            }else{
                                GilesBeta= sqrt(AvgMach - 1.0);
                            }
                            if(nDim == 2){
                                cOutjs_Re= (2.0 * AverageTurboMach[0])/(GilesBeta - AverageTurboMach[1])*cj[1] - (GilesBeta + AverageTurboMach[1])/(GilesBeta - AverageTurboMach[1])*cj[2];
                            }
                            else{
                                cOutjs_Re= (2.0 * AverageTurboMach[0])/(GilesBeta - AverageTurboMach[1])*cj[1] - (GilesBeta + AverageTurboMach[1])/(GilesBeta - AverageTurboMach[1])*cj[3];
                            }
                            if (nDim == 2){
                                dcjs[3] = cOutjs_Re - cj[3];
                            }
                            else{
                                dcjs[4] = cOutjs_Re - cj[4];
                            }
                        }else{
                            
                            /* --- subsonic Giles implementation ---*/
                            Beta_inf= I*complex<su2double>(sqrt(1.0  - AvgMach));
                            cOutjs  = complex<su2double>(0.0,0.0);
                            for(k=0; k < 2*kend_max+1; k++){
                                freq = k - kend_max;
                                if(freq >= (long)(-kend) && freq <= (long)(kend) && AverageTurboMach[0] > config->GetAverageMachLimit()){
                                    TwoPiThetaFreq_Pitch = 2*PI_NUMBER*freq*theta/pitch;
                                    cOutks  = complex<su2double>(2.0 * AverageTurboMach[0])/complex<su2double>(Beta_inf - AverageTurboMach[1])*CkOutflow1[val_marker][iSpan][k];
                                    cOutks -= complex<su2double>(Beta_inf + AverageTurboMach[1])/complex<su2double>(Beta_inf - AverageTurboMach[1])*CkOutflow2[val_marker][iSpan][k];
                                    
                                    cOutjs += cOutks*(complex<su2double>(cos(TwoPiThetaFreq_Pitch)) + I*complex<su2double>(sin(TwoPiThetaFreq_Pitch)));
                                }
                                else{
                                    cOutjs +=complex<su2double>(0.0,0.0);
                                }
                            }
                            cOutjs_Re = cOutjs.real();
                            
                            if (nDim == 2){
                                dcjs[3] = cOutjs_Re - cj[3];
                            }
                            else{
                                dcjs[4] = cOutjs_Re - cj[4];
                            }
                        }
                    }
                    else{
                        if (nDim == 2){
                            dcjs[3] = 0.0;
                        }
                        else{
                            dcjs[4] = 0.0;
                        }
                    }
                    /* --- Impose Outlet BC Non-Reflecting  --- */
                    delta_c[0] = cj[0];
                    delta_c[1] = cj[1];
                    delta_c[2] = cj[2];
                    if (nDim == 2){
                        delta_c[3] = relfacAvg*c_avg[3] + relfacFou*dcjs[3];
                    }
                    else{
                        delta_c[3] = cj[3];
                        delta_c[4] = relfacAvg*c_avg[4] + relfacFou*dcjs[4];
                    }
                    
                    
                    /*--- Automatically impose supersonic autoflow ---*/
                    if (abs(AverageTurboMach[0]) > 1.0000){
                        delta_c[0] = 0.0;
                        delta_c[1] = 0.0;
                        delta_c[2] = 0.0;
                        delta_c[2] = 0.0;
                        if (nDim == 3)delta_c[4] = 0.0;
                    }
                    
                    break;
                    
                default:
                    SU2_MPI::Error("Invalid Giles input!", CURRENT_FUNCTION);
                    break;
            }
            
            /*--- Compute primitive jump from characteristic variables  ---*/
            for (iVar = 0; iVar < nVar; iVar++)
            {
                deltaprim[iVar]=0.0;
                for (jVar = 0; jVar < nVar; jVar++)
                {
                    deltaprim[iVar] +=  R_Matrix[iVar][jVar]*delta_c[jVar];
                }
            }
            
            /*--- retrieve boundary variables ---*/
            Density_b = AverageDensity[val_marker][iSpan] + deltaprim[0];
            turboVelocity[0] = AverageTurboVelocity[val_marker][iSpan][0] + deltaprim[1];
            turboVelocity[1] = AverageTurboVelocity[val_marker][iSpan][1] + deltaprim[2];
            if(nDim == 2){
                Pressure_b = AveragePressure[val_marker][iSpan] + deltaprim[3];
            }
            else{
                turboVelocity[2] = AverageTurboVelocity[val_marker][iSpan][2] + deltaprim[3];
                Pressure_b = AveragePressure[val_marker][iSpan] + deltaprim[4];
            }
            
            
            ComputeBackVelocity(turboVelocity, turboNormal, Velocity_b, config->GetMarker_All_TurbomachineryFlag(val_marker), config->GetKind_TurboMachinery(iZone));
            Velocity2_b = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
                Velocity2_b+= Velocity_b[iDim]*Velocity_b[iDim];
            }
            
            if(Pressure_b <= 0.0 || Density_b <= 0.0 ){
                Pressure_b = Pressure_i;
                Density_b = Density_i;
                Velocity2_b = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Velocity_b[iDim] = Velocity_i[iDim];
                    Velocity2_b+= Velocity_b[iDim]*Velocity_b[iDim];
                }
            }
            
            GetFluidModel()->SetTDState_Prho(Pressure_b, Density_b);
            Energy_b = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_b;
            Temperature_b= GetFluidModel()->GetTemperature();
            Enthalpy_b = Energy_b + Pressure_b/Density_b;
            
            /*--- Primitive variables, using the derived quantities ---*/
            V_boundary[0] = Temperature_b;
            for (iDim = 0; iDim < nDim; iDim++)
                V_boundary[iDim+1] = Velocity_b[iDim];
            V_boundary[nDim+1] = Pressure_b;
            V_boundary[nDim+2] = Density_b;
            V_boundary[nDim+3] = Enthalpy_b;
            
            S_boundary[0]= GetFluidModel()->GetdPdrho_e();
            S_boundary[1]= GetFluidModel()->GetdPde_rho();
            
            
            
            /*--- Set various quantities in the solver class ---*/
            
            conv_numerics->SetPrimitive(V_domain, V_boundary);
            conv_numerics->SetSecondary(S_domain, S_boundary);
            
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            /*--- Update residual value ---*/
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration ---*/
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            /*--- Viscous contribution ---*/
            
            if (viscous) {
                
                /*--- Set laminar and eddy viscosity at the infinity ---*/
                
                V_boundary[nDim+5] = GetFluidModel()->GetLaminarViscosity();
                V_boundary[nDim+6] = nodes->GetEddyViscosity(iPoint);
                V_boundary[nDim+7] = GetFluidModel()->GetThermalConductivity();
                V_boundary[nDim+8] = GetFluidModel()->GetCp();
                
                /*--- Set the normal vector and the coordinates ---*/
                
                visc_numerics->SetNormal(Normal);
                visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
                
                /*--- Primitive variables, and gradient ---*/
                
                visc_numerics->SetPrimitive(V_domain, V_boundary);
                visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
                
                
                /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/
                
                S_boundary[0]= GetFluidModel()->GetdPdrho_e();
                S_boundary[1]= GetFluidModel()->GetdPde_rho();
                
                S_boundary[2]= GetFluidModel()->GetdTdrho_e();
                S_boundary[3]= GetFluidModel()->GetdTde_rho();
                
                /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/
                
                S_boundary[4]= GetFluidModel()->Getdmudrho_T();
                S_boundary[5]= GetFluidModel()->GetdmudT_rho();
                
                S_boundary[6]= GetFluidModel()->Getdktdrho_T();
                S_boundary[7]= GetFluidModel()->GetdktdT_rho();
                
                visc_numerics->SetSecondary(S_domain, S_boundary);
                
                /*--- Turbulent kinetic energy ---*/
                
                if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
                    visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                                        solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
                
                /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
                
                visc_numerics->SetTauWall(-1.0, -1.0);
                
                /*--- Compute and update residual ---*/
                
                auto residual = visc_numerics->ComputeResidual(config);
                LinSysRes.SubtractBlock(iPoint, residual);
                
                /*--- Jacobian contribution for implicit integration ---*/
                
                if (implicit)
                    Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
                
            }
            
        }
    }
    
    /*--- Free locally allocated memory ---*/
    delete [] Normal;
    
    delete [] Velocity_b;
    delete [] Velocity_i;
    
    delete [] S_boundary;
    delete [] delta_c;
    delete [] deltaprim;
    delete [] cj;
    for (iVar = 0; iVar < nVar; iVar++)
    {
        delete [] R_Matrix[iVar];
    }
    for (iVar = 0; iVar < nVar-1; iVar++)
    {
        delete [] R_c_inv[iVar];
        delete [] R_c[iVar];
    }
    delete [] R_Matrix;
    delete [] R_c;
    delete [] R_c_inv;
    delete [] R;
    delete [] c_avg;
    delete [] dcjs;
    
    delete [] AverageTurboMach;
    delete [] UnitNormal;
    delete [] turboNormal;
    delete [] turboVelocity;
}

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics,
                            CConfig *config, unsigned short val_marker) {
    unsigned short iDim;
    unsigned long iVertex, iPoint;
    su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
    Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
    alpha, aa, bb, cc, dd, Area, UnitNormal[3];
    su2double *V_inlet, *V_domain;
    
    bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    su2double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
    su2double Gas_Constant       = config->GetGas_ConstantND();
    unsigned short Kind_Inlet = config->GetKind_Inlet();
    string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    su2double *Normal = new su2double[nDim];
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        /*--- Allocate the value at the inlet ---*/
        
        V_inlet = GetCharacPrimVar(val_marker, iVertex);
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
        /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
        
        if (geometry->node[iPoint]->GetDomain()) {
            
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
            
            V_domain = nodes->GetPrimitive(iPoint);
            
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
                    
                    P_Total  = Inlet_Ptotal[val_marker][iVertex];
                    T_Total  = Inlet_Ttotal[val_marker][iVertex];
                    Flow_Dir = Inlet_FlowDir[val_marker][iVertex];
                    
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
                    
                    Density  = Inlet_Ttotal[val_marker][iVertex];
                    Vel_Mag  = Inlet_Ptotal[val_marker][iVertex];
                    Flow_Dir = Inlet_FlowDir[val_marker][iVertex];
                    
                    /*--- Non-dim. the inputs if necessary. ---*/
                    
                    Density /= config->GetDensity_Ref();
                    Vel_Mag /= config->GetVelocity_Ref();
                    
                    /*--- Get primitives from current inlet state. ---*/
                    
                    for (iDim = 0; iDim < nDim; iDim++)
                        Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);
                    Pressure    = nodes->GetPressure(iPoint);
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
            
            /*--- Set various quantities in the solver class ---*/
            
            conv_numerics->SetPrimitive(V_domain, V_inlet);
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            /*--- Update residual value ---*/
            
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration ---*/
            
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
            //
            //      if (viscous) {
            //
            //        /*--- Set laminar and eddy viscosity at the infinity ---*/
            //
            //        V_inlet[nDim+5] = nodes->GetLaminarViscosity(iPoint);
            //        V_inlet[nDim+6] = nodes->GetEddyViscosity(iPoint);
            //
            //        /*--- Set the normal vector and the coordinates ---*/
            //
            //        visc_numerics->SetNormal(Normal);
            //        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
            //
            //        /*--- Primitive variables, and gradient ---*/
            //
            //        visc_numerics->SetPrimitive(V_domain, V_inlet);
            //        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
            //
            //        /*--- Turbulent kinetic energy ---*/
            //
            //        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
            //          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
            //                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
            //
            //        /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
            //
            //        visc_numerics->SetTauWall(-1.0, -1.0);
            //
            //        /*--- Compute and update residual ---*/
            //
            //        auto residual = visc_numerics->ComputeResidual(config);
            //        LinSysRes.SubtractBlock(iPoint, residual);
            //
            //        /*--- Jacobian contribution for implicit integration ---*/
            //
            //        if (implicit)
            //          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
            //
            //      }
            
        }
    }
    
    /*--- Free locally allocated memory ---*/
    
    delete [] Normal;
    
}

void CEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics,
                             CConfig *config, unsigned short val_marker) {
    unsigned short iVar, iDim;
    unsigned long iVertex, iPoint;
    su2double Pressure, P_Exit, Velocity[3],
    Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
    Area, UnitNormal[3];
    su2double *V_outlet, *V_domain;
    bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    su2double Gas_Constant     = config->GetGas_ConstantND();
    string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
    bool gravity = (config->GetGravityForce());
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    su2double *Normal = new su2double[nDim];
  
    unsigned short Kind_Outlet = config->GetKind_Comp_Outlet(Marker_Tag);

    /*--- Loop over all the vertices on this boundary marker ---*/
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        /*--- Allocate the value at the outlet ---*/
        V_outlet = GetCharacPrimVar(val_marker, iVertex);
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
        /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
        if (geometry->node[iPoint]->GetDomain()) {
            
            /*--- Normal vector for this vertex (negate for outward convention) ---*/
            geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            conv_numerics->SetNormal(Normal);
            
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);
            
            for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = Normal[iDim]/Area;
            
            /*--- Current solution at this boundary node ---*/
            V_domain = nodes->GetPrimitive(iPoint);
            
            /*--- Compute a boundary value for the pressure depending on whether
                  we are prescribing a back pressure or a mass flow target. ---*/

            switch (Kind_Outlet) {

                case PRESSURE_OUTLET:

                    /*--- Build the fictitious intlet state based on characteristics ---*/
          
                    /*--- Retrieve the specified back pressure for this outlet. ---*/
                    if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->node[iPoint]->GetCoord(nDim-1)*STANDARD_GRAVITY;
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
                        for (iDim = 0; iDim < nDim; iDim++) V_outlet[iDim+1] = Velocity[iDim];
                        V_outlet[nDim+1] = Pressure;
                        V_outlet[nDim+2] = Density;
                        V_outlet[nDim+3] = Energy + Pressure/Density;
                    }
                    break;

                case MASS_FLOW_OUTLET:

                    /*--- Check whether the flow is supersonic at the exit. The type
                          of boundary update depends on this. ---*/

                    auto P_Inf  = config->GetPressure_FreeStream();
		    Density = V_domain[nDim + 2];
                    cout <<Density<<endl;
		    Velocity2 = 0.0; Vn = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Velocity[iDim] = V_domain[iDim + 1];
                        Velocity2 += Velocity[iDim]*Velocity[iDim];
                        Vn += Velocity[iDim]*UnitNormal[iDim];
                    }
                    Pressure   = V_domain[nDim + 1];
                    cout <<Pressure<<endl;
		    SoundSpeed = sqrt(Gamma*Pressure/Density);
                    Mach_Exit  = sqrt(Velocity2)/SoundSpeed;
            
                    /*--- Retrieve the specficied target mass flow at the outlet. ---*/
                    
		    su2double mDot_Target = config->GetOutlet_Pressure(Marker_Tag);
            
                    /*--- Compute the mass flow increment based on the difference
                          between the target mass flow and current mass flow. ---*/
            
                    //TODO: The same thing occurs here......Outlet_Area is not set yet.  See GetOutlet_Properties.
		    auto Area_Outlet = config->GetOutlet_Area(Marker_Tag);

                    auto mDot_Current = Density * Area_Outlet * Vn;
            
                    auto delta_mDot = mDot_Target - mDot_Current;
                    cout <<"Target: "<<mDot_Target<<endl;
		    cout <<"Area: "<<Area_Outlet<<endl;
		    cout <<"Curr: "<<mDot_Current<<endl;
		    
		    
		    
		    
		    su2double Damping, P_Outlet;

                    /*--- Update the new outlet mass flow. Note that we use damping
                          here to improve stability and convergence. ---*/
            
                    if (Mach_Exit >= 1.0) {
                
                        /*--- Supersonic exit flow: there are no incoming characteristics,
                              so no boundary condition is necessary. Set outlet state to current
                              state so that upwinding handles the direction of propagation. ---*/
                
                        Damping = 0.5;
                
                        if (mDot_Current < 0) { P_Outlet = Pressure * (1 - Damping * delta_mDot / mDot_Target); }
                        else                  { P_Outlet = P_Inf * 1.1;                                         }            
                    } else {
                
                        /*--- Subsonic exit flow: there is one incoming characteristic,
                              therefore one variable can be specified (back pressure). ---*/
                
                        Damping = 0.2;
                
                        if (mDot_Current < 0) { P_Outlet = Pressure * (1 - Damping * delta_mDot / mDot_Target); }
                        else                  { P_Outlet = P_Inf;                                               }
                    }
            
                    /*--- The pressure is prescribed at the outlet ---*/
                    V_outlet = V_domain; 
                    V_outlet[nDim + 1] = max(P_Inf, P_Outlet);
                    
                    break; 
            }

            /*--- Set various quantities in the solver class ---*/
            conv_numerics->SetPrimitive(V_domain, V_outlet);
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, residual);
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
        }
    }
    /*--- Free locally allocated memory ---*/
    delete [] Normal;
}

void CEulerSolver::BC_Massflow_Outlet(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics *conv_numerics, CNumerics *visc_numerics,
                                      CConfig *config, unsigned short val_marker) {
    unsigned short iVar, iDim;
    unsigned long iVertex, iPoint;
    su2double Pressure, Pt_Inf, P_Inf, P_Outlet, Velocity[3], Damping,
    Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
    Area, UnitNormal[3], Area_Outlet, mDot_Current, mDot_Target, delta_mDot;
    su2double *V_outlet, *V_domain;
    
    bool implicit          = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    su2double Gas_Constant = config->GetGas_ConstantND();
    su2double Mach         = config->GetMach();
    string Marker_Tag      = config->GetMarker_All_TagBound(val_marker);
    su2double *Normal      = new su2double[nDim];
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
        /*--- Allocate the value at the outlet ---*/
        V_outlet = GetCharacPrimVar(val_marker, iVertex);
        
        /*--- Get freestream pressure ---*/
        Pt_Inf = config->GetPressure_FreeStreamND();
        P_Inf  = config->GetPressure_FreeStream();
        
        /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
        if (geometry->node[iPoint]->GetDomain()) {
            
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
            V_domain = nodes->GetPrimitive(iPoint);
            
            /*--- Check whether the flow is supersonic at the exit. The type
             of boundary update depends on this. ---*/
            
            Density = V_domain[nDim + 2];
            Velocity2 = 0.0; Vn = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
                Velocity[iDim] = V_domain[iDim + 1];
                Velocity2 += Velocity[iDim]*Velocity[iDim];
                Vn += Velocity[iDim]*UnitNormal[iDim];
            }
            
            Pressure   = V_domain[nDim + 1];
            SoundSpeed = sqrt(Gamma*Pressure/Density);
            Mach_Exit  = sqrt(Velocity2)/SoundSpeed;
            
            /*--- Retrieve the specficied target mass flow at the outlet. ---*/
            
            mDot_Target = config->GetOutlet_MassFlow(Marker_Tag);
            
            /*--- Compute the mass flow increment based on the difference
             between the target mass flow and current mass flow. ---*/
            
            Area_Outlet = config->GetOutlet_Area(Marker_Tag);
            mDot_Current = Density * Area_Outlet * Vn;
            
            delta_mDot = mDot_Target - mDot_Current;
            
            /*--- Update the new outlet mass flow. Note that we use damping
             here to improve stability and convergence. ---*/
            
            if (Mach_Exit >= 1.0) {
                
                /*--- Supersonic exit flow: there are no incoming characteristics,
                 so no boundary condition is necessary. Set outlet state to current
                 state so that upwinding handles the direction of propagation. ---*/
                
                Damping = 0.5;
                
                if (mDot_Current < 0) {
                    P_Outlet = Pressure * (1 - Damping * delta_mDot / mDot_Target);
                } else {
                    P_Outlet = P_Inf * 1.1;
                }
                
            } else {
                
                /*--- Subsonic exit flow: there is one incoming characteristic,
                 therefore one variable can be specified (back pressure). ---*/
                
                Damping = 0.2;
                
                if (mDot_Current < 0) {
                    P_Outlet = Pressure * (1 - Damping * delta_mDot / mDot_Target);
                } else {
                    P_Outlet = P_Inf;
                }
            }
            
            /*--- The pressure is prescribed at the outlet ---*/
            
            V_outlet[nDim + 1] = max(P_Inf, P_Outlet);
            
            /*--- Set various quantities in the solver class ---*/
            conv_numerics->SetPrimitive(V_domain, V_outlet);
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, residual);
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
        }
    }
    
    /*--- Free locally allocated memory ---*/
    delete [] Normal;
    
}

void CEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                       CNumerics *conv_numerics, CNumerics *visc_numerics,
                                       CConfig *config, unsigned short val_marker) {
    unsigned short iDim;
    unsigned long iVertex, iPoint;
    su2double *V_inlet, *V_domain;
    
    su2double Density, Pressure, Temperature, Energy, *Vel, Velocity2;
    su2double Gas_Constant = config->GetGas_ConstantND();
    
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    su2double *Normal = new su2double[nDim];
    su2double *Velocity = new su2double[nDim];
    
    /*--- Supersonic inlet flow: there are no outgoing characteristics,
     so all flow variables can be imposed at the inlet.
     First, retrieve the specified values for the primitive variables. ---*/
    
    Temperature = config->GetInlet_Temperature(Marker_Tag);
    Pressure    = config->GetInlet_Pressure(Marker_Tag);
    Vel         = config->GetInlet_Velocity(Marker_Tag);
    
    /*--- Non-dim. the inputs if necessary. ---*/
    
    Temperature /= config->GetTemperature_Ref();
    Pressure    /= config->GetPressure_Ref();
    for (iDim = 0; iDim < nDim; iDim++)
        Velocity[iDim] = Vel[iDim] / config->GetVelocity_Ref();
    
    /*--- Density at the inlet from the gas law ---*/
    
    Density = Pressure/(Gas_Constant*Temperature);
    
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
            
            /*--- Current solution at this boundary node ---*/
            
            V_domain = nodes->GetPrimitive(iPoint);
            
            /*--- Normal vector for this vertex (negate for outward convention) ---*/
            
            geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            
            /*--- Set various quantities in the solver class ---*/
            
            conv_numerics->SetNormal(Normal);
            conv_numerics->SetPrimitive(V_domain, V_inlet);
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                          geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration ---*/
            
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
            //
            //      if (viscous) {
            //
            //        /*--- Set laminar and eddy viscosity at the infinity ---*/
            //
            //        V_inlet[nDim+5] = nodes->GetLaminarViscosity(iPoint);
            //        V_inlet[nDim+6] = nodes->GetEddyViscosity(iPoint);
            //
            //        /*--- Set the normal vector and the coordinates ---*/
            //
            //        visc_numerics->SetNormal(Normal);
            //        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
            //
            //        /*--- Primitive variables, and gradient ---*/
            //
            //        visc_numerics->SetPrimitive(V_domain, V_inlet);
            //        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
            //
            //        /*--- Turbulent kinetic energy ---*/
            //
            //        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
            //          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
            //                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
            //
            //        /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
            //
            //        visc_numerics->SetTauWall(-1.0, -1.0);
            //
            //        /*--- Compute and update residual ---*/
            //
            //        auto residual = visc_numerics->ComputeResidual(config);
            //        LinSysRes.SubtractBlock(iPoint, residual);
            //
            //        /*--- Jacobian contribution for implicit integration ---*/
            //
            //        if (implicit)
            //          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
            //      }
            
        }
    }
    
    /*--- Free locally allocated memory ---*/
    
    delete [] Normal;
    delete [] Velocity;
    
}

void CEulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                        CNumerics *conv_numerics, CNumerics *visc_numerics,
                                        CConfig *config, unsigned short val_marker) {
    unsigned short iDim;
    unsigned long iVertex, iPoint;
    su2double *V_outlet, *V_domain;
    
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
    
    su2double *Normal = new su2double[nDim];
    
    /*--- Supersonic outlet flow: there are no ingoing characteristics,
     so all flow variables can should be interpolated from the domain. ---*/
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
        /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
        
        if (geometry->node[iPoint]->GetDomain()) {
            
            /*--- Current solution at this boundary node ---*/
            
            V_domain = nodes->GetPrimitive(iPoint);
            
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
            
            V_domain = nodes->GetPrimitive(iPoint);
            
            /*--- Normal vector for this vertex (negate for outward convention) ---*/
            
            geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            
            /*--- Set various quantities in the solver class ---*/
            
            conv_numerics->SetNormal(Normal);
            conv_numerics->SetPrimitive(V_domain, V_outlet);
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                          geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration ---*/
            
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
            //
            //      if (viscous) {
            //
            //        /*--- Set laminar and eddy viscosity at the infinity ---*/
            //
            //        V_outlet[nDim+5] = nodes->GetLaminarViscosity(iPoint);
            //        V_outlet[nDim+6] = nodes->GetEddyViscosity(iPoint);
            //
            //        /*--- Set the normal vector and the coordinates ---*/
            //
            //        visc_numerics->SetNormal(Normal);
            //        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
            //
            //        /*--- Primitive variables, and gradient ---*/
            //
            //        visc_numerics->SetPrimitive(V_domain, V_outlet);
            //        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
            //
            //        /*--- Turbulent kinetic energy ---*/
            //
            //        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
            //          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
            //                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
            //
            //        /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
            //
            //        visc_numerics->SetTauWall(-1.0, -1.0);
            //
            //        /*--- Compute and update residual ---*/
            //
            //        auto residual = visc_numerics->ComputeResidual(config);
            //        LinSysRes.SubtractBlock(iPoint, residual);
            //
            //        /*--- Jacobian contribution for implicit integration ---*/
            //
            //        if (implicit)
            //          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
            //      }
            
        }
    }
    
    /*--- Free locally allocated memory ---*/
    
    delete [] Normal;
    
}

void CEulerSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
    
    unsigned short iDim;
    unsigned long iVertex, iPoint;
    su2double Pressure, Inflow_Pressure = 0.0, Velocity[3], Velocity2, Entropy, Target_Inflow_MassFlow = 0.0, Target_Inflow_Mach = 0.0, Density, Energy,
    Riemann, Area, UnitNormal[3], Vn, SoundSpeed, Vn_Exit, Inflow_Pressure_inc, Inflow_Pressure_old, Inflow_Mach_old, Inflow_MassFlow_old;
    su2double *V_inflow, *V_domain;
    
    su2double DampingFactor = config->GetDamp_Engine_Inflow();
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    unsigned short Kind_Engine_Inflow = config->GetKind_Engine_Inflow();
    su2double Gas_Constant = config->GetGas_ConstantND();
    string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
    bool Engine_HalfModel = config->GetEngine_HalfModel();
    
    su2double *Normal = new su2double[nDim];
    
    
    if (Kind_Engine_Inflow == FAN_FACE_MACH) {
        
        /*--- Retrieve the specified target fan face mach at the nacelle. ---*/
        
        Target_Inflow_Mach = config->GetEngineInflow_Target(Marker_Tag);
        
        /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
        
        Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
        Inflow_Mach_old = config->GetInflow_Mach(Marker_Tag);
        
        /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/
        
        Inflow_Pressure_inc = - (1.0 - (Inflow_Mach_old/Target_Inflow_Mach)) * Baseline_Press;
        
        /*--- Estimate the new fan face pressure ---*/
        
        Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);
        
    }
    
    if (Kind_Engine_Inflow == FAN_FACE_MDOT) {
        
        /*--- Retrieve the specified target mass flow (non-dimensional) at the nacelle. ---*/
        
        Target_Inflow_MassFlow = config->GetEngineInflow_Target(Marker_Tag) / (config->GetDensity_Ref() * config->GetVelocity_Ref());
        
        if (config->GetSystemMeasurements() == US) Target_Inflow_MassFlow /= 32.174;
        
        if (Engine_HalfModel) Target_Inflow_MassFlow /= 2.0;
        
        /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/
        
        Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
        Inflow_MassFlow_old = config->GetInflow_MassFlow(Marker_Tag);  // same here... it is a non dimensional value
        
        /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/
        
        Inflow_Pressure_inc = - (1.0 - (Inflow_MassFlow_old/Target_Inflow_MassFlow)) * Baseline_Press;
        
        /*--- Estimate the new fan face pressure ---*/
        
        Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);
        
    }
    
    /*--- No iterative scheme if we provide the static pressure ---*/
    
    if (Kind_Engine_Inflow == FAN_FACE_PRESSURE) {
        
        /*--- Retrieve the specified pressure (non-dimensional) at the nacelle. ---*/
        
        Inflow_Pressure = config->GetEngineInflow_Target(Marker_Tag) / config->GetPressure_Ref();
        
    }
    
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        /*--- Allocate the value at the outlet ---*/
        
        V_inflow = GetCharacPrimVar(val_marker, iVertex);
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
        /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
        
        if (geometry->node[iPoint]->GetDomain()) {
            
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
            
            V_domain = nodes->GetPrimitive(iPoint);
            
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
            V_inflow[nDim+4] = SoundSpeed;
            
            /*--- Set various quantities in the solver class ---*/
            
            conv_numerics->SetNormal(Normal);
            conv_numerics->SetPrimitive(V_domain, V_inflow);
            
            /*--- Set grid movement ---*/
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration ---*/
            
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
            //
            //      if (viscous) {
            //
            //        /*--- Set laminar and eddy viscosity at the infinity ---*/
            //
            //        V_inflow[nDim+5] = nodes->GetLaminarViscosity(iPoint);
            //        V_inflow[nDim+6] = nodes->GetEddyViscosity(iPoint);
            //
            //        /*--- Set the normal vector and the coordinates ---*/
            //
            //        visc_numerics->SetNormal(Normal);
            //        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
            //
            //        /*--- Primitive variables, and gradient ---*/
            //
            //        visc_numerics->SetPrimitive(V_domain, V_inflow);
            //        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
            //
            //        /*--- Turbulent kinetic energy ---*/
            //
            //        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
            //          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
            //                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
            //
            //        /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
            //
            //        visc_numerics->SetTauWall(-1.0, -1.0);
            //
            //        /*--- Compute and update residual ---*/
            //
            //        auto residual = visc_numerics->ComputeResidual(config);
            //        LinSysRes.SubtractBlock(iPoint, residual);
            //
            //        /*--- Jacobian contribution for implicit integration ---*/
            //
            //        if (implicit)
            //          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
            //
            //      }
            
        }
    }
    
    delete [] Normal;
    
}


void CEulerSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
    
    unsigned short iDim;
    unsigned long iVertex, iPoint;
    su2double Exhaust_Pressure, Exhaust_Temperature, Velocity[3], Velocity2, H_Exhaust, Temperature, Riemann, Area, UnitNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Exhaust2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
    su2double *V_exhaust, *V_domain, Target_Exhaust_Pressure, Exhaust_Pressure_old, Exhaust_Pressure_inc;
    
    su2double Gas_Constant = config->GetGas_ConstantND();
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
    bool tkeNeeded = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    su2double DampingFactor = config->GetDamp_Engine_Exhaust();
    su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
    
    su2double *Normal = new su2double[nDim];
    
    /*--- Retrieve the specified exhaust pressure in the engine (non-dimensional). ---*/
    
    Target_Exhaust_Pressure = config->GetExhaust_Pressure_Target(Marker_Tag) / config->GetPressure_Ref();
    
    /*--- Retrieve the old exhaust pressure in the engine exhaust (this has been computed in a preprocessing). ---*/
    
    Exhaust_Pressure_old = config->GetExhaust_Pressure(Marker_Tag);
    
    /*--- Compute the Pressure increment ---*/
    
    Exhaust_Pressure_inc = (1.0 - (Exhaust_Pressure_old/Target_Exhaust_Pressure)) * Baseline_Press;
    
    /*--- Estimate the new exhaust pressure ---*/
    
    Exhaust_Pressure = (1.0 - DampingFactor) * Exhaust_Pressure_old + DampingFactor * (Exhaust_Pressure_old + Exhaust_Pressure_inc);
    
    /*--- The temperature is given (no iteration is required) ---*/
    
    Exhaust_Temperature  = config->GetExhaust_Temperature_Target(Marker_Tag);
    Exhaust_Temperature /= config->GetTemperature_Ref();
    
    /*--- The pressure is given (no iteration is required) ---*/
    
    Exhaust_Pressure  = config->GetExhaust_Pressure_Target(Marker_Tag);
    Exhaust_Pressure /= config->GetPressure_Ref();
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        /*--- Allocate the value at the exhaust ---*/
        
        V_exhaust = GetCharacPrimVar(val_marker, iVertex);
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        
        /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
        
        if (geometry->node[iPoint]->GetDomain()) {
            
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
            
            V_domain = nodes->GetPrimitive(iPoint);
            
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
                V_exhaust[nDim+4] = sqrt(SoundSpeed2);
                
            }
            /*--- The flow goes in the wrong direction ---*/
            
            else {
                
                V_exhaust[0] = V_domain[0];
                for (iDim = 0; iDim < nDim; iDim++)
                    V_exhaust[iDim+1] = V_domain[iDim+1];
                V_exhaust[nDim+1] = V_domain[nDim+1];
                V_exhaust[nDim+2] = V_domain[nDim+2];
                V_exhaust[nDim+3] = V_domain[nDim+3];
                V_exhaust[nDim+4] = V_domain[nDim+4];
                
            }
            
            /*--- Set various quantities in the solver class ---*/
            
            conv_numerics->SetNormal(Normal);
            conv_numerics->SetPrimitive(V_domain, V_exhaust);
            
            /*--- Set grid movement ---*/
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration ---*/
            
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
            //
            //      if (viscous) {
            //
            //        /*--- Set laminar and eddy viscosity at the infinity ---*/
            //
            //        V_exhaust[nDim+5] = nodes->GetLaminarViscosity(iPoint);
            //        V_exhaust[nDim+6] = nodes->GetEddyViscosity(iPoint);
            //
            //        /*--- Set the normal vector and the coordinates ---*/
            //
            //        visc_numerics->SetNormal(Normal);
            //        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
            //
            //        /*--- Primitive variables, and gradient ---*/
            //
            //        visc_numerics->SetPrimitive(V_domain, V_exhaust);
            //        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
            //
            //        /*--- Turbulent kinetic energy ---*/
            //
            //        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
            //          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
            //                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
            //
            //        /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
            //
            //        visc_numerics->SetTauWall(-1.0, -1.0);
            //
            //        /*--- Compute and update residual ---*/
            //
            //        auto residual = visc_numerics->ComputeResidual(config)
            //        LinSysRes.SubtractBlock(iPoint, residual);
            //
            //        /*--- Jacobian contribution for implicit integration ---*/
            //
            //        if (implicit)
            //          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
            //
            //      }
            
        }
    }
    
    delete [] Normal;
    
}

void CEulerSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config) {
    
    unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
    unsigned short iDim, iVar, jVar, iMarker, nDonorVertex;
    
    bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    bool viscous       = config->GetViscous();
    
    su2double Normal[MAXNDIM] = {0.0};
    su2double PrimVar_i[MAXNVAR] = {0.0};
    su2double PrimVar_j[MAXNVAR] = {0.0};
    su2double Secondary_j[MAXNVAR] = {0.0};
    
    su2double weight;
    su2double P_static, rho_static;
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
        if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
            
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                
                if (geometry->node[iPoint]->GetDomain()) {
                    
                    nDonorVertex = GetnSlidingStates(iMarker, iVertex);
                    
                    /*--- Initialize Residual, this will serve to accumulate the average ---*/
                    
                    for (iVar = 0; iVar < nVar; iVar++) {
                        Residual[iVar] = 0.0;
                        for (jVar = 0; jVar < nVar; jVar++)
                            Jacobian_i[iVar][jVar] = 0.0;
                    }
                    
                    /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/
                    
                    for (jVertex = 0; jVertex < nDonorVertex; jVertex++) {
                        
                        Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
                        
                        for (iVar = 0; iVar < nPrimVar; iVar++) {
                            PrimVar_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
                            PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);
                        }
                        
                        /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/
                        
                        weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);
                        
                        /*--- Set primitive variables ---*/
                        
                        conv_numerics->SetPrimitive( PrimVar_i, PrimVar_j );
                        
                        if( !( config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS ) ) {
                            auto Secondary_i = nodes->GetSecondary(iPoint);
                            
                            P_static   = PrimVar_j[nDim+1];
                            rho_static = PrimVar_j[nDim+2];
                            GetFluidModel()->SetTDState_Prho(P_static, rho_static);
                            
                            Secondary_j[0] = GetFluidModel()->GetdPdrho_e();
                            Secondary_j[1] = GetFluidModel()->GetdPde_rho();
                            
                            conv_numerics->SetSecondary(Secondary_i, Secondary_j);
                        }
                        
                        /*--- Set the normal vector ---*/
                        
                        geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
                        for (iDim = 0; iDim < nDim; iDim++)
                            Normal[iDim] = -Normal[iDim];
                        
                        conv_numerics->SetNormal(Normal);
                        
                        if (dynamic_grid)
                            conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
                        
                        /*--- Compute the convective residual using an upwind scheme ---*/
                        
                        auto residual = conv_numerics->ComputeResidual(config);
                        
                        /*--- Accumulate the residuals to compute the average ---*/
                        
                        for (iVar = 0; iVar < nVar; iVar++) {
                            Residual[iVar] += weight*residual.residual[iVar];
                            for (jVar = 0; jVar < nVar; jVar++)
                                Jacobian_i[iVar][jVar] += weight*residual.jacobian_i[iVar][jVar];
                        }
                    }
                    
                    /*--- Add Residuals and Jacobians ---*/
                    
                    LinSysRes.AddBlock(iPoint, Residual);
                    
                    if (implicit)
                        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
                    
                    if (viscous) {
                        
                        /*--- Initialize Residual, this will serve to accumulate the average ---*/
                        
                        for (iVar = 0; iVar < nVar; iVar++) {
                            Residual[iVar] = 0.0;
                            for (jVar = 0; jVar < nVar; jVar++)
                                Jacobian_i[iVar][jVar] = 0.0;
                        }
                        
                        /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/
                        
                        for (jVertex = 0; jVertex < nDonorVertex; jVertex++) {
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
                                visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                                                    solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
                            
                            /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
                            
                            visc_numerics->SetTauWall(-1.0, -1.0);
                            
                            /*--- Compute and update residual ---*/
                            
                            auto residual = visc_numerics->ComputeResidual(config);
                            
                            /*--- Accumulate the residuals to compute the average ---*/
                            
                            for (iVar = 0; iVar < nVar; iVar++) {
                                Residual[iVar] += weight*residual.residual[iVar];
                                for (jVar = 0; jVar < nVar; jVar++)
                                    Jacobian_i[iVar][jVar] += weight*residual.jacobian_i[iVar][jVar];
                            }
                        }
                        
                        LinSysRes.SubtractBlock(iPoint, Residual);
                        
                        /*--- Jacobian contribution for implicit integration ---*/
                        
                        if (implicit)
                            Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);
                        
                    }
                }
            }
        }
    }
    
}

void CEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {
    
    unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
    unsigned short iDim, iVar;
    
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    
    su2double *Normal = new su2double[nDim];
    su2double *PrimVar_i = new su2double[nPrimVar];
    su2double *PrimVar_j = new su2double[nPrimVar];
    
    /*--- Do the send process, by the moment we are sending each
     node individually, this must be changed ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
        GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
        
        if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
                PrimVar_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
                PrimVar_j[iVar] = GetDonorPrimVar(val_marker, iVertex, iVar);
            }
            
            /*--- Set Conservative Variables ---*/
            
            numerics->SetPrimitive(PrimVar_i, PrimVar_j);
            
            /*--- Set Normal ---*/
            
            geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, residual);
            if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
        }
        
    }
    
    /*--- Free locally allocated memory ---*/
    
    delete [] Normal;
    delete [] PrimVar_i;
    delete [] PrimVar_j;
    
}

void CEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config, unsigned short val_marker) {
    
    unsigned long iVertex, iPoint, GlobalIndex_iPoint, GlobalIndex_jPoint;
    unsigned short iDim, iVar;
    
    bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    
    su2double *Normal = new su2double[nDim];
    su2double *PrimVar_i = new su2double[nPrimVar];
    su2double *PrimVar_j = new su2double[nPrimVar];
    
    /*--- Do the send process, by the moment we are sending each
     node individually, this must be changed ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        GlobalIndex_iPoint = geometry->node[iPoint]->GetGlobalIndex();
        GlobalIndex_jPoint = GetDonorGlobalIndex(val_marker, iVertex);
        
        if ((geometry->node[iPoint]->GetDomain()) && (GlobalIndex_iPoint != GlobalIndex_jPoint)) {
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++) {
                PrimVar_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
                PrimVar_j[iVar] = GetDonorPrimVar(val_marker, iVertex, iVar);
            }
            
            /*--- Set Conservative Variables ---*/
            
            numerics->SetPrimitive(PrimVar_i, PrimVar_j);
            
            /*--- Set Normal ---*/
            
            geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            
            auto residual = numerics->ComputeResidual(config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.AddBlock(iPoint, residual);
            
            if (implicit)
                Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
        }
        
    }
    
    /*--- Free locally allocated memory ---*/
    
    delete [] Normal;
    delete [] PrimVar_i;
    delete [] PrimVar_j;
    
}

void CEulerSolver::BC_ActDisk_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                    CConfig *config, unsigned short val_marker) {
    
    BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker, true);
    
}

void CEulerSolver::BC_ActDisk_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                     CConfig *config, unsigned short val_marker) {
    
    BC_ActDisk(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker, false);
    
}

void CEulerSolver::BC_ActDisk(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker, bool val_inlet_surface) {
    
    unsigned short iDim;
    unsigned long iVertex, iPoint, GlobalIndex_donor, GlobalIndex;
    su2double Pressure, Velocity[3], Target_Press_Jump, Target_Temp_Jump,
    Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Vn_Inlet, Mach_Outlet,
    Area, UnitNormal[3], *V_outlet, *V_domain, *V_inlet, P_Total, T_Total, H_Total, Temperature,
    Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag, alpha, aa, bb, cc, dd;
    su2double Factor, P_static, T_static, SoS_outlet, Rho_outlet, Rho_inlet;
    su2double Vel_normal_inlet[3], Vel_tangent_inlet[3], Vel_inlet[3];
    su2double Vel_normal_outlet[3], Vel_tangent_outlet[3], Vel_outlet[3];
    su2double Vel_normal_inlet_, Vel_tangent_inlet_, Vel_inlet_;
    su2double Vel_normal_outlet_, Vel_outlet_;
    
    su2double Pressure_out, Density_out, SoundSpeed_out, Velocity2_out,
    Mach_out, Pressure_in, Density_in, SoundSpeed_in, Velocity2_in,
    Mach_in, PressureAdj, TemperatureAdj;
    
    bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    su2double Gas_Constant  = config->GetGas_ConstantND();
    bool tkeNeeded          = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);
    bool ratio              = (config->GetActDisk_Jump() == RATIO);
    su2double SecondaryFlow = config->GetSecondaryFlow_ActDisk();
    
    su2double *Normal = new su2double[nDim];
    su2double *Flow_Dir = new su2double[nDim];
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    
    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
        
        iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
        GlobalIndex = geometry->node[iPoint]->GetGlobalIndex();
        GlobalIndex_donor = GetDonorGlobalIndex(val_marker, iVertex);
        
        /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
        
        if ((geometry->node[iPoint]->GetDomain()) &&
            (GlobalIndex != GlobalIndex_donor)) {
            
            /*--- Normal vector for this vertex (negative for outward convention) ---*/
            
            geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            conv_numerics->SetNormal(Normal);
            
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);
            
            for (iDim = 0; iDim < nDim; iDim++)
                UnitNormal[iDim] = Normal[iDim]/Area;
            
            /*--- Current solution at this boundary node and jumps values ---*/
            
            V_domain = nodes->GetPrimitive(iPoint);
            Target_Press_Jump = GetActDisk_DeltaP(val_marker, iVertex);
            Target_Temp_Jump = GetActDisk_DeltaT(val_marker, iVertex);
            
            if (val_inlet_surface) {
                V_inlet  = nodes->GetPrimitive(iPoint);
                V_outlet = GetDonorPrimVar(val_marker, iVertex);
                
                Pressure_out    = V_outlet[nDim+1];
                Density_out     = V_outlet[nDim+2];
                SoundSpeed_out  = sqrt(Gamma*Pressure_out/Density_out);
                
                Pressure_in    = V_inlet[nDim+1];
                Density_in     = V_inlet[nDim+2];
                SoundSpeed_in  = sqrt(Gamma*Pressure_in/Density_in);
                
                Velocity2_out = 0.0; Velocity2_in = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Velocity2_out += V_outlet[iDim+1]*V_outlet[iDim+1];
                    Velocity2_in  += V_inlet[iDim+1]*V_inlet[iDim+1];
                }
                
                PressureAdj = 1.0; TemperatureAdj = 1.0;
                if ((Velocity2_out > 0.0) && (Velocity2_in > 0.0)) {
                    
                    Mach_out = sqrt(Velocity2_out)/SoundSpeed_out;
                    Mach_in  = sqrt(Velocity2_in)/SoundSpeed_in;
                    
                    PressureAdj    = pow( 1.0 + Mach_out * Mach_out * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0)) /
                    pow( 1.0 + Mach_in * Mach_in * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
                    TemperatureAdj = (1.0 + Mach_out * Mach_out * 0.5 * (Gamma - 1.0)) /
                    (1.0 + Mach_in * Mach_in * 0.5 * (Gamma - 1.0));
                    
                }
                
                if (ratio) {
                    P_static = V_outlet[nDim+1] / (Target_Press_Jump/PressureAdj);
                    T_static = V_outlet[0] / (Target_Temp_Jump/TemperatureAdj);
                }
                else { P_static = V_outlet[nDim+1] - Target_Press_Jump; T_static = V_outlet[0] - Target_Temp_Jump; }
            }
            else {
                V_outlet = nodes->GetPrimitive(iPoint);
                V_inlet  = GetDonorPrimVar(val_marker, iVertex);
                
                Pressure_out    = V_outlet[nDim+1];
                Density_out     = V_outlet[nDim+2];
                SoundSpeed_out  = sqrt(Gamma*Pressure_out/Density_out);
                
                Pressure_in    = V_inlet[nDim+1];
                Density_in     = V_inlet[nDim+2];
                SoundSpeed_in  = sqrt(Gamma*Pressure_in/Density_in);
                
                Velocity2_out = 0.0; Velocity2_in = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Velocity2_out += V_outlet[iDim+1]*V_outlet[iDim+1];
                    Velocity2_in  += V_inlet[iDim+1]*V_inlet[iDim+1];
                }
                
                PressureAdj = 1.0; TemperatureAdj = 1.0;
                if ((Velocity2_out > 0.0) && (Velocity2_in > 0.0)) {
                    
                    Mach_out = sqrt(Velocity2_out)/SoundSpeed_out;
                    Mach_in  = sqrt(Velocity2_in)/SoundSpeed_in;
                    
                    PressureAdj    = pow( 1.0 + Mach_out * Mach_out * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0)) /
                    pow( 1.0 + Mach_in * Mach_in * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
                    TemperatureAdj = (1.0 + Mach_out * Mach_out * 0.5 * (Gamma - 1.0)) /
                    (1.0 + Mach_in * Mach_in * 0.5 * (Gamma - 1.0));
                }
                
                if (ratio) {
                    P_static = V_inlet[nDim+1] * (Target_Press_Jump/PressureAdj);
                    T_static = V_inlet[0] * (Target_Temp_Jump/TemperatureAdj);
                }
                else       { P_static = V_inlet[nDim+1] + Target_Press_Jump; T_static = V_inlet[0] + Target_Temp_Jump; }
            }
            
            /*--- Subsonic inlet ---*/
            
            if (val_inlet_surface) {
                
                /*--- Build the fictitious intlet state based on characteristics.
                 Retrieve the specified back pressure for this inlet ---*/
                
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
                
                Pressure   = P_static;
                Density    = pow(Pressure/Entropy,1.0/Gamma);
                SoundSpeed = sqrt(Gamma*Pressure/Density);
                Vn_Inlet    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
                
                Velocity2  = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Velocity[iDim] = Velocity[iDim] + (Vn_Inlet-Vn)*UnitNormal[iDim];
                    Velocity2 += Velocity[iDim]*Velocity[iDim];
                }
                Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
                if (tkeNeeded) Energy += GetTke_Inf();
                
                /*--- Conservative variables, using the derived quantities ---*/
                
                V_inlet[0] = Pressure / ( Gas_Constant * Density);
                for (iDim = 0; iDim < nDim; iDim++)
                    V_inlet[iDim+1] = Velocity[iDim];
                V_inlet[nDim+1] = Pressure;
                V_inlet[nDim+2] = Density;
                V_inlet[nDim+3] = Energy + Pressure/Density;
                V_inlet[nDim+4] = SoundSpeed;
                conv_numerics->SetPrimitive(V_domain, V_inlet);
                
            }
            
            /*--- Subsonic outlet ---*/
            
            else {
                
                GetFluidModel()->SetTDState_PT(P_static, T_static);
                SoS_outlet = GetFluidModel()->GetSoundSpeed();
                Rho_outlet = GetFluidModel()->GetDensity();
                
                /*--- We use the velocity and the density from the flow inlet
                 to evaluate flow direction and mass flow ---*/
                
                Rho_inlet = V_inlet[nDim+2];
                for (iDim = 0; iDim < nDim; iDim++)
                    Vel_inlet[iDim] = V_inlet[iDim+1];
                
                Vel_normal_inlet_ = 0.0; Vel_inlet_ = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Vel_normal_inlet[iDim] = -Vel_inlet[iDim]*UnitNormal[iDim];
                    Vel_normal_inlet_ += Vel_normal_inlet[iDim]*Vel_normal_inlet[iDim];
                    Vel_inlet_+= Vel_inlet[iDim]*Vel_inlet[iDim];
                }
                Vel_inlet_ = sqrt(Vel_inlet_);
                Vel_normal_inlet_ = sqrt(Vel_normal_inlet_);
                
                Vel_tangent_inlet_ = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Vel_tangent_inlet[iDim] = Vel_inlet[iDim] - Vel_normal_inlet[iDim];
                    Vel_tangent_inlet_ += Vel_tangent_inlet[iDim]*Vel_tangent_inlet[iDim];
                }
                Vel_tangent_inlet_ = sqrt(Vel_tangent_inlet_);
                
                /*--- Mass flow conservation (normal direction) and
                 no jump in the tangential velocity ---*/
                
                Vel_normal_outlet_ = (1.0-SecondaryFlow/100.0)*(Rho_inlet*Vel_normal_inlet_)/Rho_outlet;
                
                Vel_outlet_ = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Vel_normal_outlet[iDim] = -Vel_normal_outlet_*UnitNormal[iDim];
                    Vel_tangent_outlet[iDim] = Vel_tangent_inlet[iDim];
                    Vel_outlet[iDim] = Vel_normal_outlet[iDim] + Vel_tangent_outlet[iDim];
                    Vel_outlet_ += Vel_outlet[iDim]*Vel_outlet[iDim];
                }
                Vel_outlet_ = sqrt(Vel_outlet_);
                
                Mach_Outlet = min(Vel_outlet_/SoS_outlet, 1.0);
                
                /*--- Reevaluate the Total Pressure and Total Temperature using the
                 Fan Face Mach number and the static values from the jum condition ---*/
                
                Factor = 1.0 + 0.5*Mach_Outlet*Mach_Outlet*Gamma_Minus_One;
                P_Total = P_static * pow(Factor, Gamma/Gamma_Minus_One);
                T_Total = T_static * Factor;
                
                /*--- Flow direction using the velocity direction at the outlet  ---*/
                
                if (Vel_outlet_ != 0.0) {
                    for (iDim = 0; iDim < nDim; iDim++) Flow_Dir[iDim] = Vel_outlet[iDim]/Vel_outlet_;
                }
                else {
                    for (iDim = 0; iDim < nDim; iDim++) Flow_Dir[iDim] = 0.0;
                }
                
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
                cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Total2/Gamma_Minus_One;
                
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
                
                Mach2 = min(1.0, Velocity2/SoundSpeed2);
                Velocity2   = Mach2*SoundSpeed2;
                Vel_Mag     = sqrt(Velocity2);
                SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
                
                /*--- Compute new velocity vector at the exit ---*/
                
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
                
                V_outlet[0] = Temperature;
                for (iDim = 0; iDim < nDim; iDim++)
                    V_outlet[iDim+1] = Velocity[iDim];
                V_outlet[nDim+1] = Pressure;
                V_outlet[nDim+2] = Density;
                V_outlet[nDim+3] = Energy + Pressure/Density;
                V_outlet[nDim+4] = sqrt(SoundSpeed2);
                conv_numerics->SetPrimitive(V_domain, V_outlet);
                
            }
            
            /*--- Grid Movement ---*/
            
            if (dynamic_grid)
                conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
            
            /*--- Compute the residual using an upwind scheme ---*/
            
            auto residual = conv_numerics->ComputeResidual(config);
            
            /*--- Update residual value ---*/
            
            LinSysRes.AddBlock(iPoint, residual);
            
            /*--- Jacobian contribution for implicit integration ---*/
            
            if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
            
            //      /*--- Viscous contribution, commented out because serious convergence problems ---*/
            //
            //      if (viscous) {
            //
            //        /*--- Set laminar and eddy viscosity at the infinity ---*/
            //
            //        if (val_inlet_surface) {
            //          V_inlet[nDim+5] = nodes->GetLaminarViscosity(iPoint);
            //          V_inlet[nDim+6] = nodes->GetEddyViscosity(iPoint);
            //        }
            //        else {
            //          V_outlet[nDim+5] = nodes->GetLaminarViscosity(iPoint);
            //          V_outlet[nDim+6] = nodes->GetEddyViscosity(iPoint);
            //        }
            //
            //        /*--- Set the normal vector and the coordinates ---*/
            //
            //        visc_numerics->SetNormal(Normal);
            //        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint_Normal]->GetCoord());
            //
            //        /*--- Primitive variables, and gradient ---*/
            //
            //        if (val_inlet_surface) visc_numerics->SetPrimitive(V_domain, V_inlet);
            //        else visc_numerics->SetPrimitive(V_domain, V_outlet);
            //
            //        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nodes->GetGradient_Primitive(iPoint));
            //
            //        /*--- Turbulent kinetic energy ---*/
            //
            //        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
            //          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
            //                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));
            //
            //        /*--- Set the wall shear stress values (wall functions) to -1 (no evaluation using wall functions) ---*/
            //
            //        visc_numerics->SetTauWall(-1.0, -1.0);
            //
            //        /*--- Compute and update residual ---*/
            //
            //        auto residual = visc_numerics->ComputeResidual(config);
            //        LinSysRes.SubtractBlock(iPoint, residual);
            //
            //        /*--- Jacobian contribution for implicit integration ---*/
            //
            //        if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);
            //
            //      }
            
        }
        
    }
    
    /*--- Free locally allocated memory ---*/
    
    delete [] Normal;
    delete [] Flow_Dir;
    
}

void CEulerSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
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

void CEulerSolver::BC_Custom(CGeometry      *geometry,
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

void CEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
    
    /*--- Local variables ---*/
    
    unsigned short iVar, iMarker, iDim, iNeigh;
    unsigned long iPoint, jPoint, iEdge, iVertex;
    
    const su2double *U_time_nM1 = nullptr, *U_time_n = nullptr, *U_time_nP1 = nullptr;
    su2double Volume_nM1, Volume_nP1, TimeStep;
    const su2double *Normal = nullptr, *GridVel_i = nullptr, *GridVel_j = nullptr;
    su2double Residual_GCL;
    
    const bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
    const bool first_order = (config->GetTime_Marching() == DT_STEPPING_1ST);
    const bool second_order = (config->GetTime_Marching() == DT_STEPPING_2ND);
    
    /*--- Store the physical time step ---*/
    
    TimeStep = config->GetDelta_UnstTimeND();
    
    /*--- Compute the dual time-stepping source term for static meshes ---*/
    
    if (!dynamic_grid) {
        
        /*--- Loop over all nodes (excluding halos) ---*/
        
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
             we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
             previous solutions that are stored in memory. ---*/
            
            U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
            U_time_n   = nodes->GetSolution_time_n(iPoint);
            U_time_nP1 = nodes->GetSolution(iPoint);
            
            /*--- CV volume at time n+1. As we are on a static mesh, the volume
             of the CV will remained fixed for all time steps. ---*/
            
            Volume_nP1 = geometry->node[iPoint]->GetVolume();
            
            /*--- Compute the dual time-stepping source term based on the chosen
             time discretization scheme (1st- or 2nd-order).---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
                if (first_order)
                    LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
                if (second_order)
                    LinSysRes(iPoint,iVar) += ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                                               +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
            }
            
            /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
            if (implicit) {
                if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
                if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
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
        
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
            
            GridVel_i = geometry->node[iPoint]->GetGridVel();
            U_time_n = nodes->GetSolution_time_n(iPoint);
            
            for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnNeighbor(); iNeigh++) {
                
                iEdge = geometry->node[iPoint]->GetEdge(iNeigh);
                Normal = geometry->edges->GetNormal(iEdge);
                
                jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
                GridVel_j = geometry->node[jPoint]->GetGridVel();
                
                /*--- Determine whether to consider the normal outward or inward. ---*/
                su2double dir = (geometry->edges->GetNode(iEdge,0) == iPoint)? 0.5 : -0.5;
                
                Residual_GCL = 0.0;
                for (iDim = 0; iDim < nDim; iDim++)
                    Residual_GCL += dir*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
                
                for (iVar = 0; iVar < nVar; iVar++)
                    LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
            }
        }
        
        /*--- Loop over the boundary edges ---*/
        
        for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
            if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)  &&
                (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
                
                SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
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
                    
                    U_time_n = nodes->GetSolution_time_n(iPoint);
                    for (iVar = 0; iVar < nVar; iVar++)
                        LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
                }
            }
        }
        
        /*--- Loop over all nodes (excluding halos) to compute the remainder
         of the dual time-stepping source term. ---*/
        
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
            
            /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
             we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
             previous solutions that are stored in memory. ---*/
            
            U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
            U_time_n   = nodes->GetSolution_time_n(iPoint);
            U_time_nP1 = nodes->GetSolution(iPoint);
            
            /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
             grids, the volumes will change. On rigidly transforming grids, the
             volumes will remain constant. ---*/
            
            Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
            Volume_nP1 = geometry->node[iPoint]->GetVolume();
            
            /*--- Compute the dual time-stepping source residual. Due to the
             introduction of the GCL term above, the remainder of the source residual
             due to the time discretization has a new form.---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
                if (first_order)
                    LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
                if (second_order)
                    LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
                    + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
            }
            
            /*--- Compute the Jacobian contribution due to the dual time source term. ---*/
            if (implicit) {
                if (first_order) Jacobian.AddVal2Diag(iPoint, Volume_nP1/TimeStep);
                if (second_order) Jacobian.AddVal2Diag(iPoint, (Volume_nP1*3.0)/(2.0*TimeStep));
            }
        }
    }
    
}

void CEulerSolver::GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {

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
    if (config->GetKind_Comp_Outlet(Outlet_TagBound) == MASS_FLOW_OUTLET)
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

            Density = V_outlet[nDim+2];

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

void CEulerSolver::ComputeVerificationError(CGeometry *geometry,
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
                
                cout.precision(5);
                cout.setf(ios::scientific, ios::floatfield);
                
                if (!config->GetDiscrete_Adjoint()) {
                    
                    cout << endl   << "------------------------ Global Error Analysis --------------------------" << endl;
                    
                    cout << setw(20) << "RMS Error  [Rho]: " << setw(12) << VerificationSolution->GetError_RMS(0) << "     | ";
                    cout << setw(20) << "Max Error  [Rho]: " << setw(12) << VerificationSolution->GetError_Max(0);
                    cout << endl;
                    
                    cout << setw(20) << "RMS Error [RhoU]: " << setw(12) << VerificationSolution->GetError_RMS(1) << "     | ";
                    cout << setw(20) << "Max Error [RhoU]: " << setw(12) << VerificationSolution->GetError_Max(1);
                    cout << endl;
                    
                    cout << setw(20) << "RMS Error [RhoV]: " << setw(12) << VerificationSolution->GetError_RMS(2) << "     | ";
                    cout << setw(20) << "Max Error [RhoV]: " << setw(12) << VerificationSolution->GetError_Max(2);
                    cout << endl;
                    
                    if (nDim == 3) {
                        cout << setw(20) << "RMS Error [RhoW]: " << setw(12) << VerificationSolution->GetError_RMS(3) << "     | ";
                        cout << setw(20) << "Max Error [RhoW]: " << setw(12) << VerificationSolution->GetError_Max(3);
                        cout << endl;
                    }
                    
                    cout << setw(20) << "RMS Error [RhoE]: " << setw(12) << VerificationSolution->GetError_RMS(nDim+1) << "     | ";
                    cout << setw(20) << "Max Error [RhoE]: " << setw(12) << VerificationSolution->GetError_Max(nDim+1);
                    cout << endl;
                    
                    cout << "-------------------------------------------------------------------------" << endl << endl;
                    cout.unsetf(ios_base::floatfield);
                    
                }
            }
        }
    }
    
}

void CEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {
    
    /*--- Restart the solution from file information ---*/
    
    unsigned short iDim, iVar, iMesh, iMeshFine;
    unsigned long iPoint, index, iChildren, Point_Fine;
    unsigned short turb_model = config->GetKind_Turb_Model();
    su2double Area_Children, Area_Parent;
    const su2double* Solution_Fine = nullptr;
    const passivedouble* Coord = nullptr;
    bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                      (config->GetTime_Marching() == DT_STEPPING_2ND));
    bool static_fsi = ((config->GetTime_Marching() == STEADY) && config->GetFSI_Simulation());
    bool steady_restart = config->GetSteadyRestart();
    bool turbulent = (config->GetKind_Turb_Model() != NONE);
    
    string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);
    
    /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
    SU2_OMP_MASTER
    {
        
        /*--- Skip coordinates ---*/
        
        unsigned short skipVars = geometry[MESH_0]->GetnDim();
        
        /*--- Store the number of variables for the turbulence model
         (that could appear in the restart file before the grid velocities). ---*/
        unsigned short turbVars = 0;
        if (turbulent){
            if ((turb_model == SST) || (turb_model == SST_SUST)) turbVars = 2;
            else turbVars = 1;
        }
        
        /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/
        
        if (config->GetRead_Binary_Restart()) {
            Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
        } else {
            Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
        }
        
        /*--- Load data from the restart into correct containers. ---*/
        
        unsigned long counter = 0, iPoint_Global = 0;
        for (; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {
            
            /*--- Retrieve local index. If this node from the restart file lives
             on the current processor, we will load and instantiate the vars. ---*/
            
            auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);
            
            if (iPoint_Local > -1) {
                
                /*--- We need to store this point's data, so jump to the correct
                 offset in the buffer of data from the restart file and load it. ---*/
                
                index = counter*Restart_Vars[1] + skipVars;
                for (iVar = 0; iVar < nVar; ++iVar)
                    nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index+iVar]);
                
                /*--- For dynamic meshes, read in and store the
                 grid coordinates and grid velocities for each node. ---*/
                
                if (dynamic_grid && val_update_geo) {
                    
                    /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
                    /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
                    /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/
                    
                    /*--- Rewind the index to retrieve the Coords. ---*/
                    index = counter*Restart_Vars[1];
                    Coord = &Restart_Data[index];
                    
                    su2double GridVel[MAXNDIM] = {0.0};
                    if (!steady_restart) {
                        /*--- Move the index forward to get the grid velocities. ---*/
                        index += skipVars + nVar + turbVars;
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
                    Coord = &Restart_Data[index];
                    
                    for (iDim = 0; iDim < nDim; iDim++) {
                        geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
                    }
                }
                
                /*--- Increment the overall counter for how many points have been loaded. ---*/
                counter++;
            }
            
        }
        
        /*--- Detect a wrong solution file ---*/
        
        if (counter != nPointDomain) {
            SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                           string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
        }
        
        /*--- Update the geometry for flows on deforming meshes ---*/
        
        if ((dynamic_grid || static_fsi) && val_update_geo) {
            
            /*--- Communicate the new coordinates and grid velocities at the halos ---*/
            
            geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
            geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);
            
            if (dynamic_grid) {
                geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, GRID_VELOCITY);
                geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, GRID_VELOCITY);
            }
            
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
                if (dynamic_grid) {
                    geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine], config);
                }
                geometry[iMesh]->SetMaxLength(config);
            }
        }
        
        /*--- Communicate the loaded solution on the fine grid before we transfer
         it down to the coarse levels. We also call the preprocessing routine
         on the fine level in order to have all necessary quantities updated,
         especially if this is a turbulent simulation (eddy viscosity). ---*/
        
        solver[MESH_0][FLOW_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
        solver[MESH_0][FLOW_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);
        
    } // end SU2_OMP_MASTER, preprocessing is thread-safe.
    SU2_OMP_BARRIER
    
    /*--- For turbulent simulations the flow preprocessing is done by the turbulence solver
     *    after it loads its variables (they are needed to compute flow primitives). ---*/
    if (!turbulent) {
        solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
    
    /*--- Interpolate the solution down to the coarse multigrid levels ---*/
    
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
        SU2_OMP_FOR_STAT(omp_chunk_size)
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
            Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
            su2double Solution_Coarse[MAXNVAR] = {0.0};
            for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
                Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
                Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
                Solution_Fine = solver[iMesh-1][FLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
                for (iVar = 0; iVar < nVar; iVar++) {
                    Solution_Coarse[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
                }
            }
            solver[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint,Solution_Coarse);
        }
        
        SU2_OMP_MASTER
        {
            solver[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
            solver[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
        }
        SU2_OMP_BARRIER
        
        if (!turbulent) {
            solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        }
    }
    
    /*--- Go back to single threaded execution. ---*/
    SU2_OMP_MASTER
    {
        
        /*--- Update the old geometry (coordinates n and n-1) in dual time-stepping strategy. ---*/
        if (dual_time && config->GetGrid_Movement() && !config->GetDeform_Mesh() &&
            (config->GetKind_GridMovement() != RIGID_MOTION)) {
            Restart_OldGeometry(geometry[MESH_0], config);
        }
        
        /*--- Delete the class memory that is used to load the restart. ---*/
        
        delete [] Restart_Vars; Restart_Vars = nullptr;
        delete [] Restart_Data; Restart_Data = nullptr;
        
    } // end SU2_OMP_MASTER
    SU2_OMP_BARRIER
    
}

void CEulerSolver::SetFreeStream_Solution(CConfig *config) {
    
    unsigned long iPoint;
    unsigned short iDim;
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
        nodes->SetSolution(iPoint,0, Density_Inf);
        for (iDim = 0; iDim < nDim; iDim++) {
            nodes->SetSolution(iPoint,iDim+1, Density_Inf*Velocity_Inf[iDim]);
        }
        nodes->SetSolution(iPoint,nVar-1, Density_Inf*Energy_Inf);
    }
}

void CEulerSolver::SetFreeStream_TurboSolution(CConfig *config) {
    
    unsigned long iPoint;
    unsigned short iDim;
    unsigned short iZone  =  config->GetiZone();
    su2double *turboVelocity, *cartVelocity, *turboNormal;
    
    su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
    su2double Mach             = config->GetMach();
    su2double SoundSpeed;
    
    turboVelocity   = new su2double[nDim];
    cartVelocity    = new su2double[nDim];
    
    turboNormal     = config->GetFreeStreamTurboNormal();
    
    GetFluidModel()->SetTDState_Prho(Pressure_Inf, Density_Inf);
    SoundSpeed = GetFluidModel()->GetSoundSpeed();
    
    /*--- Compute the Free Stream velocity, using the Mach number ---*/
    turboVelocity[0] = cos(Alpha)*Mach*SoundSpeed;
    turboVelocity[1] = sin(Alpha)*Mach*SoundSpeed;
    
    
    if (nDim == 3) {
        turboVelocity[2] = 0.0;
    }
    
    ComputeBackVelocity(turboVelocity, turboNormal, cartVelocity, INFLOW, config->GetKind_TurboMachinery(iZone));
    
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
        nodes->SetSolution(iPoint,0, Density_Inf);
        for (iDim = 0; iDim < nDim; iDim++) {
            nodes->SetSolution(iPoint,iDim+1, Density_Inf*cartVelocity[iDim]);
        }
        nodes->SetSolution(iPoint,nVar-1, Density_Inf*Energy_Inf);
        
        nodes->SetPrimVar(iPoint, GetFluidModel());
        nodes->SetSecondaryVar(iPoint, GetFluidModel());
    }
    
    delete [] turboVelocity;
    delete [] cartVelocity;
    
    
}



void CEulerSolver::PreprocessAverage(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag) {
    
    unsigned long iVertex, iPoint;
    unsigned short iDim, iMarker, iMarkerTP, iSpan;
    su2double Pressure = 0.0, Density = 0.0, *Velocity = NULL, *TurboVelocity,
    Area, TotalArea, TotalAreaPressure, TotalAreaDensity, *TotalAreaVelocity, *UnitNormal, *TurboNormal;
    string Marker_Tag, Monitoring_Tag;
    unsigned short  iZone     = config->GetiZone();
    const su2double  *AverageTurboNormal;
    su2double VelSq;
    
    /*-- Variables declaration and allocation ---*/
    Velocity           = new su2double[nDim];
    UnitNormal         = new su2double[nDim];
    TurboNormal        = new su2double[nDim];
    TurboVelocity      = new su2double[nDim];
    TotalAreaVelocity  = new su2double[nDim];
    
    
    for (iSpan= 0; iSpan < nSpanWiseSections; iSpan++){
        
        for (iDim=0; iDim<nDim; iDim++) {
            TotalAreaVelocity[iDim] = 0.0;
        }
        
        TotalAreaPressure = 0.0;
        TotalAreaDensity  = 0.0;
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
            for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
                if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                    if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                        
                        /*--- Retrieve Old Solution ---*/
                        
                        /*--- Loop over the vertices to sum all the quantities pithc-wise ---*/
                        for (iVertex = 0; iVertex < geometry->GetnVertexSpan(iMarker,iSpan); iVertex++) {
                            iPoint = geometry->turbovertex[iMarker][iSpan][iVertex]->GetNode();
                            if (geometry->node[iPoint]->GetDomain()){
                                /*--- Compute the integral fluxes for the boundaries ---*/
                                
                                Pressure = nodes->GetPressure(iPoint);
                                Density = nodes->GetDensity(iPoint);
                                
                                /*--- Normal vector for this vertex (negate for outward convention) ---*/
                                geometry->turbovertex[iMarker][iSpan][iVertex]->GetNormal(UnitNormal);
                                geometry->turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(TurboNormal);
                                Area = geometry->turbovertex[iMarker][iSpan][iVertex]->GetArea();
                                
                                VelSq = 0.0;
                                for (iDim = 0; iDim < nDim; iDim++) {
                                    Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);
                                    VelSq += Velocity[iDim]*Velocity[iDim];
                                }
                                
                                ComputeTurboVelocity(Velocity, TurboNormal , TurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));
                                
                                /*--- Compute different integral quantities for the boundary of interest ---*/
                                
                                TotalAreaPressure += Area*Pressure;
                                TotalAreaDensity  += Area*Density;
                                for (iDim = 0; iDim < nDim; iDim++)
                                    TotalAreaVelocity[iDim] += Area*Velocity[iDim];
                            }
                        }
                    }
                }
            }
        }
        
        
#ifdef HAVE_MPI
        
        /*--- Add information using all the nodes ---*/
        
        su2double MyTotalAreaDensity = TotalAreaDensity;
        su2double MyTotalAreaPressure  = TotalAreaPressure;
        
        SU2_MPI::Allreduce(&MyTotalAreaDensity, &TotalAreaDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&MyTotalAreaPressure, &TotalAreaPressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        su2double* MyTotalAreaVelocity = new su2double[nDim];
        
        for (iDim = 0; iDim < nDim; iDim++) {
            MyTotalAreaVelocity[iDim] = TotalAreaVelocity[iDim];
        }
        
        SU2_MPI::Allreduce(MyTotalAreaVelocity, TotalAreaVelocity, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        delete [] MyTotalAreaVelocity;
        
#endif
        
        /*--- initialize spanwise average quantities ---*/
        
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
            for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
                if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                    if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                        
                        TotalArea           = geometry->GetSpanArea(iMarker,iSpan);
                        AverageTurboNormal  = geometry->GetAverageTurboNormal(iMarker,iSpan);
                        
                        /*--- Compute the averaged value for the boundary of interest for the span of interest ---*/
                        
                        AverageDensity[iMarker][iSpan]           = TotalAreaDensity / TotalArea;
                        AveragePressure[iMarker][iSpan]          = TotalAreaPressure / TotalArea;
                        for (iDim = 0; iDim < nDim; iDim++)
                            AverageVelocity[iMarker][iSpan][iDim]  = TotalAreaVelocity[iDim] / TotalArea;
                        
                        /* --- compute static averaged quantities ---*/
                        ComputeTurboVelocity(AverageVelocity[iMarker][iSpan], AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                        
                        OldAverageDensity[iMarker][iSpan]               = AverageDensity[iMarker][iSpan];
                        OldAveragePressure[iMarker][iSpan]              = AveragePressure[iMarker][iSpan];
                        for(iDim = 0; iDim < nDim;iDim++)
                            OldAverageTurboVelocity[iMarker][iSpan][iDim] = AverageTurboVelocity[iMarker][iSpan][iDim];
                        
                    }
                }
            }
        }
    }
    
    /*--- initialize 1D average quantities ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
            if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                    
                    AverageTurboNormal  = geometry->GetAverageTurboNormal(iMarker,nSpanWiseSections);
                    
                    /*--- Compute the averaged value for the boundary of interest for the span of interest ---*/
                    
                    AverageDensity[iMarker][nSpanWiseSections]          = AverageDensity[iMarker][nSpanWiseSections/2];
                    AveragePressure[iMarker][nSpanWiseSections]         = AveragePressure[iMarker][nSpanWiseSections/2];
                    for (iDim = 0; iDim < nDim; iDim++)
                        AverageVelocity[iMarker][nSpanWiseSections][iDim] = AverageVelocity[iMarker][nSpanWiseSections/2][iDim];
                    
                    /* --- compute static averaged quantities ---*/
                    ComputeTurboVelocity(AverageVelocity[iMarker][nSpanWiseSections], AverageTurboNormal , AverageTurboVelocity[iMarker][nSpanWiseSections], marker_flag, config->GetKind_TurboMachinery(iZone));
                    
                    OldAverageDensity[iMarker][nSpanWiseSections]               = AverageDensity[iMarker][nSpanWiseSections];
                    OldAveragePressure[iMarker][nSpanWiseSections]              = AveragePressure[iMarker][nSpanWiseSections];
                    for(iDim = 0; iDim < nDim;iDim++)
                        OldAverageTurboVelocity[iMarker][nSpanWiseSections][iDim] = AverageTurboVelocity[iMarker][nSpanWiseSections][iDim];
                    
                }
            }
        }
    }
    
    
    /*--- Free locally allocated memory ---*/
    delete [] Velocity;
    delete [] UnitNormal;
    delete [] TurboNormal;
    delete [] TurboVelocity;
    delete [] TotalAreaVelocity;
    
}


void CEulerSolver::TurboAverageProcess(CSolver **solver, CGeometry *geometry, CConfig *config, unsigned short marker_flag) {
    
    unsigned long iVertex, iPoint, nVert;
    unsigned short iDim, iVar, iMarker, iMarkerTP, iSpan, jSpan;
    unsigned short average_process = config->GetKind_AverageProcess();
    unsigned short performance_average_process = config->GetKind_PerformanceAverageProcess();
    su2double Pressure = 0.0, Density = 0.0, Enthalpy = 0.0,  *Velocity = NULL, *TurboVelocity,
    Area, TotalArea, Radius1, Radius2, Vt2, TotalAreaPressure, TotalAreaDensity, *TotalAreaVelocity, *UnitNormal, *TurboNormal,
    TotalMassPressure, TotalMassDensity, *TotalMassVelocity;
    string Marker_Tag, Monitoring_Tag;
    su2double val_init_pressure;
    unsigned short  iZone     = config->GetiZone();
    su2double TotalDensity, TotalPressure, *TotalVelocity, *TotalFluxes;
    const su2double *AverageTurboNormal;
    su2double TotalNu, TotalOmega, TotalKine, TotalMassNu, TotalMassOmega, TotalMassKine, TotalAreaNu, TotalAreaOmega, TotalAreaKine;
    su2double Nu, Kine, Omega;
    su2double MachTest, soundSpeed;
    bool turbulent = (config->GetKind_Turb_Model() != NONE);
    bool spalart_allmaras = (config->GetKind_Turb_Model() == SA);
    bool menter_sst       = ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST));
    
    /*-- Variables declaration and allocation ---*/
    Velocity            = new su2double[nDim];
    UnitNormal          = new su2double[nDim];
    TurboNormal         = new su2double[nDim];
    TurboVelocity       = new su2double[nDim];
    TotalVelocity       = new su2double[nDim];
    TotalAreaVelocity   = new su2double[nDim];
    TotalMassVelocity   = new su2double[nDim];
    TotalFluxes         = new su2double[nVar];
    
    su2double avgDensity, *avgVelocity, avgPressure, avgKine, avgOmega, avgNu, avgAreaDensity, *avgAreaVelocity, avgAreaPressure,
    avgAreaKine, avgAreaOmega, avgAreaNu, avgMassDensity, *avgMassVelocity, avgMassPressure, avgMassKine, avgMassOmega, avgMassNu,
    avgMixDensity, *avgMixVelocity, *avgMixTurboVelocity, avgMixPressure, avgMixKine, avgMixOmega, avgMixNu;
    
    avgVelocity         = new su2double[nDim];
    avgAreaVelocity     = new su2double[nDim];
    avgMassVelocity     = new su2double[nDim];
    avgMixVelocity      = new su2double[nDim];
    avgMixTurboVelocity = new su2double[nDim];
    
    
    for (iSpan= 0; iSpan < nSpanWiseSections + 1; iSpan++){
        
        /*--- Forces initialization for contenitors ---*/
        for (iVar=0;iVar<nVar;iVar++)
            TotalFluxes[iVar]= 0.0;
        for (iDim=0; iDim<nDim; iDim++) {
            TotalVelocity[iDim]     = 0.0;
            TotalAreaVelocity[iDim] = 0.0;
            TotalMassVelocity[iDim] = 0.0;
        }
        
        TotalDensity      = 0.0;
        TotalPressure     = 0.0;
        TotalAreaPressure = 0.0;
        TotalAreaDensity  = 0.0;
        TotalMassPressure = 0.0;
        TotalMassDensity  = 0.0;
        TotalNu           = 0.0;
        TotalOmega        = 0.0;
        TotalKine         = 0.0;
        TotalMassNu       = 0.0;
        TotalMassOmega    = 0.0;
        TotalMassKine     = 0.0;
        TotalAreaNu       = 0.0;
        TotalAreaOmega    = 0.0;
        TotalAreaKine     = 0.0;
        
        Nu    = 0.0;
        Omega = 0.0;
        Kine  = 0.0;
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
            for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
                if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                    if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                        
                        /*--- Retrieve Old Solution ---*/
                        
                        /*--- Loop over the vertices to sum all the quantities pithc-wise ---*/
                        if(iSpan < nSpanWiseSections){
                            for (iVertex = 0; iVertex < geometry->GetnVertexSpan(iMarker,iSpan); iVertex++) {
                                iPoint = geometry->turbovertex[iMarker][iSpan][iVertex]->GetNode();
                                
                                /*--- Compute the integral fluxes for the boundaries ---*/
                                Pressure = nodes->GetPressure(iPoint);
                                Density  = nodes->GetDensity(iPoint);
                                Enthalpy = nodes->GetEnthalpy(iPoint);
                                
                                /*--- Normal vector for this vertex (negate for outward convention) ---*/
                                geometry->turbovertex[iMarker][iSpan][iVertex]->GetNormal(UnitNormal);
                                geometry->turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(TurboNormal);
                                Area = geometry->turbovertex[iMarker][iSpan][iVertex]->GetArea();
                                su2double VelNormal = 0.0, VelSq = 0.0;
                                
                                for (iDim = 0; iDim < nDim; iDim++) {
                                    Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);
                                    VelNormal += UnitNormal[iDim]*Velocity[iDim];
                                    VelSq += Velocity[iDim]*Velocity[iDim];
                                }
                                
                                ComputeTurboVelocity(Velocity, TurboNormal , TurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));
                                
                                /*--- Compute different integral quantities for the boundary of interest ---*/
                                
                                TotalDensity          += Density;
                                TotalPressure         += Pressure;
                                for (iDim = 0; iDim < nDim; iDim++)
                                    TotalVelocity[iDim] += Velocity[iDim];
                                
                                TotalAreaPressure         += Area*Pressure;
                                TotalAreaDensity          += Area*Density;
                                for (iDim = 0; iDim < nDim; iDim++)
                                    TotalAreaVelocity[iDim] += Area*Velocity[iDim];
                                
                                TotalMassPressure         += Area*(Density*TurboVelocity[0] )*Pressure;
                                TotalMassDensity          += Area*(Density*TurboVelocity[0] )*Density;
                                for (iDim = 0; iDim < nDim; iDim++)
                                    TotalMassVelocity[iDim] += Area*(Density*TurboVelocity[0] )*Velocity[iDim];
                                
                                TotalFluxes[0]      += Area*(Density*TurboVelocity[0]);
                                TotalFluxes[1]      += Area*(Density*TurboVelocity[0]*TurboVelocity[0] + Pressure);
                                for (iDim = 2; iDim < nDim+1; iDim++)
                                    TotalFluxes[iDim] += Area*(Density*TurboVelocity[0]*TurboVelocity[iDim -1]);
                                TotalFluxes[nDim+1] += Area*(Density*TurboVelocity[0]*Enthalpy);
                                
                                
                                /*--- Compute turbulent integral quantities for the boundary of interest ---*/
                                
                                if(turbulent){
                                    if(menter_sst){
                                        Kine = solver[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
                                        Omega = solver[TURB_SOL]->GetNodes()->GetSolution(iPoint,1);
                                    }
                                    if(spalart_allmaras){
                                        Nu = solver[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
                                    }
                                    
                                    TotalKine   += Kine;
                                    TotalOmega  += Omega;
                                    TotalNu     += Nu;
                                    
                                    TotalAreaKine    += Area*Kine;
                                    TotalAreaOmega   += Area*Omega;
                                    TotalAreaNu      += Area*Nu;
                                    
                                    TotalMassKine    += Area*(Density*TurboVelocity[0] )*Kine;
                                    TotalMassOmega   += Area*(Density*TurboVelocity[0] )*Omega;
                                    TotalMassNu      += Area*(Density*TurboVelocity[0] )*Nu;
                                    
                                    
                                }
                            }
                        }
                        else{
                            for (jSpan= 0; jSpan < nSpanWiseSections; jSpan++){
                                for (iVertex = 0; iVertex < geometry->GetnVertexSpan(iMarker,jSpan); iVertex++) {
                                    iPoint = geometry->turbovertex[iMarker][jSpan][iVertex]->GetNode();
                                    
                                    /*--- Compute the integral fluxes for the boundaries ---*/
                                    Pressure = nodes->GetPressure(iPoint);
                                    Density  = nodes->GetDensity(iPoint);
                                    Enthalpy = nodes->GetEnthalpy(iPoint);
                                    
                                    /*--- Normal vector for this vertex (negate for outward convention) ---*/
                                    geometry->turbovertex[iMarker][jSpan][iVertex]->GetNormal(UnitNormal);
                                    geometry->turbovertex[iMarker][jSpan][iVertex]->GetTurboNormal(TurboNormal);
                                    Area = geometry->turbovertex[iMarker][jSpan][iVertex]->GetArea();
                                    su2double VelNormal = 0.0, VelSq = 0.0;
                                    
                                    for (iDim = 0; iDim < nDim; iDim++) {
                                        Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);
                                        VelNormal += UnitNormal[iDim]*Velocity[iDim];
                                        VelSq += Velocity[iDim]*Velocity[iDim];
                                    }
                                    
                                    ComputeTurboVelocity(Velocity, TurboNormal , TurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));
                                    
                                    /*--- Compute different integral quantities for the boundary of interest ---*/
                                    
                                    TotalDensity          += Density;
                                    TotalPressure         += Pressure;
                                    for (iDim = 0; iDim < nDim; iDim++)
                                        TotalVelocity[iDim] += Velocity[iDim];
                                    
                                    TotalAreaPressure         += Area*Pressure;
                                    TotalAreaDensity          += Area*Density;
                                    for (iDim = 0; iDim < nDim; iDim++)
                                        TotalAreaVelocity[iDim] += Area*Velocity[iDim];
                                    
                                    TotalMassPressure         += Area*(Density*TurboVelocity[0] )*Pressure;
                                    TotalMassDensity          += Area*(Density*TurboVelocity[0] )*Density;
                                    for (iDim = 0; iDim < nDim; iDim++)
                                        TotalMassVelocity[iDim] += Area*(Density*TurboVelocity[0] )*Velocity[iDim];
                                    
                                    TotalFluxes[0]      += Area*(Density*TurboVelocity[0]);
                                    TotalFluxes[1]      += Area*(Density*TurboVelocity[0]*TurboVelocity[0] + Pressure);
                                    for (iDim = 2; iDim < nDim+1; iDim++)
                                        TotalFluxes[iDim] += Area*(Density*TurboVelocity[0]*TurboVelocity[iDim -1]);
                                    TotalFluxes[nDim+1] += Area*(Density*TurboVelocity[0]*Enthalpy);
                                    
                                    
                                    /*--- Compute turbulent integral quantities for the boundary of interest ---*/
                                    
                                    if(turbulent){
                                        if(menter_sst){
                                            Kine  = solver[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
                                            Omega = solver[TURB_SOL]->GetNodes()->GetSolution(iPoint,1);
                                        }
                                        if(spalart_allmaras){
                                            Nu    = solver[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);
                                        }
                                        
                                        TotalKine   += Kine;
                                        TotalOmega  += Omega;
                                        TotalNu     += Nu;
                                        
                                        TotalAreaKine   += Area*Kine;
                                        TotalAreaOmega  += Area*Omega;
                                        TotalAreaNu     += Area*Nu;
                                        
                                        TotalMassKine    += Area*(Density*TurboVelocity[0] )*Kine;
                                        TotalMassOmega   += Area*(Density*TurboVelocity[0] )*Omega;
                                        TotalMassNu      += Area*(Density*TurboVelocity[0] )*Nu;
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
#ifdef HAVE_MPI
        
        /*--- Add information using all the nodes ---*/
        
        auto Allreduce = [](su2double x) {
            su2double tmp = x; x = 0.0;
            SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return x;
        };
        
        TotalDensity = Allreduce(TotalDensity);
        TotalPressure = Allreduce(TotalPressure);
        TotalAreaDensity = Allreduce(TotalAreaDensity);
        TotalAreaPressure = Allreduce(TotalAreaPressure);
        TotalMassDensity = Allreduce(TotalMassDensity);
        TotalMassPressure = Allreduce(TotalMassPressure);
        
        TotalNu = Allreduce(TotalNu);
        TotalKine = Allreduce(TotalKine);
        TotalOmega = Allreduce(TotalOmega);
        TotalAreaNu = Allreduce(TotalAreaNu);
        TotalAreaKine = Allreduce(TotalAreaKine);
        TotalAreaOmega = Allreduce(TotalAreaOmega);
        
        TotalMassNu = Allreduce(TotalMassNu);
        TotalMassKine = Allreduce(TotalMassKine);
        TotalMassOmega = Allreduce(TotalMassOmega);
        
        su2double* buffer = new su2double[max(nVar,nDim)];
        
        auto Allreduce_inplace = [buffer](int size, su2double* x) {
            SU2_MPI::Allreduce(x, buffer, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for(int i=0; i<size; ++i) x[i] = buffer[i];
        };
        
        Allreduce_inplace(nVar, TotalFluxes);
        Allreduce_inplace(nDim, TotalVelocity);
        Allreduce_inplace(nDim, TotalAreaVelocity);
        Allreduce_inplace(nDim, TotalMassVelocity);
        
        delete [] buffer;
        
#endif
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
            for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
                if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                    if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                        
                        TotalArea           = geometry->GetSpanArea(iMarker,iSpan);
                        AverageTurboNormal  = geometry->GetAverageTurboNormal(iMarker,iSpan);
                        nVert               = geometry->GetnTotVertexSpan(iMarker,iSpan);
                        
                        /*--- compute normal Mach number as a check for massflow average and mixedout average ---*/
                        GetFluidModel()->SetTDState_Prho(TotalAreaPressure/TotalArea, TotalAreaDensity / TotalArea);
                        soundSpeed = GetFluidModel()->GetSoundSpeed();
                        MachTest   = TotalFluxes[0]/(TotalAreaDensity*soundSpeed);
                        
                        /*--- Compute the averaged value for the boundary of interest for the span of interest ---*/
                        
                        /*--- compute algebraic average ---*/
                        avgDensity        = TotalDensity / nVert;
                        avgPressure       = TotalPressure / nVert;
                        for (iDim = 0; iDim < nDim; iDim++) avgVelocity[iDim] = TotalVelocity[iDim] / nVert;
                        avgKine           = TotalKine/nVert;
                        avgOmega          = TotalOmega/nVert;
                        avgNu             = TotalNu/nVert;
                        
                        /*--- compute area average ---*/
                        avgAreaDensity     = TotalAreaDensity / TotalArea;
                        avgAreaPressure    = TotalAreaPressure / TotalArea;
                        for (iDim = 0; iDim < nDim; iDim++) avgAreaVelocity[iDim] = TotalAreaVelocity[iDim] / TotalArea;
                        avgAreaKine        = TotalAreaKine / TotalArea;
                        avgAreaOmega       = TotalAreaOmega / TotalArea;
                        avgAreaNu          = TotalAreaNu / TotalArea;
                        
                        /*--- compute mass-flow average ---*/
                        if (abs(MachTest)< config->GetAverageMachLimit()) {
                            avgMassDensity   = avgAreaDensity;
                            avgMassPressure  = avgAreaPressure;
                            for (iDim = 0; iDim < nDim; iDim++) avgMassVelocity[iDim] = avgAreaVelocity[iDim];
                            avgMassKine      = avgAreaKine;
                            avgMassOmega     = avgAreaOmega;
                            avgMassNu        = avgAreaNu;
                        }else{
                            avgMassDensity     = TotalMassDensity / TotalFluxes[0];
                            avgMassPressure    = TotalMassPressure / TotalFluxes[0];
                            for (iDim = 0; iDim < nDim; iDim++) avgMassVelocity[iDim] = TotalMassVelocity[iDim] / TotalFluxes[0];
                            avgMassKine        = TotalMassKine / TotalFluxes[0];
                            avgMassOmega       = TotalMassOmega / TotalFluxes[0];
                            avgMassNu          = TotalMassNu / TotalFluxes[0];
                        }
                        /*--- compute mixed-out average ---*/
                        for (iVar = 0; iVar<nVar; iVar++){
                            AverageFlux[iMarker][iSpan][iVar]   = TotalFluxes[iVar]/TotalArea;
                            SpanTotalFlux[iMarker][iSpan][iVar] = TotalFluxes[iVar];
                        }
                        val_init_pressure = OldAveragePressure[iMarker][iSpan];
                        
                        if (abs(MachTest)< config->GetAverageMachLimit()) {
                            avgMixDensity    = avgAreaDensity;
                            avgMixPressure   = avgAreaPressure;
                            for (iDim = 0; iDim < nDim; iDim++)
                                avgMixVelocity[iDim] = avgAreaVelocity[iDim];
                            ComputeTurboVelocity(avgMixVelocity, AverageTurboNormal , avgMixTurboVelocity, marker_flag, config->GetKind_TurboMachinery(iZone));
                            avgMixKine       = avgAreaKine;
                            avgMixOmega      = avgAreaOmega;
                            avgMixNu         = avgAreaNu;
                        }else {
                            MixedOut_Average (config, val_init_pressure, AverageFlux[iMarker][iSpan], AverageTurboNormal, avgMixPressure, avgMixDensity);
                            avgMixTurboVelocity[0]         = ( AverageFlux[iMarker][iSpan][1] - avgMixPressure) / AverageFlux[iMarker][iSpan][0];
                            for (iDim = 2; iDim < nDim +1;iDim++)
                                avgMixTurboVelocity[iDim-1]  = AverageFlux[iMarker][iSpan][iDim] / AverageFlux[iMarker][iSpan][0];
                            
                            if (avgMixDensity!= avgMixDensity || avgMixPressure!= avgMixPressure || avgMixPressure < 0.0 || avgMixDensity < 0.0 ){
                                val_init_pressure = avgAreaPressure;
                                MixedOut_Average (config, val_init_pressure, AverageFlux[iMarker][iSpan], AverageTurboNormal, avgMixPressure, avgMixDensity);
                                avgMixTurboVelocity[0]          = ( AverageFlux[iMarker][iSpan][1] - avgMixPressure) / AverageFlux[iMarker][iSpan][0];
                                for (iDim = 2; iDim < nDim +1;iDim++)
                                    avgMixTurboVelocity[iDim-1]   = AverageFlux[iMarker][iSpan][iDim] / AverageFlux[iMarker][iSpan][0];
                            }
                            avgMixKine       = avgMassKine;
                            avgMixOmega      = avgMassOmega;
                            avgMixNu         = avgMassNu;
                        }
                        
                        /*--- Store averaged value for the selected average method ---*/
                        switch(average_process){
                            case ALGEBRAIC:
                                AverageDensity[iMarker][iSpan]  = avgDensity;
                                AveragePressure[iMarker][iSpan] = avgPressure;
                                ComputeTurboVelocity(avgVelocity, AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                AverageKine[iMarker][iSpan]     = avgKine;
                                AverageOmega[iMarker][iSpan]    = avgOmega;
                                AverageNu[iMarker][iSpan]       = avgNu;
                                break;
                                
                            case AREA:
                                AverageDensity[iMarker][iSpan]  = avgAreaDensity;
                                AveragePressure[iMarker][iSpan] = avgAreaPressure;
                                ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                AverageKine[iMarker][iSpan]     = avgAreaKine;
                                AverageOmega[iMarker][iSpan]    = avgAreaOmega;
                                AverageNu[iMarker][iSpan]       = avgAreaNu;
                                break;
                                
                            case MASSFLUX:
                                AverageDensity[iMarker][iSpan]  = avgMassDensity;
                                AveragePressure[iMarker][iSpan] = avgMassPressure;
                                ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , AverageTurboVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                AverageKine[iMarker][iSpan]     = avgMassKine;
                                AverageOmega[iMarker][iSpan]    = avgMassOmega;
                                AverageNu[iMarker][iSpan]       = avgMassNu;
                                break;
                                
                            case MIXEDOUT:
                                AverageDensity[iMarker][iSpan] = avgMixDensity;
                                AveragePressure[iMarker][iSpan] = avgMixPressure;
                                for (iDim = 0; iDim < nDim; iDim++) AverageTurboVelocity[iMarker][iSpan][iDim] = avgMixTurboVelocity[iDim];
                                AverageKine[iMarker][iSpan]     = avgMixKine;
                                AverageOmega[iMarker][iSpan]    = avgMixOmega;
                                AverageNu[iMarker][iSpan]       = avgMixNu;
                                break;
                                
                            default:
                                SU2_MPI::Error(" Invalid AVERAGE PROCESS input!", CURRENT_FUNCTION);
                                break;
                        }
                        
                        /* --- check if averaged quantities are correct otherwise reset the old quantities ---*/
                        if (AverageDensity[iMarker][iSpan]!= AverageDensity[iMarker][iSpan] || AveragePressure[iMarker][iSpan]!= AveragePressure[iMarker][iSpan]){
                            cout<<"nan in mixing process routine for iSpan: " << iSpan<< " in marker " << config->GetMarker_All_TagBound(iMarker)<< endl;
                            AverageDensity[iMarker][iSpan]               = OldAverageDensity[iMarker][iSpan];
                            AveragePressure[iMarker][iSpan]              = OldAveragePressure[iMarker][iSpan];
                            for(iDim = 0; iDim < nDim;iDim++)
                                AverageTurboVelocity[iMarker][iSpan][iDim] = OldAverageTurboVelocity[iMarker][iSpan][iDim];
                        }
                        
                        
                        if (AverageDensity[iMarker][iSpan] < 0.0 || AveragePressure[iMarker][iSpan] < 0.0){
                            cout << " negative density or pressure in mixing process routine for iSpan: " << iSpan<< " in marker " << config->GetMarker_All_TagBound(iMarker)<< endl;
                            AverageDensity[iMarker][iSpan]               = OldAverageDensity[iMarker][iSpan];
                            AveragePressure[iMarker][iSpan]              = OldAveragePressure[iMarker][iSpan];
                            for(iDim = 0; iDim < nDim;iDim++)
                                AverageTurboVelocity[iMarker][iSpan][iDim] = OldAverageTurboVelocity[iMarker][iSpan][iDim];
                        }
                        
                        /* --- update old average solution ---*/
                        OldAverageDensity[iMarker][iSpan]               = AverageDensity[iMarker][iSpan];
                        OldAveragePressure[iMarker][iSpan]              = AveragePressure[iMarker][iSpan];
                        for(iDim = 0; iDim < nDim;iDim++)
                            OldAverageTurboVelocity[iMarker][iSpan][iDim] = AverageTurboVelocity[iMarker][iSpan][iDim];
                        
                        /*--- to avoid back flow ---*/
                        if (AverageTurboVelocity[iMarker][iSpan][0] < 0.0){
                            AverageTurboVelocity[iMarker][iSpan][0]       = soundSpeed*config->GetAverageMachLimit();
                        }
                        
                        /*--- compute cartesian average Velocity ---*/
                        ComputeBackVelocity(AverageTurboVelocity[iMarker][iSpan], AverageTurboNormal , AverageVelocity[iMarker][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                        
                        
                        /*--- Store averaged performance value for the selected average method ---*/
                        switch(performance_average_process){
                            case ALGEBRAIC:
                                if(marker_flag == INFLOW){
                                    DensityIn[iMarkerTP - 1][iSpan]   = avgDensity;
                                    PressureIn[iMarkerTP - 1][iSpan]  = avgPressure;
                                    ComputeTurboVelocity(avgVelocity, AverageTurboNormal , TurboVelocityIn[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                    KineIn[iMarkerTP - 1][iSpan]      = avgKine;
                                    OmegaIn[iMarkerTP - 1][iSpan]     = avgOmega;
                                    NuIn[iMarkerTP - 1][iSpan]        = avgNu;
                                }
                                else{
                                    DensityOut[iMarkerTP - 1][iSpan]  = avgDensity;
                                    PressureOut[iMarkerTP - 1][iSpan] = avgPressure;
                                    ComputeTurboVelocity(avgVelocity, AverageTurboNormal , TurboVelocityOut[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                    KineOut[iMarkerTP - 1][iSpan]     = avgKine;
                                    OmegaOut[iMarkerTP - 1][iSpan]    = avgOmega;
                                    NuOut[iMarkerTP - 1][iSpan]       = avgNu;
                                }
                                
                                break;
                            case AREA:
                                if(marker_flag == INFLOW){
                                    DensityIn[iMarkerTP - 1][iSpan]   = avgAreaDensity;
                                    PressureIn[iMarkerTP - 1][iSpan]  = avgAreaPressure;
                                    ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , TurboVelocityIn[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                    KineIn[iMarkerTP - 1][iSpan]      = avgAreaKine;
                                    OmegaIn[iMarkerTP - 1][iSpan]     = avgAreaOmega;
                                    NuIn[iMarkerTP - 1][iSpan]        = avgAreaNu;
                                }
                                else{
                                    DensityOut[iMarkerTP - 1][iSpan]  = avgAreaDensity;
                                    PressureOut[iMarkerTP - 1][iSpan] = avgAreaPressure;
                                    ComputeTurboVelocity(avgAreaVelocity, AverageTurboNormal , TurboVelocityOut[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                    KineOut[iMarkerTP - 1][iSpan]     = avgAreaKine;
                                    OmegaOut[iMarkerTP - 1][iSpan]    = avgAreaOmega;
                                    NuOut[iMarkerTP - 1][iSpan]       = avgAreaNu/TotalArea;
                                }
                                break;
                                
                            case MASSFLUX:
                                if(marker_flag == INFLOW){
                                    DensityIn[iMarkerTP - 1][iSpan]   = avgMassDensity;
                                    PressureIn[iMarkerTP - 1][iSpan]  = avgMassPressure;
                                    ComputeTurboVelocity(avgMassVelocity, AverageTurboNormal , TurboVelocityIn[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                    KineIn[iMarkerTP - 1][iSpan]      = avgMassKine;
                                    OmegaIn[iMarkerTP - 1][iSpan]     = avgMassOmega;
                                    NuIn[iMarkerTP - 1][iSpan]        = avgMassNu;
                                }
                                else{
                                    DensityOut[iMarkerTP - 1][iSpan]  = avgMassDensity;
                                    PressureOut[iMarkerTP - 1][iSpan] = avgMassPressure;
                                    ComputeTurboVelocity(avgMassVelocity, AverageTurboNormal , TurboVelocityOut[iMarkerTP -1][iSpan], marker_flag, config->GetKind_TurboMachinery(iZone));
                                    KineOut[iMarkerTP - 1][iSpan]     = avgMassKine;
                                    OmegaOut[iMarkerTP - 1][iSpan]    = avgMassOmega;
                                    NuOut[iMarkerTP - 1][iSpan]       = avgMassNu;
                                }
                                
                                break;
                                
                            case MIXEDOUT:
                                if (marker_flag == INFLOW){
                                    DensityIn[iMarkerTP - 1][iSpan]  = avgMixDensity;
                                    PressureIn[iMarkerTP - 1][iSpan] = avgMixPressure;
                                    for (iDim = 0; iDim < nDim; iDim++) TurboVelocityIn[iMarkerTP -1][iSpan][iDim] = avgMixTurboVelocity[iDim];
                                    KineIn[iMarkerTP - 1][iSpan]     = avgMixKine;
                                    OmegaIn[iMarkerTP - 1][iSpan]    = avgMixOmega;
                                    NuIn[iMarkerTP - 1][iSpan]       = avgMixNu;
                                }
                                else{
                                    DensityOut[iMarkerTP - 1][iSpan]  = avgMixDensity;
                                    PressureOut[iMarkerTP - 1][iSpan] = avgMixPressure;
                                    for (iDim = 0; iDim < nDim; iDim++) TurboVelocityOut[iMarkerTP -1][iSpan][iDim] = avgMixTurboVelocity[iDim];
                                    KineOut[iMarkerTP - 1][iSpan]     = avgMixKine;
                                    OmegaOut[iMarkerTP - 1][iSpan]    = avgMixOmega;
                                    NuOut[iMarkerTP - 1][iSpan]       = avgMixNu;
                                }
                                break;
                                
                            default:
                                SU2_MPI::Error(" Invalid MIXING_PROCESS input!", CURRENT_FUNCTION);
                                break;
                        }
                    }
                }
            }
        }
    }
    
    /*--- Compute Outlet Static Pressure if Radial equilibrium is imposed ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
            if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
                    Marker_Tag         = config->GetMarker_All_TagBound(iMarker);
                    if(config->GetBoolGiles() || config->GetBoolRiemann()){
                        if(config->GetBoolRiemann()){
                            if(config->GetKind_Data_Riemann(Marker_Tag) == RADIAL_EQUILIBRIUM){
                                RadialEquilibriumPressure[iMarker][nSpanWiseSections/2] = config->GetRiemann_Var1(Marker_Tag)/config->GetPressure_Ref();
                                for (iSpan= nSpanWiseSections/2; iSpan < nSpanWiseSections-1; iSpan++){
                                    Radius2    = geometry->GetTurboRadius(iMarker,iSpan+1);
                                    Radius1    = geometry->GetTurboRadius(iMarker,iSpan);
                                    Vt2        = AverageTurboVelocity[iMarker][iSpan +1][1]*AverageTurboVelocity[iMarker][iSpan +1][1];
                                    RadialEquilibriumPressure[iMarker][iSpan +1] =  RadialEquilibriumPressure[iMarker][iSpan] + AverageDensity[iMarker][iSpan +1]*Vt2/Radius2*(Radius2 - Radius1);
                                }
                                for (iSpan= nSpanWiseSections/2; iSpan > 0; iSpan--){
                                    Radius2    = geometry->GetTurboRadius(iMarker,iSpan);
                                    Radius1    = geometry->GetTurboRadius(iMarker,iSpan-1);
                                    Vt2        = AverageTurboVelocity[iMarker][iSpan - 1][1]*AverageTurboVelocity[iMarker][iSpan - 1][1];
                                    RadialEquilibriumPressure[iMarker][iSpan -1] =  RadialEquilibriumPressure[iMarker][iSpan] - AverageDensity[iMarker][iSpan -1]*Vt2/Radius2*(Radius2 - Radius1);
                                }
                            }
                        }
                        else{
                            if(config->GetKind_Data_Giles(Marker_Tag) == RADIAL_EQUILIBRIUM){
                                RadialEquilibriumPressure[iMarker][nSpanWiseSections/2] = config->GetGiles_Var1(Marker_Tag)/config->GetPressure_Ref();
                                for (iSpan= nSpanWiseSections/2; iSpan < nSpanWiseSections-1; iSpan++){
                                    Radius2    = geometry->GetTurboRadius(iMarker,iSpan+1);
                                    Radius1    = geometry->GetTurboRadius(iMarker,iSpan);
                                    Vt2        = AverageTurboVelocity[iMarker][iSpan +1][1]*AverageTurboVelocity[iMarker][iSpan +1][1];
                                    RadialEquilibriumPressure[iMarker][iSpan +1] =  RadialEquilibriumPressure[iMarker][iSpan] + AverageDensity[iMarker][iSpan +1]*Vt2/Radius2*(Radius2 - Radius1);
                                }
                                for (iSpan= nSpanWiseSections/2; iSpan > 0; iSpan--){
                                    Radius2    = geometry->GetTurboRadius(iMarker,iSpan);
                                    Radius1    = geometry->GetTurboRadius(iMarker,iSpan-1);
                                    Vt2        = AverageTurboVelocity[iMarker][iSpan -1][1]*AverageTurboVelocity[iMarker][iSpan - 1][1];
                                    RadialEquilibriumPressure[iMarker][iSpan -1] =  RadialEquilibriumPressure[iMarker][iSpan] - AverageDensity[iMarker][iSpan -1]*Vt2/Radius2*(Radius2 - Radius1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /*--- Free locally allocated memory ---*/
    delete [] Velocity;
    delete [] UnitNormal;
    delete [] TurboNormal;
    delete [] TurboVelocity;
    delete [] TotalVelocity;
    delete [] TotalAreaVelocity;
    delete [] TotalFluxes;
    delete [] TotalMassVelocity;
    delete [] avgVelocity;
    delete [] avgAreaVelocity;
    delete [] avgMassVelocity;
    delete [] avgMixVelocity;
    delete [] avgMixTurboVelocity;
    
}

void CEulerSolver::MixedOut_Average (CConfig *config, su2double val_init_pressure, const su2double *val_Averaged_Flux,
                                     const su2double *val_normal, su2double& pressure_mix, su2double& density_mix) {
    
    su2double dx, f, df, resdl = 1.0E+05;
    unsigned short iter = 0, iDim;
    su2double relax_factor = config->GetMixedout_Coeff(0);
    su2double toll = config->GetMixedout_Coeff(1);
    unsigned short maxiter = SU2_TYPE::Int(config->GetMixedout_Coeff(2));
    su2double dhdP, dhdrho, enthalpy_mix, velsq, *vel;
    
    vel = new su2double[nDim];
    
    pressure_mix = val_init_pressure;
    
    /*--- Newton-Raphson's method with central difference formula ---*/
    
    while ( iter <= maxiter ) {
        
        density_mix = val_Averaged_Flux[0]*val_Averaged_Flux[0]/(val_Averaged_Flux[1] - pressure_mix);
        GetFluidModel()->SetTDState_Prho(pressure_mix, density_mix);
        enthalpy_mix = GetFluidModel()->GetStaticEnergy() + (pressure_mix)/(density_mix);
        
        GetFluidModel()->ComputeDerivativeNRBC_Prho(pressure_mix, density_mix);
        dhdP   = GetFluidModel()->GetdhdP_rho();
        dhdrho = GetFluidModel()->Getdhdrho_P();
        
        vel[0]  = (val_Averaged_Flux[1] - pressure_mix) / val_Averaged_Flux[0];
        for (iDim = 1; iDim < nDim; iDim++) {
            vel[iDim]  = val_Averaged_Flux[iDim+1] / val_Averaged_Flux[0];
        }
        
        velsq = 0.0;
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
            velsq += vel[iDim]*vel[iDim];
        }
        f = val_Averaged_Flux[nDim+1] - val_Averaged_Flux[0]*(enthalpy_mix + velsq/2);
        df = -val_Averaged_Flux[0]*(dhdP - 1/density_mix) - dhdrho*density_mix*density_mix/val_Averaged_Flux[0];
        dx = -f/df;
        resdl = dx/val_init_pressure;
        pressure_mix += relax_factor*dx;
        
        iter += 1;
        if ( abs(resdl) <= toll ) {
            break;
        }
        
    }
    
    density_mix = val_Averaged_Flux[0]*val_Averaged_Flux[0]/(val_Averaged_Flux[1] - pressure_mix);
    
    delete [] vel;
    
}

void CEulerSolver::GatherInOutAverageValues(CConfig *config, CGeometry *geometry){
    
    unsigned short iMarker, iMarkerTP;
    unsigned short iSpan;
    int markerTP;
    su2double     densityIn, pressureIn, normalVelocityIn, tangVelocityIn, radialVelocityIn;
    su2double     densityOut, pressureOut, normalVelocityOut, tangVelocityOut, radialVelocityOut;
    su2double     kineIn, omegaIn, nuIn, kineOut, omegaOut, nuOut;
    //TODO (turbo) implement interpolation so that Inflow and Outflow spanwise section can be different
    
    for (iSpan= 0; iSpan < nSpanWiseSections + 1 ; iSpan++) {
        
#ifdef HAVE_MPI
        unsigned short i, n1, n2, n1t,n2t;
        su2double *TurbPerfIn= NULL,*TurbPerfOut= NULL;
        su2double *TotTurbPerfIn = NULL,*TotTurbPerfOut = NULL;
        int *TotMarkerTP = NULL;
        
        n1          = 8;
        n2          = 8;
        n1t         = n1*size;
        n2t         = n2*size;
        TurbPerfIn  = new su2double[n1];
        TurbPerfOut = new su2double[n2];
        
        for (i=0;i<n1;i++)
            TurbPerfIn[i]    = -1.0;
        for (i=0;i<n2;i++)
            TurbPerfOut[i]   = -1.0;
#endif
        
        densityIn            = -1.0;
        pressureIn           = -1.0;
        normalVelocityIn     = -1.0;
        tangVelocityIn       = -1.0;
        radialVelocityIn     = -1.0;
        densityOut           = -1.0;
        pressureOut          = -1.0;
        normalVelocityOut    = -1.0;
        tangVelocityOut      = -1.0;
        radialVelocityOut    = -1.0;
        kineIn               = -1.0;
        omegaIn              = -1.0;
        nuIn                 = -1.0;
        kineOut              = -1.0;
        omegaOut             = -1.0;
        nuOut                = -1.0;
        
        markerTP             = -1;
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
            for (iMarkerTP = 1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
                if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
                    if (config->GetMarker_All_TurbomachineryFlag(iMarker) == INFLOW){
                        markerTP            = iMarkerTP;
                        densityIn           = DensityIn[iMarkerTP -1][iSpan];
                        pressureIn          = PressureIn[iMarkerTP -1][iSpan];
                        normalVelocityIn    = TurboVelocityIn[iMarkerTP -1][iSpan][0];
                        tangVelocityIn      = TurboVelocityIn[iMarkerTP -1][iSpan][1];
                        if (nDim ==3){
                            radialVelocityIn  = TurboVelocityIn[iMarkerTP -1][iSpan][2];
                        }
                        kineIn              = KineIn[iMarkerTP -1][iSpan];
                        omegaIn             = OmegaIn[iMarkerTP -1][iSpan];
                        nuIn                = NuIn[iMarkerTP -1][iSpan];
                        
#ifdef HAVE_MPI
                        TurbPerfIn[0]  = densityIn;
                        TurbPerfIn[1]  = pressureIn;
                        TurbPerfIn[2]  = normalVelocityIn;
                        TurbPerfIn[3]  = tangVelocityIn;
                        TurbPerfIn[4]  = radialVelocityIn;
                        TurbPerfIn[5]  = kineIn;
                        TurbPerfIn[6]  = omegaIn;
                        TurbPerfIn[7]  = nuIn;
#endif
                    }
                    
                    /*--- retrieve outlet information ---*/
                    if (config->GetMarker_All_TurbomachineryFlag(iMarker) == OUTFLOW){
                        densityOut           = DensityOut[iMarkerTP -1][iSpan];
                        pressureOut          = PressureOut[iMarkerTP -1][iSpan];
                        normalVelocityOut    = TurboVelocityOut[iMarkerTP -1][iSpan][0];
                        tangVelocityOut      = TurboVelocityOut[iMarkerTP -1][iSpan][1];
                        if (nDim ==3){
                            radialVelocityOut  = TurboVelocityOut[iMarkerTP -1][iSpan][2];
                        }
                        kineOut              = KineOut[iMarkerTP -1][iSpan];
                        omegaOut             = OmegaOut[iMarkerTP -1][iSpan];
                        nuOut                = NuOut[iMarkerTP -1][iSpan];
                        
#ifdef HAVE_MPI
                        TurbPerfOut[0]  = densityOut;
                        TurbPerfOut[1]  = pressureOut;
                        TurbPerfOut[2]  = normalVelocityOut;
                        TurbPerfOut[3]  = tangVelocityOut;
                        TurbPerfOut[4]  = radialVelocityOut;
                        TurbPerfOut[5]  = kineOut;
                        TurbPerfOut[6]  = omegaOut;
                        TurbPerfOut[7]  = nuOut;
#endif
                    }
                }
            }
        }
        
#ifdef HAVE_MPI
        if (rank == MASTER_NODE){
            TotTurbPerfIn       = new su2double[n1t];
            TotTurbPerfOut      = new su2double[n2t];
            for (i=0;i<n1t;i++)
                TotTurbPerfIn[i]  = -1.0;
            for (i=0;i<n2t;i++)
                TotTurbPerfOut[i] = -1.0;
            TotMarkerTP = new int[size];
            for(i=0; i<size; i++){
                TotMarkerTP[i]    = -1;
            }
        }
        SU2_MPI::Gather(TurbPerfIn, n1, MPI_DOUBLE, TotTurbPerfIn, n1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Gather(TurbPerfOut, n2, MPI_DOUBLE,TotTurbPerfOut, n2, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Gather(&markerTP, 1, MPI_INT,TotMarkerTP, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
        if (rank == MASTER_NODE){
            delete [] TurbPerfIn, delete [] TurbPerfOut;
        }
        
        if (rank == MASTER_NODE){
            for (i=0;i<size;i++){
                if(TotTurbPerfIn[n1*i] > 0.0){
                    densityIn        = TotTurbPerfIn[n1*i];
                    pressureIn       = TotTurbPerfIn[n1*i+1];
                    normalVelocityIn = TotTurbPerfIn[n1*i+2];
                    tangVelocityIn   = TotTurbPerfIn[n1*i+3];
                    radialVelocityIn = TotTurbPerfIn[n1*i+4];
                    kineIn           = TotTurbPerfIn[n1*i+5];
                    omegaIn          = TotTurbPerfIn[n1*i+6];
                    nuIn             = TotTurbPerfIn[n1*i+7];
                    markerTP         = TotMarkerTP[i];
                }
                
                if(TotTurbPerfOut[n2*i] > 0.0){
                    densityOut        = TotTurbPerfOut[n1*i];
                    pressureOut       = TotTurbPerfOut[n1*i+1];
                    normalVelocityOut = TotTurbPerfOut[n1*i+2];
                    tangVelocityOut   = TotTurbPerfOut[n1*i+3];
                    radialVelocityOut = TotTurbPerfOut[n1*i+4];
                    kineOut           = TotTurbPerfOut[n1*i+5];
                    omegaOut          = TotTurbPerfOut[n1*i+6];
                    nuOut             = TotTurbPerfOut[n1*i+7];
                }
            }
            
            delete [] TotTurbPerfIn, delete [] TotTurbPerfOut; delete [] TotMarkerTP;
        }
        
#endif
        
        if (rank == MASTER_NODE && markerTP > -1){
            /*----Quantities needed for computing the turbomachinery performance -----*/
            DensityIn[markerTP -1][iSpan]              = densityIn;
            PressureIn[markerTP -1][iSpan]             = pressureIn;
            TurboVelocityIn[markerTP -1][iSpan][0]     = normalVelocityIn;
            TurboVelocityIn[markerTP -1][iSpan][1]     = tangVelocityIn;
            if (nDim == 3)
                TurboVelocityIn[markerTP -1][iSpan][2]   = radialVelocityIn;
            KineIn[markerTP -1][iSpan]                 = kineIn;
            OmegaIn[markerTP -1][iSpan]                = omegaIn;
            NuIn[markerTP -1][iSpan]                   = nuIn;
            
            DensityOut[markerTP -1][iSpan]             = densityOut;
            PressureOut[markerTP -1][iSpan]            = pressureOut;
            TurboVelocityOut[markerTP -1][iSpan][0]    = normalVelocityOut;
            TurboVelocityOut[markerTP -1][iSpan][1]    = tangVelocityOut;
            if (nDim == 3)
                TurboVelocityOut[markerTP -1][iSpan][2]  = radialVelocityOut;
            KineOut[markerTP -1][iSpan]                = kineOut;
            OmegaOut[markerTP -1][iSpan]               = omegaOut;
            NuOut[markerTP -1][iSpan]                  = nuOut;
        }
    }
}
